#!/usr/bin/env python3
"""Train and evaluate a SILVA 138.2 V4 Naive Bayes classifier.

The Article 15 contract separates three claims:

1. a deterministic, known-lineage holdout measures interpolation inside one
   SILVA release;
2. a classifier refitted on every region-matched reference sequence is the
   production deliverable; and
3. concordance on real ASVs compares methods, not accuracy against truth.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import tempfile
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Sequence


QIIME_ENV = "microbiome-qiime2-2026.4"
QIIME_VERSION = "2026.4.0"
Q2_FEATURE_CLASSIFIER_VERSION = "2026.4.0"
RESCRIPT_VERSION = "2026.4.0"
SKLEARN_VERSION = "1.7.1"
JOBLIB_VERSION = "1.5.3"
NUMPY_VERSION = "2.4.2"
SCIPY_VERSION = "1.17.1"

REFERENCE_RECORDS = 262_752
TRAINING_RECORDS = 236_516
HOLDOUT_RECORDS = 26_236
QUERY_ASVS = 366
QUERY_READS = 5_213
SPLIT_SEED = "article15-silva1382-v4-known-lineage-v1"

REFERENCE_FASTA_SHA256 = (
    "3a251d263bd4e12c96023c84b3f2e255b8bcb0d5865dc6f299f9918e6be3a808"
)
REFERENCE_TAXONOMY_SHA256 = (
    "f49320aaa32aa70af5c5a48649786e3bb4177868da5da4d2143d0a54c5fa47b0"
)
QUERY_FASTA_SHA256 = (
    "fa71915be8007d36ec70f85d28401f5a94d0d33eca323782cc937cfc4bcf1af4"
)
QUERY_ABUNDANCE_SHA256 = (
    "a80e3b94d8632c48fb6c01e788e842c022c3b1ee64b79dc94456618328895e6e"
)

CONFIDENCES = (0.50, 0.70, 0.90)
CONFIDENCE_KEYS = {0.50: "c050", 0.70: "c070", 0.90: "c090"}
CLASSIFIER_JOBS = 4
READS_PER_BATCH = 20_000
RANKS = (
    ("domain", "d__"),
    ("phylum", "p__"),
    ("class", "c__"),
    ("order", "o__"),
    ("family", "f__"),
    ("genus", "g__"),
)
COLORS = {"c050": "#0072B2", "c070": "#009E73", "c090": "#D55E00"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument(
        "--sequence-search-taxonomy",
        type=Path,
        help=(
            "Optional existing SILVA sequence-search taxonomy Artifact. "
            "When omitted, the comparator is recomputed from the bundled "
            "reference and query inside this run."
        ),
    )
    parser.add_argument("--qiime-env", default=QIIME_ENV)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(parts: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for part in parts:
        digest.update(part.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_tsv(
    path: Path,
    rows: Iterable[dict[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_tsv(
    path: Path,
    rows: Iterable[dict[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(
                    text,
                    fieldnames=list(fields),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    identifier: str | None = None
    pieces: list[str] = []
    with opener(path, "rt", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(pieces).upper()
                identifier = line[1:].split(maxsplit=1)[0]
                pieces = []
            else:
                if identifier is None:
                    raise ValueError(f"sequence before FASTA header: {path}")
                pieces.append(line)
    if identifier is not None:
        yield identifier, "".join(pieces).upper()


def read_taxonomy(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    result: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"taxonomy has no header: {path}")
        id_field = next(
            (
                field
                for field in reader.fieldnames
                if field.strip().lower() in {"feature id", "feature-id", "featureid"}
            ),
            None,
        )
        taxon_field = next(
            (
                field
                for field in reader.fieldnames
                if field.strip().lower() == "taxon"
            ),
            None,
        )
        if id_field is None or taxon_field is None:
            raise ValueError(f"unexpected taxonomy columns: {reader.fieldnames}")
        for row in reader:
            feature_id = str(row[id_field]).strip()
            if feature_id in result:
                raise ValueError(f"duplicate taxonomy feature ID: {feature_id}")
            result[feature_id] = str(row[taxon_field]).strip()
    return result


def read_query_abundance(path: Path) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header != ["Feature ID", "Frequency", "No. of Samples Observed In"]:
            raise ValueError(f"unexpected query abundance header: {header}")
        for row in reader:
            if row and row[0] == "#q2:types":
                continue
            if len(row) != 3:
                raise ValueError(f"unexpected query abundance row: {row}")
            result[row[0]] = {
                "frequency": int(float(row[1])),
                "samples_observed": int(float(row[2])),
            }
    return result


def rank_map(taxon: str) -> dict[str, str]:
    if not taxon or taxon.strip().lower() == "unassigned":
        return {}
    result: dict[str, str] = {}
    for item in taxon.split(";"):
        label = item.strip()
        for rank, prefix in RANKS:
            if label.startswith(prefix):
                value = label[len(prefix) :].strip()
                if value and value.lower() not in {"unassigned", "uncultured"}:
                    result[rank] = value
                break
    return result


def genus_lineage(taxon: str) -> str:
    parts = [part.strip() for part in taxon.split(";")]
    parts = (parts + [""] * 6)[:6]
    return "; ".join(parts)


def parse_classification(path: Path) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"classification has no header: {path}")
        for row in reader:
            feature_id = str(
                row.get("Feature ID", row.get("FeatureID", ""))
            ).strip()
            taxon = str(row.get("Taxon", "")).strip()
            confidence_raw = str(
                row.get("Confidence", row.get("Consensus", ""))
            ).strip()
            confidence = (
                float(confidence_raw)
                if confidence_raw and confidence_raw.lower() != "nan"
                else 0.0
            )
            result[feature_id] = {
                "taxon": taxon,
                "confidence": confidence,
                "ranks": rank_map(taxon),
            }
    return result


def unpack_gzip(source: Path, destination: Path) -> None:
    with gzip.open(source, "rb") as input_handle, destination.open("wb") as output:
        shutil.copyfileobj(input_handle, output, length=1024 * 1024)


def write_fasta(
    path: Path,
    identifiers: Sequence[str],
    sequences: dict[str, str],
) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for feature_id in identifiers:
            handle.write(f">{feature_id}\n{sequences[feature_id]}\n")


def write_taxonomy(
    path: Path,
    identifiers: Sequence[str],
    taxonomy: dict[str, str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Feature ID", "Taxon"))
        for feature_id in identifiers:
            writer.writerow((feature_id, taxonomy[feature_id]))


def deterministic_split(
    sequences: dict[str, str],
    taxonomy: dict[str, str],
) -> tuple[list[str], list[str], list[dict[str, object]], dict[str, object]]:
    groups: defaultdict[str, list[str]] = defaultdict(list)
    for feature_id in sequences:
        groups[genus_lineage(taxonomy[feature_id])].append(feature_id)

    training: list[str] = []
    holdout: list[str] = []
    membership: list[dict[str, object]] = []
    singleton_lineages = 0

    for lineage in sorted(groups):
        identifiers = sorted(
            groups[lineage],
            key=lambda feature_id: hashlib.sha256(
                (
                    f"{SPLIT_SEED}\t{feature_id}\t"
                    f"{sequences[feature_id]}"
                ).encode("utf-8")
            ).hexdigest(),
        )
        if len(identifiers) == 1:
            holdout_count = 0
            singleton_lineages += 1
        else:
            holdout_count = min(
                len(identifiers) - 1,
                max(1, math.floor(len(identifiers) * 0.10)),
            )
        holdout_ids = set(identifiers[:holdout_count])
        for feature_id in identifiers:
            split = "holdout" if feature_id in holdout_ids else "training"
            split_hash = hashlib.sha256(
                (
                    f"{SPLIT_SEED}\t{feature_id}\t"
                    f"{sequences[feature_id]}"
                ).encode("utf-8")
            ).hexdigest()
            membership.append(
                {
                    "feature_id": feature_id,
                    "split": split,
                    "genus_lineage": lineage,
                    "lineage_records": len(identifiers),
                    "split_hash": split_hash,
                }
            )
            if split == "holdout":
                holdout.append(feature_id)
            else:
                training.append(feature_id)

    training.sort()
    holdout.sort()
    membership.sort(key=lambda row: str(row["feature_id"]))
    audit = {
        "lineages": len(groups),
        "singleton_lineages": singleton_lineages,
        "training_lineages": len(
            {genus_lineage(taxonomy[feature_id]) for feature_id in training}
        ),
        "holdout_lineages": len(
            {genus_lineage(taxonomy[feature_id]) for feature_id in holdout}
        ),
    }
    return training, holdout, membership, audit


def conda_command(environment: str, executable: str, *args: str) -> list[str]:
    return ["conda", "run", "-n", environment, executable, *args]


def isolated_environment(root: Path) -> dict[str, str]:
    environment = os.environ.copy()
    for name in ("xdg", "numba", "matplotlib"):
        (root / "cache" / name).mkdir(parents=True, exist_ok=True)
    environment.update(
        {
            "XDG_CACHE_HOME": str(root / "cache" / "xdg"),
            "NUMBA_CACHE_DIR": str(root / "cache" / "numba"),
            "MPLCONFIGDIR": str(root / "cache" / "matplotlib"),
            "PYTHONHASHSEED": "0",
            "LC_ALL": "C.UTF-8",
            "LANG": "C.UTF-8",
        }
    )
    return environment


class CommandRunner:
    def __init__(
        self,
        environment: dict[str, str],
        log_path: Path,
        replacements: dict[str, str],
    ) -> None:
        self.environment = environment
        self.log_path = log_path
        self.replacements = replacements
        self.records: list[dict[str, object]] = []

    def sanitize(self, value: str) -> str:
        result = value
        for source, replacement in sorted(
            self.replacements.items(), key=lambda item: len(item[0]), reverse=True
        ):
            result = result.replace(source, replacement)
        return result

    def flush(self) -> None:
        pieces: list[str] = []
        for record in self.records:
            command = " ".join(
                self.sanitize(str(part)) for part in record["command"]
            )
            pieces.extend(
                [
                    f"[{record['label']}]",
                    f"command: {command}",
                    f"returncode: {record['returncode']}",
                    f"elapsed_seconds: {record['elapsed_seconds']}",
                    "stdout:",
                    self.sanitize(str(record["stdout"])),
                    "stderr:",
                    self.sanitize(str(record["stderr"])),
                    "",
                ]
            )
        self.log_path.write_text("\n".join(pieces), encoding="utf-8")

    def run(
        self,
        label: str,
        command: Sequence[str],
        timeout: int = 14_400,
    ) -> dict[str, object]:
        started = time.monotonic()
        completed = subprocess.run(
            list(command),
            check=False,
            capture_output=True,
            text=True,
            env=self.environment,
            timeout=timeout,
        )
        record = {
            "label": label,
            "command": list(command),
            "returncode": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
            "elapsed_seconds": round(time.monotonic() - started, 3),
        }
        self.records.append(record)
        self.flush()
        if completed.returncode != 0:
            raise RuntimeError(
                f"{label} failed\nSTDOUT:\n{completed.stdout}\n"
                f"STDERR:\n{completed.stderr}"
            )
        return record

    def qiime(
        self,
        label: str,
        environment: str,
        *arguments: str,
        timeout: int = 14_400,
    ) -> dict[str, object]:
        return self.run(
            label,
            conda_command(environment, "qiime", *arguments),
            timeout=timeout,
        )


def artifact_metadata(path: Path) -> dict[str, str]:
    with zipfile.ZipFile(path) as archive:
        matches = [
            name
            for name in archive.namelist()
            if name.count("/") == 1 and name.endswith("/metadata.yaml")
        ]
        if len(matches) != 1:
            raise ValueError(f"{path.name}: root metadata.yaml is ambiguous")
        text = archive.read(matches[0]).decode("utf-8", errors="replace")
    result: dict[str, str] = {}
    for key in ("uuid", "type", "format"):
        match = re.search(rf"^{key}:\s*(.+)$", text, flags=re.MULTILINE)
        result[key] = match.group(1).strip() if match else ""
    return result


def provenance_actions(path: Path) -> list[tuple[str, str]]:
    actions: set[tuple[str, str]] = set()
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if "/provenance/" in name and name.endswith("/action/action.yaml")
        ]
        for member in members:
            text = archive.read(member).decode("utf-8", errors="replace")
            plugin = re.search(
                r"^[ \t]+plugin:[ \t]*(.+)$", text, flags=re.MULTILINE
            )
            action = re.search(
                r"^[ \t]+action:[ \t]*(.+)$", text, flags=re.MULTILINE
            )
            if plugin and action:
                plugin_value = plugin.group(1).strip().strip("'\"")
                reference = re.search(
                    r"environment:plugins:([^'\"]+)", plugin_value
                )
                if reference:
                    plugin_value = reference.group(1)
                actions.add(
                    (plugin_value, action.group(1).strip().strip("'\""))
                )
    return sorted(actions)


def export_artifact(
    runner: CommandRunner,
    environment: str,
    artifact: Path,
    output_dir: Path,
    label: str,
) -> Path:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    runner.qiime(
        label,
        environment,
        "tools",
        "export",
        "--input-path",
        str(artifact),
        "--output-path",
        str(output_dir),
    )
    matches = sorted(output_dir.rglob("taxonomy.tsv"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"{artifact.name}: expected one exported taxonomy.tsv"
        )
    return matches[0]


def import_artifact(
    runner: CommandRunner,
    environment: str,
    input_path: Path,
    output_path: Path,
    semantic_type: str,
    label: str,
    input_format: str | None = None,
) -> None:
    command = [
        "tools",
        "import",
        "--type",
        semantic_type,
        "--input-path",
        str(input_path),
        "--output-path",
        str(output_path),
        "--validate-level",
        "max",
    ]
    if input_format is not None:
        command.extend(["--input-format", input_format])
    runner.qiime(label, environment, *command)


def classifier_fit_arguments(
    reference_reads: Path,
    reference_taxonomy: Path,
    output_classifier: Path,
    cache_dir: Path,
) -> list[str]:
    return [
        "feature-classifier",
        "fit-classifier-naive-bayes",
        "--i-reference-reads",
        str(reference_reads),
        "--i-reference-taxonomy",
        str(reference_taxonomy),
        "--p-classify--alpha",
        "0.001",
        "--p-classify--chunk-size",
        "20000",
        "--p-no-classify--fit-prior",
        "--p-no-feat-ext--alternate-sign",
        "--p-feat-ext--analyzer",
        "char_wb",
        "--p-feat-ext--n-features",
        "8192",
        "--p-feat-ext--ngram-range",
        "[7, 7]",
        "--p-feat-ext--norm",
        "l2",
        "--o-classifier",
        str(output_classifier),
        "--use-cache",
        str(cache_dir),
    ]


def classify_arguments(
    reads: Path,
    classifier: Path,
    confidence: float,
    output_taxonomy: Path,
    cache_dir: Path,
) -> list[str]:
    return [
        "feature-classifier",
        "classify-sklearn",
        "--i-reads",
        str(reads),
        "--i-classifier",
        str(classifier),
        "--p-reads-per-batch",
        str(READS_PER_BATCH),
        "--p-n-jobs",
        str(CLASSIFIER_JOBS),
        "--p-confidence",
        f"{confidence:.2f}",
        "--p-read-orientation",
        "same",
        "--o-classification",
        str(output_taxonomy),
        "--use-cache",
        str(cache_dir),
    ]


def holdout_metrics(
    expected: dict[str, str],
    observed: dict[str, dict[str, object]],
    confidence_key: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rank, _prefix in RANKS:
        expected_labels = {
            feature_id: ranks[rank]
            for feature_id, taxon in expected.items()
            if rank in (ranks := rank_map(taxon))
        }
        correct = 0
        incorrect = 0
        abstained = 0
        for feature_id, expected_label in expected_labels.items():
            predicted = observed[feature_id]["ranks"]
            predicted_label = (
                predicted.get(rank, "") if isinstance(predicted, dict) else ""
            )
            if not predicted_label:
                abstained += 1
            elif predicted_label == expected_label:
                correct += 1
            else:
                incorrect += 1
        predicted_count = correct + incorrect
        expected_count = len(expected_labels)
        precision = correct / predicted_count if predicted_count else 0.0
        recall = correct / expected_count if expected_count else 0.0
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )
        rows.append(
            {
                "confidence_key": confidence_key,
                "confidence": float(confidence_key[1:]) / 100,
                "rank": rank,
                "expected": expected_count,
                "correct": correct,
                "incorrect": incorrect,
                "abstained": abstained,
                "precision": round(precision, 6),
                "recall": round(recall, 6),
                "f1": round(f1, 6),
            }
        )
    return rows


def confidence_summary(
    observed: dict[str, dict[str, object]],
    confidence_key: str,
) -> dict[str, object]:
    depths = [
        sum(1 for rank, _prefix in RANKS if rank in record["ranks"])
        for record in observed.values()
    ]
    return {
        "confidence_key": confidence_key,
        "confidence": float(confidence_key[1:]) / 100,
        "records": len(observed),
        "assigned": sum(
            1 for record in observed.values() if "domain" in record["ranks"]
        ),
        "genus_labeled": sum(
            1 for record in observed.values() if "genus" in record["ranks"]
        ),
        "median_rank_depth": round(statistics.median(depths), 3),
        "mean_rank_depth": round(statistics.fmean(depths), 3),
    }


def query_rank_coverage(
    observed_by_confidence: dict[str, dict[str, dict[str, object]]],
    abundance: dict[str, dict[str, int]],
) -> list[dict[str, object]]:
    total_reads = sum(row["frequency"] for row in abundance.values())
    rows: list[dict[str, object]] = []
    for confidence_key, observed in observed_by_confidence.items():
        for rank, _prefix in RANKS:
            ids = {
                feature_id
                for feature_id, record in observed.items()
                if rank in record["ranks"]
            }
            reads = sum(abundance[feature_id]["frequency"] for feature_id in ids)
            rows.append(
                {
                    "confidence_key": confidence_key,
                    "confidence": float(confidence_key[1:]) / 100,
                    "rank": rank,
                    "asv_count": len(ids),
                    "asv_fraction": round(len(ids) / len(abundance), 6),
                    "read_count": reads,
                    "read_fraction": round(reads / total_reads, 6),
                }
            )
    return rows


def method_concordance(
    naive_bayes: dict[str, dict[str, object]],
    sequence_search: dict[str, dict[str, object]],
    abundance: dict[str, dict[str, int]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summary: list[dict[str, object]] = []
    crosswalk: list[dict[str, object]] = []
    total_reads = sum(row["frequency"] for row in abundance.values())

    for feature_id in sorted(abundance):
        nb_record = naive_bayes[feature_id]
        search_record = sequence_search[feature_id]
        row: dict[str, object] = {
            "feature_id": feature_id,
            "frequency": abundance[feature_id]["frequency"],
            "samples_observed": abundance[feature_id]["samples_observed"],
            "naive_bayes_taxon": nb_record["taxon"],
            "naive_bayes_confidence": nb_record["confidence"],
            "sequence_search_taxon": search_record["taxon"],
            "sequence_search_confidence": search_record["confidence"],
        }
        for rank, _prefix in RANKS:
            nb_label = nb_record["ranks"].get(rank, "")
            search_label = search_record["ranks"].get(rank, "")
            row[f"{rank}_naive_bayes"] = nb_label
            row[f"{rank}_sequence_search"] = search_label
            row[f"{rank}_exact_match"] = bool(
                nb_label and search_label and nb_label == search_label
            )
        crosswalk.append(row)

    for rank, _prefix in RANKS:
        nb_ids = {
            feature_id
            for feature_id, record in naive_bayes.items()
            if rank in record["ranks"]
        }
        search_ids = {
            feature_id
            for feature_id, record in sequence_search.items()
            if rank in record["ranks"]
        }
        comparable = nb_ids & search_ids
        agreeing = {
            feature_id
            for feature_id in comparable
            if naive_bayes[feature_id]["ranks"][rank]
            == sequence_search[feature_id]["ranks"][rank]
        }
        comparable_reads = sum(
            abundance[feature_id]["frequency"] for feature_id in comparable
        )
        agreeing_reads = sum(
            abundance[feature_id]["frequency"] for feature_id in agreeing
        )
        summary.append(
            {
                "rank": rank,
                "naive_bayes_asvs": len(nb_ids),
                "naive_bayes_reads": sum(
                    abundance[feature_id]["frequency"] for feature_id in nb_ids
                ),
                "sequence_search_asvs": len(search_ids),
                "sequence_search_reads": sum(
                    abundance[feature_id]["frequency"]
                    for feature_id in search_ids
                ),
                "comparable_asvs": len(comparable),
                "comparable_reads": comparable_reads,
                "exact_match_asvs": len(agreeing),
                "exact_match_reads": agreeing_reads,
                "exact_match_asv_fraction": round(
                    len(agreeing) / len(comparable), 6
                )
                if comparable
                else 0.0,
                "exact_match_read_fraction": round(
                    agreeing_reads / comparable_reads, 6
                )
                if comparable_reads
                else 0.0,
                "naive_bayes_read_coverage": round(
                    sum(
                        abundance[feature_id]["frequency"]
                        for feature_id in nb_ids
                    )
                    / total_reads,
                    6,
                ),
                "sequence_search_read_coverage": round(
                    sum(
                        abundance[feature_id]["frequency"]
                        for feature_id in search_ids
                    )
                    / total_reads,
                    6,
                ),
            }
        )
    return summary, crosswalk


def configure_plotting() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "figure.dpi": 150,
            "savefig.dpi": 350,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: object, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(
        base.with_suffix(".png"),
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def draw_build_contract(figure_dir: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    configure_plotting()
    fig, ax = plt.subplots(figsize=(11.4, 5.6))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 7)
    ax.axis("off")

    nodes = [
        (0.4, 4.9, 2.3, 1.2, "SILVA 138.2\nfull SSU reference", "#E5E7EB"),
        (3.1, 4.9, 2.3, 1.2, "515F/806R extraction\nexact V4 LCA", "#BFDBFE"),
        (5.8, 4.9, 2.3, 1.2, "Known-lineage split\n90% train · 10% holdout", "#D1FAE5"),
        (8.5, 4.9, 2.3, 1.2, "Fit + threshold audit\n0.50 · 0.70 · 0.90", "#FEF3C7"),
        (10.2, 2.4, 2.3, 1.2, "Refit all 262,752\nproduction classifier", "#FDE68A"),
        (6.8, 2.4, 2.3, 1.2, "366 real V4 ASVs\nconfidence sensitivity", "#DDD6FE"),
        (3.4, 2.4, 2.3, 1.2, "Sequence-search\nmethod concordance", "#FBCFE8"),
        (0.4, 2.4, 2.3, 1.2, "Archive classifier\nversions + provenance", "#C7D2FE"),
    ]
    for x, y, width, height, text, color in nodes:
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.04,rounding_size=0.08",
            linewidth=1.2,
            edgecolor="#374151",
            facecolor=color,
        )
        ax.add_patch(patch)
        ax.text(
            x + width / 2,
            y + height / 2,
            text,
            ha="center",
            va="center",
            fontsize=9.2,
            weight="bold" if "production" in text else "normal",
        )

    arrow = dict(arrowstyle="-|>", lw=1.5, color="#4B5563")
    for x1, x2 in ((2.7, 3.1), (5.4, 5.8), (8.1, 8.5)):
        ax.annotate("", xy=(x2, 5.5), xytext=(x1, 5.5), arrowprops=arrow)
    ax.annotate("", xy=(11.35, 3.6), xytext=(9.65, 4.9), arrowprops=arrow)
    ax.annotate("", xy=(9.1, 3.0), xytext=(10.2, 3.0), arrowprops=arrow)
    ax.annotate("", xy=(5.7, 3.0), xytext=(6.8, 3.0), arrowprops=arrow)
    ax.annotate("", xy=(2.7, 3.0), xytext=(3.4, 3.0), arrowprops=arrow)

    ax.text(
        6.5,
        6.65,
        "A region-specific classifier is a versioned model, not just a database file",
        ha="center",
        va="center",
        fontsize=15,
        weight="bold",
    )
    ax.text(
        6.5,
        0.85,
        "Internal holdout = interpolation evidence; real-ASV concordance = method sensitivity",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#374151",
    )
    save_figure(fig, figure_dir / "15-classifier-build-contract")
    plt.close(fig)


def draw_holdout_performance(
    rows: Sequence[dict[str, object]],
    figure_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    configure_plotting()
    ranks = [rank for rank, _prefix in RANKS]
    metrics = (("precision", "Precision"), ("recall", "Recall"), ("f1", "F1"))
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.1), sharey=True)
    for axis, (metric, title) in zip(axes, metrics):
        for confidence_key in ("c050", "c070", "c090"):
            subset = {
                str(row["rank"]): float(row[metric])
                for row in rows
                if row["confidence_key"] == confidence_key
            }
            axis.plot(
                ranks,
                [subset[rank] for rank in ranks],
                marker="o",
                linewidth=2,
                color=COLORS[confidence_key],
                label=confidence_key.replace("c0", "0."),
            )
        axis.set_title(title, weight="bold")
        axis.set_xlabel("Taxonomic rank")
        axis.set_ylim(0.70, 1.01)
        axis.grid(axis="y", alpha=0.25)
        axis.tick_params(axis="x", rotation=35)
    axes[0].set_ylabel("Holdout performance")
    axes[-1].legend(title="Confidence", frameon=False, loc="lower left")
    fig.suptitle(
        "Confidence trades recall for precision",
        fontsize=14,
        weight="bold",
        y=1.04,
    )
    fig.text(
        0.5,
        -0.03,
        "Same-release interpolation; not external mock-community accuracy",
        ha="center",
        color="#4B5563",
    )
    fig.tight_layout()
    save_figure(fig, figure_dir / "15-holdout-rank-performance")
    plt.close(fig)


def draw_query_sensitivity(
    rows: Sequence[dict[str, object]],
    figure_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    configure_plotting()
    ranks = [rank for rank, _prefix in RANKS]
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    for confidence_key in ("c050", "c070", "c090"):
        subset = {
            str(row["rank"]): 100 * float(row["read_fraction"])
            for row in rows
            if row["confidence_key"] == confidence_key
        }
        values = [subset[rank] for rank in ranks]
        ax.plot(
            ranks,
            values,
            marker="o",
            linewidth=2.2,
            color=COLORS[confidence_key],
            label=confidence_key.replace("c0", "0."),
        )
        for index, value in enumerate(values):
            if index in {4, 5}:
                offset = {
                    "c050": 1.0,
                    "c070": 0.8,
                    "c090": -2.5,
                }[confidence_key]
                ax.text(
                    index,
                    value + offset,
                    f"{value:.1f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=COLORS[confidence_key],
                )
    ax.set_ylim(65, 102)
    ax.set_ylabel("Read-weighted ASVs with a label (%)")
    ax.set_xlabel("Taxonomic rank")
    ax.set_title(
        "Confidence changes usable depth on 5,213 real V4 reads",
        weight="bold",
    )
    ax.grid(axis="y", alpha=0.25)
    ax.legend(title="Confidence", frameon=False, loc="lower left")
    fig.tight_layout()
    save_figure(fig, figure_dir / "15-query-confidence-sensitivity")
    plt.close(fig)


def draw_method_concordance(
    rows: Sequence[dict[str, object]],
    figure_dir: Path,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_plotting()
    ranks = [rank for rank, _prefix in RANKS]
    row_map = {str(row["rank"]): row for row in rows}
    x = np.arange(len(ranks))
    width = 0.24
    nb = [
        100 * float(row_map[rank]["naive_bayes_read_coverage"]) for rank in ranks
    ]
    search = [
        100 * float(row_map[rank]["sequence_search_read_coverage"])
        for rank in ranks
    ]
    agreement = [
        100 * float(row_map[rank]["exact_match_read_fraction"]) for rank in ranks
    ]

    fig, ax = plt.subplots(figsize=(9.2, 5.1))
    ax.bar(
        x - width,
        nb,
        width,
        label="Naive Bayes coverage",
        color="#0072B2",
    )
    ax.bar(
        x,
        search,
        width,
        label="Sequence-search coverage",
        color="#009E73",
    )
    ax.bar(
        x + width,
        agreement,
        width,
        label="Exact agreement among comparable reads",
        color="#CC79A7",
    )
    ax.set_xticks(x, [rank.title() for rank in ranks], rotation=30, ha="right")
    ax.set_ylim(0, 105)
    ax.set_ylabel("Read-weighted percentage (%)")
    ax.set_title(
        "Classifier concordance is a sensitivity audit, not a truth benchmark",
        weight="bold",
    )
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    save_figure(fig, figure_dir / "15-method-concordance")
    plt.close(fig)


def audit_row(
    check_id: str,
    observed: object,
    expected: object,
    passed: bool,
) -> dict[str, object]:
    return {
        "check_id": check_id,
        "observed": observed,
        "expected": expected,
        "status": "PASS" if passed else "FAIL",
    }


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    sequence_search_supplied = args.sequence_search_taxonomy is not None
    sequence_search_artifact = (
        args.sequence_search_taxonomy.resolve()
        if args.sequence_search_taxonomy is not None
        else None
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    reference_dir = project_root / "data/small/taxonomy-databases"
    reference_fasta = reference_dir / "silva-138.2-v4.fasta.gz"
    reference_taxonomy = reference_dir / "silva-138.2-v4-taxonomy.tsv.gz"
    query_fasta = reference_dir / "query-v4-asvs.fasta.gz"
    query_abundance_path = reference_dir / "query-v4-asv-abundance.tsv"
    source_summary_path = reference_dir / "silva-138.2-source-summary.json"

    required_paths = [
        reference_fasta,
        reference_taxonomy,
        query_fasta,
        query_abundance_path,
        source_summary_path,
    ]
    if sequence_search_artifact is not None:
        required_paths.append(sequence_search_artifact)
    for path in required_paths:
        if not path.exists():
            raise FileNotFoundError(path)

    command_log = output_dir / "qiime-classifier.log"
    validation_log = output_dir / "validation.log"

    with tempfile.TemporaryDirectory(
        prefix=".article15-work-",
        dir=str(output_dir.parent),
    ) as temporary:
        work_dir = Path(temporary)
        cache_dir = work_dir / "qiime-cache"
        runtime_environment = isolated_environment(work_dir)
        for variable in (
            "XDG_CACHE_HOME",
            "NUMBA_CACHE_DIR",
            "MPLCONFIGDIR",
            "PYTHONHASHSEED",
            "LC_ALL",
            "LANG",
        ):
            os.environ[variable] = runtime_environment[variable]
        runner = CommandRunner(
            runtime_environment,
            command_log,
            {
                str(work_dir): "<WORK_DIR>",
                str(output_dir): "<OUTPUT_DIR>",
                str(project_root): "<PROJECT_ROOT>",
            },
        )
        runner.qiime(
            "qiime-cache-create",
            args.qiime_env,
            "tools",
            "cache-create",
            "--cache",
            str(cache_dir),
        )

        package_result = runner.run(
            "conda-package-audit",
            ["conda", "list", "-n", args.qiime_env, "--json"],
        )
        package_rows = json.loads(str(package_result["stdout"]))
        packages = {row["name"]: row["version"] for row in package_rows}
        versions = {
            "qiime2": packages.get("qiime2", ""),
            "q2-feature-classifier": packages.get("q2-feature-classifier", ""),
            "rescript": packages.get("rescript", ""),
            "scikit-learn": packages.get("scikit-learn", ""),
            "joblib": packages.get("joblib", ""),
            "numpy": packages.get("numpy", ""),
            "scipy": packages.get("scipy", ""),
        }

        source_summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
        reference_sequences = dict(iter_fasta(reference_fasta))
        reference_taxonomy_map = read_taxonomy(reference_taxonomy)
        query_sequences = dict(iter_fasta(query_fasta))
        query_abundance = read_query_abundance(query_abundance_path)

        if set(reference_sequences) != set(reference_taxonomy_map):
            raise ValueError("reference sequence and taxonomy IDs differ")
        if set(query_sequences) != set(query_abundance):
            raise ValueError("query sequence and abundance IDs differ")

        training_ids, holdout_ids, membership, split_audit = deterministic_split(
            reference_sequences, reference_taxonomy_map
        )
        membership_path = output_dir / "reference-split-membership.tsv.gz"
        write_gzip_tsv(
            membership_path,
            membership,
            (
                "feature_id",
                "split",
                "genus_lineage",
                "lineage_records",
                "split_hash",
            ),
        )

        train_fasta = work_dir / "silva-v4-training.fasta"
        train_taxonomy = work_dir / "silva-v4-training-taxonomy.tsv"
        holdout_fasta = work_dir / "silva-v4-holdout.fasta"
        holdout_taxonomy = work_dir / "silva-v4-holdout-taxonomy.tsv"
        full_fasta = work_dir / "silva-v4-full.fasta"
        full_taxonomy = work_dir / "silva-v4-full-taxonomy.tsv"
        query_fasta_plain = work_dir / "query-v4-asvs.fasta"

        write_fasta(train_fasta, training_ids, reference_sequences)
        write_taxonomy(train_taxonomy, training_ids, reference_taxonomy_map)
        write_fasta(holdout_fasta, holdout_ids, reference_sequences)
        write_taxonomy(holdout_taxonomy, holdout_ids, reference_taxonomy_map)
        unpack_gzip(reference_fasta, full_fasta)
        unpack_gzip(reference_taxonomy, full_taxonomy)
        unpack_gzip(query_fasta, query_fasta_plain)

        train_sequences_qza = work_dir / "training-sequences.qza"
        train_taxonomy_qza = work_dir / "training-taxonomy.qza"
        holdout_sequences_qza = work_dir / "holdout-sequences.qza"
        holdout_taxonomy_qza = work_dir / "holdout-expected-taxonomy.qza"
        full_sequences_qza = work_dir / "full-sequences.qza"
        full_taxonomy_qza = work_dir / "full-taxonomy.qza"
        query_sequences_qza = work_dir / "query-sequences.qza"

        for input_path, output_path, semantic_type, label, input_format in (
            (
                train_fasta,
                train_sequences_qza,
                "FeatureData[Sequence]",
                "import-training-sequences",
                None,
            ),
            (
                train_taxonomy,
                train_taxonomy_qza,
                "FeatureData[Taxonomy]",
                "import-training-taxonomy",
                "TSVTaxonomyFormat",
            ),
            (
                holdout_fasta,
                holdout_sequences_qza,
                "FeatureData[Sequence]",
                "import-holdout-sequences",
                None,
            ),
            (
                holdout_taxonomy,
                holdout_taxonomy_qza,
                "FeatureData[Taxonomy]",
                "import-holdout-taxonomy",
                "TSVTaxonomyFormat",
            ),
            (
                full_fasta,
                full_sequences_qza,
                "FeatureData[Sequence]",
                "import-full-sequences",
                None,
            ),
            (
                full_taxonomy,
                full_taxonomy_qza,
                "FeatureData[Taxonomy]",
                "import-full-taxonomy",
                "TSVTaxonomyFormat",
            ),
            (
                query_fasta_plain,
                query_sequences_qza,
                "FeatureData[Sequence]",
                "import-query-sequences",
                None,
            ),
        ):
            import_artifact(
                runner,
                args.qiime_env,
                input_path,
                output_path,
                semantic_type,
                label,
                input_format,
            )

        holdout_classifier_work = work_dir / "silva-v4-holdout-classifier.qza"
        production_classifier_work = work_dir / "silva-v4-classifier.qza"
        fit_holdout = runner.qiime(
            "fit-holdout-classifier",
            args.qiime_env,
            *classifier_fit_arguments(
                train_sequences_qza,
                train_taxonomy_qza,
                holdout_classifier_work,
                cache_dir,
            ),
            timeout=28_800,
        )
        fit_production = runner.qiime(
            "fit-production-classifier",
            args.qiime_env,
            *classifier_fit_arguments(
                full_sequences_qza,
                full_taxonomy_qza,
                production_classifier_work,
                cache_dir,
            ),
            timeout=28_800,
        )

        holdout_taxonomy_work: dict[str, Path] = {}
        query_taxonomy_work: dict[str, Path] = {}
        primary_command_labels = [
            "fit-holdout-classifier",
            "fit-production-classifier",
        ]
        for confidence in CONFIDENCES:
            confidence_key = CONFIDENCE_KEYS[confidence]
            holdout_output = work_dir / f"holdout-taxonomy-{confidence_key}.qza"
            query_output = work_dir / f"query-taxonomy-{confidence_key}.qza"
            runner.qiime(
                f"classify-holdout-{confidence_key}",
                args.qiime_env,
                *classify_arguments(
                    holdout_sequences_qza,
                    holdout_classifier_work,
                    confidence,
                    holdout_output,
                    cache_dir,
                ),
                timeout=14_400,
            )
            runner.qiime(
                f"classify-query-{confidence_key}",
                args.qiime_env,
                *classify_arguments(
                    query_sequences_qza,
                    production_classifier_work,
                    confidence,
                    query_output,
                    cache_dir,
                ),
                timeout=14_400,
            )
            primary_command_labels.extend(
                [
                    f"classify-holdout-{confidence_key}",
                    f"classify-query-{confidence_key}",
                ]
            )
            holdout_taxonomy_work[confidence_key] = holdout_output
            query_taxonomy_work[confidence_key] = query_output

        if sequence_search_artifact is None:
            sequence_search_artifact = (
                work_dir / "silva-sequence-search-taxonomy.qza"
            )
            sequence_search_results = (
                work_dir / "silva-sequence-search-results.qza"
            )
            runner.qiime(
                "classify-sequence-search-comparator",
                args.qiime_env,
                "feature-classifier",
                "classify-consensus-vsearch",
                "--i-query",
                str(query_sequences_qza),
                "--i-reference-reads",
                str(full_sequences_qza),
                "--i-reference-taxonomy",
                str(full_taxonomy_qza),
                "--p-maxaccepts",
                "50",
                "--p-maxrejects",
                "all",
                "--p-perc-identity",
                "0.8",
                "--p-query-cov",
                "0.8",
                "--p-strand",
                "both",
                "--p-top-hits-only",
                "--p-maxhits",
                "all",
                "--p-output-no-hits",
                "--p-threads",
                str(CLASSIFIER_JOBS),
                "--p-min-consensus",
                "0.51",
                "--p-unassignable-label",
                "Unassigned",
                "--o-classification",
                str(sequence_search_artifact),
                "--o-search-results",
                str(sequence_search_results),
                "--use-cache",
                str(cache_dir),
                "--no-recycle",
                timeout=14_400,
            )

        evaluation_work = work_dir / "holdout-evaluation-c070.qzv"
        runner.qiime(
            "evaluate-holdout-c070",
            args.qiime_env,
            "quality-control",
            "evaluate-taxonomy",
            "--i-expected-taxa",
            str(holdout_taxonomy_qza),
            "--i-observed-taxa",
            str(holdout_taxonomy_work["c070"]),
            "--p-depth",
            "6",
            "--p-palette",
            "tab10",
            "--o-visualization",
            str(evaluation_work),
            "--use-cache",
            str(cache_dir),
        )
        primary_command_labels.append("evaluate-holdout-c070")

        output_artifact_sources = {
            "silva-v4-holdout-classifier.qza": holdout_classifier_work,
            "silva-v4-classifier.qza": production_classifier_work,
            "silva-v4-holdout-sequences.qza": holdout_sequences_qza,
            "silva-v4-holdout-expected-taxonomy.qza": holdout_taxonomy_qza,
            "holdout-evaluation-c070.qzv": evaluation_work,
        }
        for confidence_key, source in holdout_taxonomy_work.items():
            output_artifact_sources[
                f"holdout-taxonomy-{confidence_key}.qza"
            ] = source
        for confidence_key, source in query_taxonomy_work.items():
            output_artifact_sources[
                f"query-taxonomy-{confidence_key}.qza"
            ] = source

        output_artifacts: dict[str, Path] = {}
        for filename, source in output_artifact_sources.items():
            destination = output_dir / filename
            shutil.copy2(source, destination)
            output_artifacts[filename] = destination

        holdout_observed: dict[str, dict[str, dict[str, object]]] = {}
        query_observed: dict[str, dict[str, dict[str, object]]] = {}
        holdout_metric_rows: list[dict[str, object]] = []
        holdout_summary_rows: list[dict[str, object]] = []

        holdout_expected = {
            feature_id: reference_taxonomy_map[feature_id]
            for feature_id in holdout_ids
        }
        for confidence_key, artifact in holdout_taxonomy_work.items():
            exported = export_artifact(
                runner,
                args.qiime_env,
                artifact,
                work_dir / f"export-holdout-{confidence_key}",
                f"export-holdout-{confidence_key}",
            )
            observed = parse_classification(exported)
            if set(observed) != set(holdout_expected):
                raise ValueError(
                    f"{confidence_key}: holdout classification IDs differ"
                )
            holdout_observed[confidence_key] = observed
            holdout_metric_rows.extend(
                holdout_metrics(holdout_expected, observed, confidence_key)
            )
            holdout_summary_rows.append(
                confidence_summary(observed, confidence_key)
            )

        for confidence_key, artifact in query_taxonomy_work.items():
            exported = export_artifact(
                runner,
                args.qiime_env,
                artifact,
                work_dir / f"export-query-{confidence_key}",
                f"export-query-{confidence_key}",
            )
            observed = parse_classification(exported)
            if set(observed) != set(query_sequences):
                raise ValueError(
                    f"{confidence_key}: query classification IDs differ"
                )
            query_observed[confidence_key] = observed

        sequence_search_export = export_artifact(
            runner,
            args.qiime_env,
            sequence_search_artifact,
            work_dir / "export-sequence-search",
            "export-sequence-search-taxonomy",
        )
        sequence_search = parse_classification(sequence_search_export)
        if set(sequence_search) != set(query_sequences):
            raise ValueError("sequence-search taxonomy IDs differ from query")

        query_coverage_rows = query_rank_coverage(
            query_observed, query_abundance
        )
        method_rows, crosswalk_rows = method_concordance(
            query_observed["c070"],
            sequence_search,
            query_abundance,
        )

        write_tsv(
            output_dir / "holdout-rank-metrics.tsv",
            holdout_metric_rows,
            (
                "confidence_key",
                "confidence",
                "rank",
                "expected",
                "correct",
                "incorrect",
                "abstained",
                "precision",
                "recall",
                "f1",
            ),
        )
        write_tsv(
            output_dir / "holdout-confidence-summary.tsv",
            holdout_summary_rows,
            (
                "confidence_key",
                "confidence",
                "records",
                "assigned",
                "genus_labeled",
                "median_rank_depth",
                "mean_rank_depth",
            ),
        )
        write_tsv(
            output_dir / "query-rank-coverage.tsv",
            query_coverage_rows,
            (
                "confidence_key",
                "confidence",
                "rank",
                "asv_count",
                "asv_fraction",
                "read_count",
                "read_fraction",
            ),
        )
        write_tsv(
            output_dir / "query-method-concordance.tsv",
            method_rows,
            (
                "rank",
                "naive_bayes_asvs",
                "naive_bayes_reads",
                "sequence_search_asvs",
                "sequence_search_reads",
                "comparable_asvs",
                "comparable_reads",
                "exact_match_asvs",
                "exact_match_reads",
                "exact_match_asv_fraction",
                "exact_match_read_fraction",
                "naive_bayes_read_coverage",
                "sequence_search_read_coverage",
            ),
        )
        crosswalk_fields = [
            "feature_id",
            "frequency",
            "samples_observed",
            "naive_bayes_taxon",
            "naive_bayes_confidence",
            "sequence_search_taxon",
            "sequence_search_confidence",
        ]
        for rank, _prefix in RANKS:
            crosswalk_fields.extend(
                [
                    f"{rank}_naive_bayes",
                    f"{rank}_sequence_search",
                    f"{rank}_exact_match",
                ]
            )
        write_tsv(
            output_dir / "query-taxonomy-crosswalk.tsv",
            crosswalk_rows,
            crosswalk_fields,
        )

        split_rows = [
            {
                "split": "all",
                "records": len(reference_sequences),
                "genus_lineages": split_audit["lineages"],
                "singleton_lineages": split_audit["singleton_lineages"],
                "record_fraction": 1.0,
                "identifier_sha256": sha256_text(sorted(reference_sequences)),
            },
            {
                "split": "training",
                "records": len(training_ids),
                "genus_lineages": split_audit["training_lineages"],
                "singleton_lineages": split_audit["singleton_lineages"],
                "record_fraction": round(
                    len(training_ids) / len(reference_sequences), 6
                ),
                "identifier_sha256": sha256_text(training_ids),
            },
            {
                "split": "holdout",
                "records": len(holdout_ids),
                "genus_lineages": split_audit["holdout_lineages"],
                "singleton_lineages": 0,
                "record_fraction": round(
                    len(holdout_ids) / len(reference_sequences), 6
                ),
                "identifier_sha256": sha256_text(holdout_ids),
            },
        ]
        write_tsv(
            output_dir / "reference-split-audit.tsv",
            split_rows,
            (
                "split",
                "records",
                "genus_lineages",
                "singleton_lineages",
                "record_fraction",
                "identifier_sha256",
            ),
        )

        provenance_rows: list[dict[str, object]] = []
        expected_types = {
            "silva-v4-holdout-classifier.qza": "TaxonomicClassifier",
            "silva-v4-classifier.qza": "TaxonomicClassifier",
            "silva-v4-holdout-sequences.qza": "FeatureData[Sequence]",
            "silva-v4-holdout-expected-taxonomy.qza": "FeatureData[Taxonomy]",
            "holdout-taxonomy-c050.qza": "FeatureData[Taxonomy]",
            "holdout-taxonomy-c070.qza": "FeatureData[Taxonomy]",
            "holdout-taxonomy-c090.qza": "FeatureData[Taxonomy]",
            "query-taxonomy-c050.qza": "FeatureData[Taxonomy]",
            "query-taxonomy-c070.qza": "FeatureData[Taxonomy]",
            "query-taxonomy-c090.qza": "FeatureData[Taxonomy]",
            "holdout-evaluation-c070.qzv": "Visualization",
        }
        artifact_validation: dict[str, bool] = {}
        artifact_actions: dict[str, list[tuple[str, str]]] = {}
        for filename in expected_types:
            artifact = output_artifacts[filename]
            runner.qiime(
                f"validate-{filename}",
                args.qiime_env,
                "tools",
                "validate",
                str(artifact),
                "--level",
                "max",
            )
            metadata = artifact_metadata(artifact)
            actions = provenance_actions(artifact)
            artifact_validation[filename] = (
                metadata["type"] == expected_types[filename]
            )
            artifact_actions[filename] = actions
            provenance_rows.append(
                {
                    "artifact": filename,
                    "uuid": metadata["uuid"],
                    "semantic_type": metadata["type"],
                    "data_format": metadata["format"],
                    "maximum_validation": "PASS",
                    "sha256": sha256_file(artifact),
                    "bytes": artifact.stat().st_size,
                    "provenance_actions": "; ".join(
                        f"{plugin}:{action}" for plugin, action in actions
                    ),
                }
            )
        write_tsv(
            output_dir / "provenance-audit.tsv",
            provenance_rows,
            (
                "artifact",
                "uuid",
                "semantic_type",
                "data_format",
                "maximum_validation",
                "sha256",
                "bytes",
                "provenance_actions",
            ),
        )

        environment_rows = [
            {
                "component": component,
                "observed_version": observed,
                "expected_version": expected,
                "status": "PASS" if observed == expected else "FAIL",
            }
            for component, observed, expected in (
                ("QIIME 2", versions["qiime2"], QIIME_VERSION),
                (
                    "q2-feature-classifier",
                    versions["q2-feature-classifier"],
                    Q2_FEATURE_CLASSIFIER_VERSION,
                ),
                ("RESCRIPt", versions["rescript"], RESCRIPT_VERSION),
                ("scikit-learn", versions["scikit-learn"], SKLEARN_VERSION),
                ("joblib", versions["joblib"], JOBLIB_VERSION),
                ("NumPy", versions["numpy"], NUMPY_VERSION),
                ("SciPy", versions["scipy"], SCIPY_VERSION),
            )
        ]
        write_tsv(
            output_dir / "environment-audit.tsv",
            environment_rows,
            (
                "component",
                "observed_version",
                "expected_version",
                "status",
            ),
        )

        draw_build_contract(figure_dir)
        draw_holdout_performance(holdout_metric_rows, figure_dir)
        draw_query_sensitivity(query_coverage_rows, figure_dir)
        draw_method_concordance(method_rows, figure_dir)

        metric_lookup = {
            (str(row["confidence_key"]), str(row["rank"])): row
            for row in holdout_metric_rows
        }
        query_lookup = {
            (str(row["confidence_key"]), str(row["rank"])): row
            for row in query_coverage_rows
        }
        method_lookup = {str(row["rank"]): row for row in method_rows}
        training_lineages = {
            genus_lineage(reference_taxonomy_map[feature_id])
            for feature_id in training_ids
        }
        holdout_lineages = {
            genus_lineage(reference_taxonomy_map[feature_id])
            for feature_id in holdout_ids
        }
        sequence_overlap = len(
            {reference_sequences[feature_id] for feature_id in training_ids}
            & {reference_sequences[feature_id] for feature_id in holdout_ids}
        )

        checks: list[dict[str, object]] = []

        def add(
            check_id: str,
            observed: object,
            expected: object,
            passed: bool,
        ) -> None:
            checks.append(audit_row(check_id, observed, expected, passed))

        # Environment and immutable inputs: 20 checks.
        for check_id, observed, expected in (
            ("environment-qiime2", versions["qiime2"], QIIME_VERSION),
            (
                "environment-q2-feature-classifier",
                versions["q2-feature-classifier"],
                Q2_FEATURE_CLASSIFIER_VERSION,
            ),
            ("environment-rescript", versions["rescript"], RESCRIPT_VERSION),
            (
                "environment-scikit-learn",
                versions["scikit-learn"],
                SKLEARN_VERSION,
            ),
            ("environment-joblib", versions["joblib"], JOBLIB_VERSION),
            ("environment-numpy", versions["numpy"], NUMPY_VERSION),
            ("environment-scipy", versions["scipy"], SCIPY_VERSION),
            (
                "input-reference-fasta-sha256",
                sha256_file(reference_fasta),
                REFERENCE_FASTA_SHA256,
            ),
            (
                "input-reference-taxonomy-sha256",
                sha256_file(reference_taxonomy),
                REFERENCE_TAXONOMY_SHA256,
            ),
            (
                "input-query-fasta-sha256",
                sha256_file(query_fasta),
                QUERY_FASTA_SHA256,
            ),
            (
                "input-query-abundance-sha256",
                sha256_file(query_abundance_path),
                QUERY_ABUNDANCE_SHA256,
            ),
            (
                "input-reference-records",
                len(reference_sequences),
                REFERENCE_RECORDS,
            ),
            (
                "input-reference-unique-ids",
                len(set(reference_sequences)),
                REFERENCE_RECORDS,
            ),
            (
                "input-reference-taxonomy-records",
                len(reference_taxonomy_map),
                REFERENCE_RECORDS,
            ),
            (
                "input-reference-id-match",
                set(reference_sequences) == set(reference_taxonomy_map),
                True,
            ),
            ("input-query-asvs", len(query_sequences), QUERY_ASVS),
            (
                "input-query-unique-ids",
                len(set(query_sequences)),
                QUERY_ASVS,
            ),
            (
                "input-query-abundance-records",
                len(query_abundance),
                QUERY_ASVS,
            ),
            (
                "input-query-id-match",
                set(query_sequences) == set(query_abundance),
                True,
            ),
            (
                "input-query-reads",
                sum(row["frequency"] for row in query_abundance.values()),
                QUERY_READS,
            ),
        ):
            add(check_id, observed, expected, observed == expected)

        # Deterministic split: 14 checks.
        split_checks = (
            ("split-training-records", len(training_ids), TRAINING_RECORDS),
            ("split-holdout-records", len(holdout_ids), HOLDOUT_RECORDS),
            (
                "split-record-conservation",
                len(training_ids) + len(holdout_ids),
                REFERENCE_RECORDS,
            ),
            (
                "split-id-overlap",
                len(set(training_ids) & set(holdout_ids)),
                0,
            ),
            ("split-sequence-overlap", sequence_overlap, 0),
            (
                "split-training-subset",
                set(training_ids) <= set(reference_sequences),
                True,
            ),
            (
                "split-holdout-subset",
                set(holdout_ids) <= set(reference_sequences),
                True,
            ),
            ("split-membership-records", len(membership), REFERENCE_RECORDS),
            (
                "split-holdout-lineages-in-training",
                len(holdout_lineages - training_lineages),
                0,
            ),
            (
                "split-singletons-in-holdout",
                sum(
                    1
                    for row in membership
                    if row["split"] == "holdout"
                    and int(row["lineage_records"]) == 1
                ),
                0,
            ),
            (
                "split-holdout-fraction",
                round(len(holdout_ids) / len(reference_sequences), 6),
                0.099851,
            ),
            (
                "split-training-taxonomy-id-match",
                set(training_ids)
                == {feature_id for feature_id in training_ids if feature_id in reference_taxonomy_map},
                True,
            ),
            (
                "split-holdout-taxonomy-id-match",
                set(holdout_ids)
                == {feature_id for feature_id in holdout_ids if feature_id in reference_taxonomy_map},
                True,
            ),
            (
                "split-membership-sha256",
                bool(re.fullmatch(r"[0-9a-f]{64}", sha256_file(membership_path))),
                True,
            ),
        )
        for check_id, observed, expected in split_checks:
            add(check_id, observed, expected, observed == expected)

        # Nine primary fit/classify/evaluate commands.
        primary_records = {
            str(record["label"]): record for record in runner.records
        }
        for label in primary_command_labels:
            observed = int(primary_records[label]["returncode"])
            add(f"command-{label}", observed, 0, observed == 0)

        # Eleven retained Artifacts/Visualization objects.
        for filename, expected_type in expected_types.items():
            observed = artifact_validation[filename]
            add(
                f"artifact-{filename}-maximum-validation",
                f"{expected_type}; {observed}",
                f"{expected_type}; True",
                observed,
            )

        # Four provenance lineage checks.
        holdout_fit_action = (
            "feature-classifier",
            "fit_classifier_naive_bayes",
        ) in artifact_actions["silva-v4-holdout-classifier.qza"]
        production_fit_action = (
            "feature-classifier",
            "fit_classifier_naive_bayes",
        ) in artifact_actions["silva-v4-classifier.qza"]
        classification_action = all(
            ("feature-classifier", "classify_sklearn")
            in artifact_actions[filename]
            for filename in (
                "holdout-taxonomy-c050.qza",
                "holdout-taxonomy-c070.qza",
                "holdout-taxonomy-c090.qza",
                "query-taxonomy-c050.qza",
                "query-taxonomy-c070.qza",
                "query-taxonomy-c090.qza",
            )
        )
        evaluation_action = (
            "quality-control",
            "evaluate_taxonomy",
        ) in artifact_actions["holdout-evaluation-c070.qzv"]
        for check_id, observed in (
            ("provenance-holdout-fit-action", holdout_fit_action),
            ("provenance-production-fit-action", production_fit_action),
            ("provenance-all-classification-actions", classification_action),
            ("provenance-evaluation-action", evaluation_action),
        ):
            add(check_id, observed, True, observed)

        # Ten performance and sensitivity checks.
        add(
            "holdout-c050-genus-f1",
            metric_lookup[("c050", "genus")]["f1"],
            ">=0.65",
            float(metric_lookup[("c050", "genus")]["f1"]) >= 0.65,
        )
        add(
            "holdout-c070-genus-f1",
            metric_lookup[("c070", "genus")]["f1"],
            ">=0.65",
            float(metric_lookup[("c070", "genus")]["f1"]) >= 0.65,
        )
        add(
            "holdout-c090-genus-f1",
            metric_lookup[("c090", "genus")]["f1"],
            ">=0.50",
            float(metric_lookup[("c090", "genus")]["f1"]) >= 0.50,
        )
        domain_recalls = [
            float(metric_lookup[(key, "domain")]["recall"])
            for key in ("c050", "c070", "c090")
        ]
        genus_recalls = [
            float(metric_lookup[(key, "genus")]["recall"])
            for key in ("c050", "c070", "c090")
        ]
        add(
            "holdout-domain-recall-monotonic",
            domain_recalls,
            "non-increasing",
            domain_recalls[0] >= domain_recalls[1] >= domain_recalls[2],
        )
        add(
            "holdout-genus-recall-monotonic",
            genus_recalls,
            "non-increasing",
            genus_recalls[0] >= genus_recalls[1] >= genus_recalls[2],
        )
        add(
            "holdout-metric-row-count",
            len(holdout_metric_rows),
            18,
            len(holdout_metric_rows) == 18,
        )
        add(
            "query-all-threshold-id-count",
            [len(query_observed[key]) for key in ("c050", "c070", "c090")],
            [QUERY_ASVS] * 3,
            all(
                len(query_observed[key]) == QUERY_ASVS
                for key in ("c050", "c070", "c090")
            ),
        )
        assigned_reads = [
            int(query_lookup[(key, "domain")]["read_count"])
            for key in ("c050", "c070", "c090")
        ]
        add(
            "query-assigned-reads-monotonic",
            assigned_reads,
            "non-increasing",
            assigned_reads[0] >= assigned_reads[1] >= assigned_reads[2],
        )
        add(
            "query-c070-assigned-reads",
            assigned_reads[1],
            ">=3000",
            assigned_reads[1] >= 3000,
        )
        concordance_valid = len(method_rows) == 6 and all(
            0 <= float(row["exact_match_read_fraction"]) <= 1
            for row in method_rows
        )
        add(
            "query-method-concordance-contract",
            f"{len(method_rows)} ranks; fractions bounded",
            "6 ranks; fractions bounded",
            concordance_valid,
        )

        # Twelve publication figure files.
        for stem in (
            "15-classifier-build-contract",
            "15-holdout-rank-performance",
            "15-query-confidence-sensitivity",
            "15-method-concordance",
        ):
            for suffix in (".pdf", ".png", ".tiff"):
                path = figure_dir / f"{stem}{suffix}"
                add(
                    f"figure-{stem}{suffix}",
                    path.stat().st_size if path.exists() else 0,
                    ">0 bytes",
                    path.exists() and path.stat().st_size > 0,
                )

        if len(checks) != 80:
            raise AssertionError(f"expected 80 checks, observed {len(checks)}")
        checks_passed = sum(row["status"] == "PASS" for row in checks)
        write_tsv(
            output_dir / "validation-checks.tsv",
            checks,
            ("check_id", "observed", "expected", "status"),
        )

        nested_holdout: dict[str, dict[str, dict[str, object]]] = {}
        for row in holdout_metric_rows:
            nested_holdout.setdefault(
                str(row["confidence_key"]), {}
            )[str(row["rank"])] = {
                key: row[key]
                for key in (
                    "expected",
                    "correct",
                    "incorrect",
                    "abstained",
                    "precision",
                    "recall",
                    "f1",
                )
            }
        nested_query: dict[str, dict[str, object]] = {}
        for confidence_key in ("c050", "c070", "c090"):
            domain_row = query_lookup[(confidence_key, "domain")]
            genus_row = query_lookup[(confidence_key, "genus")]
            nested_query[confidence_key] = {
                "assigned_asvs": domain_row["asv_count"],
                "assigned_reads": domain_row["read_count"],
                "genus_asvs": genus_row["asv_count"],
                "genus_reads": genus_row["read_count"],
            }

        summary = {
            "status": "passed" if checks_passed == len(checks) else "failed",
            "validation_date": "2026-07-20",
            "reference_database": "SILVA 138.2 SSURef NR99",
            "marker_region": "16S V4 515F/806R",
            "reference_records": len(reference_sequences),
            "query_asvs": len(query_sequences),
            "query_reads": sum(
                row["frequency"] for row in query_abundance.values()
            ),
            "versions": versions,
            "split": {
                "seed": SPLIT_SEED,
                "strategy": (
                    "Known-genus-lineage stratified deterministic SHA-256 "
                    "holdout; singleton lineages remain in training."
                ),
                "training_records": len(training_ids),
                "holdout_records": len(holdout_ids),
                "training_lineages": split_audit["training_lineages"],
                "holdout_lineages": split_audit["holdout_lineages"],
                "singleton_lineages": split_audit["singleton_lineages"],
                "membership_sha256": sha256_file(membership_path),
                "interpretation_boundary": (
                    "Same-release known-lineage interpolation; not an "
                    "independent mock-community or cross-release accuracy estimate."
                ),
            },
            "classifier_contract": {
                "method": "q2-feature-classifier fit-classifier-naive-bayes",
                "feature_analyzer": "char_wb",
                "ngram_range": [7, 7],
                "hashed_features": 8192,
                "alternate_sign": False,
                "norm": "l2",
                "alpha": 0.001,
                "fit_prior": False,
                "class_weight": "not supplied",
                "classification_orientation": "same",
                "reads_per_batch": READS_PER_BATCH,
                "n_jobs": CLASSIFIER_JOBS,
                "confidence_thresholds": list(CONFIDENCES),
            },
            "training": {
                "holdout_classifier_seconds": fit_holdout["elapsed_seconds"],
                "production_classifier_seconds": fit_production[
                    "elapsed_seconds"
                ],
                "holdout_classifier_bytes": output_artifacts[
                    "silva-v4-holdout-classifier.qza"
                ].stat().st_size,
                "production_classifier_bytes": output_artifacts[
                    "silva-v4-classifier.qza"
                ].stat().st_size,
                "production_classifier_sha256": sha256_file(
                    output_artifacts["silva-v4-classifier.qza"]
                ),
            },
            "holdout_metrics": nested_holdout,
            "query_metrics": nested_query,
            "method_concordance": {
                row["rank"]: {
                    "comparable_asvs": row["comparable_asvs"],
                    "exact_match_asvs": row["exact_match_asvs"],
                    "exact_match_asv_fraction": row[
                        "exact_match_asv_fraction"
                    ],
                    "exact_match_read_fraction": row[
                        "exact_match_read_fraction"
                    ],
                }
                for row in method_rows
            },
            "sequence_search_contract": {
                "generated_within_step": not sequence_search_supplied,
                "method": (
                    "q2-feature-classifier classify-consensus-vsearch"
                ),
                "maxaccepts": 50,
                "maxrejects": "all",
                "identity": 0.8,
                "query_coverage": 0.8,
                "strand": "both",
                "top_hits_only": True,
                "maxhits": "all",
                "minimum_consensus": 0.51,
            },
            "source_release": {
                "release": source_summary["release"],
                "license": source_summary["license"],
                "reference_fasta_sha256": sha256_file(reference_fasta),
                "reference_taxonomy_sha256": sha256_file(reference_taxonomy),
            },
            "checks_total": len(checks),
            "checks_passed": checks_passed,
            "checks_failed": len(checks) - checks_passed,
        }
        write_json(output_dir / "classifier-training-summary.json", summary)

        validation_lines = [
            f"status={summary['status']}",
            f"checks={checks_passed}/{len(checks)}",
            f"reference_records={len(reference_sequences)}",
            f"training_records={len(training_ids)}",
            f"holdout_records={len(holdout_ids)}",
            f"query_asvs={len(query_sequences)}",
            f"query_reads={summary['query_reads']}",
            "",
            "FAILED CHECKS:",
        ]
        validation_lines.extend(
            row["check_id"] for row in checks if row["status"] != "PASS"
        )
        validation_lines.extend(
            [
                "",
                "INTERPRETATION BOUNDARY:",
                summary["split"]["interpretation_boundary"],
                (
                    "Naive Bayes versus sequence-search agreement is a "
                    "method-sensitivity audit, not proof that either label is true."
                ),
                "",
            ]
        )
        validation_log.write_text(
            "\n".join(str(line) for line in validation_lines),
            encoding="utf-8",
        )

    return 0 if checks_passed == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
