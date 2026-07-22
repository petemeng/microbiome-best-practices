#!/usr/bin/env python3
"""Build region-matched reference assets for Article 14.

This is a one-time preparation utility. It starts from immutable official
SILVA, GTDB, and Greengenes2 releases, derives the same 515F/806R V4 region,
collapses exact duplicate amplicons with least-common-ancestor taxonomy, and
writes deterministic gzip files plus source audits. The lightweight tutorial
validator consumes those derived files and never downloads a database.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


QIIME_ENV = "microbiome-qiime2-2026.4"
FORWARD_PRIMER = "GTGYCAGCMGCCGCGGTAA"
REVERSE_PRIMER = "GGACTACNVGGGTWTCTAAT"
EXTRACTION_IDENTITY = 0.80
MIN_LENGTH = 200
MAX_LENGTH = 400
READ_ORIENTATION = "both"
PREPARATION_DATE = "2026-07-20"

EXPECTED_MD5 = {
    "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz": (
        "c3b7638cab8bf2fe1e5cad01ce406a53"
    ),
    "ar53_ssu_reps_r232.fna.gz": "96dddfbdb5bebafd8e459c02785d5a02",
    "bac120_ssu_reps_r232.fna.gz": "3f0a53ff7fd8e3c56fc7246fe3b04342",
    "ar53_taxonomy_r232.tsv.gz": "9fa569a0804c46e4ba8dcd85a5bc5240",
    "bac120_taxonomy_r232.tsv.gz": "807e3221f5553142eb0544ce3247d079",
}

RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


@dataclass(frozen=True)
class DatabasePaths:
    key: str
    full_sequences: Path
    full_taxonomy: Path
    rank_handles: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--qiime-env", default=QIIME_ENV)
    parser.add_argument("--conda-executable", type=Path)
    parser.add_argument(
        "--silva-taxonomy-artifact",
        type=Path,
        help=(
            "Optional fixed-rank SILVA taxonomy produced by "
            "RESCRIPt get-silva-data. If omitted, RESCRIPt downloads only "
            "the taxonomy mapping files."
        ),
    )
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def deterministic_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as input_handle, destination.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_output,
            mtime=0,
            compresslevel=9,
        ) as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield identifier, description, and uppercase sequence."""
    identifier: str | None = None
    description = ""
    sequence_parts: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, description, "".join(sequence_parts).upper()
                description = line[1:]
                identifier = description.split()[0]
                sequence_parts = []
            else:
                if identifier is None:
                    raise ValueError(f"sequence before FASTA header in {path}")
                sequence_parts.append(line)
    if identifier is not None:
        yield identifier, description, "".join(sequence_parts).upper()


def quantile(sorted_values: Sequence[int], probability: float) -> float:
    if not sorted_values:
        raise ValueError("cannot calculate a quantile of an empty sequence")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = (len(sorted_values) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    fraction = position - lower
    return (
        float(sorted_values[lower]) * (1.0 - fraction)
        + float(sorted_values[upper]) * fraction
    )


def length_summary(lengths: Sequence[int]) -> dict[str, int | float]:
    ordered = sorted(lengths)
    return {
        "records": len(ordered),
        "minimum": ordered[0],
        "q1": round(quantile(ordered, 0.25), 2),
        "median": round(statistics.median(ordered), 2),
        "mean": round(statistics.fmean(ordered), 2),
        "q3": round(quantile(ordered, 0.75), 2),
        "maximum": ordered[-1],
    }


def parse_taxonomy(path: Path) -> dict[str, str]:
    taxonomy: dict[str, str] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader)
        has_header = (
            len(first) >= 2
            and first[0].strip().lower() in {"feature id", "feature-id"}
        )
        rows: Iterable[list[str]] = (
            reader if has_header else itertools.chain((first,), reader)
        )
        for row in rows:
            if len(row) < 2:
                continue
            feature_id = row[0].strip()
            taxon = row[1].strip()
            if feature_id in taxonomy:
                raise ValueError(f"duplicate taxonomy ID: {feature_id}")
            taxonomy[feature_id] = taxon
    return taxonomy


def taxonomy_summary(taxonomy: dict[str, str]) -> dict[str, object]:
    domains: Counter[str] = Counter()
    rank_nonempty: Counter[str] = Counter()
    organelle = 0
    for taxon in taxonomy.values():
        labels = [item.strip() for item in taxon.split(";")]
        for label in labels:
            for prefix in RANK_PREFIXES:
                if label.startswith(prefix) and label != prefix:
                    rank_nonempty[prefix] += 1
                    if prefix == "d__":
                        domains[label.removeprefix(prefix)] += 1
                    break
        if re.search(r"chloroplast|mitochond", taxon, flags=re.IGNORECASE):
            organelle += 1
    return {
        "records": len(taxonomy),
        "domain_counts": dict(sorted(domains.items())),
        "nonempty_rank_counts": {
            prefix: rank_nonempty[prefix] for prefix in RANK_PREFIXES
        },
        "organelle_label_records": organelle,
    }


class CommandRunner:
    def __init__(
        self,
        conda_executable: Path,
        qiime_env: str,
        environment: dict[str, str],
        log_path: Path,
    ) -> None:
        self.conda_executable = conda_executable
        self.qiime_env = qiime_env
        self.environment = environment
        self.log_path = log_path
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        self.log_path.write_text("", encoding="utf-8")

    def qiime(self, arguments: Sequence[str], timeout: int = 7200) -> None:
        command = [
            str(self.conda_executable),
            "run",
            "-n",
            self.qiime_env,
            "qiime",
            *arguments,
        ]
        started = time.monotonic()
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            env=self.environment,
            timeout=timeout,
        )
        elapsed = time.monotonic() - started
        public_command = " ".join(arguments)
        with self.log_path.open("a", encoding="utf-8") as log:
            log.write(f"$ qiime {public_command}\n")
            if completed.stdout:
                log.write(completed.stdout.rstrip() + "\n")
            if completed.stderr:
                log.write(completed.stderr.rstrip() + "\n")
            log.write(
                f"[exit={completed.returncode}; elapsed={elapsed:.2f}s]\n\n"
            )
        if completed.returncode != 0:
            raise RuntimeError(
                f"QIIME command failed ({completed.returncode}): "
                f"qiime {public_command}\n{completed.stderr}"
            )


def ensure_qiime_cache(runner: CommandRunner, cache_dir: Path) -> None:
    if (cache_dir / "VERSION").exists():
        return
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    runner.qiime(["tools", "cache-create", "--cache", str(cache_dir)])


def ensure_silva_taxonomy(
    runner: CommandRunner,
    supplied: Path | None,
    work_dir: Path,
    cache_dir: Path,
) -> Path:
    destination = work_dir / "silva-fixed-taxonomy.qza"
    if destination.exists():
        return destination
    if supplied is not None:
        shutil.copy2(supplied, destination)
        return destination
    empty_sequences = work_dir / "silva-empty-sequences.qza"
    runner.qiime(
        [
            "rescript",
            "get-silva-data",
            "--p-version",
            "138.2",
            "--p-target",
            "SSURef_NR99",
            "--p-no-include-species-labels",
            "--p-rank-propagation",
            "--p-no-download-sequences",
            "--o-silva-sequences",
            str(empty_sequences),
            "--o-silva-taxonomy",
            str(destination),
            "--use-cache",
            str(cache_dir),
        ]
    )
    return destination


def ensure_silva_sequences(
    runner: CommandRunner,
    raw_fasta_gz: Path,
    work_dir: Path,
    cache_dir: Path,
) -> Path:
    dna_artifact = work_dir / "silva-full-dna.qza"
    if dna_artifact.exists():
        return dna_artifact
    fasta_path = work_dir / "silva-full-rna.fasta"
    if not fasta_path.exists():
        with gzip.open(raw_fasta_gz, "rb") as source, fasta_path.open("wb") as dest:
            shutil.copyfileobj(source, dest, length=1024 * 1024)
    rna_artifact = work_dir / "silva-full-rna.qza"
    if not rna_artifact.exists():
        runner.qiime(
            [
                "tools",
                "import",
                "--input-path",
                str(fasta_path),
                "--output-path",
                str(rna_artifact),
                "--type",
                "FeatureData[RNASequence]",
            ]
        )
    runner.qiime(
        [
            "rescript",
            "reverse-transcribe",
            "--i-rna-sequences",
            str(rna_artifact),
            "--o-dna-sequences",
            str(dna_artifact),
            "--use-cache",
            str(cache_dir),
        ]
    )
    return dna_artifact


def prepare_gtdb_flat_files(
    raw_dir: Path,
    work_dir: Path,
) -> tuple[Path, Path, dict[str, object]]:
    fasta_output = work_dir / "gtdb-r232-full.fasta"
    taxonomy_output = work_dir / "gtdb-r232-full-taxonomy.tsv"
    sequence_sources = (
        raw_dir / "ar53_ssu_reps_r232.fna.gz",
        raw_dir / "bac120_ssu_reps_r232.fna.gz",
    )
    taxonomy_sources = (
        raw_dir / "ar53_taxonomy_r232.tsv.gz",
        raw_dir / "bac120_taxonomy_r232.tsv.gz",
    )

    header_taxonomy: dict[str, str] = {}
    lengths: list[int] = []
    sequences: list[tuple[str, str]] = []
    invalid_bases = 0
    for source in sequence_sources:
        for feature_id, description, sequence in iter_fasta(source):
            if feature_id in header_taxonomy:
                raise ValueError(f"duplicate GTDB SSU representative: {feature_id}")
            match = re.match(r"^\S+\s+(.+?)\s+\[locus_tag=", description)
            if match is None:
                raise ValueError(f"cannot parse GTDB FASTA header: {description}")
            header_taxonomy[feature_id] = match.group(1)
            sequence = sequence.replace("U", "T")
            invalid_bases += sum(base not in "ACGTRYSWKMBDHVN" for base in sequence)
            sequences.append((feature_id, sequence))
            lengths.append(len(sequence))

    canonical_taxonomy: dict[str, str] = {}
    taxonomy_rows = 0
    header_disagreements = 0
    for source in taxonomy_sources:
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            for raw_line in handle:
                taxonomy_rows += 1
                feature_id, taxon = raw_line.rstrip("\n").split("\t", 1)
                if feature_id in header_taxonomy:
                    canonical_taxonomy[feature_id] = taxon
                    if taxon != header_taxonomy[feature_id]:
                        header_disagreements += 1

    missing = sorted(set(header_taxonomy) - set(canonical_taxonomy))
    if missing:
        raise ValueError(
            f"{len(missing)} GTDB SSU representatives lack taxonomy; "
            f"first={missing[0]}"
        )
    if header_disagreements:
        raise ValueError(
            f"{header_disagreements} GTDB FASTA headers disagree with taxonomy files"
        )
    if invalid_bases:
        raise ValueError(f"GTDB references contain {invalid_bases} invalid bases")

    with fasta_output.open("w", encoding="ascii", newline="\n") as handle:
        for feature_id, sequence in sequences:
            handle.write(f">{feature_id}\n{sequence}\n")
    with taxonomy_output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Feature ID", "Taxon"))
        for feature_id, _sequence in sequences:
            writer.writerow((feature_id, canonical_taxonomy[feature_id]))

    audit = {
        "full_sequence_lengths": length_summary(lengths),
        "full_taxonomy_rows": taxonomy_rows,
        "representative_taxonomy": taxonomy_summary(canonical_taxonomy),
        "header_taxonomy_disagreements": header_disagreements,
        "invalid_sequence_bases": invalid_bases,
    }
    return fasta_output, taxonomy_output, audit


def ensure_import(
    runner: CommandRunner,
    input_path: Path,
    output_path: Path,
    semantic_type: str,
    input_format: str | None = None,
) -> Path:
    if output_path.exists():
        return output_path
    arguments = [
        "tools",
        "import",
        "--input-path",
        str(input_path),
        "--output-path",
        str(output_path),
        "--type",
        semantic_type,
    ]
    if input_format is not None:
        arguments.extend(["--input-format", input_format])
    runner.qiime(arguments)
    return output_path


def ensure_export(
    runner: CommandRunner,
    artifact: Path,
    output_dir: Path,
) -> Path:
    if output_dir.exists() and any(output_dir.iterdir()):
        return output_dir
    if output_dir.exists():
        shutil.rmtree(output_dir)
    runner.qiime(
        [
            "tools",
            "export",
            "--input-path",
            str(artifact),
            "--output-path",
            str(output_dir),
        ]
    )
    return output_dir


def locate_export(output_dir: Path, expected_names: Sequence[str]) -> Path:
    for name in expected_names:
        candidate = output_dir / name
        if candidate.exists():
            return candidate
    files = [path for path in output_dir.rglob("*") if path.is_file()]
    raise FileNotFoundError(
        f"none of {expected_names} found in {output_dir}; files={files}"
    )


def write_taxonomy_subset(
    full_taxonomy: Path,
    identifiers: set[str],
    destination: Path,
) -> dict[str, str]:
    taxonomy = parse_taxonomy(full_taxonomy)
    missing = sorted(identifiers - set(taxonomy))
    if missing:
        raise ValueError(
            f"{len(missing)} extracted sequence IDs lack taxonomy; first={missing[0]}"
        )
    selected = {feature_id: taxonomy[feature_id] for feature_id in identifiers}
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Feature ID", "Taxon"))
        for feature_id in sorted(selected):
            writer.writerow((feature_id, selected[feature_id]))
    return selected


def normalize_fasta(source: Path, destination: Path) -> list[int]:
    lengths: list[int] = []
    with destination.open("w", encoding="ascii", newline="\n") as handle:
        for feature_id, _description, sequence in iter_fasta(source):
            sequence = sequence.replace("U", "T")
            handle.write(f">{feature_id}\n{sequence}\n")
            lengths.append(len(sequence))
    return lengths


def normalize_taxonomy(source: Path, destination: Path) -> dict[str, str]:
    taxonomy = parse_taxonomy(source)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Feature ID", "Taxon"))
        for feature_id in sorted(taxonomy):
            taxon = "; ".join(
                part.strip() for part in taxonomy[feature_id].split(";")
            )
            writer.writerow((feature_id, taxon))
    return taxonomy


def derive_region_reference(
    runner: CommandRunner,
    database: DatabasePaths,
    work_dir: Path,
    output_dir: Path,
    cache_dir: Path,
    threads: int,
) -> dict[str, object]:
    prefix = database.key
    extracted_artifact = work_dir / f"{prefix}-v4-extracted.qza"
    if not extracted_artifact.exists():
        runner.qiime(
            [
                "feature-classifier",
                "extract-reads",
                "--i-sequences",
                str(database.full_sequences),
                "--p-f-primer",
                FORWARD_PRIMER,
                "--p-r-primer",
                REVERSE_PRIMER,
                "--p-identity",
                str(EXTRACTION_IDENTITY),
                "--p-min-length",
                str(MIN_LENGTH),
                "--p-max-length",
                str(MAX_LENGTH),
                "--p-n-jobs",
                str(threads),
                "--p-read-orientation",
                READ_ORIENTATION,
                "--o-reads",
                str(extracted_artifact),
                "--use-cache",
                str(cache_dir),
            ],
            timeout=14400,
        )

    extract_export_dir = ensure_export(
        runner, extracted_artifact, work_dir / f"{prefix}-v4-extracted-export"
    )
    extracted_fasta = locate_export(
        extract_export_dir, ("dna-sequences.fasta", "rna-sequences.fasta")
    )
    extracted_records = list(iter_fasta(extracted_fasta))
    extracted_ids = {record[0] for record in extracted_records}
    extracted_lengths = [len(record[2]) for record in extracted_records]

    region_taxonomy_tsv = work_dir / f"{prefix}-v4-extracted-taxonomy.tsv"
    selected_taxonomy = write_taxonomy_subset(
        database.full_taxonomy, extracted_ids, region_taxonomy_tsv
    )
    region_taxonomy_artifact = ensure_import(
        runner,
        region_taxonomy_tsv,
        work_dir / f"{prefix}-v4-extracted-taxonomy.qza",
        "FeatureData[Taxonomy]",
        "TSVTaxonomyFormat",
    )

    derep_sequences = work_dir / f"{prefix}-v4-lca-sequences.qza"
    derep_taxonomy = work_dir / f"{prefix}-v4-lca-taxonomy.qza"
    if not derep_sequences.exists() or not derep_taxonomy.exists():
        runner.qiime(
            [
                "rescript",
                "dereplicate",
                "--i-sequences",
                str(extracted_artifact),
                "--i-taxa",
                str(region_taxonomy_artifact),
                "--p-mode",
                "lca",
                "--p-perc-identity",
                "1.0",
                "--p-threads",
                str(threads),
                "--p-rank-handles",
                *database.rank_handles,
                "--o-dereplicated-sequences",
                str(derep_sequences),
                "--o-dereplicated-taxa",
                str(derep_taxonomy),
                "--use-cache",
                str(cache_dir),
            ],
            timeout=14400,
        )

    derep_seq_export = ensure_export(
        runner, derep_sequences, work_dir / f"{prefix}-v4-lca-sequences-export"
    )
    derep_tax_export = ensure_export(
        runner, derep_taxonomy, work_dir / f"{prefix}-v4-lca-taxonomy-export"
    )
    source_fasta = locate_export(
        derep_seq_export, ("dna-sequences.fasta", "rna-sequences.fasta")
    )
    source_taxonomy = locate_export(derep_tax_export, ("taxonomy.tsv",))

    normalized_fasta = work_dir / f"{prefix}-v4.fasta"
    normalized_taxonomy = work_dir / f"{prefix}-v4-taxonomy.tsv"
    final_lengths = normalize_fasta(source_fasta, normalized_fasta)
    final_taxonomy = normalize_taxonomy(source_taxonomy, normalized_taxonomy)
    final_ids = {record[0] for record in iter_fasta(normalized_fasta)}
    if final_ids != set(final_taxonomy):
        raise ValueError(f"{prefix} final sequence/taxonomy IDs do not match")

    fasta_gz = output_dir / f"{prefix}-v4.fasta.gz"
    taxonomy_gz = output_dir / f"{prefix}-v4-taxonomy.tsv.gz"
    deterministic_gzip(normalized_fasta, fasta_gz)
    deterministic_gzip(normalized_taxonomy, taxonomy_gz)

    return {
        "extraction_contract": {
            "marker": "16S V4",
            "forward_primer_5_to_3": FORWARD_PRIMER,
            "reverse_primer_5_to_3": REVERSE_PRIMER,
            "combined_primer_identity": EXTRACTION_IDENTITY,
            "minimum_amplicon_length": MIN_LENGTH,
            "maximum_amplicon_length": MAX_LENGTH,
            "read_orientation": READ_ORIENTATION,
        },
        "dereplication_contract": {
            "mode": "least common ancestor",
            "sequence_identity": 1.0,
            "reason": (
                "Prevent identical V4 copies from weighting consensus by "
                "database redundancy while preserving label ambiguity."
            ),
        },
        "raw_v4_records": len(extracted_records),
        "raw_v4_length": length_summary(extracted_lengths),
        "raw_v4_taxonomy": taxonomy_summary(selected_taxonomy),
        "lca_v4_records": len(final_ids),
        "lca_v4_length": length_summary(final_lengths),
        "lca_v4_taxonomy": taxonomy_summary(final_taxonomy),
        "derived_files": {
            "sequences": fasta_gz.name,
            "sequences_sha256": sha256_file(fasta_gz),
            "taxonomy": taxonomy_gz.name,
            "taxonomy_sha256": sha256_file(taxonomy_gz),
        },
    }


def validate_raw_checksums(raw_dir: Path) -> dict[str, dict[str, object]]:
    audit: dict[str, dict[str, object]] = {}
    for name, expected_md5 in EXPECTED_MD5.items():
        path = raw_dir / name
        observed_md5 = md5_file(path)
        if observed_md5 != expected_md5:
            raise ValueError(
                f"MD5 mismatch for {name}: {observed_md5} != {expected_md5}"
            )
        audit[name] = {
            "bytes": path.stat().st_size,
            "md5": observed_md5,
            "sha256": sha256_file(path),
        }
    for name in (
        "2024.09.backbone.full-length.fna.qza",
        "2024.09.backbone.tax.qza",
    ):
        path = raw_dir / name
        audit[name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return audit


def summarize_fasta(path: Path) -> dict[str, object]:
    identifiers: set[str] = set()
    lengths: list[int] = []
    invalid_bases = 0
    for feature_id, _description, sequence in iter_fasta(path):
        if feature_id in identifiers:
            raise ValueError(f"duplicate FASTA ID in {path}: {feature_id}")
        identifiers.add(feature_id)
        lengths.append(len(sequence))
        invalid_bases += sum(base not in "ACGTURYSWKMBDHVN" for base in sequence)
    return {
        "length": length_summary(lengths),
        "unique_identifiers": len(identifiers),
        "invalid_sequence_bases": invalid_bases,
    }


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    raw_dir = (
        args.raw_dir.resolve()
        if args.raw_dir
        else project_root / "data/raw/taxonomy-databases"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else project_root / "data/small/taxonomy-databases"
    )
    work_dir = args.work_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    conda_executable = args.conda_executable
    if conda_executable is None:
        discovered = shutil.which("conda")
        if discovered is None:
            raise FileNotFoundError(
                "conda executable not found; pass --conda-executable"
            )
        conda_executable = Path(discovered)
    conda_executable = conda_executable.resolve()

    cache_root = project_root / ".cache"
    xdg_cache = cache_root / "xdg"
    matplotlib_cache = cache_root / "matplotlib"
    numba_cache = cache_root / "numba"
    qiime_cache = cache_root / "qiime2-article14-prep"
    for directory in (xdg_cache, matplotlib_cache, numba_cache):
        directory.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(
        {
            "XDG_CACHE_HOME": str(xdg_cache),
            "MPLCONFIGDIR": str(matplotlib_cache),
            "NUMBA_CACHE_DIR": str(numba_cache),
        }
    )
    runner = CommandRunner(
        conda_executable,
        args.qiime_env,
        environment,
        work_dir / "reference-preparation.log",
    )
    ensure_qiime_cache(runner, qiime_cache)
    raw_checksums = validate_raw_checksums(raw_dir)

    silva_taxonomy_artifact = ensure_silva_taxonomy(
        runner,
        args.silva_taxonomy_artifact.resolve()
        if args.silva_taxonomy_artifact
        else None,
        work_dir,
        qiime_cache,
    )
    silva_tax_export = ensure_export(
        runner, silva_taxonomy_artifact, work_dir / "silva-taxonomy-export"
    )
    silva_taxonomy_tsv = locate_export(silva_tax_export, ("taxonomy.tsv",))
    silva_sequence_artifact = ensure_silva_sequences(
        runner,
        raw_dir / "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz",
        work_dir,
        qiime_cache,
    )

    gtdb_fasta, gtdb_taxonomy, gtdb_flat_audit = prepare_gtdb_flat_files(
        raw_dir, work_dir
    )
    gtdb_sequence_artifact = ensure_import(
        runner,
        gtdb_fasta,
        work_dir / "gtdb-r232-full-sequences.qza",
        "FeatureData[Sequence]",
        "DNAFASTAFormat",
    )
    gtdb_taxonomy_artifact = ensure_import(
        runner,
        gtdb_taxonomy,
        work_dir / "gtdb-r232-full-taxonomy.qza",
        "FeatureData[Taxonomy]",
        "TSVTaxonomyFormat",
    )

    gg2_sequences_artifact = raw_dir / "2024.09.backbone.full-length.fna.qza"
    gg2_taxonomy_artifact = raw_dir / "2024.09.backbone.tax.qza"
    runner.qiime(
        ["tools", "validate", str(gg2_sequences_artifact), "--level", "max"]
    )
    runner.qiime(
        ["tools", "validate", str(gg2_taxonomy_artifact), "--level", "max"]
    )
    gg2_seq_export = ensure_export(
        runner, gg2_sequences_artifact, work_dir / "gg2-full-sequences-export"
    )
    gg2_tax_export = ensure_export(
        runner, gg2_taxonomy_artifact, work_dir / "gg2-full-taxonomy-export"
    )
    gg2_fasta = locate_export(gg2_seq_export, ("dna-sequences.fasta",))
    gg2_taxonomy = locate_export(gg2_tax_export, ("taxonomy.tsv",))

    database_paths = {
        "silva-138.2": DatabasePaths(
            key="silva-138.2",
            full_sequences=silva_sequence_artifact,
            full_taxonomy=silva_taxonomy_tsv,
            rank_handles=("domain", "phylum", "class", "order", "family", "genus"),
        ),
        "gtdb-r232": DatabasePaths(
            key="gtdb-r232",
            full_sequences=gtdb_sequence_artifact,
            full_taxonomy=gtdb_taxonomy,
            rank_handles=(
                "domain",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ),
        ),
        "greengenes2-2024.09": DatabasePaths(
            key="greengenes2-2024.09",
            full_sequences=gg2_sequences_artifact,
            full_taxonomy=gg2_taxonomy,
            rank_handles=(
                "domain",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ),
        ),
    }

    region_audits = {
        key: derive_region_reference(
            runner,
            paths,
            work_dir,
            output_dir,
            qiime_cache,
            args.threads,
        )
        for key, paths in database_paths.items()
    }

    silva_taxonomy = parse_taxonomy(silva_taxonomy_tsv)
    gg2_taxonomy_map = parse_taxonomy(gg2_taxonomy)
    common_contract = {
        "preparation_date": PREPARATION_DATE,
        "qiime_environment": args.qiime_env,
        "qiime_release": "2026.4",
        "q2_feature_classifier": "2026.4.0",
        "rescript": "2026.4.0",
        "region_processing": region_audits["silva-138.2"][
            "extraction_contract"
        ],
        "exact_sequence_dereplication": region_audits["silva-138.2"][
            "dereplication_contract"
        ],
    }

    silva_summary = {
        "database": "SILVA",
        "release": "138.2",
        "reference_subset": "SSURef NR99",
        "release_date": "2024-07-11",
        "official_url": "https://www.arb-silva.de/documentation/release-1382/",
        "source_url": (
            "https://ftp.arb-silva.de/release_138.2/Exports/"
            "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
        ),
        "license": "CC BY 4.0",
        "taxonomy_contract": (
            "RESCRIPt fixed six ranks (domain through genus), rank propagation "
            "enabled, species labels intentionally excluded"
        ),
        "source_files": {
            name: audit
            for name, audit in raw_checksums.items()
            if name.startswith("SILVA")
        },
        "full_sequences": summarize_fasta(
            raw_dir / "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
        ),
        "full_taxonomy": taxonomy_summary(silva_taxonomy),
        **common_contract,
        **region_audits["silva-138.2"],
        "interpretation_boundary": (
            "SILVA species text in raw FASTA descriptions was not promoted to "
            "a fixed species rank; short V4 reads do not prove species."
        ),
    }
    gtdb_summary = {
        "database": "GTDB",
        "release": "R11-RS232",
        "reference_subset": "species-representative SSU",
        "release_announcement_date": "2026-04-15",
        "official_url": (
            "https://data.gtdb.ecogenomic.org/releases/release232/232.0/"
        ),
        "license": "CC BY-SA 4.0",
        "release_scope": {
            "genomes": 901341,
            "species_clusters": 199923,
        },
        "source_files": {
            name: audit
            for name, audit in raw_checksums.items()
            if "r232" in name
        },
        "full_reference_audit": gtdb_flat_audit,
        "method_documentation_boundary": {
            "release_notes": (
                "R232 changed 16S identification/extraction to mirror Prokka "
                "barnap and set a 396 bp minimum."
            ),
            "legacy_generic_files": (
                "Bundled METHODS and FILE_DESCRIPTIONS still describe nhmmer "
                "and a 200 bp threshold; actual FASTA lengths are audited."
            ),
        },
        **common_contract,
        **region_audits["gtdb-r232"],
        "interpretation_boundary": (
            "Taxonomy is assigned from the representative genome phylogeny; "
            "the extracted 16S copy may be absent, partial, contaminated, or "
            "incongruent with the genome taxonomy."
        ),
    }
    gg2_summary = {
        "database": "Greengenes2",
        "release": "2024.09",
        "reference_subset": "full-length backbone",
        "release_date": "2024-09-26",
        "official_url": (
            "https://ftp.microbio.me/greengenes_release/2024.09/"
        ),
        "license": "BSD-3-Clause",
        "taxonomy_sources": {
            "gtdb": "R220",
            "ltp": "08/2023",
            "organelle_sequences": "SILVA 138.1",
        },
        "construction_scope": (
            "Web of Life genome backbone updated with full-length 16S, with "
            "large-scale Qiita V4 ASVs placed by DEPP in the complete resource"
        ),
        "source_files": {
            name: audit
            for name, audit in raw_checksums.items()
            if name.startswith("2024.09")
        },
        "full_sequences": summarize_fasta(gg2_fasta),
        "full_taxonomy": taxonomy_summary(gg2_taxonomy_map),
        **common_contract,
        **region_audits["greengenes2-2024.09"],
        "interpretation_boundary": (
            "This comparison uses the official full-length backbone and "
            "derives the same V4 locus; it does not substitute the legacy "
            "Greengenes 13_8 database."
        ),
    }

    write_json(output_dir / "silva-138.2-source-summary.json", silva_summary)
    write_json(output_dir / "gtdb-r232-source-summary.json", gtdb_summary)
    write_json(
        output_dir / "greengenes2-2024.09-source-summary.json", gg2_summary
    )
    print(
        json.dumps(
            {
                "status": "passed",
                "output_dir": str(output_dir),
                "databases": {
                    key: audit["lca_v4_records"]
                    for key, audit in region_audits.items()
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
