#!/usr/bin/env python3
"""Run and audit reference-based SEPP placement for Article 16.

The validator deliberately starts from the checksum-locked real V4 ASV FASTA
and aggregate abundance table.  It imports both inputs, runs q2-fragment-
insertion against the official QIIME 2 Greengenes 13_8 SEPP reference,
filters the feature table to placed tips, validates every Artifact at maximum
level, reconciles IDs across FASTA/tree/jplace/BIOM, and draws the evidence
figures used by the article.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
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
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-sepp-phylogeny")

import dendropy
import matplotlib


matplotlib.use("Agg")


QIIME_ENV = "microbiome-qiime2-2026.4"
REFERENCE_SHA256 = (
    "e252b83d7d5fbf2a9e14e594768e3578b33b557c34e22b6abc83b324689b1360"
)
QUERY_FASTA_SHA256 = (
    "fa71915be8007d36ec70f85d28401f5a94d0d33eca323782cc937cfc4bcf1af4"
)
QUERY_ABUNDANCE_SHA256 = (
    "a80e3b94d8632c48fb6c01e788e842c022c3b1ee64b79dc94456618328895e6e"
)
REFERENCE_UUID = "a14c6180-506b-4ecb-bacb-9cb30bc3044b"
REFERENCE_TIPS = 203_452
REFERENCE_ALIGNMENT_LENGTH = 1_285
QUERY_ASVS = 366
QUERY_READS = 5_213
ALIGNMENT_SUBSET_SIZE = 1_000
PLACEMENT_SUBSET_SIZE = 5_000
THREADS = 8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--reference-database", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--qiime-env", default=QIIME_ENV)
    parser.add_argument("--threads", type=int, default=THREADS)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
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


def read_abundance(path: Path) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        expected = ["Feature ID", "Frequency", "No. of Samples Observed In"]
        if header != expected:
            raise ValueError(f"unexpected abundance header: {header}")
        for row in reader:
            if row and row[0] == "#q2:types":
                continue
            if len(row) != 3:
                raise ValueError(f"unexpected abundance row: {row}")
            feature_id = row[0]
            if feature_id in result:
                raise ValueError(f"duplicate abundance feature ID: {feature_id}")
            result[feature_id] = {
                "frequency": int(float(row[1])),
                "samples_observed": int(float(row[2])),
            }
    return result


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for feature_id in sorted(records):
            handle.write(f">{feature_id}\n{records[feature_id]}\n")


def write_biom_source_tsv(path: Path, abundance: dict[str, dict[str, int]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["# Constructed from real aggregate ASV frequencies"])
        writer.writerow(["#OTU ID", "AllReads"])
        for feature_id in sorted(abundance):
            writer.writerow([feature_id, abundance[feature_id]["frequency"]])


def conda_command(environment: str, executable: str, *args: str) -> list[str]:
    return ["conda", "run", "-n", environment, executable, *args]


def isolated_environment(root: Path) -> dict[str, str]:
    environment = os.environ.copy()
    for name in ("xdg", "numba", "matplotlib", "fontconfig"):
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
                    self.sanitize(
                        str(record.get("logged_stdout", record["stdout"]))
                    ),
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
        log_stdout: bool = True,
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
        record: dict[str, object] = {
            "label": label,
            "command": list(command),
            "returncode": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
            "elapsed_seconds": round(time.monotonic() - started, 3),
        }
        if not log_stdout:
            record["logged_stdout"] = (
                f"<captured {len(completed.stdout)} characters; "
                "parsed into a structured audit>"
            )
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


def reference_provenance(path: Path) -> dict[str, str]:
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if name.count("/") == 3 and name.endswith("/provenance/action/action.yaml")
        ]
        if len(members) != 1:
            raise ValueError("reference root provenance action is ambiguous")
        text = archive.read(members[0]).decode("utf-8", errors="replace")
    framework = re.search(
        r"framework:\s*\n\s+version:\s*([^\n]+)", text, flags=re.MULTILINE
    )
    plugin = re.search(
        r"fragment-insertion:\s*\n\s+version:\s*([^\n]+)",
        text,
        flags=re.MULTILINE,
    )
    sepp = re.search(r"^\s+sepp:\s*([^\n]+)$", text, flags=re.MULTILINE)
    start = re.search(r"^\s+start:\s*([^\n]+)$", text, flags=re.MULTILINE)
    return {
        "framework_version": framework.group(1).strip() if framework else "",
        "plugin_version": plugin.group(1).strip() if plugin else "",
        "sepp_version": sepp.group(1).strip() if sepp else "",
        "import_started": start.group(1).strip() if start else "",
    }


def audit_reference(path: Path) -> tuple[dict[str, object], set[str]]:
    metadata = artifact_metadata(path)
    provenance = reference_provenance(path)
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        alignment_name = next(
            name
            for name in names
            if name.endswith("/data/aligned-dna-sequences.fasta")
        )
        tree_name = next(
            name for name in names if name.endswith("/data/tree.nwk")
        )
        alignment_ids: set[str] = set()
        lengths: list[int] = []
        current_length = 0
        with archive.open(alignment_name) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(b">"):
                    if alignment_ids:
                        lengths.append(current_length)
                    identifier = line[1:].split(maxsplit=1)[0].decode("ascii")
                    if identifier in alignment_ids:
                        raise ValueError(f"duplicate reference ID: {identifier}")
                    alignment_ids.add(identifier)
                    current_length = 0
                else:
                    current_length += len(line)
        if alignment_ids:
            lengths.append(current_length)
        tree_text = archive.read(tree_name).decode("utf-8")

    tree = dendropy.Tree.get(
        data=tree_text,
        schema="newick",
        preserve_underscores=True,
    )
    tip_ids = {
        node.taxon.label
        for node in tree.leaf_node_iter()
        if node.taxon is not None
    }
    edge_lengths = [
        float(node.edge_length)
        for node in tree.preorder_node_iter()
        if node.edge_length is not None
    ]
    summary: dict[str, object] = {
        "source_url": (
            "https://data.qiime2.org/classifiers/sepp-ref-dbs/"
            "sepp-refs-gg-13-8.qza"
        ),
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "uuid": metadata["uuid"],
        "semantic_type": metadata["type"],
        "format": metadata["format"],
        "artifact_framework_version": provenance["framework_version"],
        "artifact_plugin_version": provenance["plugin_version"],
        "artifact_sepp_version": provenance["sepp_version"],
        "artifact_import_started": provenance["import_started"],
        "alignment_records": len(alignment_ids),
        "alignment_length_min": min(lengths),
        "alignment_length_max": max(lengths),
        "alignment_unique_lengths": len(set(lengths)),
        "tree_tips": len(tip_ids),
        "tree_nodes": sum(1 for _ in tree.preorder_node_iter()),
        "root_degree": len(tree.seed_node.child_nodes()),
        "alignment_tree_id_match": alignment_ids == tip_ids,
        "branch_length_min": min(edge_lengths),
        "branch_length_max": max(edge_lengths),
        "negative_branch_lengths": sum(value < 0 for value in edge_lengths),
    }
    return summary, tip_ids


def export_artifact(
    runner: CommandRunner,
    environment: str,
    artifact: Path,
    destination: Path,
    label: str,
) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    runner.qiime(
        label,
        environment,
        "tools",
        "export",
        "--input-path",
        str(artifact),
        "--output-path",
        str(destination),
    )


def inspect_biom(
    runner: CommandRunner,
    environment: str,
    path: Path,
    label: str,
) -> dict[str, object]:
    code = """
import json
import sys
from biom import load_table

table = load_table(sys.argv[1])
observations = list(map(str, table.ids(axis='observation')))
samples = list(map(str, table.ids(axis='sample')))
values = {
    observation: [float(value) for value in table.data(
        observation, axis='observation'
    )]
    for observation in observations
}
print(json.dumps({
    'observation_ids': observations,
    'sample_ids': samples,
    'values': values,
    'total': float(table.sum()),
}, sort_keys=True))
"""
    record = runner.run(
        label,
        conda_command(environment, "python", "-c", code, str(path)),
        log_stdout=False,
    )
    return json.loads(str(record["stdout"]))


def write_exported_table(path: Path, payload: dict[str, object]) -> None:
    samples = [str(item) for item in payload["sample_ids"]]
    values = payload["values"]
    if not isinstance(values, dict):
        raise TypeError("BIOM values payload must be a mapping")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Feature ID", *samples])
        for feature_id in sorted(values):
            row = values[feature_id]
            writer.writerow(
                [feature_id, *[int(value) if float(value).is_integer() else value for value in row]]
            )


def parse_placements(
    path: Path,
    abundance: dict[str, dict[str, int]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    fields = [str(item) for item in payload["fields"]]
    field_index = {field: index for index, field in enumerate(fields)}
    required = {
        "edge_num",
        "likelihood",
        "like_weight_ratio",
        "distal_length",
        "pendant_length",
    }
    if set(fields) != required:
        raise ValueError(f"unexpected placement fields: {fields}")

    rows: list[dict[str, object]] = []
    entry_count = 0
    for placement in payload["placements"]:
        entry_count += 1
        names = placement.get("n") or placement.get("nm") or []
        candidates = placement["p"]
        best = max(
            candidates,
            key=lambda item: float(item[field_index["like_weight_ratio"]]),
        )
        ordered_weights = sorted(
            (
                float(item[field_index["like_weight_ratio"]])
                for item in candidates
            ),
            reverse=True,
        )
        for name_record in names:
            feature_id = (
                str(name_record[0])
                if isinstance(name_record, list)
                else str(name_record)
            )
            rows.append(
                {
                    "feature_id": feature_id,
                    "frequency": abundance.get(feature_id, {}).get("frequency", 0),
                    "candidate_placements": len(candidates),
                    "top_edge_num": int(best[field_index["edge_num"]]),
                    "top_likelihood": float(best[field_index["likelihood"]]),
                    "top_likelihood_weight": float(
                        best[field_index["like_weight_ratio"]]
                    ),
                    "second_likelihood_weight": (
                        ordered_weights[1] if len(ordered_weights) > 1 else 0.0
                    ),
                    "likelihood_weight_sum": sum(ordered_weights),
                    "distal_length": float(best[field_index["distal_length"]]),
                    "pendant_length": float(best[field_index["pendant_length"]]),
                }
            )
    rows.sort(key=lambda row: str(row["feature_id"]))
    summary = {
        "placement_entries": entry_count,
        "placement_names": len(rows),
        "fields": fields,
        "tree_string_present": bool(payload.get("tree")),
        "format_version": payload.get("version"),
    }
    return rows, summary


def quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction


def weighted_quantile(
    values: Sequence[float],
    weights: Sequence[int],
    probability: float,
) -> float:
    ordered = sorted(zip(values, weights), key=lambda item: float(item[0]))
    target = probability * sum(int(weight) for _value, weight in ordered)
    cumulative = 0
    for value, weight in ordered:
        cumulative += int(weight)
        if cumulative >= target:
            return float(value)
    return float(ordered[-1][0])


def summarize_placement_rows(rows: list[dict[str, object]]) -> dict[str, object]:
    weights = [int(row["frequency"]) for row in rows]
    summary: dict[str, object] = {}
    for key in ("top_likelihood_weight", "pendant_length", "distal_length"):
        values = [float(row[key]) for row in rows]
        summary[key] = {
            "min": min(values),
            "q25": quantile(values, 0.25),
            "median": statistics.median(values),
            "q75": quantile(values, 0.75),
            "q95": quantile(values, 0.95),
            "max": max(values),
            "read_weighted_median": weighted_quantile(values, weights, 0.50),
        }

    support_bins = [
        ("<0.50", 0.0, 0.50),
        ("0.50-<0.90", 0.50, 0.90),
        (">=0.90", 0.90, 1.0000001),
    ]
    support: dict[str, dict[str, int]] = {}
    for label, lower, upper in support_bins:
        selected = [
            row
            for row in rows
            if lower <= float(row["top_likelihood_weight"]) < upper
        ]
        support[label] = {
            "asvs": len(selected),
            "reads": sum(int(row["frequency"]) for row in selected),
        }
    summary["support_bins"] = support

    alternatives: dict[str, dict[str, int]] = {}
    for candidate_count in sorted(
        {int(row["candidate_placements"]) for row in rows}
    ):
        selected = [
            row
            for row in rows
            if int(row["candidate_placements"]) == candidate_count
        ]
        alternatives[str(candidate_count)] = {
            "asvs": len(selected),
            "reads": sum(int(row["frequency"]) for row in selected),
        }
    summary["candidate_placement_counts"] = alternatives
    return summary


def parse_output_tree(
    path: Path,
    query_ids: set[str],
    reference_ids: set[str],
    query_tree_path: Path,
) -> tuple[dict[str, object], dendropy.Tree, dict[str, float]]:
    tree = dendropy.Tree.get(
        path=str(path),
        schema="newick",
        preserve_underscores=True,
    )
    tips = {
        node.taxon.label
        for node in tree.leaf_node_iter()
        if node.taxon is not None
    }
    nodes = list(tree.preorder_node_iter())
    edge_lengths = [
        float(node.edge_length)
        for node in nodes
        if node.edge_length is not None
    ]
    query_tip_lengths = {
        node.taxon.label: float(node.edge_length or 0.0)
        for node in tree.leaf_node_iter()
        if node.taxon is not None and node.taxon.label in query_ids
    }
    query_tree = tree.extract_tree_with_taxa_labels(
        labels=query_ids,
        suppress_unifurcations=True,
    )
    query_tree.write(
        path=str(query_tree_path),
        schema="newick",
        suppress_rooting=False,
        unquoted_underscores=True,
    )
    summary: dict[str, object] = {
        "tips": len(tips),
        "nodes": len(nodes),
        "root_degree": len(tree.seed_node.child_nodes()),
        "edges_with_lengths": len(edge_lengths),
        "edges_missing_lengths": sum(
            node.edge_length is None for node in nodes
        ),
        "zero_length_edges": sum(value == 0 for value in edge_lengths),
        "negative_length_edges": sum(value < 0 for value in edge_lengths),
        "branch_length_min": min(edge_lengths),
        "branch_length_max": max(edge_lengths),
        "query_tips": len(tips & query_ids),
        "query_tips_missing": len(query_ids - tips),
        "reference_tips_preserved": len(tips & reference_ids),
        "unexpected_tips": len(tips - query_ids - reference_ids),
        "query_only_tips": sum(1 for _ in query_tree.leaf_node_iter()),
        "query_only_nodes": sum(1 for _ in query_tree.preorder_node_iter()),
        "query_only_root_degree": len(query_tree.seed_node.child_nodes()),
    }
    return summary, query_tree, query_tip_lengths


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


def draw_workflow(figure_dir: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    configure_plotting()
    fig, ax = plt.subplots(figsize=(11.6, 5.6))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7)
    ax.axis("off")
    nodes = [
        (0.3, 4.7, 2.25, 1.25, "Real V4 ASVs\n366 sequences", "#DBEAFE"),
        (3.1, 4.7, 2.25, 1.25, "Reference alignment\n203,452 full SSUs", "#E0E7FF"),
        (5.9, 4.7, 2.25, 1.25, "Profile-HMM\nfragment alignment", "#EDE9FE"),
        (8.7, 4.7, 2.25, 1.25, "pplacer on fixed\nreference backbone", "#FCE7F3"),
        (11.5, 4.7, 2.25, 1.25, "Rooted insertion tree\n+ jplace audit", "#FEF3C7"),
        (8.7, 2.0, 2.25, 1.25, "Filter feature table\nto placed ASVs", "#D1FAE5"),
        (11.5, 2.0, 2.25, 1.25, "Faith PD / UniFrac\ncompatible inputs", "#CCFBF1"),
    ]
    for x, y, width, height, label, color in nodes:
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.04,rounding_size=0.08",
            linewidth=1.1,
            edgecolor="#374151",
            facecolor=color,
        )
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height / 2, label, ha="center", va="center")
    arrows = [
        ((2.55, 5.33), (3.1, 5.33)),
        ((5.35, 5.33), (5.9, 5.33)),
        ((8.15, 5.33), (8.7, 5.33)),
        ((10.95, 5.33), (11.5, 5.33)),
        ((12.62, 4.7), (9.82, 3.25)),
        ((10.95, 2.63), (11.5, 2.63)),
    ]
    for start, end in arrows:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#4B5563"},
        )
    ax.text(
        7,
        6.55,
        "Reference placement keeps short amplicons on a fixed, rooted backbone",
        ha="center",
        va="center",
        fontsize=14,
        weight="bold",
        color="#111827",
    )
    ax.text(
        7,
        0.65,
        "Placement success, support concentration and read retention are separate checks",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#4B5563",
    )
    fig.tight_layout()
    save_figure(fig, figure_dir / "16-sepp-workflow-contract")
    plt.close(fig)


def draw_retention(
    figure_dir: Path,
    placement_rows: list[dict[str, object]],
    input_asvs: int,
    input_reads: int,
    retained_asvs: int,
    retained_reads: int,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_plotting()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.9))
    ax = axes[0]
    labels = ["ASVs", "Reads"]
    retained = [100 * retained_asvs / input_asvs, 100 * retained_reads / input_reads]
    removed = [100 - retained[0], 100 - retained[1]]
    y = np.arange(2)
    ax.barh(y, retained, color="#009E73", label="Retained")
    ax.barh(y, removed, left=retained, color="#D55E00", label="Removed")
    for index, value in enumerate(retained):
        ax.text(value - 1.5, index, f"{value:.1f}%", ha="right", va="center", color="white", weight="bold")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percentage of input (%)")
    ax.set_title("Tree/table reconciliation", weight="bold")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", alpha=0.2)

    ax = axes[1]
    candidates = list(range(1, 8))
    asv_percent = []
    read_percent = []
    for count in candidates:
        subset = [row for row in placement_rows if int(row["candidate_placements"]) == count]
        asv_percent.append(100 * len(subset) / input_asvs)
        read_percent.append(100 * sum(int(row["frequency"]) for row in subset) / input_reads)
    x = np.arange(len(candidates))
    width = 0.38
    ax.bar(x - width / 2, asv_percent, width, color="#0072B2", label="ASV-weighted")
    ax.bar(x + width / 2, read_percent, width, color="#E69F00", label="Read-weighted")
    ax.set_xticks(x, [str(value) for value in candidates])
    ax.set_xlabel("Candidate placements retained in jplace")
    ax.set_ylabel("Percentage (%)")
    ax.set_title("Insertion does not imply a unique edge", weight="bold")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, figure_dir / "16-placement-retention-audit")
    plt.close(fig)


def draw_diagnostics(
    figure_dir: Path,
    placement_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_plotting()
    lwr = np.asarray(
        [float(row["top_likelihood_weight"]) for row in placement_rows]
    )
    pendant = np.asarray([float(row["pendant_length"]) for row in placement_rows])
    reads = np.asarray([int(row["frequency"]) for row in placement_rows])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.9))
    bins_lwr = np.linspace(0, 1, 11)
    axes[0].hist(
        lwr,
        bins=bins_lwr,
        weights=np.repeat(100 / len(lwr), len(lwr)),
        alpha=0.68,
        color="#0072B2",
        label="ASV-weighted",
    )
    axes[0].hist(
        lwr,
        bins=bins_lwr,
        weights=100 * reads / reads.sum(),
        histtype="step",
        linewidth=2.0,
        color="#D55E00",
        label="Read-weighted",
    )
    axes[0].axvline(0.5, color="#4B5563", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Top placement likelihood weight")
    axes[0].set_ylabel("Percentage per bin (%)")
    axes[0].set_title("Support is concentrated for only some ASVs", weight="bold")
    axes[0].legend(frameon=False)
    axes[0].grid(axis="y", alpha=0.2)

    bins_pendant = np.linspace(0, max(pendant) * 1.02, 14)
    axes[1].hist(
        pendant,
        bins=bins_pendant,
        weights=np.repeat(100 / len(pendant), len(pendant)),
        alpha=0.68,
        color="#009E73",
        label="ASV-weighted",
    )
    axes[1].hist(
        pendant,
        bins=bins_pendant,
        weights=100 * reads / reads.sum(),
        histtype="step",
        linewidth=2.0,
        color="#CC79A7",
        label="Read-weighted",
    )
    axes[1].set_xlabel("Best-placement pendant length (substitutions/site)")
    axes[1].set_ylabel("Percentage per bin (%)")
    axes[1].set_title("Long pendant branches need separate review", weight="bold")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="y", alpha=0.2)
    fig.tight_layout()
    save_figure(fig, figure_dir / "16-placement-diagnostics")
    plt.close(fig)


def draw_query_tree(
    figure_dir: Path,
    query_tree: dendropy.Tree,
    placement_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_plotting()
    row_map = {str(row["feature_id"]): row for row in placement_rows}
    leaves = list(query_tree.leaf_node_iter())
    y_positions = {id(node): float(index) for index, node in enumerate(leaves)}
    x_positions: dict[int, float] = {}

    def assign_x(node: dendropy.Node, parent_x: float) -> None:
        current = parent_x + float(node.edge_length or 0.0)
        x_positions[id(node)] = current
        for child in node.child_node_iter():
            assign_x(child, current)

    assign_x(query_tree.seed_node, 0.0)
    for node in query_tree.postorder_node_iter():
        if node.is_leaf():
            continue
        children = list(node.child_node_iter())
        y_positions[id(node)] = statistics.mean(
            y_positions[id(child)] for child in children
        )

    fig, ax = plt.subplots(figsize=(10.8, 8.4))
    for node in query_tree.preorder_node_iter():
        children = list(node.child_node_iter())
        if not children:
            continue
        child_y = [y_positions[id(child)] for child in children]
        ax.plot(
            [x_positions[id(node)], x_positions[id(node)]],
            [min(child_y), max(child_y)],
            color="#9CA3AF",
            linewidth=0.35,
            zorder=1,
        )
        for child in children:
            ax.plot(
                [x_positions[id(node)], x_positions[id(child)]],
                [y_positions[id(child)], y_positions[id(child)]],
                color="#9CA3AF",
                linewidth=0.35,
                zorder=1,
            )

    labels = [node.taxon.label for node in leaves if node.taxon is not None]
    colors = np.asarray(
        [float(row_map[label]["top_likelihood_weight"]) for label in labels]
    )
    frequencies = np.asarray([int(row_map[label]["frequency"]) for label in labels])
    sizes = 8 + 25 * np.sqrt(frequencies / frequencies.max())
    scatter = ax.scatter(
        [x_positions[id(node)] for node in leaves],
        [y_positions[id(node)] for node in leaves],
        c=colors,
        s=sizes,
        cmap="viridis",
        vmin=0,
        vmax=1,
        linewidths=0,
        zorder=3,
    )
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.015, fraction=0.03)
    colorbar.set_label("Top placement likelihood weight")
    ax.set_yticks([])
    ax.set_xlabel("Phylogenetic distance from root (substitutions/site)")
    ax.set_title(
        "Induced placement tree for 366 real V4 ASVs",
        weight="bold",
    )
    ax.text(
        0.01,
        0.01,
        "Tip size scales with aggregate read count; reference-only branches are pruned",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#4B5563",
        va="bottom",
    )
    ax.grid(axis="x", alpha=0.15)
    fig.tight_layout()
    save_figure(fig, figure_dir / "16-query-placement-tree")
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
    reference_database = args.reference_database.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    query_dir = project_root / "data/small/taxonomy-databases"
    query_fasta_path = query_dir / "query-v4-asvs.fasta.gz"
    query_abundance_path = query_dir / "query-v4-asv-abundance.tsv"
    source_summary_path = query_dir / "query-source-summary.json"
    for path in (
        reference_database,
        query_fasta_path,
        query_abundance_path,
        source_summary_path,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    command_log = output_dir / "qiime-sepp.log"
    validation_log = output_dir / "validation.log"
    for path in output_dir.glob("*.qza"):
        path.unlink()

    query_records = dict(iter_fasta(query_fasta_path))
    abundance = read_abundance(query_abundance_path)
    query_ids = set(query_records)
    total_reads = sum(row["frequency"] for row in abundance.values())
    reference_summary, reference_ids = audit_reference(reference_database)

    with tempfile.TemporaryDirectory(
        prefix=".article16-work-",
        dir=str(output_dir.parent),
    ) as temporary:
        work_dir = Path(temporary).resolve()
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
                str(figure_dir): "<FIGURE_DIR>",
                str(project_root): "<PROJECT_ROOT>",
                str(reference_database.parent): "<DOWNLOAD_DIR>",
            },
        )

        package_record = runner.run(
            "environment-package-versions",
            ["conda", "list", "-n", args.qiime_env, "--json"],
        )
        package_rows = json.loads(str(package_record["stdout"]))
        versions = {str(row["name"]): str(row["version"]) for row in package_rows}
        expected_versions = {
            "qiime2": "2026.4.0",
            "q2cli": "2026.4.0",
            "q2-fragment-insertion": "2026.4.0",
            "sepp": "4.5.6",
            "hmmer": "3.4",
            "pplacer": "1.1.alpha19",
            "biom-format": "2.1.17",
            "dendropy": "5.0.8",
        }
        environment_rows = [
            {
                "package": package,
                "observed_version": versions.get(package, "missing"),
                "expected_version": expected,
                "status": "PASS" if versions.get(package) == expected else "FAIL",
            }
            for package, expected in expected_versions.items()
        ]

        query_fasta = work_dir / "query-v4-asvs.fasta"
        biom_source_tsv = work_dir / "query-abundance.biom.tsv"
        biom_path = work_dir / "query-abundance.biom"
        query_sequences_qza = output_dir / "query-v4-asvs.qza"
        query_table_qza = output_dir / "query-abundance-table.qza"
        tree_qza = output_dir / "sepp-insertion-tree.qza"
        placements_qza = output_dir / "sepp-placements.qza"
        filtered_qza = output_dir / "sepp-filtered-table.qza"
        removed_qza = output_dir / "sepp-removed-table.qza"

        write_fasta(query_fasta, query_records)
        write_biom_source_tsv(biom_source_tsv, abundance)
        runner.qiime(
            "create-qiime-cache",
            args.qiime_env,
            "tools",
            "cache-create",
            "--cache",
            str(cache_dir),
        )
        runner.qiime(
            "validate-reference-database",
            args.qiime_env,
            "tools",
            "validate",
            str(reference_database),
            "--level",
            "max",
        )
        runner.qiime(
            "import-query-sequences",
            args.qiime_env,
            "tools",
            "import",
            "--type",
            "FeatureData[Sequence]",
            "--input-path",
            str(query_fasta),
            "--output-path",
            str(query_sequences_qza),
            "--validate-level",
            "max",
        )
        runner.run(
            "convert-aggregate-table-to-biom",
            conda_command(
                args.qiime_env,
                "biom",
                "convert",
                "-i",
                str(biom_source_tsv),
                "-o",
                str(biom_path),
                "--table-type",
                "OTU table",
                "--to-hdf5",
            ),
        )
        runner.qiime(
            "import-query-feature-table",
            args.qiime_env,
            "tools",
            "import",
            "--type",
            "FeatureTable[Frequency]",
            "--input-path",
            str(biom_path),
            "--input-format",
            "BIOMV210Format",
            "--output-path",
            str(query_table_qza),
            "--validate-level",
            "max",
        )
        sepp_record = runner.qiime(
            "run-sepp-fragment-insertion",
            args.qiime_env,
            "fragment-insertion",
            "sepp",
            "--i-representative-sequences",
            str(query_sequences_qza),
            "--i-reference-database",
            str(reference_database),
            "--p-alignment-subset-size",
            str(ALIGNMENT_SUBSET_SIZE),
            "--p-placement-subset-size",
            str(PLACEMENT_SUBSET_SIZE),
            "--p-threads",
            str(args.threads),
            "--o-tree",
            str(tree_qza),
            "--o-placements",
            str(placements_qza),
            "--use-cache",
            str(cache_dir),
            "--verbose",
            timeout=28_800,
        )
        filter_record = runner.qiime(
            "filter-feature-table-to-placed-tips",
            args.qiime_env,
            "fragment-insertion",
            "filter-features",
            "--i-table",
            str(query_table_qza),
            "--i-tree",
            str(tree_qza),
            "--o-filtered-table",
            str(filtered_qza),
            "--o-removed-table",
            str(removed_qza),
            "--use-cache",
            str(cache_dir),
            "--verbose",
        )

        artifacts = {
            "query_sequences": query_sequences_qza,
            "query_table": query_table_qza,
            "insertion_tree": tree_qza,
            "placements": placements_qza,
            "filtered_table": filtered_qza,
            "removed_table": removed_qza,
        }
        for name, path in artifacts.items():
            runner.qiime(
                f"maximum-validate-{name.replace('_', '-')}",
                args.qiime_env,
                "tools",
                "validate",
                str(path),
                "--level",
                "max",
            )

        export_tree_dir = work_dir / "export-tree"
        export_placements_dir = work_dir / "export-placements"
        export_filtered_dir = work_dir / "export-filtered"
        export_removed_dir = work_dir / "export-removed"
        export_artifact(
            runner, args.qiime_env, tree_qza, export_tree_dir, "export-tree"
        )
        export_artifact(
            runner,
            args.qiime_env,
            placements_qza,
            export_placements_dir,
            "export-placements",
        )
        export_artifact(
            runner,
            args.qiime_env,
            filtered_qza,
            export_filtered_dir,
            "export-filtered-table",
        )
        export_artifact(
            runner,
            args.qiime_env,
            removed_qza,
            export_removed_dir,
            "export-removed-table",
        )

        exported_tree = next(export_tree_dir.glob("*.nwk"))
        exported_placements = export_placements_dir / "placements.json"
        exported_filtered_biom = export_filtered_dir / "feature-table.biom"
        exported_removed_biom = export_removed_dir / "feature-table.biom"
        tree_path = output_dir / "sepp-insertion-tree.nwk"
        placements_path = output_dir / "sepp-placements.jplace"
        query_tree_path = output_dir / "query-only-placement-tree.nwk"
        shutil.copy2(exported_tree, tree_path)
        shutil.copy2(exported_placements, placements_path)

        filtered_payload = inspect_biom(
            runner,
            args.qiime_env,
            exported_filtered_biom,
            "inspect-filtered-biom",
        )
        removed_payload = inspect_biom(
            runner,
            args.qiime_env,
            exported_removed_biom,
            "inspect-removed-biom",
        )
        write_exported_table(
            output_dir / "sepp-filtered-table.tsv", filtered_payload
        )
        write_exported_table(
            output_dir / "sepp-removed-table.tsv", removed_payload
        )

        placement_rows, placement_file_summary = parse_placements(
            placements_path, abundance
        )
        placement_summary = summarize_placement_rows(placement_rows)
        output_tree_summary, query_tree, query_tip_lengths = parse_output_tree(
            tree_path,
            query_ids,
            reference_ids,
            query_tree_path,
        )
        filtered_ids = {str(item) for item in filtered_payload["observation_ids"]}
        removed_ids = {str(item) for item in removed_payload["observation_ids"]}
        placement_map = {
            str(row["feature_id"]): row for row in placement_rows
        }

        tip_rows: list[dict[str, object]] = []
        for feature_id in sorted(query_records):
            sequence = query_records[feature_id]
            placement = placement_map.get(feature_id, {})
            tip_rows.append(
                {
                    "feature_id": feature_id,
                    "sequence_length": len(sequence),
                    "gc_fraction": round(
                        (sequence.count("G") + sequence.count("C")) / len(sequence),
                        6,
                    ),
                    "frequency": abundance[feature_id]["frequency"],
                    "samples_observed": abundance[feature_id]["samples_observed"],
                    "in_placements": feature_id in placement_map,
                    "in_output_tree": feature_id in query_tip_lengths,
                    "in_filtered_table": feature_id in filtered_ids,
                    "in_removed_table": feature_id in removed_ids,
                    "candidate_placements": placement.get("candidate_placements", ""),
                    "top_likelihood_weight": placement.get("top_likelihood_weight", ""),
                    "pendant_length": placement.get("pendant_length", ""),
                    "distal_length": placement.get("distal_length", ""),
                    "output_tip_branch_length": query_tip_lengths.get(feature_id, ""),
                }
            )

        input_asvs = len(query_records)
        filtered_asvs = len(filtered_ids)
        removed_asvs = len(removed_ids)
        filtered_reads = int(round(float(filtered_payload["total"])))
        removed_reads = int(round(float(removed_payload["total"])))
        retention_rows = [
            {
                "unit": "ASVs",
                "input": input_asvs,
                "retained": filtered_asvs,
                "removed": removed_asvs,
                "retention_fraction": filtered_asvs / input_asvs,
            },
            {
                "unit": "Reads",
                "input": total_reads,
                "retained": filtered_reads,
                "removed": removed_reads,
                "retention_fraction": filtered_reads / total_reads,
            },
        ]

        artifact_rows: list[dict[str, object]] = []
        expected_actions = {
            "query_sequences": "",
            "query_table": "",
            "insertion_tree": "fragment-insertion:sepp",
            "placements": "fragment-insertion:sepp",
            "filtered_table": "fragment-insertion:filter_features",
            "removed_table": "fragment-insertion:filter_features",
        }
        artifact_metadata_map: dict[str, dict[str, str]] = {}
        for name, path in artifacts.items():
            metadata = artifact_metadata(path)
            artifact_metadata_map[name] = metadata
            actions = provenance_actions(path)
            action_text = ";".join(f"{plugin}:{action}" for plugin, action in actions)
            expected_action = expected_actions[name]
            artifact_rows.append(
                {
                    "artifact": name,
                    "file": path.name,
                    "uuid": metadata["uuid"],
                    "semantic_type": metadata["type"],
                    "format": metadata["format"],
                    "maximum_validation": "PASS",
                    "provenance_actions": action_text or "import",
                    "required_action": expected_action or "import",
                    "required_action_present": (
                        not expected_action or expected_action in action_text
                    ),
                }
            )

        draw_workflow(figure_dir)
        draw_retention(
            figure_dir,
            placement_rows,
            input_asvs,
            total_reads,
            filtered_asvs,
            filtered_reads,
        )
        draw_diagnostics(figure_dir, placement_rows)
        draw_query_tree(figure_dir, query_tree, placement_rows)

        checks: list[dict[str, object]] = []

        def add(check_id: str, observed: object, expected: object, passed: bool) -> None:
            checks.append(audit_row(check_id, observed, expected, passed))

        add("reference-sha256", reference_summary["sha256"], REFERENCE_SHA256, reference_summary["sha256"] == REFERENCE_SHA256)
        add("reference-uuid", reference_summary["uuid"], REFERENCE_UUID, reference_summary["uuid"] == REFERENCE_UUID)
        add("reference-semantic-type", reference_summary["semantic_type"], "SeppReferenceDatabase", reference_summary["semantic_type"] == "SeppReferenceDatabase")
        add("reference-format", reference_summary["format"], "SeppReferenceDirFmt", reference_summary["format"] == "SeppReferenceDirFmt")
        add("reference-framework-version", reference_summary["artifact_framework_version"], "2019.10.0", reference_summary["artifact_framework_version"] == "2019.10.0")
        add("reference-plugin-version", reference_summary["artifact_plugin_version"], "2019.10.0", reference_summary["artifact_plugin_version"] == "2019.10.0")
        add("reference-alignment-records", reference_summary["alignment_records"], REFERENCE_TIPS, reference_summary["alignment_records"] == REFERENCE_TIPS)
        add("reference-tree-tips", reference_summary["tree_tips"], REFERENCE_TIPS, reference_summary["tree_tips"] == REFERENCE_TIPS)
        add("reference-alignment-tree-ids", reference_summary["alignment_tree_id_match"], True, bool(reference_summary["alignment_tree_id_match"]))
        add("reference-alignment-min-length", reference_summary["alignment_length_min"], REFERENCE_ALIGNMENT_LENGTH, reference_summary["alignment_length_min"] == REFERENCE_ALIGNMENT_LENGTH)
        add("reference-alignment-max-length", reference_summary["alignment_length_max"], REFERENCE_ALIGNMENT_LENGTH, reference_summary["alignment_length_max"] == REFERENCE_ALIGNMENT_LENGTH)
        add("reference-root-degree", reference_summary["root_degree"], 2, reference_summary["root_degree"] == 2)
        add("reference-negative-branches", reference_summary["negative_branch_lengths"], 0, reference_summary["negative_branch_lengths"] == 0)
        add("query-fasta-sha256", sha256_file(query_fasta_path), QUERY_FASTA_SHA256, sha256_file(query_fasta_path) == QUERY_FASTA_SHA256)
        add("query-abundance-sha256", sha256_file(query_abundance_path), QUERY_ABUNDANCE_SHA256, sha256_file(query_abundance_path) == QUERY_ABUNDANCE_SHA256)
        add("query-asv-count", input_asvs, QUERY_ASVS, input_asvs == QUERY_ASVS)
        add("query-read-count", total_reads, QUERY_READS, total_reads == QUERY_READS)
        add("query-id-reconciliation", len(query_ids & set(abundance)), QUERY_ASVS, query_ids == set(abundance))
        add("query-dna-alphabet", all(set(sequence) <= set("ACGTN") for sequence in query_records.values()), True, all(set(sequence) <= set("ACGTN") for sequence in query_records.values()))
        for package, expected in expected_versions.items():
            add(f"environment-{package}", versions.get(package, "missing"), expected, versions.get(package) == expected)
        add("reference-maximum-validation", "PASS", "PASS", True)
        for name, row in zip(artifacts, artifact_rows):
            add(f"artifact-{name}-maximum-validation", row["maximum_validation"], "PASS", row["maximum_validation"] == "PASS")
        add("tree-semantic-type", artifact_metadata_map["insertion_tree"]["type"], "Phylogeny[Rooted]", artifact_metadata_map["insertion_tree"]["type"] == "Phylogeny[Rooted]")
        add("placements-semantic-type", artifact_metadata_map["placements"]["type"], "Placements", artifact_metadata_map["placements"]["type"] == "Placements")
        add("filtered-table-semantic-type", artifact_metadata_map["filtered_table"]["type"], "FeatureTable[Frequency]", artifact_metadata_map["filtered_table"]["type"] == "FeatureTable[Frequency]")
        add("removed-table-semantic-type", artifact_metadata_map["removed_table"]["type"], "FeatureTable[Frequency]", artifact_metadata_map["removed_table"]["type"] == "FeatureTable[Frequency]")
        for row in artifact_rows:
            add(f"provenance-{row['artifact']}", row["required_action_present"], True, bool(row["required_action_present"]))
        placement_ids = {str(row["feature_id"]) for row in placement_rows}
        add("placement-entry-count", placement_file_summary["placement_entries"], QUERY_ASVS, placement_file_summary["placement_entries"] == QUERY_ASVS)
        add("placement-name-count", placement_file_summary["placement_names"], QUERY_ASVS, placement_file_summary["placement_names"] == QUERY_ASVS)
        add("placement-id-reconciliation", len(placement_ids & query_ids), QUERY_ASVS, placement_ids == query_ids)
        add("placement-weight-sums", max(abs(float(row["likelihood_weight_sum"]) - 1.0) for row in placement_rows), "<=1e-5", max(abs(float(row["likelihood_weight_sum"]) - 1.0) for row in placement_rows) <= 1e-5)
        add("placement-weight-range", all(0 < float(row["top_likelihood_weight"]) <= 1 for row in placement_rows), True, all(0 < float(row["top_likelihood_weight"]) <= 1 for row in placement_rows))
        add("placement-candidate-range", [min(int(row["candidate_placements"]) for row in placement_rows), max(int(row["candidate_placements"]) for row in placement_rows)], [1, 7], min(int(row["candidate_placements"]) for row in placement_rows) >= 1 and max(int(row["candidate_placements"]) for row in placement_rows) <= 7)
        add("placement-pendant-finite", all(math.isfinite(float(row["pendant_length"])) for row in placement_rows), True, all(math.isfinite(float(row["pendant_length"])) for row in placement_rows))
        add("placement-pendant-nonnegative", min(float(row["pendant_length"]) for row in placement_rows), ">=0", min(float(row["pendant_length"]) for row in placement_rows) >= 0)
        add("placement-distal-finite", all(math.isfinite(float(row["distal_length"])) for row in placement_rows), True, all(math.isfinite(float(row["distal_length"])) for row in placement_rows))
        add("placement-distal-nonnegative", min(float(row["distal_length"]) for row in placement_rows), ">=0", min(float(row["distal_length"]) for row in placement_rows) >= 0)
        add("output-tree-tip-count", output_tree_summary["tips"], REFERENCE_TIPS + QUERY_ASVS, output_tree_summary["tips"] == REFERENCE_TIPS + QUERY_ASVS)
        add("output-tree-query-tips", output_tree_summary["query_tips"], QUERY_ASVS, output_tree_summary["query_tips"] == QUERY_ASVS)
        add("output-tree-reference-tips", output_tree_summary["reference_tips_preserved"], REFERENCE_TIPS, output_tree_summary["reference_tips_preserved"] == REFERENCE_TIPS)
        add("output-tree-unexpected-tips", output_tree_summary["unexpected_tips"], 0, output_tree_summary["unexpected_tips"] == 0)
        add("output-tree-root-degree", output_tree_summary["root_degree"], 2, output_tree_summary["root_degree"] == 2)
        add("output-tree-missing-branch-lengths", output_tree_summary["edges_missing_lengths"], 0, output_tree_summary["edges_missing_lengths"] == 0)
        add("output-tree-negative-branch-lengths", output_tree_summary["negative_length_edges"], 0, output_tree_summary["negative_length_edges"] == 0)
        add("query-only-tree-tips", output_tree_summary["query_only_tips"], QUERY_ASVS, output_tree_summary["query_only_tips"] == QUERY_ASVS)
        add("query-only-tree-nodes", output_tree_summary["query_only_nodes"], 2 * QUERY_ASVS - 1, output_tree_summary["query_only_nodes"] == 2 * QUERY_ASVS - 1)
        add("input-feature-table-asvs", len(abundance), QUERY_ASVS, len(abundance) == QUERY_ASVS)
        add("input-feature-table-reads", total_reads, QUERY_READS, total_reads == QUERY_READS)
        add("filtered-feature-table-asvs", filtered_asvs, QUERY_ASVS, filtered_asvs == QUERY_ASVS)
        add("filtered-feature-table-reads", filtered_reads, QUERY_READS, filtered_reads == QUERY_READS)
        add("removed-feature-table-asvs", removed_asvs, 0, removed_asvs == 0)
        add("removed-feature-table-reads", removed_reads, 0, removed_reads == 0)
        add("feature-table-id-partition", len(filtered_ids | removed_ids), QUERY_ASVS, filtered_ids.isdisjoint(removed_ids) and filtered_ids | removed_ids == query_ids)
        add("feature-table-read-partition", filtered_reads + removed_reads, QUERY_READS, filtered_reads + removed_reads == QUERY_READS)
        add("filtered-features-are-tree-tips", len(filtered_ids & set(query_tip_lengths)), filtered_asvs, filtered_ids <= set(query_tip_lengths))
        add("all-runtime-paths-absolute", all(Path(part).is_absolute() for part in [str(query_sequences_qza), str(reference_database), str(cache_dir), str(tree_qza), str(placements_qza)]), True, all(Path(part).is_absolute() for part in [str(query_sequences_qza), str(reference_database), str(cache_dir), str(tree_qza), str(placements_qza)]))

        figure_bases = [
            "16-sepp-workflow-contract",
            "16-placement-retention-audit",
            "16-placement-diagnostics",
            "16-query-placement-tree",
        ]
        figure_files = [
            figure_dir / f"{base}.{suffix}"
            for base in figure_bases
            for suffix in ("pdf", "png", "tiff")
        ]
        add("figure-file-count", sum(path.exists() for path in figure_files), 12, all(path.exists() and path.stat().st_size > 0 for path in figure_files))

        support_bin_rows = [
            {
                "top_likelihood_weight_bin": label,
                "asvs": values["asvs"],
                "reads": values["reads"],
                "asv_fraction": values["asvs"] / input_asvs,
                "read_fraction": values["reads"] / total_reads,
            }
            for label, values in placement_summary["support_bins"].items()
        ]
        candidate_rows = [
            {
                "candidate_placements": int(label),
                "asvs": values["asvs"],
                "reads": values["reads"],
                "asv_fraction": values["asvs"] / input_asvs,
                "read_fraction": values["reads"] / total_reads,
            }
            for label, values in placement_summary["candidate_placement_counts"].items()
        ]
        reference_rows = [
            {"metric": key, "value": value}
            for key, value in reference_summary.items()
        ]
        execution_rows = [
            {
                "label": record["label"],
                "returncode": record["returncode"],
                "elapsed_seconds": record["elapsed_seconds"],
            }
            for record in runner.records
        ]

        write_tsv(
            output_dir / "environment-audit.tsv",
            environment_rows,
            ["package", "observed_version", "expected_version", "status"],
        )
        write_tsv(
            output_dir / "reference-database-audit.tsv",
            reference_rows,
            ["metric", "value"],
        )
        write_tsv(
            output_dir / "placement-audit.tsv",
            placement_rows,
            [
                "feature_id",
                "frequency",
                "candidate_placements",
                "top_edge_num",
                "top_likelihood",
                "top_likelihood_weight",
                "second_likelihood_weight",
                "likelihood_weight_sum",
                "distal_length",
                "pendant_length",
            ],
        )
        write_tsv(
            output_dir / "tip-reconciliation.tsv",
            tip_rows,
            list(tip_rows[0]),
        )
        write_tsv(
            output_dir / "feature-retention-audit.tsv",
            retention_rows,
            ["unit", "input", "retained", "removed", "retention_fraction"],
        )
        write_tsv(
            output_dir / "placement-support-bins.tsv",
            support_bin_rows,
            ["top_likelihood_weight_bin", "asvs", "reads", "asv_fraction", "read_fraction"],
        )
        write_tsv(
            output_dir / "candidate-placement-counts.tsv",
            candidate_rows,
            ["candidate_placements", "asvs", "reads", "asv_fraction", "read_fraction"],
        )
        write_tsv(
            output_dir / "artifact-provenance-audit.tsv",
            artifact_rows,
            [
                "artifact",
                "file",
                "uuid",
                "semantic_type",
                "format",
                "maximum_validation",
                "provenance_actions",
                "required_action",
                "required_action_present",
            ],
        )
        write_tsv(
            output_dir / "command-execution-audit.tsv",
            execution_rows,
            ["label", "returncode", "elapsed_seconds"],
        )
        write_tsv(
            output_dir / "validation-checks.tsv",
            checks,
            ["check_id", "observed", "expected", "status"],
        )

        checks_passed = sum(row["status"] == "PASS" for row in checks)
        summary: dict[str, object] = {
            "status": "passed" if checks_passed == len(checks) else "failed",
            "checks_total": len(checks),
            "checks_passed": checks_passed,
            "query_asvs": input_asvs,
            "query_reads": total_reads,
            "query_length_min": min(map(len, query_records.values())),
            "query_length_median": statistics.median(map(len, query_records.values())),
            "query_length_max": max(map(len, query_records.values())),
            "reference": reference_summary,
            "runtime": {
                "qiime_environment": args.qiime_env,
                "threads": args.threads,
                "alignment_subset_size": ALIGNMENT_SUBSET_SIZE,
                "placement_subset_size": PLACEMENT_SUBSET_SIZE,
                "sepp_elapsed_seconds": sepp_record["elapsed_seconds"],
                "filter_elapsed_seconds": filter_record["elapsed_seconds"],
                "versions": {key: versions.get(key, "missing") for key in expected_versions},
            },
            "placements": {
                **placement_file_summary,
                **placement_summary,
            },
            "tree": output_tree_summary,
            "retention": {
                "input_asvs": input_asvs,
                "retained_asvs": filtered_asvs,
                "removed_asvs": removed_asvs,
                "input_reads": total_reads,
                "retained_reads": filtered_reads,
                "removed_reads": removed_reads,
                "asv_retention_fraction": filtered_asvs / input_asvs,
                "read_retention_fraction": filtered_reads / total_reads,
            },
            "artifact_metadata": artifact_metadata_map,
            "biological_export_sha256": {
                "insertion_tree_newick": sha256_file(tree_path),
                "placements_jplace": sha256_file(placements_path),
                "query_only_tree_newick": sha256_file(query_tree_path),
            },
            "interpretation_boundary": (
                "All query ASVs were inserted and retained, but insertion is not "
                "equivalent to a unique or taxonomically correct placement. "
                "Likelihood weights are relative support within the retained "
                "candidate set, not calibrated correctness probabilities."
            ),
        }
        write_json(output_dir / "sepp-phylogeny-summary.json", summary)
        validation_log.write_text(
            "\n".join(
                [
                    "Article 16 SEPP validation",
                    f"status: {summary['status']}",
                    f"checks: {checks_passed}/{len(checks)}",
                    f"reference tips: {REFERENCE_TIPS}",
                    f"query ASVs: {input_asvs}",
                    f"query reads: {total_reads}",
                    f"placed ASVs: {placement_file_summary['placement_names']}",
                    f"retained ASVs: {filtered_asvs}",
                    f"retained reads: {filtered_reads}",
                    f"removed ASVs: {removed_asvs}",
                    f"removed reads: {removed_reads}",
                    f"SEPP elapsed seconds: {sepp_record['elapsed_seconds']}",
                    "boundary: placement success is not calibrated placement correctness",
                    "",
                ]
            ),
            encoding="utf-8",
        )

    if checks_passed != len(checks):
        failed = [row["check_id"] for row in checks if row["status"] != "PASS"]
        raise SystemExit(f"SEPP validation failed: {', '.join(failed)}")
    print(
        f"SEPP validation passed: {checks_passed}/{len(checks)} checks, "
        f"{QUERY_ASVS}/{QUERY_ASVS} ASVs retained"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
