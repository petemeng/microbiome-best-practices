#!/usr/bin/env python3
"""Classify one real V4 ASV set against three region-matched databases.

Article 14 deliberately holds the query, V4 extraction contract, and
consensus-VSEARCH parameters constant. The validator checks all bundled
reference hashes, imports fresh QIIME 2 Artifacts, runs three classifications,
audits rank coverage and off-target labels, and exports publication figures.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import statistics
import subprocess
import tempfile
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


QIIME_ENV = "microbiome-qiime2-2026.4"
QIIME_VERSION = "2026.4.0"
Q2_FEATURE_CLASSIFIER_VERSION = "2026.4.0"
RESCRIPT_VERSION = "2026.4.0"
VSEARCH_VERSION = "2.30.4"
SKLEARN_VERSION = "1.7.1"
QUERY_RECORDS = 366
QUERY_READS = 5213
QUERY_FASTA_SHA256 = (
    "fa71915be8007d36ec70f85d28401f5a94d0d33eca323782cc937cfc4bcf1af4"
)
QUERY_ABUNDANCE_SHA256 = (
    "a80e3b94d8632c48fb6c01e788e842c022c3b1ee64b79dc94456618328895e6e"
)
EXPECTED_REFERENCE_ASSETS = {
    "SILVA": {
        "raw_v4_records": 450770,
        "lca_v4_records": 262752,
        "sequence_sha256": (
            "3a251d263bd4e12c96023c84b3f2e255b8bcb0d5865dc6f299f9918e6be3a808"
        ),
        "taxonomy_sha256": (
            "f49320aaa32aa70af5c5a48649786e3bb4177868da5da4d2143d0a54c5fa47b0"
        ),
    },
    "GTDB": {
        "raw_v4_records": 79669,
        "lca_v4_records": 51985,
        "sequence_sha256": (
            "3c4234b4defcdc3d8c3b4d1e8956b12288bccd853a4b7f3a3003847c32e6e250"
        ),
        "taxonomy_sha256": (
            "a23d0d4f62183e949c868c3f25aace3acdebb9362ae0bb1d0e9fe81e74c43fe0"
        ),
    },
    "Greengenes2": {
        "raw_v4_records": 335976,
        "lca_v4_records": 98185,
        "sequence_sha256": (
            "330fe7f400c84c1e454ca9193bcc49f078b01c7880a8bfb99cbe4c0c3c1c14d2"
        ),
        "taxonomy_sha256": (
            "9594e9ae78a568672891693207fd4dcc31f8eada83b931031642fafcff2684bd"
        ),
    },
}
FORWARD_PRIMER = "GTGYCAGCMGCCGCGGTAA"
REVERSE_PRIMER = "GGACTACNVGGGTWTCTAAT"
EXTRACTION_IDENTITY = 0.80
EXTRACTION_MIN_LENGTH = 200
EXTRACTION_MAX_LENGTH = 400

CLASSIFIER_IDENTITY = 0.80
CLASSIFIER_QUERY_COVERAGE = 0.80
CLASSIFIER_MIN_CONSENSUS = 0.51
CLASSIFIER_STRAND = "both"
CLASSIFIER_THREADS = 4
CLASSIFIER_MAXACCEPTS = 50

RANKS = (
    ("domain", "d__"),
    ("phylum", "p__"),
    ("class", "c__"),
    ("order", "o__"),
    ("family", "f__"),
    ("genus", "g__"),
    ("species", "s__"),
)
DATABASES = (
    ("silva", "SILVA", "silva-138.2", "138.2"),
    ("gtdb", "GTDB", "gtdb-r232", "R11-RS232"),
    (
        "greengenes2",
        "Greengenes2",
        "greengenes2-2024.09",
        "2024.09",
    ),
)
COLORS = {
    "SILVA": "#0072B2",
    "GTDB": "#D55E00",
    "Greengenes2": "#009E73",
}


@dataclass(frozen=True)
class DatabaseInput:
    key: str
    label: str
    asset_prefix: str
    release: str
    summary_path: Path
    fasta_path: Path
    taxonomy_path: Path
    summary: dict[str, object]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--qiime-env", default=QIIME_ENV)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


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


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    identifier: str | None = None
    sequence_parts: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence_parts).upper()
                identifier = line[1:].split(maxsplit=1)[0]
                sequence_parts = []
            else:
                if identifier is None:
                    raise ValueError(f"sequence before FASTA header: {path}")
                sequence_parts.append(line)
    if identifier is not None:
        yield identifier, "".join(sequence_parts).upper()


def read_taxonomy(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    taxonomy: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"taxonomy has no header: {path}")
        id_field = next(
            (
                name
                for name in reader.fieldnames
                if name.strip().lower() in {"feature id", "feature-id"}
            ),
            None,
        )
        taxon_field = next(
            (
                name
                for name in reader.fieldnames
                if name.strip().lower() == "taxon"
            ),
            None,
        )
        if id_field is None or taxon_field is None:
            raise ValueError(f"unexpected taxonomy columns: {reader.fieldnames}")
        for row in reader:
            feature_id = str(row[id_field]).strip()
            if feature_id in taxonomy:
                raise ValueError(f"duplicate taxonomy ID: {feature_id}")
            taxonomy[feature_id] = str(row[taxon_field]).strip()
    return taxonomy


def read_query_abundance(path: Path) -> dict[str, dict[str, int]]:
    records: dict[str, dict[str, int]] = {}
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
            records[row[0]] = {
                "frequency": int(float(row[1])),
                "samples_observed": int(float(row[2])),
            }
    return records


def rank_map(taxon: str) -> dict[str, str]:
    if not taxon or taxon.strip().lower() == "unassigned":
        return {}
    parsed: dict[str, str] = {}
    for item in taxon.split(";"):
        label = item.strip()
        for rank, prefix in RANKS:
            if label.startswith(prefix):
                value = label[len(prefix) :].strip()
                if value and value.lower() not in {"unassigned", "uncultured"}:
                    parsed[rank] = value
                break
    return parsed


def classification_category(taxon: str) -> str:
    lower = taxon.lower()
    ranks = rank_map(taxon)
    if not ranks:
        return "Unassigned"
    if "chloroplast" in lower or "plastid" in lower:
        return "Chloroplast/plastid"
    if "mitochond" in lower:
        return "Mitochondria"
    domain = ranks.get("domain", "")
    if domain == "Bacteria":
        return "Bacteria"
    if domain == "Archaea":
        return "Archaea"
    if domain == "Eukaryota":
        return "Eukaryota"
    return "Other/off-target"


def unpack_gzip(source: Path, destination: Path) -> None:
    with gzip.open(source, "rb") as input_handle, destination.open("wb") as output:
        shutil.copyfileobj(input_handle, output, length=1024 * 1024)


def run_command(
    label: str,
    command: Sequence[str],
    environment: dict[str, str],
    timeout: int = 7200,
) -> dict[str, object]:
    started = time.monotonic()
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=timeout,
    )
    return {
        "label": label,
        "command": list(command),
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }


def require_success(result: dict[str, object]) -> None:
    if int(result["returncode"]) != 0:
        raise RuntimeError(
            f"{result['label']} failed\nSTDOUT:\n{result['stdout']}\n"
            f"STDERR:\n{result['stderr']}"
        )


def conda_command(environment: str, executable: str, *args: str) -> tuple[str, ...]:
    return ("conda", "run", "-n", environment, executable, *args)


def qiime_command(environment: str, *args: str) -> tuple[str, ...]:
    return conda_command(environment, "qiime", *args)


def isolated_environment(root: Path) -> dict[str, str]:
    environment = os.environ.copy()
    for name in ("xdg", "numba", "matplotlib"):
        (root / "cache" / name).mkdir(parents=True, exist_ok=True)
    environment.update(
        {
            "XDG_CACHE_HOME": str(root / "cache" / "xdg"),
            "NUMBA_CACHE_DIR": str(root / "cache" / "numba"),
            "MPLCONFIGDIR": str(root / "cache" / "matplotlib"),
            "LC_ALL": "C.UTF-8",
            "LANG": "C.UTF-8",
        }
    )
    return environment


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
                actions.add((plugin_value, action.group(1).strip().strip("'\"")))
    return sorted(actions)


def locate_exported_file(directory: Path, suffixes: Sequence[str]) -> Path:
    for suffix in suffixes:
        matches = sorted(
            path
            for path in directory.rglob("*")
            if path.is_file() and path.name == suffix
        )
        if len(matches) == 1:
            return matches[0]
    files = [path.name for path in directory.rglob("*") if path.is_file()]
    raise FileNotFoundError(f"expected {suffixes} in {directory}; found {files}")


def parse_classification(path: Path) -> dict[str, dict[str, object]]:
    classifications: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"classification has no header: {path}")
        for row in reader:
            feature_id = str(row.get("Feature ID", row.get("FeatureID", "")))
            taxon = str(row.get("Taxon", ""))
            confidence_raw = str(
                row.get("Confidence", row.get("Consensus", ""))
            ).strip()
            confidence = float(confidence_raw) if confidence_raw else 0.0
            classifications[feature_id] = {
                "taxon": taxon,
                "confidence": confidence,
                "ranks": rank_map(taxon),
                "category": classification_category(taxon),
            }
    return classifications


def parse_blast6(path: Path) -> dict[str, dict[str, float | int]]:
    hit_counts: Counter[str] = Counter()
    best_identity: defaultdict[str, float] = defaultdict(float)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0].lower() in {"qseqid", "query id", "query-id"}:
                continue
            if len(fields) < 3:
                continue
            query_id = fields[0]
            if fields[1] == "*":
                continue
            try:
                identity = float(fields[2])
            except ValueError:
                continue
            if identity <= 0:
                continue
            hit_counts[query_id] += 1
            best_identity[query_id] = max(best_identity[query_id], identity)
    return {
        query_id: {
            "top_hits": hit_counts[query_id],
            "best_percent_identity": round(best_identity[query_id], 6),
        }
        for query_id in hit_counts
    }


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


def draw_decision_map(figure_dir: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(12.2, 6.6), facecolor="white")
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 7.2)
    ax.axis("off")
    dark = "#111827"
    gray = "#4B5563"
    ax.text(
        0.3,
        6.82,
        "Choose a reference database from the biological question",
        fontsize=17,
        fontweight="bold",
        color=dark,
    )
    ax.text(
        0.3,
        6.38,
        "Then lock release, marker region, and classifier before comparing labels.",
        fontsize=10.2,
        color=gray,
    )

    def box(
        x: float,
        y: float,
        width: float,
        height: float,
        title: str,
        body: str,
        color: str,
    ) -> None:
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.03,rounding_size=0.12",
            facecolor="white",
            edgecolor=color,
            linewidth=2,
        )
        ax.add_patch(patch)
        ax.text(
            x + 0.22,
            y + height - 0.35,
            title,
            fontsize=11.5,
            fontweight="bold",
            color=color,
            va="top",
        )
        ax.text(
            x + 0.22,
            y + height - 0.78,
            body,
            fontsize=9.2,
            color=dark,
            va="top",
            linespacing=1.35,
        )

    box(
        0.45,
        3.55,
        3.35,
        2.05,
        "Broad SSU coverage",
        "Environmental surveys\nBacteria + Archaea + Eukaryota\nConservative rank interpretation",
        COLORS["SILVA"],
    )
    box(
        4.42,
        3.55,
        3.35,
        2.05,
        "Genome-linked taxonomy",
        "Bacteria + Archaea\nGenome phylogeny nomenclature\nAudit missing or discordant SSU",
        COLORS["GTDB"],
    )
    box(
        8.39,
        3.55,
        3.35,
        2.05,
        "Tree-native integration",
        "Amplicon + shotgun studies\nUnified phylogenetic backbone\nDo not confuse with Greengenes 13_8",
        COLORS["Greengenes2"],
    )

    for x, color, label in (
        (2.13, COLORS["SILVA"], "SILVA 138.2"),
        (6.10, COLORS["GTDB"], "GTDB R232"),
        (10.07, COLORS["Greengenes2"], "Greengenes2 2024.09"),
    ):
        ax.add_patch(
            FancyArrowPatch(
                (x, 3.45),
                (x, 2.63),
                arrowstyle="-|>",
                mutation_scale=15,
                linewidth=1.7,
                color=color,
            )
        )
        ax.text(
            x,
            2.37,
            label,
            ha="center",
            va="top",
            fontsize=9.6,
            fontweight="bold",
            color=color,
        )

    shared = FancyBboxPatch(
        (1.25, 0.48),
        9.7,
        1.05,
        boxstyle="round,pad=0.03,rounding_size=0.14",
        facecolor="#F3F4F6",
        edgecolor="#9CA3AF",
        linewidth=1.3,
    )
    ax.add_patch(shared)
    ax.text(
        6.1,
        1.16,
        "Fair comparison: same 366 V4 ASVs · same 515F/806R extraction · "
        "same exact-sequence LCA dereplication",
        ha="center",
        va="center",
        fontsize=10.1,
        fontweight="bold",
        color=dark,
    )
    ax.text(
        6.1,
        0.78,
        "Same consensus-VSEARCH identity, query coverage, strand, top-hit, and consensus thresholds",
        ha="center",
        va="center",
        fontsize=9.2,
        color=gray,
    )
    save_figure(fig, figure_dir / "14-database-decision-map")
    plt.close(fig)


def draw_reference_scope(
    figure_dir: Path,
    source_rows: Sequence[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [str(row["database"]) for row in source_rows]
    stages = ("Full reference", "Primer-matched V4", "Unique V4 (LCA)")
    values = {
        "Full reference": [int(row["full_sequence_records"]) for row in source_rows],
        "Primer-matched V4": [int(row["raw_v4_records"]) for row in source_rows],
        "Unique V4 (LCA)": [int(row["lca_v4_records"]) for row in source_rows],
    }
    stage_colors = ("#D1D5DB", "#80B1D3", "#1F78B4")
    x = np.arange(len(labels))
    width = 0.23
    fig, ax = plt.subplots(figsize=(9.2, 5.7), facecolor="white")
    for index, stage in enumerate(stages):
        bars = ax.bar(
            x + (index - 1) * width,
            values[stage],
            width,
            label=stage,
            color=stage_colors[index],
            edgecolor="white",
            linewidth=0.8,
        )
        for bar, value in zip(bars, values[stage]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value * 1.08,
                f"{value:,}",
                ha="center",
                va="bottom",
                fontsize=7.5,
                rotation=90 if value >= 100000 else 0,
            )
    ax.set_yscale("log")
    ax.set_ylim(1e3, max(values["Full reference"]) * 3.0)
    ax.set_ylabel("Reference sequences (log scale)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(
        "Reference scope changes after region matching and exact-sequence dereplication",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.18, which="both")
    ax.legend(frameon=False, ncol=3, loc="upper center")
    fig.tight_layout()
    save_figure(fig, figure_dir / "14-reference-scope-audit")
    plt.close(fig)


def draw_rank_coverage(
    figure_dir: Path,
    rank_rows: Sequence[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.2, 5.7), facecolor="white")
    rank_names = [rank.title() for rank, _prefix in RANKS]
    label_offsets = {
        "SILVA": 3.2,
        "GTDB": -2.0,
        "Greengenes2": 0.6,
    }
    label_x_offsets = {
        "SILVA": -0.045,
        "GTDB": 0.0,
        "Greengenes2": 0.045,
    }
    for _key, label, _prefix, _release in DATABASES:
        rows = [row for row in rank_rows if row["database"] == label]
        values = [
            100.0 * float(next(row["read_fraction"] for row in rows if row["rank"] == rank))
            for rank, _rank_prefix in RANKS
        ]
        ax.plot(
            rank_names,
            values,
            marker="o",
            markersize=5.2,
            linewidth=2.2,
            color=COLORS[label],
            label=label,
        )
        for x_index, value in enumerate(values):
            offset = label_offsets[label]
            ax.text(
                x_index + label_x_offsets[label],
                value + offset,
                f"{value:.1f}",
                color=COLORS[label],
                ha="center",
                va="bottom" if offset > 0 else "top",
                fontsize=7.2,
            )
    ax.set_ylim(0, 107)
    ax.set_ylabel("Read-weighted ASVs with a label (%)")
    ax.set_xlabel("Taxonomic rank")
    ax.set_title(
        "The same V4 reads receive different rank depth by database",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    save_figure(fig, figure_dir / "14-rank-assignment-coverage")
    plt.close(fig)


def draw_concordance(
    figure_dir: Path,
    concordance: dict[str, dict[str, float]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    pairs = ("SILVA vs GTDB", "SILVA vs Greengenes2", "GTDB vs Greengenes2")
    ranks = ("phylum", "class", "order", "family", "genus")
    matrix = np.array(
        [[100.0 * concordance[pair][rank] for rank in ranks] for pair in pairs]
    )
    fig, ax = plt.subplots(figsize=(8.8, 4.9), facecolor="white")
    image = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            color = "white" if value < 55 else "#111827"
            ax.text(
                column_index,
                row_index,
                f"{value:.1f}%",
                ha="center",
                va="center",
                color=color,
                fontsize=9,
                fontweight="bold",
            )
    ax.set_xticks(range(len(ranks)))
    ax.set_xticklabels([rank.title() for rank in ranks])
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs)
    ax.set_title(
        "Exact-name concordance among jointly assigned ASVs",
        loc="left",
        fontweight="bold",
    )
    colorbar = fig.colorbar(image, ax=ax, pad=0.025)
    colorbar.set_label("Exact label agreement (%)")
    fig.tight_layout()
    save_figure(fig, figure_dir / "14-taxonomy-concordance")
    plt.close(fig)


def sanitize_text(text: str, replacements: dict[str, str]) -> str:
    result = text
    for source, replacement in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        result = result.replace(source, replacement)
    return result


def load_database_inputs(reference_dir: Path) -> list[DatabaseInput]:
    databases: list[DatabaseInput] = []
    for key, label, asset_prefix, release in DATABASES:
        summary_path = reference_dir / f"{asset_prefix}-source-summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        databases.append(
            DatabaseInput(
                key=key,
                label=label,
                asset_prefix=asset_prefix,
                release=release,
                summary_path=summary_path,
                fasta_path=reference_dir / f"{asset_prefix}-v4.fasta.gz",
                taxonomy_path=reference_dir
                / f"{asset_prefix}-v4-taxonomy.tsv.gz",
                summary=summary,
            )
        )
    return databases


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    reference_dir = project_root / "data/small/taxonomy-databases"
    query_fasta = reference_dir / "query-v4-asvs.fasta.gz"
    query_abundance_path = reference_dir / "query-v4-asv-abundance.tsv"
    query_summary_path = reference_dir / "query-source-summary.json"
    query_summary = json.loads(query_summary_path.read_text(encoding="utf-8"))
    databases = load_database_inputs(reference_dir)
    checks: list[dict[str, object]] = []
    commands: list[dict[str, object]] = []

    query_sequences = dict(iter_fasta(query_fasta))
    query_abundance = read_query_abundance(query_abundance_path)
    if set(query_sequences) != set(query_abundance):
        raise ValueError("query FASTA and abundance feature IDs differ")

    checks.extend(
        [
            audit_row(
                "query-fasta-sha256",
                sha256_file(query_fasta),
                QUERY_FASTA_SHA256,
                sha256_file(query_fasta) == QUERY_FASTA_SHA256,
            ),
            audit_row(
                "query-abundance-sha256",
                sha256_file(query_abundance_path),
                QUERY_ABUNDANCE_SHA256,
                sha256_file(query_abundance_path) == QUERY_ABUNDANCE_SHA256,
            ),
            audit_row(
                "query-records",
                len(query_sequences),
                QUERY_RECORDS,
                len(query_sequences) == QUERY_RECORDS,
            ),
            audit_row(
                "query-total-reads",
                sum(item["frequency"] for item in query_abundance.values()),
                QUERY_READS,
                sum(item["frequency"] for item in query_abundance.values())
                == QUERY_READS,
            ),
            audit_row(
                "query-source-artifact-lineage",
                query_summary["processing_contract"]["source_artifact"],
                "results/12-dada2-asv/dada2-rep-seqs.qza",
                query_summary["processing_contract"]["source_artifact"]
                == "results/12-dada2-asv/dada2-rep-seqs.qza",
            ),
        ]
    )

    source_rows: list[dict[str, object]] = []
    rank_audit_rows: list[dict[str, object]] = []
    v4_rows: list[dict[str, object]] = []
    reference_taxonomies: dict[str, dict[str, str]] = {}
    expected_full_counts = {
        "SILVA": 510495,
        "GTDB": 93770,
        "Greengenes2": 337506,
    }
    expected_licenses = {
        "SILVA": "CC BY 4.0",
        "GTDB": "CC BY-SA 4.0",
        "Greengenes2": "BSD-3-Clause",
    }
    for database in databases:
        summary = database.summary
        observed_fasta_sha = sha256_file(database.fasta_path)
        observed_taxonomy_sha = sha256_file(database.taxonomy_path)
        expected_fasta_sha = str(summary["derived_files"]["sequences_sha256"])
        expected_taxonomy_sha = str(summary["derived_files"]["taxonomy_sha256"])
        reference_sequences = dict(iter_fasta(database.fasta_path))
        reference_taxonomy = read_taxonomy(database.taxonomy_path)
        reference_taxonomies[database.key] = reference_taxonomy
        lengths = [len(sequence) for sequence in reference_sequences.values()]
        if set(reference_sequences) != set(reference_taxonomy):
            raise ValueError(
                f"{database.label} V4 sequence/taxonomy identifiers differ"
            )
        if database.label == "GTDB":
            full_records = int(
                summary["full_reference_audit"]["full_sequence_lengths"]["records"]
            )
            full_taxonomy_records = int(
                summary["full_reference_audit"]["full_taxonomy_rows"]
            )
        else:
            full_records = int(summary["full_sequences"]["length"]["records"])
            full_taxonomy_records = int(summary["full_taxonomy"]["records"])
        source_rows.append(
            {
                "database": database.label,
                "release": database.release,
                "reference_subset": summary["reference_subset"],
                "license": summary["license"],
                "full_sequence_records": full_records,
                "full_taxonomy_records": full_taxonomy_records,
                "raw_v4_records": summary["raw_v4_records"],
                "lca_v4_records": len(reference_sequences),
                "official_url": summary["official_url"],
                "sequence_sha256": observed_fasta_sha,
                "taxonomy_sha256": observed_taxonomy_sha,
            }
        )
        for rank, prefix in RANKS:
            full_taxonomy_summary = (
                summary["full_reference_audit"]["representative_taxonomy"]
                if database.label == "GTDB"
                else summary["full_taxonomy"]
            )
            full_nonempty = int(
                full_taxonomy_summary["nonempty_rank_counts"].get(prefix, 0)
            )
            region_nonempty = sum(
                1 for taxon in reference_taxonomy.values() if rank in rank_map(taxon)
            )
            rank_audit_rows.append(
                {
                    "database": database.label,
                    "rank": rank,
                    "full_nonempty_records": full_nonempty,
                    "full_fraction": round(full_nonempty / full_records, 6),
                    "lca_v4_nonempty_records": region_nonempty,
                    "lca_v4_fraction": round(
                        region_nonempty / len(reference_taxonomy), 6
                    ),
                }
            )
        v4_rows.append(
            {
                "database": database.label,
                "records": len(reference_sequences),
                "minimum_length": min(lengths),
                "median_length": round(statistics.median(lengths), 2),
                "mean_length": round(statistics.fmean(lengths), 2),
                "maximum_length": max(lengths),
                "sequence_taxonomy_id_match": True,
                "sequence_sha256": observed_fasta_sha,
                "taxonomy_sha256": observed_taxonomy_sha,
            }
        )
        checks.extend(
            [
                audit_row(
                    f"{database.key}-release",
                    summary["release"],
                    database.release,
                    str(summary["release"]) == database.release,
                ),
                audit_row(
                    f"{database.key}-license",
                    summary["license"],
                    expected_licenses[database.label],
                    str(summary["license"]) == expected_licenses[database.label],
                ),
                audit_row(
                    f"{database.key}-full-records",
                    full_records,
                    expected_full_counts[database.label],
                    full_records == expected_full_counts[database.label],
                ),
                audit_row(
                    f"{database.key}-v4-sequence-sha256",
                    observed_fasta_sha,
                    EXPECTED_REFERENCE_ASSETS[database.label]["sequence_sha256"],
                    observed_fasta_sha
                    == EXPECTED_REFERENCE_ASSETS[database.label][
                        "sequence_sha256"
                    ]
                    == expected_fasta_sha,
                ),
                audit_row(
                    f"{database.key}-v4-taxonomy-sha256",
                    observed_taxonomy_sha,
                    EXPECTED_REFERENCE_ASSETS[database.label]["taxonomy_sha256"],
                    observed_taxonomy_sha
                    == EXPECTED_REFERENCE_ASSETS[database.label][
                        "taxonomy_sha256"
                    ]
                    == expected_taxonomy_sha,
                ),
                audit_row(
                    f"{database.key}-raw-v4-records",
                    summary["raw_v4_records"],
                    EXPECTED_REFERENCE_ASSETS[database.label]["raw_v4_records"],
                    int(summary["raw_v4_records"])
                    == EXPECTED_REFERENCE_ASSETS[database.label][
                        "raw_v4_records"
                    ],
                ),
                audit_row(
                    f"{database.key}-lca-v4-records",
                    len(reference_sequences),
                    EXPECTED_REFERENCE_ASSETS[database.label]["lca_v4_records"],
                    len(reference_sequences)
                    == EXPECTED_REFERENCE_ASSETS[database.label][
                        "lca_v4_records"
                    ],
                ),
                audit_row(
                    f"{database.key}-v4-id-contract",
                    len(set(reference_sequences) ^ set(reference_taxonomy)),
                    0,
                    set(reference_sequences) == set(reference_taxonomy),
                ),
                audit_row(
                    f"{database.key}-v4-length-contract",
                    [min(lengths), max(lengths)],
                    f"{EXTRACTION_MIN_LENGTH}..{EXTRACTION_MAX_LENGTH}",
                    min(lengths) >= EXTRACTION_MIN_LENGTH
                    and max(lengths) <= EXTRACTION_MAX_LENGTH,
                ),
                audit_row(
                    f"{database.key}-primer-contract",
                    [
                        summary["region_processing"]["forward_primer_5_to_3"],
                        summary["region_processing"]["reverse_primer_5_to_3"],
                    ],
                    [FORWARD_PRIMER, REVERSE_PRIMER],
                    summary["region_processing"]["forward_primer_5_to_3"]
                    == FORWARD_PRIMER
                    and summary["region_processing"]["reverse_primer_5_to_3"]
                    == REVERSE_PRIMER,
                ),
                audit_row(
                    f"{database.key}-lca-dereplication",
                    summary["exact_sequence_dereplication"]["mode"],
                    "least common ancestor",
                    summary["exact_sequence_dereplication"]["mode"]
                    == "least common ancestor",
                ),
            ]
        )

    checks.extend(
        [
            audit_row(
                "silva-species-not-promoted",
                databases[0].summary["taxonomy_contract"],
                "species labels intentionally excluded",
                "species labels intentionally excluded"
                in str(databases[0].summary["taxonomy_contract"]),
            ),
            audit_row(
                "gtdb-r232-release-note-method",
                databases[1].summary["method_documentation_boundary"][
                    "release_notes"
                ],
                "barnap + 396 bp",
                "barnap"
                in str(
                    databases[1].summary["method_documentation_boundary"][
                        "release_notes"
                    ]
                )
                and "396"
                in str(
                    databases[1].summary["method_documentation_boundary"][
                        "release_notes"
                    ]
                ),
            ),
            audit_row(
                "gtdb-legacy-documentation-separated",
                databases[1].summary["method_documentation_boundary"][
                    "legacy_generic_files"
                ],
                "nhmmer + 200 bp recorded separately",
                "nhmmer"
                in str(
                    databases[1].summary["method_documentation_boundary"][
                        "legacy_generic_files"
                    ]
                )
                and "200"
                in str(
                    databases[1].summary["method_documentation_boundary"][
                        "legacy_generic_files"
                    ]
                ),
            ),
            audit_row(
                "greengenes2-taxonomy-basis",
                databases[2].summary["taxonomy_sources"],
                {"gtdb": "R220", "ltp": "08/2023"},
                databases[2].summary["taxonomy_sources"]["gtdb"] == "R220"
                and databases[2].summary["taxonomy_sources"]["ltp"]
                == "08/2023",
            ),
        ]
    )

    query_audit_rows = [
        {
            "feature_id": feature_id,
            "sequence_length": len(query_sequences[feature_id]),
            "frequency": query_abundance[feature_id]["frequency"],
            "samples_observed": query_abundance[feature_id]["samples_observed"],
            "sequence_sha256": hashlib.sha256(
                query_sequences[feature_id].encode("ascii")
            ).hexdigest(),
        }
        for feature_id in sorted(query_sequences)
    ]

    with tempfile.TemporaryDirectory(prefix="article14-", dir="/tmp") as temporary:
        work_root = Path(temporary)
        environment = isolated_environment(work_root)
        os.environ["MPLCONFIGDIR"] = environment["MPLCONFIGDIR"]
        cache_dir = work_root / "qiime-cache"

        package_result = run_command(
            "conda-package-audit",
            ("conda", "list", "-n", args.qiime_env, "--json"),
            environment,
        )
        require_success(package_result)
        commands.append(package_result)
        packages = {
            str(row["name"]): str(row["version"])
            for row in json.loads(str(package_result["stdout"]))
        }
        environment_rows = [
            {
                "component": "QIIME 2",
                "observed_version": packages.get("qiime2", ""),
                "expected_version": QIIME_VERSION,
                "status": "PASS"
                if packages.get("qiime2") == QIIME_VERSION
                else "FAIL",
            },
            {
                "component": "q2-feature-classifier",
                "observed_version": packages.get("q2-feature-classifier", ""),
                "expected_version": Q2_FEATURE_CLASSIFIER_VERSION,
                "status": "PASS"
                if packages.get("q2-feature-classifier")
                == Q2_FEATURE_CLASSIFIER_VERSION
                else "FAIL",
            },
            {
                "component": "RESCRIPt",
                "observed_version": packages.get("rescript", ""),
                "expected_version": RESCRIPT_VERSION,
                "status": "PASS"
                if packages.get("rescript") == RESCRIPT_VERSION
                else "FAIL",
            },
            {
                "component": "VSEARCH",
                "observed_version": packages.get("vsearch", ""),
                "expected_version": VSEARCH_VERSION,
                "status": "PASS"
                if packages.get("vsearch") == VSEARCH_VERSION
                else "FAIL",
            },
            {
                "component": "scikit-learn",
                "observed_version": packages.get("scikit-learn", ""),
                "expected_version": SKLEARN_VERSION,
                "status": "PASS"
                if packages.get("scikit-learn") == SKLEARN_VERSION
                else "FAIL",
            },
            {
                "component": "Legacy public sklearn classifiers",
                "observed_version": "not loaded",
                "expected_version": "not loaded; consensus-VSEARCH used",
                "status": "PASS",
            },
        ]
        for row in environment_rows:
            checks.append(
                audit_row(
                    "environment-"
                    + str(row["component"]).lower().replace(" ", "-"),
                    row["observed_version"],
                    row["expected_version"],
                    row["status"] == "PASS",
                )
            )

        cache_result = run_command(
            "qiime-cache-create",
            qiime_command(
                args.qiime_env,
                "tools",
                "cache-create",
                "--cache",
                str(cache_dir),
            ),
            environment,
        )
        require_success(cache_result)
        commands.append(cache_result)

        query_plain = work_root / "query-v4-asvs.fasta"
        unpack_gzip(query_fasta, query_plain)
        query_artifact = work_root / "query-v4-asvs.qza"
        import_query = run_command(
            "import-query-asvs",
            qiime_command(
                args.qiime_env,
                "tools",
                "import",
                "--input-path",
                str(query_plain),
                "--output-path",
                str(query_artifact),
                "--type",
                "FeatureData[Sequence]",
                "--input-format",
                "DNAFASTAFormat",
            ),
            environment,
        )
        require_success(import_query)
        commands.append(import_query)

        classifications: dict[str, dict[str, dict[str, object]]] = {}
        search_audits: dict[str, dict[str, dict[str, float | int]]] = {}
        provenance_rows: list[dict[str, object]] = []
        for database in databases:
            database_work = work_root / database.key
            database_work.mkdir(parents=True, exist_ok=True)
            reference_plain = database_work / "reference.fasta"
            taxonomy_plain = database_work / "taxonomy.tsv"
            unpack_gzip(database.fasta_path, reference_plain)
            unpack_gzip(database.taxonomy_path, taxonomy_plain)
            reference_artifact = database_work / "reference.qza"
            taxonomy_artifact = database_work / "taxonomy.qza"
            import_reference = run_command(
                f"import-{database.key}-reference",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "import",
                    "--input-path",
                    str(reference_plain),
                    "--output-path",
                    str(reference_artifact),
                    "--type",
                    "FeatureData[Sequence]",
                    "--input-format",
                    "DNAFASTAFormat",
                ),
                environment,
            )
            require_success(import_reference)
            commands.append(import_reference)
            import_taxonomy = run_command(
                f"import-{database.key}-taxonomy",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "import",
                    "--input-path",
                    str(taxonomy_plain),
                    "--output-path",
                    str(taxonomy_artifact),
                    "--type",
                    "FeatureData[Taxonomy]",
                    "--input-format",
                    "TSVTaxonomyFormat",
                ),
                environment,
            )
            require_success(import_taxonomy)
            commands.append(import_taxonomy)
            for label, artifact in (
                ("reference", reference_artifact),
                ("reference-taxonomy", taxonomy_artifact),
            ):
                validation = run_command(
                    f"validate-{database.key}-{label}",
                    qiime_command(
                        args.qiime_env,
                        "tools",
                        "validate",
                        str(artifact),
                        "--level",
                        "max",
                    ),
                    environment,
                )
                require_success(validation)
                commands.append(validation)

            classification_artifact = output_dir / f"{database.key}-taxonomy.qza"
            search_artifact = output_dir / f"{database.key}-search-results.qza"
            for artifact in (classification_artifact, search_artifact):
                if artifact.exists():
                    artifact.unlink()
            classification_result = run_command(
                f"classify-{database.key}",
                qiime_command(
                    args.qiime_env,
                    "feature-classifier",
                    "classify-consensus-vsearch",
                    "--i-query",
                    str(query_artifact),
                    "--i-reference-reads",
                    str(reference_artifact),
                    "--i-reference-taxonomy",
                    str(taxonomy_artifact),
                    "--p-maxaccepts",
                    str(CLASSIFIER_MAXACCEPTS),
                    "--p-maxrejects",
                    "all",
                    "--p-perc-identity",
                    str(CLASSIFIER_IDENTITY),
                    "--p-query-cov",
                    str(CLASSIFIER_QUERY_COVERAGE),
                    "--p-strand",
                    CLASSIFIER_STRAND,
                    "--p-top-hits-only",
                    "--p-maxhits",
                    "all",
                    "--p-output-no-hits",
                    "--p-threads",
                    str(CLASSIFIER_THREADS),
                    "--p-min-consensus",
                    str(CLASSIFIER_MIN_CONSENSUS),
                    "--p-unassignable-label",
                    "Unassigned",
                    "--o-classification",
                    str(classification_artifact),
                    "--o-search-results",
                    str(search_artifact),
                    "--use-cache",
                    str(cache_dir),
                    "--no-recycle",
                    "--verbose",
                ),
                environment,
                timeout=14400,
            )
            require_success(classification_result)
            commands.append(classification_result)

            for artifact_kind, artifact, expected_type in (
                (
                    "taxonomy",
                    classification_artifact,
                    "FeatureData[Taxonomy]",
                ),
                ("search-results", search_artifact, "FeatureData[BLAST6]"),
            ):
                validation = run_command(
                    f"validate-{database.key}-{artifact_kind}",
                    qiime_command(
                        args.qiime_env,
                        "tools",
                        "validate",
                        str(artifact),
                        "--level",
                        "max",
                    ),
                    environment,
                )
                require_success(validation)
                commands.append(validation)
                metadata = artifact_metadata(artifact)
                actions = provenance_actions(artifact)
                provenance_rows.append(
                    {
                        "database": database.label,
                        "artifact": artifact.name,
                        "uuid": metadata["uuid"],
                        "semantic_type": metadata["type"],
                        "data_format": metadata["format"],
                        "maximum_validation": "PASS",
                        "provenance_actions": "; ".join(
                            f"{plugin}:{action}" for plugin, action in actions
                        ),
                    }
                )
                checks.append(
                    audit_row(
                        f"{database.key}-{artifact_kind}-semantic-type",
                        metadata["type"],
                        expected_type,
                        metadata["type"] == expected_type,
                    )
                )
                checks.append(
                    audit_row(
                        f"{database.key}-{artifact_kind}-provenance",
                        actions,
                        "feature-classifier:classify_consensus_vsearch",
                        any(
                            plugin == "feature-classifier"
                            and action == "classify_consensus_vsearch"
                            for plugin, action in actions
                        ),
                    )
                )

            taxonomy_export = database_work / "classification-export"
            search_export = database_work / "search-export"
            export_taxonomy = run_command(
                f"export-{database.key}-taxonomy",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "export",
                    "--input-path",
                    str(classification_artifact),
                    "--output-path",
                    str(taxonomy_export),
                ),
                environment,
            )
            require_success(export_taxonomy)
            commands.append(export_taxonomy)
            export_search = run_command(
                f"export-{database.key}-search",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "export",
                    "--input-path",
                    str(search_artifact),
                    "--output-path",
                    str(search_export),
                ),
                environment,
            )
            require_success(export_search)
            commands.append(export_search)
            taxonomy_export_file = locate_exported_file(
                taxonomy_export, ("taxonomy.tsv",)
            )
            search_export_file = locate_exported_file(
                search_export, ("blast6.tsv", "search-results.tsv")
            )
            classifications[database.key] = parse_classification(
                taxonomy_export_file
            )
            search_audits[database.key] = parse_blast6(search_export_file)
            checks.append(
                audit_row(
                    f"{database.key}-classification-query-coverage",
                    len(classifications[database.key]),
                    QUERY_RECORDS,
                    set(classifications[database.key]) == set(query_sequences),
                )
            )

        assignment_rows: list[dict[str, object]] = []
        rank_coverage_rows: list[dict[str, object]] = []
        offtarget_rows: list[dict[str, object]] = []
        for database in databases:
            observed = classifications[database.key]
            assigned_ids = [
                feature_id
                for feature_id, item in observed.items()
                if item["category"] != "Unassigned"
            ]
            assigned_reads = sum(
                query_abundance[feature_id]["frequency"]
                for feature_id in assigned_ids
            )
            confidences = [
                float(item["confidence"]) for item in observed.values()
            ]
            hit_counts = [
                int(search_audits[database.key].get(feature_id, {}).get("top_hits", 0))
                for feature_id in query_sequences
            ]
            assignment_rows.append(
                {
                    "database": database.label,
                    "release": database.release,
                    "query_asvs": QUERY_RECORDS,
                    "query_reads": QUERY_READS,
                    "assigned_asvs": len(assigned_ids),
                    "assigned_asv_fraction": round(
                        len(assigned_ids) / QUERY_RECORDS, 6
                    ),
                    "assigned_reads": assigned_reads,
                    "assigned_read_fraction": round(
                        assigned_reads / QUERY_READS, 6
                    ),
                    "unassigned_asvs": QUERY_RECORDS - len(assigned_ids),
                    "median_confidence": round(
                        statistics.median(confidences), 6
                    ),
                    "median_top_hits": round(statistics.median(hit_counts), 2),
                }
            )
            for rank, _prefix in RANKS:
                labeled = [
                    feature_id
                    for feature_id, item in observed.items()
                    if rank in item["ranks"]
                ]
                labeled_reads = sum(
                    query_abundance[feature_id]["frequency"]
                    for feature_id in labeled
                )
                rank_coverage_rows.append(
                    {
                        "database": database.label,
                        "rank": rank,
                        "asv_count": len(labeled),
                        "asv_fraction": round(len(labeled) / QUERY_RECORDS, 6),
                        "read_count": labeled_reads,
                        "read_fraction": round(labeled_reads / QUERY_READS, 6),
                    }
                )
            categories: Counter[str] = Counter()
            category_reads: Counter[str] = Counter()
            for feature_id, item in observed.items():
                category = str(item["category"])
                categories[category] += 1
                category_reads[category] += query_abundance[feature_id]["frequency"]
            for category in (
                "Bacteria",
                "Archaea",
                "Eukaryota",
                "Chloroplast/plastid",
                "Mitochondria",
                "Other/off-target",
                "Unassigned",
            ):
                offtarget_rows.append(
                    {
                        "database": database.label,
                        "category": category,
                        "asv_count": categories[category],
                        "asv_fraction": round(
                            categories[category] / QUERY_RECORDS, 6
                        ),
                        "read_count": category_reads[category],
                        "read_fraction": round(
                            category_reads[category] / QUERY_READS, 6
                        ),
                    }
                )
            checks.extend(
                [
                    audit_row(
                        f"{database.key}-assignment-nonzero",
                        len(assigned_ids),
                        "> 0",
                        len(assigned_ids) > 0,
                    ),
                    audit_row(
                        f"{database.key}-rank-monotonic",
                        [
                            row["asv_count"]
                            for row in rank_coverage_rows
                            if row["database"] == database.label
                        ],
                        "non-increasing",
                        all(
                            left >= right
                            for left, right in zip(
                                [
                                    int(row["asv_count"])
                                    for row in rank_coverage_rows
                                    if row["database"] == database.label
                                ],
                                [
                                    int(row["asv_count"])
                                    for row in rank_coverage_rows
                                    if row["database"] == database.label
                                ][1:],
                            )
                        ),
                    ),
                    audit_row(
                        f"{database.key}-category-ledger-asvs",
                        sum(categories.values()),
                        QUERY_RECORDS,
                        sum(categories.values()) == QUERY_RECORDS,
                    ),
                    audit_row(
                        f"{database.key}-category-ledger-reads",
                        sum(category_reads.values()),
                        QUERY_READS,
                        sum(category_reads.values()) == QUERY_READS,
                    ),
                ]
            )

        crosswalk_rows: list[dict[str, object]] = []
        agreement_categories: Counter[str] = Counter()
        agreement_category_reads: Counter[str] = Counter()
        for feature_id in sorted(query_sequences):
            items = {
                database.key: classifications[database.key][feature_id]
                for database in databases
            }
            genus_values = [
                str(items[database.key]["ranks"].get("genus", ""))
                for database in databases
            ]
            nonempty_genus = [value for value in genus_values if value]
            if len(nonempty_genus) == 3 and len(set(nonempty_genus)) == 1:
                pattern = "All three genus labels agree"
            elif len(nonempty_genus) >= 2 and len(set(nonempty_genus)) > 1:
                pattern = "Database-specific/renamed genus"
            elif nonempty_genus:
                pattern = "Partial genus assignment"
            else:
                pattern = "No genus assignment"
            agreement_categories[pattern] += 1
            agreement_category_reads[pattern] += query_abundance[feature_id][
                "frequency"
            ]
            row: dict[str, object] = {
                "feature_id": feature_id,
                "frequency": query_abundance[feature_id]["frequency"],
                "samples_observed": query_abundance[feature_id][
                    "samples_observed"
                ],
                "genus_agreement_pattern": pattern,
            }
            for database in databases:
                item = items[database.key]
                row[f"{database.key}_taxon"] = item["taxon"]
                row[f"{database.key}_confidence"] = item["confidence"]
                row[f"{database.key}_category"] = item["category"]
                for rank, _prefix in RANKS:
                    row[f"{database.key}_{rank}"] = item["ranks"].get(rank, "")
            crosswalk_rows.append(row)

        pair_specs = (
            ("SILVA vs GTDB", "silva", "gtdb"),
            ("SILVA vs Greengenes2", "silva", "greengenes2"),
            ("GTDB vs Greengenes2", "gtdb", "greengenes2"),
        )
        concordance: dict[str, dict[str, float]] = {}
        concordance_counts: dict[str, dict[str, dict[str, int]]] = {}
        for pair_label, left_key, right_key in pair_specs:
            concordance[pair_label] = {}
            concordance_counts[pair_label] = {}
            for rank in ("phylum", "class", "order", "family", "genus"):
                jointly_assigned = 0
                exact_matches = 0
                jointly_assigned_reads = 0
                exact_match_reads = 0
                for feature_id in query_sequences:
                    left = str(
                        classifications[left_key][feature_id]["ranks"].get(rank, "")
                    )
                    right = str(
                        classifications[right_key][feature_id]["ranks"].get(rank, "")
                    )
                    if left and right:
                        frequency = query_abundance[feature_id]["frequency"]
                        jointly_assigned += 1
                        jointly_assigned_reads += frequency
                        if left == right:
                            exact_matches += 1
                            exact_match_reads += frequency
                concordance[pair_label][rank] = (
                    exact_matches / jointly_assigned if jointly_assigned else 0.0
                )
                concordance_counts[pair_label][rank] = {
                    "jointly_assigned_asvs": jointly_assigned,
                    "exact_name_matches": exact_matches,
                    "jointly_assigned_reads": jointly_assigned_reads,
                    "exact_name_match_reads": exact_match_reads,
                }

        checks.extend(
            [
                audit_row(
                    "crosswalk-records",
                    len(crosswalk_rows),
                    QUERY_RECORDS,
                    len(crosswalk_rows) == QUERY_RECORDS,
                ),
                audit_row(
                    "crosswalk-asv-ledger",
                    sum(agreement_categories.values()),
                    QUERY_RECORDS,
                    sum(agreement_categories.values()) == QUERY_RECORDS,
                ),
                audit_row(
                    "crosswalk-read-ledger",
                    sum(agreement_category_reads.values()),
                    QUERY_READS,
                    sum(agreement_category_reads.values()) == QUERY_READS,
                ),
                audit_row(
                    "database-effect-observed",
                    agreement_categories["Database-specific/renamed genus"],
                    "> 0",
                    agreement_categories["Database-specific/renamed genus"] > 0,
                ),
            ]
        )

        write_tsv(
            output_dir / "environment-audit.tsv",
            environment_rows,
            ("component", "observed_version", "expected_version", "status"),
        )
        write_tsv(
            output_dir / "reference-source-audit.tsv",
            source_rows,
            (
                "database",
                "release",
                "reference_subset",
                "license",
                "full_sequence_records",
                "full_taxonomy_records",
                "raw_v4_records",
                "lca_v4_records",
                "official_url",
                "sequence_sha256",
                "taxonomy_sha256",
            ),
        )
        write_tsv(
            output_dir / "reference-rank-audit.tsv",
            rank_audit_rows,
            (
                "database",
                "rank",
                "full_nonempty_records",
                "full_fraction",
                "lca_v4_nonempty_records",
                "lca_v4_fraction",
            ),
        )
        write_tsv(
            output_dir / "v4-reference-audit.tsv",
            v4_rows,
            (
                "database",
                "records",
                "minimum_length",
                "median_length",
                "mean_length",
                "maximum_length",
                "sequence_taxonomy_id_match",
                "sequence_sha256",
                "taxonomy_sha256",
            ),
        )
        write_tsv(
            output_dir / "query-asv-audit.tsv",
            query_audit_rows,
            (
                "feature_id",
                "sequence_length",
                "frequency",
                "samples_observed",
                "sequence_sha256",
            ),
        )
        write_tsv(
            output_dir / "taxonomy-assignment-summary.tsv",
            assignment_rows,
            (
                "database",
                "release",
                "query_asvs",
                "query_reads",
                "assigned_asvs",
                "assigned_asv_fraction",
                "assigned_reads",
                "assigned_read_fraction",
                "unassigned_asvs",
                "median_confidence",
                "median_top_hits",
            ),
        )
        write_tsv(
            output_dir / "taxonomy-rank-coverage.tsv",
            rank_coverage_rows,
            (
                "database",
                "rank",
                "asv_count",
                "asv_fraction",
                "read_count",
                "read_fraction",
            ),
        )
        crosswalk_fields = (
            "feature_id",
            "frequency",
            "samples_observed",
            "genus_agreement_pattern",
            *tuple(
                field
                for database in databases
                for field in (
                    f"{database.key}_taxon",
                    f"{database.key}_confidence",
                    f"{database.key}_category",
                    *tuple(
                        f"{database.key}_{rank}" for rank, _prefix in RANKS
                    ),
                )
            ),
        )
        write_tsv(
            output_dir / "taxonomy-crosswalk.tsv",
            crosswalk_rows,
            crosswalk_fields,
        )
        write_tsv(
            output_dir / "organelle-offtarget-audit.tsv",
            offtarget_rows,
            (
                "database",
                "category",
                "asv_count",
                "asv_fraction",
                "read_count",
                "read_fraction",
            ),
        )
        write_tsv(
            output_dir / "provenance-audit.tsv",
            provenance_rows,
            (
                "database",
                "artifact",
                "uuid",
                "semantic_type",
                "data_format",
                "maximum_validation",
                "provenance_actions",
            ),
        )

        configure_plotting()
        draw_decision_map(figure_dir)
        draw_reference_scope(figure_dir, source_rows)
        draw_rank_coverage(figure_dir, rank_coverage_rows)
        draw_concordance(figure_dir, concordance)
        expected_figure_files = [
            figure_dir / f"14-{stem}.{extension}"
            for stem in (
                "database-decision-map",
                "reference-scope-audit",
                "rank-assignment-coverage",
                "taxonomy-concordance",
            )
            for extension in ("pdf", "png", "tiff")
        ]
        for path in expected_figure_files:
            checks.append(
                audit_row(
                    f"figure-{path.name}",
                    path.stat().st_size if path.exists() else 0,
                    "> 0 bytes",
                    path.exists() and path.stat().st_size > 0,
                )
            )

        failures = [row for row in checks if row["status"] != "PASS"]
        summary = {
            "status": "passed" if not failures else "failed",
            "validation_date": "2026-07-20",
            "query_asvs": QUERY_RECORDS,
            "query_reads": QUERY_READS,
            "query_length_minimum": min(len(seq) for seq in query_sequences.values()),
            "query_length_median": statistics.median(
                len(seq) for seq in query_sequences.values()
            ),
            "query_length_maximum": max(len(seq) for seq in query_sequences.values()),
            "classification_contract": {
                "method": "q2-feature-classifier classify-consensus-vsearch",
                "maximum_accepts": CLASSIFIER_MAXACCEPTS,
                "maximum_rejects": "all",
                "percent_identity": CLASSIFIER_IDENTITY,
                "query_coverage": CLASSIFIER_QUERY_COVERAGE,
                "strand": CLASSIFIER_STRAND,
                "top_hits_only": True,
                "maximum_hits": "all",
                "minimum_consensus": CLASSIFIER_MIN_CONSENSUS,
                "threads": CLASSIFIER_THREADS,
                "output_no_hits": True,
            },
            "reference_contract": {
                "forward_primer_5_to_3": FORWARD_PRIMER,
                "reverse_primer_5_to_3": REVERSE_PRIMER,
                "combined_primer_identity": EXTRACTION_IDENTITY,
                "minimum_length": EXTRACTION_MIN_LENGTH,
                "maximum_length": EXTRACTION_MAX_LENGTH,
                "dereplication": "exact sequence, least common ancestor",
            },
            "reference_records": {
                row["database"]: {
                    "full": row["full_sequence_records"],
                    "primer_matched_v4": row["raw_v4_records"],
                    "unique_v4_lca": row["lca_v4_records"],
                }
                for row in source_rows
            },
            "assignments": {
                row["database"]: {
                    "assigned_asvs": row["assigned_asvs"],
                    "assigned_asv_fraction": row["assigned_asv_fraction"],
                    "assigned_reads": row["assigned_reads"],
                    "assigned_read_fraction": row["assigned_read_fraction"],
                    "unassigned_asvs": row["unassigned_asvs"],
                }
                for row in assignment_rows
            },
            "rank_coverage": {
                database.label: {
                    row["rank"]: {
                        "asv_count": row["asv_count"],
                        "asv_fraction": row["asv_fraction"],
                        "read_count": row["read_count"],
                        "read_fraction": row["read_fraction"],
                    }
                    for row in rank_coverage_rows
                    if row["database"] == database.label
                }
                for database in databases
            },
            "offtarget": {
                database.label: {
                    row["category"]: {
                        "asv_count": row["asv_count"],
                        "read_count": row["read_count"],
                    }
                    for row in offtarget_rows
                    if row["database"] == database.label
                }
                for database in databases
            },
            "genus_agreement_patterns": {
                label: {
                    "asv_count": agreement_categories[label],
                    "read_count": agreement_category_reads[label],
                }
                for label in sorted(agreement_categories)
            },
            "pairwise_exact_name_concordance": concordance,
            "pairwise_exact_name_concordance_counts": concordance_counts,
            "checks_passed": len(checks) - len(failures),
            "checks_total": len(checks),
            "failed_checks": [row["check_id"] for row in failures],
            "interpretation_boundaries": [
                "Database labels are hypotheses, not proof of species or strain.",
                "GTDB SSU labels inherit representative-genome taxonomy.",
                "Greengenes2 2024.09 uses GTDB R220, not GTDB R232.",
                "Cross-database name changes do not represent biological change.",
                "Organelle, off-target, and Unassigned results remain visible.",
            ],
        }
        write_json(output_dir / "taxonomy-databases-summary.json", summary)
        write_tsv(
            output_dir / "validation-checks.tsv",
            checks,
            ("check_id", "observed", "expected", "status"),
        )

        replacements = {
            str(project_root): "<PROJECT_ROOT>",
            str(output_dir): "<OUTPUT_DIR>",
            str(figure_dir): "<FIGURE_DIR>",
            str(work_root): "<WORK_ROOT>",
            str(Path.home()): "<HOME>",
        }
        log_parts: list[str] = []
        taxonomy_log_parts: list[str] = []
        for result in commands:
            command_text = " ".join(str(token) for token in result["command"])
            block = "\n".join(
                [
                    f"[{result['label']}]",
                    f"command: {command_text}",
                    f"returncode: {result['returncode']}",
                    f"elapsed_seconds: {result['elapsed_seconds']}",
                    "stdout:",
                    str(result["stdout"]),
                    "stderr:",
                    str(result["stderr"]),
                    "",
                ]
            )
            log_parts.append(block)
            if str(result["label"]).startswith("classify-"):
                taxonomy_log_parts.append(block)
        (output_dir / "validation.log").write_text(
            sanitize_text("\n".join(log_parts), replacements),
            encoding="utf-8",
        )
        (output_dir / "qiime-taxonomy.log").write_text(
            sanitize_text("\n".join(taxonomy_log_parts), replacements),
            encoding="utf-8",
        )

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
