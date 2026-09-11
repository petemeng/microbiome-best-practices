#!/usr/bin/env python3
"""One-time import, provenance, FastQC, and MultiQC audit for Article 10."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time
import zipfile
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Iterator, Sequence


EXPECTED_QIIME_ENV = "microbiome-qiime2-2026.4"
EXPECTED_QC_ENV = "microbiome-fastqc-multiqc"
EXPECTED_QIIME_RELEASE = "2026.4"
EXPECTED_QIIME_VERSION = "2026.4.0"
EXPECTED_Q2CLI_VERSION = "2026.4.0"
EXPECTED_DEMUX_VERSION = "2026.4.0"
EXPECTED_FASTQC_VERSION = "0.12.1"
EXPECTED_MULTIQC_VERSION = "1.33"
EXPECTED_QC_ENV_SHA256 = (
    "191d516070ab0c683c684ee4bc6a1584c705698b309f5be32df4aa67460d2e2d"
)
EXPECTED_RECORDS = 2000
EXPECTED_METADATA_SAMPLES = 75
EXPECTED_ASSIGNED = 692
EXPECTED_DEMULTIPLEXED_SAMPLES = 56
EXPECTED_FASTQ_SHA256 = {
    "forward.fastq.gz": (
        "e6fabdbd31db519c719447f1caf2d1cd334c6ab932353e534ea2433ead88d254"
    ),
    "reverse.fastq.gz": (
        "d2fa5841d150fead5a73067a267ec531221de11c512a38ec32ada03f1621bf65"
    ),
    "barcodes.fastq.gz": (
        "9cce1bb9a57975d14b54a5cf07c89b2b7c2095da4a8f0d1bd2720f1349e74057"
    ),
    "metadata.tsv": (
        "8d5342f0a80abf4197088bb1ce92f4903c1bd1212f6a9182e0e1c87ebc322bcb"
    ),
}
ROLE_BY_STEM = {
    "forward": "Forward",
    "reverse": "Reverse",
    "barcodes": "Barcode",
}
EXPECTED_LENGTH = {"Forward": "151", "Reverse": "151", "Barcode": "12"}


@dataclass(frozen=True)
class CommandResult:
    label: str
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    elapsed_seconds: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--qiime-env", default=EXPECTED_QIIME_ENV)
    parser.add_argument("--qc-env", default=EXPECTED_QC_ENV)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(
    path: Path,
    rows: Iterable[dict[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_command(
    label: str,
    command: Sequence[str],
    environment: dict[str, str],
    timeout: int = 600,
) -> CommandResult:
    started = time.monotonic()
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=timeout,
    )
    return CommandResult(
        label=label,
        command=tuple(command),
        returncode=completed.returncode,
        stdout=completed.stdout.strip(),
        stderr=completed.stderr.strip(),
        elapsed_seconds=round(time.monotonic() - started, 3),
    )


def conda_command(environment: str, executable: str, *arguments: str) -> tuple[str, ...]:
    return ("conda", "run", "-n", environment, executable, *arguments)


def qiime_command(environment: str, *arguments: str) -> tuple[str, ...]:
    return conda_command(environment, "qiime", *arguments)


def parse_qiime_info(text: str) -> dict[str, object]:
    fields: dict[str, object] = {}
    for key, pattern in (
        ("python_version", r"^Python version:\s*(.+)$"),
        ("qiime_release", r"^QIIME 2 release:\s*(.+)$"),
        ("qiime_version", r"^QIIME 2 version:\s*(.+)$"),
        ("q2cli_version", r"^q2cli version:\s*(.+)$"),
    ):
        match = re.search(pattern, text, flags=re.MULTILINE)
        if not match:
            raise ValueError(f"qiime info is missing {key}")
        fields[key] = match.group(1).strip()

    plugins: dict[str, str] = {}
    in_plugins = False
    for line in text.splitlines():
        if line.strip() == "Installed plugins":
            in_plugins = True
            continue
        if in_plugins and not line.strip():
            break
        if in_plugins and ":" in line:
            plugin, version = line.split(":", 1)
            plugins[plugin.strip()] = version.strip()
    fields["plugins"] = plugins
    return fields


def parse_peek_tsv(text: str) -> dict[str, str]:
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("qiime tools peek --tsv did not return one row")
    return {str(key): str(value) for key, value in rows[0].items()}


def fastq_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            index += 1
            if not sequence or not separator or not quality:
                raise ValueError(f"{path.name}: truncated record {index}")
            record = tuple(
                line.rstrip("\r\n")
                for line in (header, sequence, separator, quality)
            )
            if not record[0].startswith("@"):
                raise ValueError(f"{path.name}: invalid header at record {index}")
            if not record[2].startswith("+"):
                raise ValueError(f"{path.name}: invalid separator at record {index}")
            if len(record[1]) != len(record[3]):
                raise ValueError(
                    f"{path.name}: sequence/quality length mismatch at record {index}"
                )
            yield record


def summarize_fastq(path: Path) -> dict[str, object]:
    identifiers: list[str] = []
    lengths: list[int] = []
    for header, sequence, _, _ in fastq_records(path):
        identifiers.append(header.split(maxsplit=1)[0].removeprefix("@"))
        lengths.append(len(sequence))
    if not lengths:
        raise ValueError(f"{path.name}: no records")
    return {
        "records": len(lengths),
        "identifiers": identifiers,
        "read_length_min": min(lengths),
        "read_length_max": max(lengths),
    }


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 3:
        raise ValueError("metadata has no data rows")
    header = rows[0]
    records: list[dict[str, str]] = []
    for row in rows[1:]:
        if not row or row[0].startswith("#"):
            continue
        records.append(dict(zip(header, row)))
    return records


def extract_zip_text(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"{path.name}: expected one {suffix}, found {len(matches)}")
        return archive.read(matches[0]).decode("utf-8", errors="replace")


def extract_qzv_table(path: Path, filename: str) -> list[list[str]]:
    text = extract_zip_text(path, f"/data/{filename}")
    return list(csv.reader(text.splitlines(), delimiter="\t"))


def parse_qzv_counts(path: Path) -> list[dict[str, object]]:
    rows = extract_qzv_table(path, "per-sample-fastq-counts.tsv")
    header = rows[0]
    parsed: list[dict[str, object]] = []
    for row in rows[1:]:
        record = dict(zip(header, row))
        parsed.append(
            {
                "sample_id": record["sample ID"],
                "forward_sequence_count": int(record["forward sequence count"]),
                "reverse_sequence_count": int(record["reverse sequence count"]),
            }
        )
    return parsed


def parse_qzv_quality(path: Path, direction: str) -> dict[str, list[float]]:
    rows = extract_qzv_table(path, f"{direction}-seven-number-summaries.tsv")
    positions = [int(value) for value in rows[0][1:]]
    values = {row[0]: [float(value) for value in row[1:]] for row in rows[1:]}
    return {
        "position": [float(value) for value in positions],
        "count": values["count"],
        "median": values["50%"],
        "q25": values["25%"],
        "q75": values["75%"],
    }


def parse_details(path: Path) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    grouped: Counter[tuple[str, int]] = Counter()
    for row in rows:
        mapped = "Mapped" if row["sample"] not in {"", "None"} else "Unmapped"
        grouped[(mapped, int(row["barcode-errors"]))] += 1
    audit_rows = [
        {
            "mapping_status": mapping_status,
            "barcode_errors": errors,
            "records": count,
            "percent_total": round(100 * count / len(rows), 3),
        }
        for (mapping_status, errors), count in sorted(grouped.items())
    ]
    return rows, audit_rows


def parse_fastqc_zip(path: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    summary = extract_zip_text(path, "/summary.txt")
    module_rows = []
    for line in summary.splitlines():
        if not line.strip():
            continue
        status, module, filename = line.split("\t")
        module_rows.append(
            {"status": status, "module": module, "filename": filename}
        )

    data = extract_zip_text(path, "/fastqc_data.txt")
    basic: dict[str, str] = {}
    in_basic = False
    for line in data.splitlines():
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic and line.startswith(">>END_MODULE"):
            break
        if in_basic and line and not line.startswith("#"):
            key, value = line.split("\t", 1)
            basic[key] = value
    return module_rows, basic


def parse_multiqc_general_stats(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    parsed = []
    for row in rows:
        parsed.append(
            {
                "sample": row["Sample"],
                "total_sequences": int(
                    round(float(row["fastqc-total_sequences"]) * 1_000_000)
                ),
                "percent_gc": float(row["fastqc-percent_gc"]),
                "percent_duplicates": float(row["fastqc-percent_duplicates"]),
                "average_sequence_length": float(
                    row["fastqc-avg_sequence_length"]
                ),
                "median_sequence_length": float(
                    row["fastqc-median_sequence_length"]
                ),
                "percent_failed_modules": float(row["fastqc-percent_fails"]),
            }
        )
    return parsed


def provenance_actions(path: Path) -> list[dict[str, str]]:
    actions: list[dict[str, str]] = []
    with zipfile.ZipFile(path) as archive:
        names = [
            name
            for name in archive.namelist()
            if "/provenance/" in name and name.endswith("/action/action.yaml")
        ]
        for name in names:
            text = archive.read(name).decode("utf-8", errors="replace")
            block_match = re.search(
                r"^action:\n(?P<block>.*?)(?=^[a-z][a-z_-]*:\n|\Z)",
                text,
                flags=re.MULTILINE | re.DOTALL,
            )
            block = block_match.group("block") if block_match else ""
            kind = re.search(r"^    type:\s*(.+)$", block, flags=re.MULTILINE)
            action = re.search(r"^    action:\s*(.+)$", block, flags=re.MULTILINE)
            actions.append(
                {
                    "member": name,
                    "kind": kind.group(1).strip() if kind else "unknown",
                    "action": action.group(1).strip() if action else "import",
                    "text": text,
                }
            )
    return actions


def artifact_audit_row(
    stage: str,
    path: Path,
    peek: dict[str, str],
    action: str,
    validation: str,
    action_count: int,
) -> dict[str, object]:
    return {
        "stage": stage,
        "file": path.name,
        "uuid": peek.get("UUID", ""),
        "semantic_type": peek.get("Type", ""),
        "data_format": peek.get("Data Format", ""),
        "terminal_action": action,
        "provenance_actions": action_count,
        "maximum_validation": validation,
    }


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


def sanitize_text(text: str, replacements: dict[str, str]) -> str:
    sanitized = text
    for source, target in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        if source:
            sanitized = sanitized.replace(source, target)
    return sanitized


def save_figure(fig: object, file_base: Path) -> None:
    file_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(file_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(
        file_base.with_suffix(".png"),
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        file_base.with_suffix(".tiff"),
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def draw_provenance_map(
    figure_dir: Path,
    assigned: int,
    samples: int,
    action_count: int,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    dark = "#111827"
    gray = "#4B5563"
    blue = "#0072B2"
    green = "#009E73"
    orange = "#E69F00"
    purple = "#CC79A7"
    red = "#D55E00"

    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.text(
        0.3,
        6.62,
        "Raw-read acceptance and provenance map",
        fontsize=19,
        fontweight="bold",
        color=dark,
        va="top",
    )
    ax.text(
        0.3,
        6.20,
        "File-level QC and QIIME 2 sample-level QC answer different questions.",
        fontsize=10.6,
        color=gray,
        va="top",
    )

    nodes = [
        (0.35, 4.25, 2.25, 1.2, "Multiplexed input", "3 × 2,000 FASTQ records\n75 mapping barcodes", blue),
        (3.05, 4.25, 2.25, 1.2, "QIIME 2 import", "EMPPairedEndSequences\nmaximum validation", purple),
        (5.75, 4.25, 2.25, 1.2, "EMP demultiplex", f"{assigned:,} mapped pairs\n{samples} observed samples", green),
        (8.45, 4.25, 2.25, 1.2, "Quality summary", f"All {assigned:,} mapped pairs\n{action_count} provenance actions", orange),
        (1.72, 1.55, 2.25, 1.2, "FastQC", "3 independent reports\n33 module states", red),
        (5.25, 1.55, 2.25, 1.2, "MultiQC", "3 reports aggregated\none portable HTML", blue),
    ]
    for x, y, width, height, title, detail, color in nodes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.035,rounding_size=0.12",
                linewidth=1.3,
                edgecolor=color,
                facecolor=color + "12",
            )
        )
        ax.text(
            x + width / 2,
            y + 0.84,
            title,
            fontsize=11.0,
            fontweight="bold",
            color=dark,
            ha="center",
        )
        ax.text(
            x + width / 2,
            y + 0.38,
            detail,
            fontsize=9.1,
            color=gray,
            ha="center",
            va="center",
            linespacing=1.25,
        )

    for start, end, color in (
        ((2.62, 4.85), (3.00, 4.85), purple),
        ((5.32, 4.85), (5.70, 4.85), green),
        ((8.02, 4.85), (8.40, 4.85), orange),
        ((1.45, 4.20), (2.75, 2.82), red),
        ((4.00, 2.15), (5.20, 2.15), blue),
    ):
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.25,
                color=color,
            )
        )

    ax.add_patch(
        FancyBboxPatch(
            (0.38, 0.45),
            10.30,
            0.58,
            boxstyle="round,pad=0.025,rounding_size=0.08",
            linewidth=0,
            facecolor="#F3F4F6",
        )
    )
    ax.text(
        5.53,
        0.74,
        "Stable evidence: input hashes · semantic types · parameters · mapped counts · module tables",
        fontsize=9.5,
        fontweight="bold",
        color=dark,
        ha="center",
        va="center",
    )
    fig.subplots_adjust(left=0.01, right=0.995, top=0.99, bottom=0.01)
    save_figure(fig, figure_dir / "10-import-provenance-map")
    plt.close(fig)


def draw_fastqc_heatmap(
    figure_dir: Path,
    module_rows: list[dict[str, str]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    roles = ["Forward", "Reverse", "Barcode"]
    modules = []
    for row in module_rows:
        if row["module"] not in modules:
            modules.append(row["module"])
    status_value = {"FAIL": 0, "WARN": 1, "PASS": 2}
    lookup = {(row["role"], row["module"]): row["status"] for row in module_rows}
    matrix = np.array(
        [[status_value[lookup[(role, module)]] for module in modules] for role in roles]
    )

    fig, ax = plt.subplots(figsize=(12.2, 6.0))
    fig.patch.set_facecolor("white")
    cmap = ListedColormap(["#D55E00", "#E69F00", "#009E73"])
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=2.5)
    ax.set_xticks(range(len(modules)))
    ax.set_xticklabels(modules, rotation=35, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(roles)))
    ax.set_yticklabels(roles, fontsize=10.5)
    fig.text(
        0.09,
        0.97,
        "FastQC module states require amplicon-aware interpretation",
        fontsize=16,
        fontweight="bold",
        color="#111827",
        ha="left",
        va="top",
    )
    fig.text(
        0.09,
        0.905,
        "Per-base quality passes for every stream; composition, GC, duplication, and overrepresentation are not independent exclusion rules.",
        fontsize=9.4,
        color="#4B5563",
        ha="left",
        va="top",
    )
    for row_index, role in enumerate(roles):
        for column_index, module in enumerate(modules):
            status = lookup[(role, module)]
            ax.text(
                column_index,
                row_index,
                status[0],
                ha="center",
                va="center",
                fontsize=9.0,
                fontweight="bold",
                color="white" if status != "WARN" else "#111827",
            )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    fig.legend(
        handles=[
            Patch(facecolor="#009E73", label="PASS"),
            Patch(facecolor="#E69F00", label="WARN"),
            Patch(facecolor="#D55E00", label="FAIL"),
        ],
        frameon=False,
        ncol=3,
        loc="lower left",
        bbox_to_anchor=(0.09, 0.02),
    )
    fig.subplots_adjust(left=0.09, right=0.99, top=0.80, bottom=0.39)
    save_figure(fig, figure_dir / "10-fastqc-module-audit")
    plt.close(fig)


def draw_demultiplexing_qc(
    figure_dir: Path,
    counts: list[dict[str, object]],
    forward_quality: dict[str, list[float]],
    reverse_quality: dict[str, list[float]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    ordered_counts = sorted(
        (int(row["forward_sequence_count"]) for row in counts), reverse=True
    )
    median_count = float(np.median(ordered_counts))
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), gridspec_kw={"wspace": 0.28})
    fig.patch.set_facecolor("white")
    fig.suptitle(
        "Demultiplexing yield and positional quality",
        fontsize=17,
        fontweight="bold",
        color="#111827",
        x=0.055,
        ha="left",
        y=0.98,
    )

    ax = axes[0]
    ax.bar(
        range(1, len(ordered_counts) + 1),
        ordered_counts,
        color="#0072B2",
        width=0.86,
    )
    ax.axhline(
        median_count,
        color="#D55E00",
        linestyle="--",
        linewidth=1.3,
        label=f"Median = {median_count:g}",
    )
    ax.set_xlabel("Observed samples (ranked)")
    ax.set_ylabel("Mapped read pairs")
    ax.set_title("A  Deterministic 2,000-pair excerpt", loc="left", fontweight="bold")
    ax.text(
        0.98,
        0.94,
        f"{sum(ordered_counts):,} mapped\n{len(ordered_counts)} samples\nrange {min(ordered_counts)}–{max(ordered_counts)}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.2,
        color="#4B5563",
    )
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)

    ax = axes[1]
    for quality, label, color in (
        (forward_quality, "Forward median", "#0072B2"),
        (reverse_quality, "Reverse median", "#D55E00"),
    ):
        positions = quality["position"]
        median = quality["median"]
        ax.plot(positions, median, color=color, linewidth=1.8, label=label)
        ax.fill_between(
            positions,
            quality["q25"],
            quality["q75"],
            color=color,
            alpha=0.10,
            linewidth=0,
        )
    ax.axhline(30, color="#6B7280", linestyle="--", linewidth=1.1, label="Q30")
    ax.set_xlabel("Read position (nt)")
    ax.set_ylabel("Phred score")
    ax.set_ylim(0, 42)
    ax.set_title("B  QIIME 2 summary uses all mapped pairs", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="lower left")
    ax.grid(color="#E5E7EB", linewidth=0.8)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.text(
        0.055,
        0.01,
        "Excerpt counts validate the workflow; they are not full-study sample-retention thresholds.",
        fontsize=9.0,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.07, right=0.985, top=0.84, bottom=0.16)
    save_figure(fig, figure_dir / "10-demultiplexing-qc")
    plt.close(fig)


def scan_for_paths(paths: Iterable[Path], forbidden: Sequence[str]) -> list[str]:
    leaks: list[str] = []
    needles = [value.encode() for value in forbidden if value]
    for path in paths:
        if not path.is_file():
            continue
        payloads: list[tuple[str, bytes]] = [(path.name, path.read_bytes())]
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                payloads.extend(
                    (f"{path.name}:{name}", archive.read(name))
                    for name in archive.namelist()
                    if not name.endswith("/")
                )
        for label, payload in payloads:
            for needle in needles:
                if needle in payload:
                    leaks.append(label)
                    break
    return sorted(set(leaks))


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    source_dir = root / "data/small/fastq"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "fastqc").mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    final_artifacts = {
        "raw": output_dir / "atacama-emp-paired-end.qza",
        "demux": output_dir / "atacama-demux-paired.qza",
        "details": output_dir / "atacama-demux-details.qza",
        "summary": output_dir / "atacama-demux-summary.qzv",
    }
    command_results: list[CommandResult] = []
    checks: list[dict[str, object]] = []

    def check(
        check_id: str,
        observed: object,
        expected: object,
        passed: bool,
    ) -> None:
        checks.append(audit_row(check_id, observed, expected, bool(passed)))

    with tempfile.TemporaryDirectory(prefix="article10-", dir=output_dir) as temporary:
        work = Path(temporary)
        input_dir = work / "emp-paired-end-sequences"
        fastqc_dir = work / "fastqc"
        multiqc_dir = work / "multiqc"
        details_export = work / "details-export"
        cache_root = work / "cache"
        home_dir = work / "home"
        for path in (
            input_dir,
            fastqc_dir,
            multiqc_dir,
            details_export,
            cache_root / "xdg",
            cache_root / "matplotlib",
            cache_root / "numba",
            cache_root / "fontconfig",
            home_dir / ".cache/fontconfig",
        ):
            path.mkdir(parents=True, exist_ok=True)
        for filename in ("forward.fastq.gz", "reverse.fastq.gz", "barcodes.fastq.gz"):
            shutil.copy2(source_dir / filename, input_dir / filename)

        environment = {
            **os.environ,
            "LC_ALL": "C",
            "LANG": "C",
            "TZ": "Asia/Shanghai",
            "HOME": str(home_dir),
            "PYTHONNOUSERSITE": "1",
            "MPLBACKEND": "Agg",
            "XDG_CACHE_HOME": str(cache_root / "xdg"),
            "MPLCONFIGDIR": str(cache_root / "matplotlib"),
            "NUMBA_CACHE_DIR": str(cache_root / "numba"),
            "FC_CACHEDIR": str(cache_root / "fontconfig"),
        }
        if Path("/etc/fonts/fonts.conf").exists():
            environment["FONTCONFIG_FILE"] = "/etc/fonts/fonts.conf"

        work_artifacts = {
            "raw": work / final_artifacts["raw"].name,
            "demux": work / final_artifacts["demux"].name,
            "details": work / final_artifacts["details"].name,
            "summary": work / final_artifacts["summary"].name,
        }

        info_result = run_command(
            "qiime-info", qiime_command(args.qiime_env, "info"), environment
        )
        command_results.append(info_result)
        if info_result.returncode != 0:
            raise RuntimeError(f"qiime info failed: {info_result.stderr}")
        qiime_info = parse_qiime_info(info_result.stdout)

        fastqc_version_result = run_command(
            "fastqc-version",
            conda_command(args.qc_env, "fastqc", "--version"),
            environment,
        )
        multiqc_version_result = run_command(
            "multiqc-version",
            conda_command(args.qc_env, "python", "-m", "multiqc", "--version"),
            environment,
        )
        command_results.extend([fastqc_version_result, multiqc_version_result])
        fastqc_version = re.search(
            r"FastQC v([0-9.]+)",
            fastqc_version_result.stdout + fastqc_version_result.stderr,
        )
        multiqc_version = re.search(
            r"(?:multiqc, version|multiqc=)\s*([0-9.]+)",
            multiqc_version_result.stdout + multiqc_version_result.stderr,
            flags=re.IGNORECASE,
        )
        observed_fastqc_version = fastqc_version.group(1) if fastqc_version else ""
        observed_multiqc_version = multiqc_version.group(1) if multiqc_version else ""

        source_summary_path = source_dir / "source_summary.json"
        source_summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
        input_summaries = {
            filename: summarize_fastq(source_dir / filename)
            for filename in ("forward.fastq.gz", "reverse.fastq.gz", "barcodes.fastq.gz")
        }
        synchronized = (
            input_summaries["forward.fastq.gz"]["identifiers"]
            == input_summaries["reverse.fastq.gz"]["identifiers"]
            == input_summaries["barcodes.fastq.gz"]["identifiers"]
        )
        metadata_rows = read_metadata(source_dir / "metadata.tsv")
        sample_ids = [row["sample-id"] for row in metadata_rows]
        barcodes = [row["barcode-sequence"] for row in metadata_rows]

        check("qiime-info-command", info_result.returncode, 0, info_result.returncode == 0)
        check(
            "qiime-release",
            qiime_info["qiime_release"],
            EXPECTED_QIIME_RELEASE,
            qiime_info["qiime_release"] == EXPECTED_QIIME_RELEASE,
        )
        check(
            "qiime-version",
            qiime_info["qiime_version"],
            EXPECTED_QIIME_VERSION,
            qiime_info["qiime_version"] == EXPECTED_QIIME_VERSION,
        )
        check(
            "q2cli-version",
            qiime_info["q2cli_version"],
            EXPECTED_Q2CLI_VERSION,
            qiime_info["q2cli_version"] == EXPECTED_Q2CLI_VERSION,
        )
        check(
            "demux-plugin-version",
            qiime_info["plugins"].get("demux"),
            EXPECTED_DEMUX_VERSION,
            qiime_info["plugins"].get("demux") == EXPECTED_DEMUX_VERSION,
        )
        check(
            "fastqc-version",
            observed_fastqc_version,
            EXPECTED_FASTQC_VERSION,
            fastqc_version_result.returncode == 0
            and observed_fastqc_version == EXPECTED_FASTQC_VERSION,
        )
        check(
            "multiqc-version",
            observed_multiqc_version,
            EXPECTED_MULTIQC_VERSION,
            multiqc_version_result.returncode == 0
            and observed_multiqc_version == EXPECTED_MULTIQC_VERSION,
        )
        qc_env_hash = sha256(root / "env/raw-qc.yml")
        check(
            "raw-qc-environment-hash",
            qc_env_hash,
            EXPECTED_QC_ENV_SHA256,
            qc_env_hash == EXPECTED_QC_ENV_SHA256,
        )
        check(
            "source-summary-readable",
            source_summary["selection"]["selected_records"],
            EXPECTED_RECORDS,
            source_summary["selection"]["selected_records"] == EXPECTED_RECORDS,
        )
        for check_id, filename in (
            ("forward-sha256", "forward.fastq.gz"),
            ("reverse-sha256", "reverse.fastq.gz"),
            ("barcode-sha256", "barcodes.fastq.gz"),
            ("metadata-sha256", "metadata.tsv"),
        ):
            observed_hash = sha256(source_dir / filename)
            check(
                check_id,
                observed_hash,
                EXPECTED_FASTQ_SHA256[filename],
                observed_hash == EXPECTED_FASTQ_SHA256[filename],
            )
        for check_id, filename in (
            ("forward-record-count", "forward.fastq.gz"),
            ("reverse-record-count", "reverse.fastq.gz"),
            ("barcode-record-count", "barcodes.fastq.gz"),
        ):
            observed_records = input_summaries[filename]["records"]
            check(
                check_id,
                observed_records,
                EXPECTED_RECORDS,
                observed_records == EXPECTED_RECORDS,
            )
        check("fastq-id-synchronization", synchronized, True, synchronized)
        metadata_valid = (
            len(metadata_rows) == EXPECTED_METADATA_SAMPLES
            and len(set(sample_ids)) == EXPECTED_METADATA_SAMPLES
        )
        check(
            "metadata-sample-count",
            len(metadata_rows),
            EXPECTED_METADATA_SAMPLES,
            metadata_valid,
        )
        barcode_valid = (
            len(set(barcodes)) == EXPECTED_METADATA_SAMPLES
            and all(re.fullmatch(r"[ACGT]{12}", barcode) for barcode in barcodes)
        )
        check(
            "metadata-barcode-contract",
            f"{len(set(barcodes))} unique; lengths {sorted(set(map(len, barcodes)))}",
            "75 unique 12-nt DNA barcodes",
            barcode_valid,
        )

        import_result = run_command(
            "qiime-import",
            qiime_command(
                args.qiime_env,
                "tools",
                "import",
                "--type",
                "EMPPairedEndSequences",
                "--input-path",
                str(input_dir),
                "--input-format",
                "EMPPairedEndDirFmt",
                "--output-path",
                str(work_artifacts["raw"]),
            ),
            environment,
        )
        command_results.append(import_result)
        if import_result.returncode != 0:
            raise RuntimeError(f"QIIME 2 import failed: {import_result.stderr}")

        raw_validate = run_command(
            "raw-maximum-validation",
            qiime_command(
                args.qiime_env,
                "tools",
                "validate",
                str(work_artifacts["raw"]),
                "--level",
                "max",
            ),
            environment,
        )
        raw_peek_result = run_command(
            "raw-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(work_artifacts["raw"]),
            ),
            environment,
        )
        command_results.extend([raw_validate, raw_peek_result])
        if raw_validate.returncode != 0 or raw_peek_result.returncode != 0:
            raise RuntimeError("raw Artifact validation or peek failed")
        raw_peek = parse_peek_tsv(raw_peek_result.stdout)
        raw_actions = provenance_actions(work_artifacts["raw"])

        demux_result = run_command(
            "emp-paired-demultiplex",
            qiime_command(
                args.qiime_env,
                "demux",
                "emp-paired",
                "--i-seqs",
                str(work_artifacts["raw"]),
                "--m-barcodes-file",
                str(source_dir / "metadata.tsv"),
                "--m-barcodes-column",
                "barcode-sequence",
                "--p-rev-comp-mapping-barcodes",
                "--o-per-sample-sequences",
                str(work_artifacts["demux"]),
                "--o-error-correction-details",
                str(work_artifacts["details"]),
            ),
            environment,
        )
        command_results.append(demux_result)
        if demux_result.returncode != 0:
            raise RuntimeError(f"EMP paired demultiplex failed: {demux_result.stderr}")

        demux_validate = run_command(
            "demux-maximum-validation",
            qiime_command(
                args.qiime_env,
                "tools",
                "validate",
                str(work_artifacts["demux"]),
                "--level",
                "max",
            ),
            environment,
        )
        demux_peek_result = run_command(
            "demux-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(work_artifacts["demux"]),
            ),
            environment,
        )
        details_peek_result = run_command(
            "demux-details-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(work_artifacts["details"]),
            ),
            environment,
        )
        command_results.extend(
            [demux_validate, demux_peek_result, details_peek_result]
        )
        if any(
            result.returncode != 0
            for result in (demux_validate, demux_peek_result, details_peek_result)
        ):
            raise RuntimeError("demultiplexed Artifact validation or peek failed")
        demux_peek = parse_peek_tsv(demux_peek_result.stdout)
        details_peek = parse_peek_tsv(details_peek_result.stdout)

        summarize_result = run_command(
            "demux-summarize-all-reads",
            qiime_command(
                args.qiime_env,
                "demux",
                "summarize",
                "--i-data",
                str(work_artifacts["demux"]),
                "--p-n",
                "10000",
                "--o-visualization",
                str(work_artifacts["summary"]),
            ),
            environment,
        )
        summary_validate = run_command(
            "summary-maximum-validation",
            qiime_command(
                args.qiime_env,
                "tools",
                "validate",
                str(work_artifacts["summary"]),
                "--level",
                "max",
            ),
            environment,
        )
        summary_peek_result = run_command(
            "summary-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(work_artifacts["summary"]),
            ),
            environment,
        )
        command_results.extend(
            [summarize_result, summary_validate, summary_peek_result]
        )
        if any(
            result.returncode != 0
            for result in (summarize_result, summary_validate, summary_peek_result)
        ):
            raise RuntimeError("demux summarize, validation, or peek failed")
        summary_peek = parse_peek_tsv(summary_peek_result.stdout)
        final_actions = provenance_actions(work_artifacts["summary"])
        rev_comp_recorded = any(
            re.search(r"rev_comp_mapping_barcodes:\s*true", row["text"])
            for row in final_actions
        )

        details_export_result = run_command(
            "export-demux-details",
            qiime_command(
                args.qiime_env,
                "tools",
                "export",
                "--input-path",
                str(work_artifacts["details"]),
                "--output-path",
                str(details_export),
            ),
            environment,
        )
        command_results.append(details_export_result)
        if details_export_result.returncode != 0:
            raise RuntimeError(f"details export failed: {details_export_result.stderr}")

        count_rows = parse_qzv_counts(work_artifacts["summary"])
        forward_quality = parse_qzv_quality(work_artifacts["summary"], "forward")
        reverse_quality = parse_qzv_quality(work_artifacts["summary"], "reverse")
        detail_rows, correction_rows = parse_details(details_export / "details.tsv")
        assigned = sum(row["sample"] not in {"", "None"} for row in detail_rows)
        unmapped = len(detail_rows) - assigned
        paired_counts_match = all(
            row["forward_sequence_count"] == row["reverse_sequence_count"]
            for row in count_rows
        )
        quality_count = {
            int(round(value))
            for value in forward_quality["count"] + reverse_quality["count"]
        }

        check("qiime-import-exit", import_result.returncode, 0, import_result.returncode == 0)
        check(
            "raw-semantic-type",
            raw_peek.get("Type"),
            "EMPPairedEndSequences",
            raw_peek.get("Type") == "EMPPairedEndSequences",
        )
        check(
            "raw-data-format",
            raw_peek.get("Data Format"),
            "EMPPairedEndDirFmt",
            raw_peek.get("Data Format") == "EMPPairedEndDirFmt",
        )
        check(
            "raw-maximum-validation",
            raw_validate.returncode,
            0,
            raw_validate.returncode == 0,
        )
        check(
            "raw-import-provenance",
            [row["action"] for row in raw_actions],
            ["import"],
            len(raw_actions) == 1 and raw_actions[0]["action"] == "import",
        )
        check(
            "emp-paired-exit",
            demux_result.returncode,
            0,
            demux_result.returncode == 0,
        )
        check(
            "demux-semantic-type",
            demux_peek.get("Type"),
            "SampleData[PairedEndSequencesWithQuality]",
            demux_peek.get("Type") == "SampleData[PairedEndSequencesWithQuality]",
        )
        check(
            "demux-data-format",
            demux_peek.get("Data Format"),
            "SingleLanePerSamplePairedEndFastqDirFmt",
            demux_peek.get("Data Format")
            == "SingleLanePerSamplePairedEndFastqDirFmt",
        )
        check(
            "demux-maximum-validation",
            demux_validate.returncode,
            0,
            demux_validate.returncode == 0,
        )
        check(
            "details-semantic-type",
            details_peek.get("Type"),
            "ErrorCorrectionDetails",
            details_peek.get("Type") == "ErrorCorrectionDetails",
        )
        check(
            "summary-result-type",
            summary_peek.get("Type"),
            "Visualization",
            summary_peek.get("Type") == "Visualization",
        )
        check(
            "summary-maximum-validation",
            summary_validate.returncode,
            0,
            summary_validate.returncode == 0,
        )
        check(
            "reverse-complement-mapping-barcodes",
            rev_comp_recorded,
            True,
            rev_comp_recorded,
        )
        check(
            "provenance-action-count",
            len(final_actions),
            3,
            len(final_actions) == 3,
        )
        check(
            "barcode-detail-row-count",
            len(detail_rows),
            EXPECTED_RECORDS,
            len(detail_rows) == EXPECTED_RECORDS,
        )
        check(
            "assigned-read-pairs",
            assigned,
            EXPECTED_ASSIGNED,
            assigned == EXPECTED_ASSIGNED,
        )
        check(
            "unmapped-read-pairs",
            unmapped,
            EXPECTED_RECORDS - EXPECTED_ASSIGNED,
            unmapped == EXPECTED_RECORDS - EXPECTED_ASSIGNED,
        )
        check(
            "demultiplexed-sample-count",
            len(count_rows),
            EXPECTED_DEMULTIPLEXED_SAMPLES,
            len(count_rows) == EXPECTED_DEMULTIPLEXED_SAMPLES,
        )
        check(
            "paired-count-identity",
            paired_counts_match,
            True,
            paired_counts_match
            and sum(int(row["forward_sequence_count"]) for row in count_rows)
            == EXPECTED_ASSIGNED,
        )
        check(
            "quality-summary-read-count",
            sorted(quality_count),
            [EXPECTED_ASSIGNED],
            quality_count == {EXPECTED_ASSIGNED},
        )

        fastqc_result = run_command(
            "fastqc-three-streams",
            conda_command(
                args.qc_env,
                "fastqc",
                "--threads",
                "3",
                "--noextract",
                "--outdir",
                str(fastqc_dir),
                str(source_dir / "forward.fastq.gz"),
                str(source_dir / "reverse.fastq.gz"),
                str(source_dir / "barcodes.fastq.gz"),
            ),
            environment,
        )
        command_results.append(fastqc_result)
        if fastqc_result.returncode != 0:
            raise RuntimeError(f"FastQC failed: {fastqc_result.stderr}")

        fastqc_rows: list[dict[str, str]] = []
        fastqc_basic: list[dict[str, str]] = []
        for path in sorted(fastqc_dir.glob("*_fastqc.zip")):
            stem = path.name.removesuffix("_fastqc.zip")
            role = ROLE_BY_STEM[stem]
            modules, basic = parse_fastqc_zip(path)
            for row in modules:
                fastqc_rows.append({"role": role, **row})
            fastqc_basic.append({"role": role, **basic})

        multiqc_result = run_command(
            "multiqc-fastqc-aggregation",
            conda_command(
                args.qc_env,
                "python",
                "-m",
                "multiqc",
                str(fastqc_dir),
                "--module",
                "fastqc",
                "--require-logs",
                "--force",
                "--no-version-check",
                "--no-ansi",
                "--data-dir",
                "--outdir",
                str(multiqc_dir),
                "--filename",
                "multiqc_report.html",
            ),
            environment,
        )
        command_results.append(multiqc_result)
        if multiqc_result.returncode != 0:
            raise RuntimeError(f"MultiQC failed: {multiqc_result.stderr}")

        multiqc_data_dirs = list(multiqc_dir.glob("*_data"))
        if len(multiqc_data_dirs) != 1:
            raise RuntimeError("MultiQC data directory was not created exactly once")
        multiqc_general = parse_multiqc_general_stats(
            multiqc_data_dirs[0] / "multiqc_general_stats.txt"
        )
        expected_lengths = {
            row["role"]: EXPECTED_LENGTH[row["role"]] for row in fastqc_basic
        }
        fastqc_reports_valid = (
            len(fastqc_basic) == 3
            and all(int(row["Total Sequences"]) == EXPECTED_RECORDS for row in fastqc_basic)
            and all(
                row["Sequence length"] == expected_lengths[row["role"]]
                for row in fastqc_basic
            )
        )
        per_base_pass = sum(
            row["module"] == "Per base sequence quality" and row["status"] == "PASS"
            for row in fastqc_rows
        )
        multiqc_valid = (
            len(multiqc_general) == 3
            and all(row["total_sequences"] == EXPECTED_RECORDS for row in multiqc_general)
        )
        check(
            "fastqc-report-contract",
            f"{len(fastqc_basic)} reports; "
            f"{[row['Total Sequences'] for row in fastqc_basic]} sequences",
            "3 reports; 2,000 sequences each; lengths 151/151/12",
            fastqc_reports_valid,
        )
        check(
            "fastqc-module-count",
            len(fastqc_rows),
            33,
            len(fastqc_rows) == 33,
        )
        check(
            "fastqc-per-base-quality-pass",
            per_base_pass,
            3,
            per_base_pass == 3,
        )
        check(
            "multiqc-fastqc-report-count",
            len(multiqc_general),
            3,
            multiqc_result.returncode == 0 and multiqc_valid,
        )

        for key, source in work_artifacts.items():
            shutil.copy2(source, final_artifacts[key])
        for source in fastqc_dir.glob("*_fastqc.*"):
            target = output_dir / "fastqc" / source.name
            if source.suffix == ".html":
                replacements = {
                    str(root): "<PROJECT_ROOT>",
                    str(work): "<WORK_DIR>",
                    str(Path.home()): "<HOME>",
                }
                target.write_text(
                    sanitize_text(source.read_text(encoding="utf-8"), replacements),
                    encoding="utf-8",
                )
            else:
                shutil.copy2(source, target)

        replacements = {
            str(root): "<PROJECT_ROOT>",
            str(work): "<WORK_DIR>",
            str(Path.home()): "<HOME>",
        }
        multiqc_html_source = multiqc_dir / "multiqc_report.html"
        (output_dir / "multiqc_report.html").write_text(
            sanitize_text(
                multiqc_html_source.read_text(encoding="utf-8"), replacements
            ),
            encoding="utf-8",
        )

        artifact_rows = [
            artifact_audit_row(
                "Raw import",
                final_artifacts["raw"],
                raw_peek,
                "import",
                "PASS",
                len(raw_actions),
            ),
            artifact_audit_row(
                "Demultiplexed reads",
                final_artifacts["demux"],
                demux_peek,
                "demux.emp_paired",
                "PASS",
                2,
            ),
            artifact_audit_row(
                "Barcode correction",
                final_artifacts["details"],
                details_peek,
                "demux.emp_paired",
                "PASS",
                2,
            ),
            artifact_audit_row(
                "Quality visualization",
                final_artifacts["summary"],
                summary_peek,
                "demux.summarize",
                "PASS",
                len(final_actions),
            ),
        ]
        write_tsv(
            output_dir / "import-provenance-audit.tsv",
            artifact_rows,
            [
                "stage",
                "file",
                "uuid",
                "semantic_type",
                "data_format",
                "terminal_action",
                "provenance_actions",
                "maximum_validation",
            ],
        )
        write_tsv(
            output_dir / "demux-sample-counts.tsv",
            count_rows,
            ["sample_id", "forward_sequence_count", "reverse_sequence_count"],
        )
        write_tsv(
            output_dir / "barcode-correction-audit.tsv",
            correction_rows,
            ["mapping_status", "barcode_errors", "records", "percent_total"],
        )
        write_tsv(
            output_dir / "fastqc-module-audit.tsv",
            fastqc_rows,
            ["role", "filename", "module", "status"],
        )
        write_tsv(
            output_dir / "multiqc-general-stats.tsv",
            multiqc_general,
            [
                "sample",
                "total_sequences",
                "percent_gc",
                "percent_duplicates",
                "average_sequence_length",
                "median_sequence_length",
                "percent_failed_modules",
            ],
        )

        os.environ["MPLBACKEND"] = "Agg"
        os.environ["MPLCONFIGDIR"] = str(cache_root / "matplotlib")
        draw_provenance_map(
            figure_dir, assigned, len(count_rows), len(final_actions)
        )
        draw_fastqc_heatmap(figure_dir, fastqc_rows)
        draw_demultiplexing_qc(
            figure_dir, count_rows, forward_quality, reverse_quality
        )
        figure_paths = [
            figure_dir / f"{stem}.{extension}"
            for stem in (
                "10-import-provenance-map",
                "10-fastqc-module-audit",
                "10-demultiplexing-qc",
            )
            for extension in ("pdf", "png", "tiff")
        ]
        check(
            "publication-figure-contract",
            sum(path.is_file() and path.stat().st_size > 0 for path in figure_paths),
            9,
            all(path.is_file() and path.stat().st_size > 0 for path in figure_paths),
        )

        log_lines = [
            "Article 10 import, provenance, and raw-QC validation",
            f"Validation date: {date.today().isoformat()}",
            "",
        ]
        for result in command_results:
            log_lines.extend(
                [
                    f"[{result.label}]",
                    f"command: {shlex.join(result.command)}",
                    f"returncode: {result.returncode}",
                    f"elapsed_seconds: {result.elapsed_seconds}",
                    "stdout:",
                    result.stdout or "(empty)",
                    "stderr:",
                    result.stderr or "(empty)",
                    "",
                ]
            )
        log_text = sanitize_text("\n".join(log_lines), replacements)
        (output_dir / "validation.log").write_text(log_text, encoding="utf-8")

        public_paths = [
            *final_artifacts.values(),
            *(output_dir / "fastqc").glob("*"),
            output_dir / "multiqc_report.html",
            output_dir / "import-provenance-audit.tsv",
            output_dir / "demux-sample-counts.tsv",
            output_dir / "barcode-correction-audit.tsv",
            output_dir / "fastqc-module-audit.tsv",
            output_dir / "multiqc-general-stats.tsv",
            output_dir / "validation.log",
            *figure_paths,
        ]
        leaks = scan_for_paths(
            public_paths,
            [str(root), str(work), str(Path.home())],
        )
        check(
            "public-path-redaction",
            leaks,
            [],
            not leaks,
        )

    if len(checks) != 45:
        raise RuntimeError(f"validator defined {len(checks)} checks; expected 45")
    passed = sum(row["status"] == "PASS" for row in checks)
    failed = len(checks) - passed
    summary = {
        "status": "passed" if failed == 0 else "failed",
        "validation_date": date.today().isoformat(),
        "qiime_release": EXPECTED_QIIME_RELEASE,
        "qiime_version": EXPECTED_QIIME_VERSION,
        "q2cli_version": EXPECTED_Q2CLI_VERSION,
        "demux_plugin_version": EXPECTED_DEMUX_VERSION,
        "fastqc_version": EXPECTED_FASTQC_VERSION,
        "multiqc_version": EXPECTED_MULTIQC_VERSION,
        "raw_qc_environment_sha256": sha256(root / "env/raw-qc.yml"),
        "input_read_pairs": EXPECTED_RECORDS,
        "metadata_samples": EXPECTED_METADATA_SAMPLES,
        "assigned_read_pairs": EXPECTED_ASSIGNED,
        "unmapped_read_pairs": EXPECTED_RECORDS - EXPECTED_ASSIGNED,
        "assignment_rate": round(EXPECTED_ASSIGNED / EXPECTED_RECORDS, 4),
        "demultiplexed_samples": EXPECTED_DEMULTIPLEXED_SAMPLES,
        "fastqc_reports": 3,
        "fastqc_module_states": 33,
        "multiqc_reports": 3,
        "provenance_actions": 3,
        "checks_total": len(checks),
        "checks_passed": passed,
        "checks_failed": failed,
        "figures": [
            "10-import-provenance-map",
            "10-fastqc-module-audit",
            "10-demultiplexing-qc",
        ],
    }
    (output_dir / "qc-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    if failed:
        failed_ids = [row["check_id"] for row in checks if row["status"] == "FAIL"]
        raise RuntimeError(f"Article 10 validation failed: {', '.join(failed_ids)}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
