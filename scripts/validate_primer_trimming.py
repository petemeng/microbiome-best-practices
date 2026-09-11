#!/usr/bin/env python3
"""One-time primer-presence and q2-cutadapt audit for Article 11."""

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
import statistics
import subprocess
import tempfile
import textwrap
import time
import zipfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Iterator, Sequence


EXPECTED_QIIME_ENV = "microbiome-qiime2-2026.4"
EXPECTED_QIIME_RELEASE = "2026.4"
EXPECTED_QIIME_VERSION = "2026.4.0"
EXPECTED_Q2CLI_VERSION = "2026.4.0"
EXPECTED_Q2_CUTADAPT_VERSION = "2026.4.0"
EXPECTED_CUTADAPT_VERSION = "5.2"
EXPECTED_QIIME_ENV_SHA256 = (
    "68b0e9b974937525b71cca3dec48b7dc5b6135f588a941d1313c118c2a41f4b5"
)
EXPECTED_POSITIVE_SUMMARY_SHA256 = (
    "1efaf83d7670a8c9a08cf55670e9536b593e54897960fe5c83e3f60d56856be7"
)
EXPECTED_ATACAMA_SUMMARY_SHA256 = (
    "daa3e140b0bf7e91c94392a48085da8e1e154920a0272eef93c15a26956e9983"
)
EXPECTED_SOURCE_COMMIT = "f282ad2bafb3ca90190ec87290066ee1985df655"
EXPECTED_FORWARD_PRIMER = "GTGYCAGCMGCCGCGGTAA"
EXPECTED_REVERSE_PRIMER = "GGACTACNVGGGTWTCTAAT"
EXPECTED_ATACAMA_FORWARD = "GTGCCAGCMGCCGCGGTAA"
EXPECTED_ATACAMA_REVERSE = "GGACTACHVGGGTWTCTAAT"
WRONG_FORWARD_PRIMER = "CCTACGGGNGGCWGCAG"
WRONG_REVERSE_PRIMER = "GACTACHVGGGTATCTAATCC"
EXPECTED_PER_SAMPLE = {"1": 1810, "1a": 2219, "2": 2347, "2a": 2238}
EXPECTED_INPUT_PAIRS = 10000
EXPECTED_RETAINED_PAIRS = 8614
EXPECTED_UNANCHORED_PAIRS = 8627
EXPECTED_ANCHOR_DELTA = 13
EXPECTED_ANCHORED_R1_MATCHES = 8841
EXPECTED_ANCHORED_R2_MATCHES = 8708
EXPECTED_ATACAMA_PAIRS = 2000
EXPECTED_CHECKS = 50


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
    timeout: int = 900,
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


def normalized_read_id(header: str) -> str:
    identifier = header.split(maxsplit=1)[0].removeprefix("@")
    # SRA-derived pairs commonly encode the mate as a final ".1"/".2";
    # Illumina exports also use "/1"/"/2".  Remove only that terminal token.
    return re.sub(r"(?:/[12]|\.[12])$", "", identifier)


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
            if any(ord(character) < 33 or ord(character) > 126 for character in record[3]):
                raise ValueError(f"{path.name}: invalid FASTQ quality at record {index}")
            yield record


def summarize_fastq(path: Path) -> dict[str, object]:
    identifiers: list[str] = []
    lengths: list[int] = []
    for header, sequence, _, _ in fastq_records(path):
        identifiers.append(normalized_read_id(header))
        lengths.append(len(sequence))
    if not lengths:
        raise ValueError(f"{path.name}: no records")
    return {
        "records": len(lengths),
        "identifiers": identifiers,
        "lengths": lengths,
        "read_length_min": min(lengths),
        "read_length_median": float(statistics.median(lengths)),
        "read_length_max": max(lengths),
    }


def extract_zip_text(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"{path.name}: expected one {suffix}, found {len(matches)}")
        return archive.read(matches[0]).decode("utf-8", errors="replace")


def parse_qzv_counts(path: Path) -> list[dict[str, object]]:
    text = extract_zip_text(path, "/data/per-sample-fastq-counts.tsv")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    parsed: list[dict[str, object]] = []
    for row in rows:
        parsed.append(
            {
                "sample_id": row["sample ID"],
                "forward_sequence_count": int(row["forward sequence count"]),
                "reverse_sequence_count": int(row["reverse sequence count"]),
            }
        )
    return parsed


def provenance_actions(path: Path) -> list[dict[str, str]]:
    actions: list[dict[str, str]] = []
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if "/provenance/" in name and name.endswith("/action/action.yaml")
        ]
        for member in members:
            text = archive.read(member).decode("utf-8", errors="replace")
            action_block = re.search(
                r"^action:\n(?P<block>.*?)(?=^[a-z][a-z_-]*:\n|\Z)",
                text,
                flags=re.MULTILINE | re.DOTALL,
            )
            block = action_block.group("block") if action_block else ""
            kind = re.search(r"^    type:\s*(.+)$", block, flags=re.MULTILINE)
            plugin = re.search(r"^    plugin:\s*(.+)$", block, flags=re.MULTILINE)
            action = re.search(r"^    action:\s*(.+)$", block, flags=re.MULTILINE)
            uuid = re.search(r"^    uuid:\s*(.+)$", text, flags=re.MULTILINE)
            actions.append(
                {
                    "member": member,
                    "execution_uuid": uuid.group(1).strip() if uuid else "",
                    "kind": kind.group(1).strip() if kind else "unknown",
                    "plugin": plugin.group(1).strip() if plugin else "",
                    "action": action.group(1).strip() if action else "import",
                    "text": text,
                }
            )
    return actions


def run_cutadapt_probe(
    *,
    label: str,
    environment_name: str,
    environment: dict[str, str],
    input_r1: Path,
    input_r2: Path,
    output_root: Path,
    front_f: str,
    front_r: str,
    anchored: bool,
) -> tuple[CommandResult, dict[str, object]]:
    output_root.parent.mkdir(parents=True, exist_ok=True)
    adapter_f = f"^{front_f}" if anchored else front_f
    adapter_r = f"^{front_r}" if anchored else front_r
    json_path = output_root.with_suffix(".json")
    output_r1 = output_root.with_suffix(".R1.fastq.gz")
    output_r2 = output_root.with_suffix(".R2.fastq.gz")
    command = conda_command(
        environment_name,
        "cutadapt",
        "--cores",
        "1",
        "-e",
        "0.1",
        "-O",
        "15",
        "-g",
        adapter_f,
        "-G",
        adapter_r,
        "--discard-untrimmed",
        "--json",
        str(json_path),
        "-o",
        str(output_r1),
        "-p",
        str(output_r2),
        str(input_r1),
        str(input_r2),
    )
    result = run_command(label, command, environment)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed: {result.stderr}")
    report = json.loads(json_path.read_text(encoding="utf-8"))
    return result, report


def report_counts(report: dict[str, object]) -> dict[str, int]:
    counts = report["read_counts"]
    if not isinstance(counts, dict):
        raise ValueError("Cutadapt report has no read_counts object")
    return {
        "input": int(counts["input"]),
        "output": int(counts["output"]),
        "read1_with_adapter": int(counts["read1_with_adapter"]),
        "read2_with_adapter": int(counts["read2_with_adapter"]),
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


def draw_primer_presence_gate(
    figure_dir: Path,
    positive_r1: int,
    positive_r2: int,
    positive_paired: int,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    dark = "#111827"
    gray = "#4B5563"
    blue = "#0072B2"
    green = "#009E73"
    orange = "#E69F00"
    red = "#D55E00"

    fig = plt.figure(figsize=(12.2, 6.6), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.35], wspace=0.22)
    ax = fig.add_subplot(grid[0, 0])
    labels = [
        "Primer-retained V4 · R1",
        "Primer-retained V4 · R2",
        "Primer-retained V4 · pair",
        "Atacama V4 · pair",
        "Wrong V3–V4 pair · pair",
    ]
    values = [
        100 * positive_r1 / EXPECTED_INPUT_PAIRS,
        100 * positive_r2 / EXPECTED_INPUT_PAIRS,
        100 * positive_paired / EXPECTED_INPUT_PAIRS,
        0.0,
        0.0,
    ]
    colors = [blue, blue, green, orange, red]
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.3)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("5′ primer match / retained pairs (%)")
    ax.set_title("A  Sequence evidence", loc="left", fontsize=13, fontweight="bold")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    for index, value in enumerate(values):
        ax.text(
            value + 1.3,
            index,
            f"{value:.2f}%" if value else "0%",
            va="center",
            fontsize=9.0,
            color=dark,
            fontweight="bold",
        )
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)

    flow = fig.add_subplot(grid[0, 1])
    flow.set_xlim(0, 10)
    flow.set_ylim(0, 10)
    flow.axis("off")
    flow.set_title("B  Primer-presence gate", loc="left", fontsize=13, fontweight="bold")
    nodes = [
        (0.5, 7.8, 2.6, 1.05, "Document", "Exact biological primers", blue),
        (3.7, 7.8, 2.6, 1.05, "Probe", "R1 and R2 5′ matches", orange),
        (6.9, 7.8, 2.6, 1.05, "Compare", "Expected location + length", green),
        (1.2, 4.7, 3.2, 1.25, "Coherent 5′ peak", "Trim both ends\nthen re-summarize", green),
        (5.6, 4.7, 3.2, 1.25, "No 5′ peak", "Skip trimming\nverify library protocol", orange),
        (3.4, 1.55, 3.2, 1.25, "Ambiguous pattern", "Stop and resolve\norientation / spacer / primer", red),
    ]
    for x, y0, width, height, title, detail, color in nodes:
        flow.add_patch(
            FancyBboxPatch(
                (x, y0),
                width,
                height,
                boxstyle="round,pad=0.04,rounding_size=0.12",
                linewidth=1.25,
                edgecolor=color,
                facecolor=color + "12",
            )
        )
        flow.text(
            x + width / 2,
            y0 + height * 0.68,
            title,
            ha="center",
            va="center",
            fontsize=10.2,
            fontweight="bold",
            color=dark,
        )
        flow.text(
            x + width / 2,
            y0 + height * 0.30,
            detail,
            ha="center",
            va="center",
            fontsize=8.7,
            color=gray,
            linespacing=1.15,
        )
    arrows = [
        ((3.12, 8.33), (3.65, 8.33), blue),
        ((6.32, 8.33), (6.85, 8.33), orange),
        ((7.65, 7.75), (3.7, 6.0), green),
        ((8.7, 7.75), (7.7, 6.0), orange),
        ((5.0, 7.75), (5.0, 2.85), red),
    ]
    for start, end, color in arrows:
        flow.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.2,
                color=color,
                connectionstyle="arc3,rad=0.05",
            )
        )
    fig.suptitle(
        "Primer identity is measured in the reads—not inferred from the region label",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color=dark,
    )
    fig.text(
        0.055,
        0.925,
        "The same V4 label leads to opposite actions in two public datasets.",
        fontsize=10.2,
        color=gray,
    )
    fig.subplots_adjust(left=0.08, right=0.985, top=0.84, bottom=0.12)
    save_figure(fig, figure_dir / "11-primer-presence-gate")
    plt.close(fig)


def draw_trimming_retention_audit(
    figure_dir: Path,
    retention_rows: list[dict[str, object]],
    length_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    dark = "#111827"
    gray = "#4B5563"
    blue = "#0072B2"
    green = "#009E73"
    orange = "#E69F00"

    samples = [str(row["sample_id"]) for row in retention_rows]
    retained = [int(row["retained_pairs"]) for row in retention_rows]
    retention = [float(row["retention_percent"]) for row in retention_rows]
    length_lookup = {
        (str(row["sample_id"]), str(row["direction"])): row for row in length_rows
    }
    r1_delta = [
        float(length_lookup[(sample, "R1")]["median_bases_removed"])
        for sample in samples
    ]
    r2_delta = [
        float(length_lookup[(sample, "R2")]["median_bases_removed"])
        for sample in samples
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.8), facecolor="white")
    x = np.arange(len(samples))
    axes[0].bar(x, retained, color=green, width=0.62)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(samples)
    axes[0].set_xlabel("Sample")
    axes[0].set_ylabel("Retained read pairs")
    axes[0].set_ylim(0, 2700)
    axes[0].set_title("A  Paired retention after anchored trimming", loc="left", fontweight="bold")
    axes[0].grid(axis="y", color="#E5E7EB", linewidth=0.8)
    axes[0].set_axisbelow(True)
    for index, (count, percent) in enumerate(zip(retained, retention)):
        axes[0].text(
            index,
            count + 65,
            f"{count:,}\n{percent:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8.8,
            color=dark,
            fontweight="bold",
        )

    width = 0.34
    axes[1].bar(x - width / 2, r1_delta, width, color=blue, label="R1 · 515F")
    axes[1].bar(x + width / 2, r2_delta, width, color=orange, label="R2 · 806R")
    axes[1].axhline(19, color=blue, linestyle=":", linewidth=1)
    axes[1].axhline(20, color=orange, linestyle=":", linewidth=1)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(samples)
    axes[1].set_xlabel("Sample")
    axes[1].set_ylabel("Median bases removed")
    axes[1].set_ylim(0, 23)
    axes[1].set_title("B  Expected primer-length shift", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, loc="lower right")
    axes[1].grid(axis="y", color="#E5E7EB", linewidth=0.8)
    axes[1].set_axisbelow(True)
    for index, value in enumerate(r1_delta):
        axes[1].text(index - width / 2, value + 0.5, f"{value:.0f}", ha="center", fontsize=8.8)
    for index, value in enumerate(r2_delta):
        axes[1].text(index + width / 2, value + 0.5, f"{value:.0f}", ha="center", fontsize=8.8)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(
        "A valid trim changes both abundance and read geometry",
        x=0.06,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color=dark,
    )
    fig.text(
        0.06,
        0.92,
        "Retention is sample-specific; the 19/20-nt median shift confirms that the intended 5′ primers—not arbitrary internal sequence—were removed.",
        fontsize=9.8,
        color=gray,
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.13, wspace=0.25)
    save_figure(fig, figure_dir / "11-trimming-retention-audit")
    plt.close(fig)


def draw_primer_read_architecture(figure_dir: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    dark = "#111827"
    gray = "#4B5563"
    blue = "#0072B2"
    green = "#009E73"
    orange = "#E69F00"
    red = "#D55E00"

    fig, ax = plt.subplots(figsize=(12.2, 6.2), facecolor="white")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.text(
        0.3,
        7.7,
        "Two V4 libraries can expose different sequence starts",
        fontsize=18,
        fontweight="bold",
        color=dark,
        va="top",
    )
    ax.text(
        0.3,
        7.18,
        "Trim only sequence that is physically present in the FASTQ record.",
        fontsize=10.2,
        color=gray,
        va="top",
    )

    ax.text(0.35, 6.28, "A  Primer-retained paired-end library", fontsize=12.5, fontweight="bold", color=dark)
    ax.add_patch(FancyBboxPatch((0.55, 5.10), 1.75, 0.68, boxstyle="round,pad=0.03", facecolor=blue + "20", edgecolor=blue))
    ax.text(1.425, 5.44, "515F · 19 nt", ha="center", va="center", fontsize=9.5, fontweight="bold", color=dark)
    ax.add_patch(FancyBboxPatch((2.30, 5.10), 6.45, 0.68, boxstyle="round,pad=0.03", facecolor=green + "18", edgecolor=green))
    ax.text(5.525, 5.44, "V4 biological insert", ha="center", va="center", fontsize=10.2, fontweight="bold", color=dark)
    ax.add_patch(FancyBboxPatch((8.75, 5.10), 1.85, 0.68, boxstyle="round,pad=0.03", facecolor=orange + "20", edgecolor=orange))
    ax.text(9.675, 5.44, "806R · 20 nt", ha="center", va="center", fontsize=9.5, fontweight="bold", color=dark)
    ax.add_patch(FancyArrowPatch((0.62, 4.76), (5.0, 4.76), arrowstyle="-|>", mutation_scale=13, color=blue, linewidth=1.6))
    ax.text(2.8, 4.45, "R1 starts with 515F", ha="center", fontsize=9.2, color=blue)
    ax.add_patch(FancyArrowPatch((10.52, 4.15), (6.15, 4.15), arrowstyle="-|>", mutation_scale=13, color=orange, linewidth=1.6))
    ax.text(8.32, 3.84, "R2 starts with 806R", ha="center", fontsize=9.2, color=orange)
    for xpos, color in ((2.30, blue), (8.75, orange)):
        ax.plot([xpos, xpos], [4.9, 6.02], color=red, linewidth=2.1)
        ax.text(xpos, 6.12, "CUT", ha="center", fontsize=8.5, fontweight="bold", color=red)

    ax.text(0.35, 2.95, "B  EMP custom sequencing-primer library", fontsize=12.5, fontweight="bold", color=dark)
    ax.add_patch(FancyBboxPatch((0.55, 1.76), 1.75, 0.68, boxstyle="round,pad=0.03", facecolor="#F3F4F6", edgecolor="#9CA3AF", linestyle="--"))
    ax.text(1.425, 2.10, "PCR primer\nnot sequenced", ha="center", va="center", fontsize=8.8, color=gray)
    ax.add_patch(FancyBboxPatch((2.30, 1.76), 6.45, 0.68, boxstyle="round,pad=0.03", facecolor=green + "18", edgecolor=green))
    ax.text(5.525, 2.10, "V4 biological insert", ha="center", va="center", fontsize=10.2, fontweight="bold", color=dark)
    ax.add_patch(FancyBboxPatch((8.75, 1.76), 1.85, 0.68, boxstyle="round,pad=0.03", facecolor="#F3F4F6", edgecolor="#9CA3AF", linestyle="--"))
    ax.text(9.675, 2.10, "PCR primer\nnot sequenced", ha="center", va="center", fontsize=8.8, color=gray)
    ax.add_patch(FancyArrowPatch((2.35, 1.38), (5.7, 1.38), arrowstyle="-|>", mutation_scale=13, color=green, linewidth=1.6))
    ax.text(4.0, 1.06, "R1 begins inside the amplicon", ha="center", fontsize=9.2, color=green)
    ax.add_patch(FancyArrowPatch((8.70, 0.74), (6.0, 0.74), arrowstyle="-|>", mutation_scale=13, color=green, linewidth=1.6))
    ax.text(7.35, 0.42, "R2 begins inside the amplicon", ha="center", fontsize=9.2, color=green)
    ax.text(
        10.95,
        2.10,
        "SKIP",
        fontsize=11,
        fontweight="bold",
        color=orange,
        ha="center",
        va="center",
    )
    ax.text(
        10.95,
        5.44,
        "TRIM",
        fontsize=11,
        fontweight="bold",
        color=green,
        ha="center",
        va="center",
    )
    fig.subplots_adjust(left=0.015, right=0.99, top=0.99, bottom=0.02)
    save_figure(fig, figure_dir / "11-primer-read-architecture")
    plt.close(fig)


def draw_region_branch_decision(
    figure_dir: Path,
    rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt

    dark = "#111827"
    gray = "#4B5563"
    header = "#1F4E79"
    alternating = ["#F8FAFC", "#EEF6FB"]

    columns = [
        ("assay_branch", "Assay branch", 0.18),
        ("read_layout", "Read layout", 0.16),
        ("primer_evidence", "Primer evidence", 0.22),
        ("trim_action", "Trim action", 0.22),
        ("downstream_implication", "Downstream implication", 0.22),
    ]
    wrap_widths = {
        "assay_branch": 22,
        "read_layout": 20,
        "primer_evidence": 30,
        "trim_action": 29,
        "downstream_implication": 31,
    }
    cell_text = [
        [
            textwrap.fill(
                str(row[key]),
                width=wrap_widths[key],
                break_long_words=False,
                break_on_hyphens=False,
            )
            for key, _, _ in columns
        ]
        for row in rows
    ]
    labels = [
        textwrap.fill(label, width=wrap_widths[key], break_long_words=False)
        for key, label, _ in columns
    ]
    widths = [width for _, _, width in columns]

    fig, ax = plt.subplots(figsize=(14.2, 7.6), facecolor="white")
    ax.axis("off")
    table = ax.table(
        cellText=cell_text,
        colLabels=labels,
        colWidths=widths,
        cellLoc="left",
        colLoc="left",
        loc="center",
        bbox=[0.02, 0.065, 0.96, 0.775],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.2)
    for (row_index, column_index), cell in table.get_celld().items():
        cell.set_linewidth(0.6)
        cell.set_edgecolor("#D1D5DB")
        cell.PAD = 0.08
        cell.get_text().set_verticalalignment("center")
        if row_index == 0:
            cell.set_facecolor(header)
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
            cell.get_text().set_fontsize(8.7)
        else:
            cell.set_facecolor(alternating[(row_index - 1) % 2])
            cell.get_text().set_color(dark)
            if column_index == 0:
                cell.get_text().set_weight("bold")
    fig.text(
        0.035,
        0.955,
        "Region and read layout determine the trimming branch",
        fontsize=18,
        fontweight="bold",
        color=dark,
        ha="left",
        va="top",
    )
    fig.text(
        0.035,
        0.895,
        "Region names guide the search; observed read starts and the library protocol decide the action.",
        fontsize=10.1,
        color=gray,
        ha="left",
        va="top",
    )
    fig.text(
        0.035,
        0.025,
        "Never merge runs from different primer regions before run-specific trimming and denoising.",
        fontsize=9.2,
        fontweight="bold",
        color="#D55E00",
        ha="left",
    )
    save_figure(fig, figure_dir / "11-region-branch-decision")
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
                    (f"{path.name}:{member}", archive.read(member))
                    for member in archive.namelist()
                    if not member.endswith("/")
                )
        for label, payload in payloads:
            if any(needle in payload for needle in needles):
                leaks.append(label)
    return sorted(set(leaks))


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    positive_dir = root / "data/small/primer-trimming"
    atacama_dir = root / "data/small/fastq"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    final_artifacts = {
        "input": output_dir / "nfcore-v4-demux.qza",
        "input_summary": output_dir / "nfcore-v4-demux-summary.qzv",
        "trimmed": output_dir / "nfcore-v4-trimmed.qza",
        "trimmed_summary": output_dir / "nfcore-v4-trimmed-summary.qzv",
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

    with tempfile.TemporaryDirectory(prefix="article11-", dir=output_dir) as temporary:
        work = Path(temporary)
        casava_dir = work / "casava"
        export_dir = work / "trimmed-export"
        probe_dir = work / "cutadapt-probes"
        cache_root = work / "cache"
        home_dir = work / "home"
        for path in (
            casava_dir,
            export_dir,
            probe_dir,
            cache_root / "xdg",
            cache_root / "matplotlib",
            cache_root / "numba",
            cache_root / "fontconfig",
            home_dir / ".cache/fontconfig",
        ):
            path.mkdir(parents=True, exist_ok=True)

        positive_summary_path = positive_dir / "source_summary.json"
        positive_source = json.loads(positive_summary_path.read_text(encoding="utf-8"))
        atacama_summary_path = atacama_dir / "source_summary.json"
        atacama_source = json.loads(atacama_summary_path.read_text(encoding="utf-8"))
        positive_files = sorted(positive_dir.glob("*.fastq.gz"))
        for path in positive_files:
            shutil.copy2(path, casava_dir / path.name)

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
            key: work / path.name for key, path in final_artifacts.items()
        }

        info_result = run_command(
            "qiime-info",
            qiime_command(args.qiime_env, "info"),
            environment,
        )
        command_results.append(info_result)
        if info_result.returncode != 0:
            raise RuntimeError(f"qiime info failed: {info_result.stderr}")
        qiime_info = parse_qiime_info(info_result.stdout)

        cutadapt_version_result = run_command(
            "cutadapt-version",
            conda_command(args.qiime_env, "cutadapt", "--version"),
            environment,
        )
        command_results.append(cutadapt_version_result)
        observed_cutadapt_version = (
            cutadapt_version_result.stdout or cutadapt_version_result.stderr
        ).strip()

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
            "q2-cutadapt-plugin-version",
            qiime_info["plugins"].get("cutadapt"),
            EXPECTED_Q2_CUTADAPT_VERSION,
            qiime_info["plugins"].get("cutadapt") == EXPECTED_Q2_CUTADAPT_VERSION,
        )
        check(
            "cutadapt-cli-version",
            observed_cutadapt_version,
            EXPECTED_CUTADAPT_VERSION,
            cutadapt_version_result.returncode == 0
            and observed_cutadapt_version == EXPECTED_CUTADAPT_VERSION,
        )
        observed_env_hash = sha256(root / "env/qiime2.yml")
        check(
            "qiime-environment-hash",
            observed_env_hash,
            EXPECTED_QIIME_ENV_SHA256,
            observed_env_hash == EXPECTED_QIIME_ENV_SHA256,
        )

        source_contract = (
            sha256(positive_summary_path) == EXPECTED_POSITIVE_SUMMARY_SHA256
            and positive_source["source_commit"] == EXPECTED_SOURCE_COMMIT
            and positive_source["assay"]["forward_primer_5to3"]
            == EXPECTED_FORWARD_PRIMER
            and positive_source["assay"]["reverse_primer_5to3"]
            == EXPECTED_REVERSE_PRIMER
        )
        check(
            "positive-source-contract",
            {
                "summary_sha256": sha256(positive_summary_path),
                "source_commit": positive_source["source_commit"],
                "region": positive_source["assay"]["region"],
            },
            {
                "summary_sha256": EXPECTED_POSITIVE_SUMMARY_SHA256,
                "source_commit": EXPECTED_SOURCE_COMMIT,
                "region": "V4",
            },
            source_contract,
        )

        input_summaries: dict[str, dict[str, object]] = {}
        expected_files = positive_source["files"]
        for filename in sorted(expected_files):
            path = positive_dir / filename
            input_summaries[filename] = summarize_fastq(path)
            observed_hash = sha256(path)
            check(
                f"positive-sha256-{filename}",
                observed_hash,
                expected_files[filename]["sha256"],
                observed_hash == expected_files[filename]["sha256"],
            )

        byte_sizes_valid = all(
            (positive_dir / filename).stat().st_size == details["bytes"]
            for filename, details in expected_files.items()
        )
        check(
            "positive-byte-sizes",
            {
                filename: (positive_dir / filename).stat().st_size
                for filename in sorted(expected_files)
            },
            {
                filename: details["bytes"]
                for filename, details in sorted(expected_files.items())
            },
            byte_sizes_valid,
        )
        record_counts_valid = all(
            input_summaries[filename]["records"] == details["records"] == 2500
            for filename, details in expected_files.items()
        )
        check(
            "positive-record-counts",
            {
                filename: input_summaries[filename]["records"]
                for filename in sorted(input_summaries)
            },
            "2,500 records in each of 8 FASTQ files",
            record_counts_valid,
        )

        sample_pairs: dict[str, tuple[Path, Path]] = {}
        for sample, source_code in (
            ("1", "1_S103"),
            ("1a", "1a_S103"),
            ("2", "2_S115"),
            ("2a", "2a_S115"),
        ):
            sample_pairs[sample] = (
                positive_dir / f"{source_code}_L001_R1_001.fastq.gz",
                positive_dir / f"{source_code}_L001_R2_001.fastq.gz",
            )
        pair_sync = all(
            input_summaries[r1.name]["identifiers"]
            == input_summaries[r2.name]["identifiers"]
            for r1, r2 in sample_pairs.values()
        )
        check("positive-pair-synchronization", pair_sync, True, pair_sync)
        observed_input_pairs = sum(
            int(input_summaries[r1.name]["records"])
            for r1, _ in sample_pairs.values()
        )
        check(
            "positive-input-contract",
            {"samples": len(sample_pairs), "read_pairs": observed_input_pairs},
            {"samples": 4, "read_pairs": EXPECTED_INPUT_PAIRS},
            len(sample_pairs) == 4 and observed_input_pairs == EXPECTED_INPUT_PAIRS,
        )
        input_fastq_valid = all(
            int(summary["read_length_min"]) > 0
            and int(summary["read_length_max"]) <= 250
            for summary in input_summaries.values()
        )
        check(
            "positive-fastq-structure",
            {
                filename: [
                    summary["read_length_min"],
                    summary["read_length_max"],
                ]
                for filename, summary in sorted(input_summaries.items())
            },
            "valid Phred+33 FASTQ records with read lengths 1–250 nt",
            input_fastq_valid,
        )

        atacama_source_contract = (
            sha256(atacama_summary_path) == EXPECTED_ATACAMA_SUMMARY_SHA256
            and atacama_source["selection"]["selected_records"] == EXPECTED_ATACAMA_PAIRS
        )
        check(
            "atacama-source-contract",
            {
                "summary_sha256": sha256(atacama_summary_path),
                "selected_records": atacama_source["selection"]["selected_records"],
            },
            {
                "summary_sha256": EXPECTED_ATACAMA_SUMMARY_SHA256,
                "selected_records": EXPECTED_ATACAMA_PAIRS,
            },
            atacama_source_contract,
        )
        atacama_hashes_valid = all(
            sha256(atacama_dir / filename)
            == atacama_source["output_files"][filename]["sha256"]
            for filename in ("forward.fastq.gz", "reverse.fastq.gz")
        )
        check(
            "atacama-fastq-hashes",
            {
                filename: sha256(atacama_dir / filename)
                for filename in ("forward.fastq.gz", "reverse.fastq.gz")
            },
            {
                filename: atacama_source["output_files"][filename]["sha256"]
                for filename in ("forward.fastq.gz", "reverse.fastq.gz")
            },
            atacama_hashes_valid,
        )
        atacama_forward_summary = summarize_fastq(atacama_dir / "forward.fastq.gz")
        atacama_reverse_summary = summarize_fastq(atacama_dir / "reverse.fastq.gz")
        atacama_sync = (
            atacama_forward_summary["identifiers"]
            == atacama_reverse_summary["identifiers"]
        )
        check(
            "atacama-fastq-contract",
            {
                "forward_records": atacama_forward_summary["records"],
                "reverse_records": atacama_reverse_summary["records"],
                "synchronized": atacama_sync,
            },
            {
                "forward_records": EXPECTED_ATACAMA_PAIRS,
                "reverse_records": EXPECTED_ATACAMA_PAIRS,
                "synchronized": True,
            },
            atacama_forward_summary["records"] == EXPECTED_ATACAMA_PAIRS
            and atacama_reverse_summary["records"] == EXPECTED_ATACAMA_PAIRS
            and atacama_sync,
        )

        manifest_rows = []
        for sample, (r1, r2) in sample_pairs.items():
            manifest_rows.append(
                {
                    "sample-id": sample,
                    "forward-absolute-filepath": (
                        f"file://<PROJECT_ROOT>/data/small/primer-trimming/{r1.name}"
                    ),
                    "reverse-absolute-filepath": (
                        f"file://<PROJECT_ROOT>/data/small/primer-trimming/{r2.name}"
                    ),
                }
            )
        write_tsv(
            output_dir / "nfcore-v4-manifest.tsv",
            manifest_rows,
            [
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            ],
        )
        portable_manifest_valid = (
            len(manifest_rows) == 4
            and {str(row["sample-id"]) for row in manifest_rows}
            == {"1", "1a", "2", "2a"}
            and all(
                "<PROJECT_ROOT>" in str(row["forward-absolute-filepath"])
                and "<PROJECT_ROOT>" in str(row["reverse-absolute-filepath"])
                and str(root) not in str(row)
                for row in manifest_rows
            )
        )
        check(
            "portable-manifest-contract",
            {"rows": len(manifest_rows), "placeholder": "<PROJECT_ROOT>"},
            {"rows": 4, "placeholder": "<PROJECT_ROOT>"},
            portable_manifest_valid,
        )

        import_result = run_command(
            "qiime-import-casava",
            qiime_command(
                args.qiime_env,
                "tools",
                "import",
                "--type",
                "SampleData[PairedEndSequencesWithQuality]",
                "--input-path",
                str(casava_dir),
                "--input-format",
                "CasavaOneEightSingleLanePerSampleDirFmt",
                "--output-path",
                str(work_artifacts["input"]),
            ),
            environment,
        )
        command_results.append(import_result)
        if import_result.returncode != 0:
            raise RuntimeError(f"QIIME 2 import failed: {import_result.stderr}")
        check("qiime-import-exit", import_result.returncode, 0, import_result.returncode == 0)

        input_validate = run_command(
            "input-maximum-validation",
            qiime_command(
                args.qiime_env,
                "tools",
                "validate",
                str(work_artifacts["input"]),
                "--level",
                "max",
            ),
            environment,
        )
        input_peek_result = run_command(
            "input-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(work_artifacts["input"]),
            ),
            environment,
        )
        command_results.extend([input_validate, input_peek_result])
        if input_validate.returncode != 0 or input_peek_result.returncode != 0:
            raise RuntimeError("input Artifact validation or peek failed")
        input_peek = parse_peek_tsv(input_peek_result.stdout)
        input_type_valid = (
            input_peek.get("Type") == "SampleData[PairedEndSequencesWithQuality]"
            and input_peek.get("Data Format")
            == "SingleLanePerSamplePairedEndFastqDirFmt"
        )
        check(
            "input-artifact-type",
            {
                "type": input_peek.get("Type"),
                "format": input_peek.get("Data Format"),
            },
            {
                "type": "SampleData[PairedEndSequencesWithQuality]",
                "format": "SingleLanePerSamplePairedEndFastqDirFmt",
            },
            input_type_valid,
        )
        check(
            "input-maximum-validation",
            input_validate.returncode,
            0,
            input_validate.returncode == 0,
        )

        input_summary_result = run_command(
            "input-demux-summary",
            qiime_command(
                args.qiime_env,
                "demux",
                "summarize",
                "--i-data",
                str(work_artifacts["input"]),
                "--p-n",
                "10000",
                "--o-visualization",
                str(work_artifacts["input_summary"]),
            ),
            environment,
        )
        command_results.append(input_summary_result)
        if input_summary_result.returncode != 0:
            raise RuntimeError(f"input demux summarize failed: {input_summary_result.stderr}")
        input_count_rows = parse_qzv_counts(work_artifacts["input_summary"])
        input_summary_valid = (
            len(input_count_rows) == 4
            and sum(int(row["forward_sequence_count"]) for row in input_count_rows)
            == EXPECTED_INPUT_PAIRS
            and all(
                row["forward_sequence_count"] == row["reverse_sequence_count"] == 2500
                for row in input_count_rows
            )
        )
        check(
            "input-summary-counts",
            input_count_rows,
            "4 samples × 2,500 synchronized pairs",
            input_summary_valid,
        )

        trim_result = run_command(
            "q2-cutadapt-trim-paired",
            qiime_command(
                args.qiime_env,
                "cutadapt",
                "trim-paired",
                "--i-demultiplexed-sequences",
                str(work_artifacts["input"]),
                "--p-front-f",
                f"^{EXPECTED_FORWARD_PRIMER}",
                "--p-front-r",
                f"^{EXPECTED_REVERSE_PRIMER}",
                "--p-error-rate",
                "0.1",
                "--p-overlap",
                "15",
                "--p-match-adapter-wildcards",
                "--p-no-match-read-wildcards",
                "--p-indels",
                "--p-discard-untrimmed",
                "--p-minimum-length",
                "1",
                "--p-cores",
                "1",
                "--o-trimmed-sequences",
                str(work_artifacts["trimmed"]),
            ),
            environment,
        )
        command_results.append(trim_result)
        if trim_result.returncode != 0:
            raise RuntimeError(f"q2-cutadapt trim-paired failed: {trim_result.stderr}")
        check("q2-cutadapt-exit", trim_result.returncode, 0, trim_result.returncode == 0)

        trimmed_validate = run_command(
            "trimmed-maximum-validation",
            qiime_command(
                args.qiime_env,
                "tools",
                "validate",
                str(work_artifacts["trimmed"]),
                "--level",
                "max",
            ),
            environment,
        )
        trimmed_peek_result = run_command(
            "trimmed-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(work_artifacts["trimmed"]),
            ),
            environment,
        )
        command_results.extend([trimmed_validate, trimmed_peek_result])
        if trimmed_validate.returncode != 0 or trimmed_peek_result.returncode != 0:
            raise RuntimeError("trimmed Artifact validation or peek failed")
        trimmed_peek = parse_peek_tsv(trimmed_peek_result.stdout)
        trimmed_type_valid = (
            trimmed_peek.get("Type") == "SampleData[PairedEndSequencesWithQuality]"
            and trimmed_peek.get("Data Format")
            == "SingleLanePerSamplePairedEndFastqDirFmt"
        )
        check(
            "trimmed-artifact-type",
            {
                "type": trimmed_peek.get("Type"),
                "format": trimmed_peek.get("Data Format"),
            },
            {
                "type": "SampleData[PairedEndSequencesWithQuality]",
                "format": "SingleLanePerSamplePairedEndFastqDirFmt",
            },
            trimmed_type_valid,
        )
        check(
            "trimmed-maximum-validation",
            trimmed_validate.returncode,
            0,
            trimmed_validate.returncode == 0,
        )

        trimmed_summary_result = run_command(
            "trimmed-demux-summary",
            qiime_command(
                args.qiime_env,
                "demux",
                "summarize",
                "--i-data",
                str(work_artifacts["trimmed"]),
                "--p-n",
                "10000",
                "--o-visualization",
                str(work_artifacts["trimmed_summary"]),
            ),
            environment,
        )
        command_results.append(trimmed_summary_result)
        if trimmed_summary_result.returncode != 0:
            raise RuntimeError(
                f"trimmed demux summarize failed: {trimmed_summary_result.stderr}"
            )
        trimmed_count_rows = parse_qzv_counts(work_artifacts["trimmed_summary"])
        check(
            "trimmed-summary-output",
            {
                "exit": trimmed_summary_result.returncode,
                "samples": len(trimmed_count_rows),
            },
            {"exit": 0, "samples": 4},
            trimmed_summary_result.returncode == 0 and len(trimmed_count_rows) == 4,
        )
        retained_total = sum(
            int(row["forward_sequence_count"]) for row in trimmed_count_rows
        )
        check(
            "trimmed-retained-total",
            retained_total,
            EXPECTED_RETAINED_PAIRS,
            retained_total == EXPECTED_RETAINED_PAIRS,
        )
        observed_per_sample = {
            str(row["sample_id"]): int(row["forward_sequence_count"])
            for row in trimmed_count_rows
        }
        check(
            "trimmed-per-sample-counts",
            observed_per_sample,
            EXPECTED_PER_SAMPLE,
            observed_per_sample == EXPECTED_PER_SAMPLE,
        )
        paired_counts_match = all(
            row["forward_sequence_count"] == row["reverse_sequence_count"]
            for row in trimmed_count_rows
        )
        check(
            "trimmed-pair-count-identity",
            paired_counts_match,
            True,
            paired_counts_match,
        )

        export_result = run_command(
            "export-trimmed-artifact",
            qiime_command(
                args.qiime_env,
                "tools",
                "export",
                "--input-path",
                str(work_artifacts["trimmed"]),
                "--output-path",
                str(export_dir),
            ),
            environment,
        )
        command_results.append(export_result)
        if export_result.returncode != 0:
            raise RuntimeError(f"trimmed Artifact export failed: {export_result.stderr}")
        exported_fastqs = sorted(export_dir.glob("*.fastq.gz"))
        exported_summaries = {
            path.name: summarize_fastq(path) for path in exported_fastqs
        }
        exported_sync = all(
            exported_summaries[r1.name]["identifiers"]
            == exported_summaries[r2.name]["identifiers"]
            for r1, r2 in sample_pairs.values()
        )
        check(
            "trimmed-export-contract",
            {
                "fastq_files": len(exported_fastqs),
                "paired_ids_synchronized": exported_sync,
            },
            {"fastq_files": 8, "paired_ids_synchronized": True},
            len(exported_fastqs) == 8 and exported_sync,
        )

        retention_rows: list[dict[str, object]] = []
        length_rows: list[dict[str, object]] = []
        for sample, (input_r1, input_r2) in sample_pairs.items():
            retained = observed_per_sample[sample]
            retention_rows.append(
                {
                    "sample_id": sample,
                    "input_pairs": 2500,
                    "retained_pairs": retained,
                    "discarded_pairs": 2500 - retained,
                    "retention_percent": round(100 * retained / 2500, 2),
                }
            )
            for direction, input_path in (("R1", input_r1), ("R2", input_r2)):
                output_path = export_dir / input_path.name
                before = input_summaries[input_path.name]
                after = exported_summaries[output_path.name]
                median_removed = float(before["read_length_median"]) - float(
                    after["read_length_median"]
                )
                length_rows.append(
                    {
                        "sample_id": sample,
                        "direction": direction,
                        "input_records": before["records"],
                        "retained_records": after["records"],
                        "median_length_before": before["read_length_median"],
                        "median_length_after": after["read_length_median"],
                        "median_bases_removed": median_removed,
                    }
                )
        r1_shift_valid = all(
            float(row["median_bases_removed"]) == 19.0
            for row in length_rows
            if row["direction"] == "R1"
        )
        r2_shift_valid = all(
            float(row["median_bases_removed"]) == 20.0
            for row in length_rows
            if row["direction"] == "R2"
        )
        check(
            "r1-median-length-shift",
            [
                row["median_bases_removed"]
                for row in length_rows
                if row["direction"] == "R1"
            ],
            [19.0, 19.0, 19.0, 19.0],
            r1_shift_valid,
        )
        check(
            "r2-median-length-shift",
            [
                row["median_bases_removed"]
                for row in length_rows
                if row["direction"] == "R2"
            ],
            [20.0, 20.0, 20.0, 20.0],
            r2_shift_valid,
        )

        probe_reports: dict[tuple[str, str], dict[str, object]] = {}
        probe_results: list[CommandResult] = []
        for mode, front_f, front_r, anchored in (
            (
                "correct_anchored",
                EXPECTED_FORWARD_PRIMER,
                EXPECTED_REVERSE_PRIMER,
                True,
            ),
            (
                "correct_unanchored",
                EXPECTED_FORWARD_PRIMER,
                EXPECTED_REVERSE_PRIMER,
                False,
            ),
            (
                "wrong_v3v4_anchored",
                WRONG_FORWARD_PRIMER,
                WRONG_REVERSE_PRIMER,
                True,
            ),
        ):
            for sample, (input_r1, input_r2) in sample_pairs.items():
                result, report = run_cutadapt_probe(
                    label=f"cutadapt-{mode}-{sample}",
                    environment_name=args.qiime_env,
                    environment=environment,
                    input_r1=input_r1,
                    input_r2=input_r2,
                    output_root=probe_dir / f"{mode}-{sample}",
                    front_f=front_f,
                    front_r=front_r,
                    anchored=anchored,
                )
                probe_results.append(result)
                probe_reports[(mode, sample)] = report
        command_results.extend(probe_results)

        anchored_counts = [
            report_counts(probe_reports[("correct_anchored", sample)])
            for sample in sample_pairs
        ]
        anchored_r1 = sum(row["read1_with_adapter"] for row in anchored_counts)
        anchored_r2 = sum(row["read2_with_adapter"] for row in anchored_counts)
        anchored_output = sum(row["output"] for row in anchored_counts)
        check(
            "anchored-r1-primer-matches",
            anchored_r1,
            EXPECTED_ANCHORED_R1_MATCHES,
            anchored_r1 == EXPECTED_ANCHORED_R1_MATCHES,
        )
        check(
            "anchored-r2-primer-matches",
            anchored_r2,
            EXPECTED_ANCHORED_R2_MATCHES,
            anchored_r2 == EXPECTED_ANCHORED_R2_MATCHES,
        )
        check(
            "anchored-cli-q2-parity",
            {"cutadapt_cli": anchored_output, "q2_cutadapt": retained_total},
            {
                "cutadapt_cli": EXPECTED_RETAINED_PAIRS,
                "q2_cutadapt": EXPECTED_RETAINED_PAIRS,
            },
            anchored_output == retained_total == EXPECTED_RETAINED_PAIRS,
        )
        unanchored_output = sum(
            report_counts(probe_reports[("correct_unanchored", sample)])["output"]
            for sample in sample_pairs
        )
        check(
            "unanchored-retained-pairs",
            unanchored_output,
            EXPECTED_UNANCHORED_PAIRS,
            unanchored_output == EXPECTED_UNANCHORED_PAIRS,
        )
        anchor_delta = unanchored_output - anchored_output
        check(
            "anchor-delta-pairs",
            anchor_delta,
            EXPECTED_ANCHOR_DELTA,
            anchor_delta == EXPECTED_ANCHOR_DELTA,
        )
        wrong_output = sum(
            report_counts(probe_reports[("wrong_v3v4_anchored", sample)])["output"]
            for sample in sample_pairs
        )
        check(
            "wrong-region-retained-pairs",
            wrong_output,
            0,
            wrong_output == 0,
        )

        atacama_result, atacama_report = run_cutadapt_probe(
            label="cutadapt-atacama-original-v4-anchored",
            environment_name=args.qiime_env,
            environment=environment,
            input_r1=atacama_dir / "forward.fastq.gz",
            input_r2=atacama_dir / "reverse.fastq.gz",
            output_root=probe_dir / "atacama-original-v4-anchored",
            front_f=EXPECTED_ATACAMA_FORWARD,
            front_r=EXPECTED_ATACAMA_REVERSE,
            anchored=True,
        )
        command_results.append(atacama_result)
        atacama_counts = report_counts(atacama_report)
        check(
            "atacama-no-primer-result",
            atacama_counts,
            {
                "input": EXPECTED_ATACAMA_PAIRS,
                "output": 0,
                "read1_with_adapter": 0,
                "read2_with_adapter": 0,
            },
            atacama_counts
            == {
                "input": EXPECTED_ATACAMA_PAIRS,
                "output": 0,
                "read1_with_adapter": 0,
                "read2_with_adapter": 0,
            },
        )

        actions = provenance_actions(work_artifacts["trimmed"])
        check(
            "trimmed-provenance-actions",
            [row["action"] for row in actions],
            ["import", "trim_paired"],
            len(actions) == 2
            and sorted(row["action"] for row in actions) == ["import", "trim_paired"],
        )
        trim_action = next(
            (row for row in actions if row["action"] == "trim_paired"),
            None,
        )
        trim_text = trim_action["text"] if trim_action else ""
        required_parameters = [
            f"- ^{EXPECTED_FORWARD_PRIMER}",
            f"- ^{EXPECTED_REVERSE_PRIMER}",
            "error_rate: 0.1",
            "indels: true",
            "overlap: 15",
            "match_read_wildcards: false",
            "match_adapter_wildcards: true",
            "minimum_length: 1",
            "discard_untrimmed: true",
            "cores: 1",
        ]
        provenance_parameters_valid = all(
            parameter in trim_text for parameter in required_parameters
        )
        check(
            "trimmed-provenance-parameters",
            [parameter for parameter in required_parameters if parameter in trim_text],
            required_parameters,
            provenance_parameters_valid,
        )

        primer_source_rows = [
            {
                "dataset": "nf-core/ampliseq MiSeq v2",
                "region": "V4",
                "primer_version": "515F (Parada) / 806R (Apprill)",
                "forward_primer_5to3": EXPECTED_FORWARD_PRIMER,
                "reverse_primer_5to3": EXPECTED_REVERSE_PRIMER,
                "input_pairs": EXPECTED_INPUT_PAIRS,
                "r1_5prime_matches": anchored_r1,
                "r2_5prime_matches": anchored_r2,
                "paired_retained": anchored_output,
                "decision": "TRIM",
                "reason": "Both read directions show coherent anchored 5-prime primer matches.",
            },
            {
                "dataset": "QIIME 2 Atacama EMP excerpt",
                "region": "V4",
                "primer_version": "Original 515F / 806R",
                "forward_primer_5to3": EXPECTED_ATACAMA_FORWARD,
                "reverse_primer_5to3": EXPECTED_ATACAMA_REVERSE,
                "input_pairs": EXPECTED_ATACAMA_PAIRS,
                "r1_5prime_matches": 0,
                "r2_5prime_matches": 0,
                "paired_retained": 0,
                "decision": "SKIP",
                "reason": "Custom EMP sequencing primers start reads inside the amplicon.",
            },
        ]
        match_rows: list[dict[str, object]] = []
        for sample in sample_pairs:
            counts = report_counts(probe_reports[("correct_anchored", sample)])
            for direction, key in (
                ("R1", "read1_with_adapter"),
                ("R2", "read2_with_adapter"),
            ):
                match_rows.append(
                    {
                        "dataset": "nf-core/ampliseq MiSeq v2",
                        "sample_id": sample,
                        "mode": "Correct V4 primers · anchored",
                        "direction": direction,
                        "input_reads": counts["input"],
                        "matched_reads": counts[key],
                        "match_percent": round(100 * counts[key] / counts["input"], 2),
                        "paired_retained": counts["output"],
                    }
                )
        for direction, key in (
            ("R1", "read1_with_adapter"),
            ("R2", "read2_with_adapter"),
        ):
            match_rows.append(
                {
                    "dataset": "QIIME 2 Atacama EMP excerpt",
                    "sample_id": "multiplexed-excerpt",
                    "mode": "Original V4 primers · anchored",
                    "direction": direction,
                    "input_reads": atacama_counts["input"],
                    "matched_reads": atacama_counts[key],
                    "match_percent": 0.0,
                    "paired_retained": atacama_counts["output"],
                }
            )
        sensitivity_rows = [
            {
                "dataset": "nf-core/ampliseq MiSeq v2",
                "primer_pair": "515F (Parada) / 806R (Apprill)",
                "search_mode": "Anchored at 5-prime",
                "input_pairs": EXPECTED_INPUT_PAIRS,
                "retained_pairs": anchored_output,
                "retention_percent": round(100 * anchored_output / EXPECTED_INPUT_PAIRS, 2),
                "interpretation": "Primary analysis",
            },
            {
                "dataset": "nf-core/ampliseq MiSeq v2",
                "primer_pair": "515F (Parada) / 806R (Apprill)",
                "search_mode": "Unanchored front search",
                "input_pairs": EXPECTED_INPUT_PAIRS,
                "retained_pairs": unanchored_output,
                "retention_percent": round(100 * unanchored_output / EXPECTED_INPUT_PAIRS, 2),
                "interpretation": "Sensitivity only; 13 extra pairs require protocol review",
            },
            {
                "dataset": "nf-core/ampliseq MiSeq v2",
                "primer_pair": "341F / 805R (wrong region)",
                "search_mode": "Anchored at 5-prime",
                "input_pairs": EXPECTED_INPUT_PAIRS,
                "retained_pairs": wrong_output,
                "retention_percent": 0.0,
                "interpretation": "Reject primer hypothesis",
            },
            {
                "dataset": "QIIME 2 Atacama EMP excerpt",
                "primer_pair": "Original 515F / 806R",
                "search_mode": "Anchored at 5-prime",
                "input_pairs": EXPECTED_ATACAMA_PAIRS,
                "retained_pairs": atacama_counts["output"],
                "retention_percent": 0.0,
                "interpretation": "Skip trimming",
            },
        ]
        region_rows = [
            {
                "assay_branch": "V4 · primer retained",
                "read_layout": "Paired short reads",
                "primer_evidence": "Both R1 and R2 start with the documented pair",
                "trim_action": "Trim anchored front-f and front-r; audit paired retention",
                "downstream_implication": "Recalculate quality and overlap; DADA2 trim-left = 0 for primers",
            },
            {
                "assay_branch": "V4 · primer absent",
                "read_layout": "Paired short reads",
                "primer_evidence": "No coherent 5′ match; protocol uses custom sequencing primers",
                "trim_action": "Skip primer trimming",
                "downstream_implication": "Proceed to quality-based truncation; do not discard all reads",
            },
            {
                "assay_branch": "V3–V4",
                "read_layout": "Paired short reads",
                "primer_evidence": "Confirm the exact 341F/805R-family sequences and orientation",
                "trim_action": "Trim run-specific primers before denoising",
                "downstream_implication": "Longer insert makes overlap chemistry- and quality-dependent",
            },
            {
                "assay_branch": "Full-length 16S",
                "read_layout": "Single long read",
                "primer_evidence": "Forward and reverse boundaries occur in one molecule",
                "trim_action": "Use linked/orientation-aware primer logic",
                "downstream_implication": "No paired merge; use platform-appropriate long-read denoising",
            },
            {
                "assay_branch": "Short-read single-end",
                "read_layout": "One acquired direction",
                "primer_evidence": "Probe only the sequenced 5′ end",
                "trim_action": "Trim the documented front primer if present",
                "downstream_implication": "No paired merge and no evidence about the missing mate",
            },
        ]
        provenance_rows = []
        for row in sorted(actions, key=lambda item: item["action"]):
            provenance_rows.append(
                {
                    "execution_uuid": row["execution_uuid"],
                    "action_type": row["kind"],
                    "plugin": row["plugin"] or "framework",
                    "action": row["action"],
                    "parameter_contract": (
                        "PASS"
                        if row["action"] == "import" or provenance_parameters_valid
                        else "FAIL"
                    ),
                }
            )

        write_tsv(
            output_dir / "primer-source-audit.tsv",
            primer_source_rows,
            [
                "dataset",
                "region",
                "primer_version",
                "forward_primer_5to3",
                "reverse_primer_5to3",
                "input_pairs",
                "r1_5prime_matches",
                "r2_5prime_matches",
                "paired_retained",
                "decision",
                "reason",
            ],
        )
        write_tsv(
            output_dir / "primer-match-audit.tsv",
            match_rows,
            [
                "dataset",
                "sample_id",
                "mode",
                "direction",
                "input_reads",
                "matched_reads",
                "match_percent",
                "paired_retained",
            ],
        )
        write_tsv(
            output_dir / "trim-retention-audit.tsv",
            retention_rows,
            [
                "sample_id",
                "input_pairs",
                "retained_pairs",
                "discarded_pairs",
                "retention_percent",
            ],
        )
        write_tsv(
            output_dir / "length-shift-audit.tsv",
            length_rows,
            [
                "sample_id",
                "direction",
                "input_records",
                "retained_records",
                "median_length_before",
                "median_length_after",
                "median_bases_removed",
            ],
        )
        write_tsv(
            output_dir / "parameter-sensitivity.tsv",
            sensitivity_rows,
            [
                "dataset",
                "primer_pair",
                "search_mode",
                "input_pairs",
                "retained_pairs",
                "retention_percent",
                "interpretation",
            ],
        )
        write_tsv(
            output_dir / "region-branch-decision.tsv",
            region_rows,
            [
                "assay_branch",
                "read_layout",
                "primer_evidence",
                "trim_action",
                "downstream_implication",
            ],
        )
        write_tsv(
            output_dir / "primer-provenance-audit.tsv",
            provenance_rows,
            [
                "execution_uuid",
                "action_type",
                "plugin",
                "action",
                "parameter_contract",
            ],
        )

        for key, source in work_artifacts.items():
            shutil.copy2(source, final_artifacts[key])

        replacements = {
            str(root): "<PROJECT_ROOT>",
            str(work): "<WORK_DIR>",
            str(Path.home()): "<HOME>",
        }
        atacama_log = "\n".join(
            [
                "Article 11 Atacama no-primer probe",
                f"command: {shlex.join(atacama_result.command)}",
                f"returncode: {atacama_result.returncode}",
                "stdout:",
                atacama_result.stdout or "(empty)",
                "stderr:",
                atacama_result.stderr or "(empty)",
                "parsed_counts:",
                json.dumps(atacama_counts, ensure_ascii=False, sort_keys=True),
            ]
        )
        (output_dir / "atacama-no-primer-probe.log").write_text(
            sanitize_text(atacama_log, replacements) + "\n",
            encoding="utf-8",
        )
        qiime_cutadapt_log = "\n".join(
            [
                "Article 11 q2-cutadapt primary run",
                f"command: {shlex.join(trim_result.command)}",
                f"returncode: {trim_result.returncode}",
                f"elapsed_seconds: {trim_result.elapsed_seconds}",
                "stdout:",
                trim_result.stdout or "(empty)",
                "stderr:",
                trim_result.stderr or "(empty)",
                f"retained_pairs: {retained_total}",
                f"retention_percent: {100 * retained_total / EXPECTED_INPUT_PAIRS:.2f}",
            ]
        )
        (output_dir / "qiime-cutadapt.log").write_text(
            sanitize_text(qiime_cutadapt_log, replacements) + "\n",
            encoding="utf-8",
        )

        os.environ["MPLBACKEND"] = "Agg"
        os.environ["MPLCONFIGDIR"] = str(cache_root / "matplotlib")
        draw_primer_presence_gate(
            figure_dir,
            anchored_r1,
            anchored_r2,
            anchored_output,
        )
        draw_trimming_retention_audit(figure_dir, retention_rows, length_rows)
        draw_primer_read_architecture(figure_dir)
        draw_region_branch_decision(figure_dir, region_rows)
        figure_paths = [
            figure_dir / f"{stem}.{extension}"
            for stem in (
                "11-primer-presence-gate",
                "11-trimming-retention-audit",
                "11-primer-read-architecture",
                "11-region-branch-decision",
            )
            for extension in ("pdf", "png", "tiff")
        ]

        from PIL import Image

        figure_contract = all(path.is_file() and path.stat().st_size > 0 for path in figure_paths)
        for path in figure_paths:
            if path.suffix not in {".png", ".tiff"}:
                continue
            with Image.open(path) as image:
                dpi = image.info.get("dpi", (0, 0))
                figure_contract = figure_contract and all(
                    abs(float(value) - 350.0) <= 1.0 for value in dpi
                )
                if path.suffix == ".tiff":
                    figure_contract = figure_contract and image.info.get(
                        "compression"
                    ) in {"tiff_lzw", "lzma"}
        check(
            "publication-figure-contract",
            sum(path.is_file() and path.stat().st_size > 0 for path in figure_paths),
            12,
            figure_contract,
        )

        log_lines = [
            "Article 11 primer trimming and region-branch validation",
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
        (output_dir / "validation.log").write_text(
            sanitize_text("\n".join(log_lines), replacements) + "\n",
            encoding="utf-8",
        )

        public_paths = [
            *final_artifacts.values(),
            output_dir / "nfcore-v4-manifest.tsv",
            output_dir / "primer-source-audit.tsv",
            output_dir / "primer-match-audit.tsv",
            output_dir / "trim-retention-audit.tsv",
            output_dir / "length-shift-audit.tsv",
            output_dir / "parameter-sensitivity.tsv",
            output_dir / "region-branch-decision.tsv",
            output_dir / "primer-provenance-audit.tsv",
            output_dir / "atacama-no-primer-probe.log",
            output_dir / "qiime-cutadapt.log",
            output_dir / "validation.log",
            *figure_paths,
        ]
        leaks = scan_for_paths(
            public_paths,
            [str(root), str(work), str(Path.home())],
        )
        check("public-path-redaction", leaks, [], not leaks)

    if len(checks) != EXPECTED_CHECKS:
        raise RuntimeError(
            f"validator defined {len(checks)} checks; expected {EXPECTED_CHECKS}"
        )
    passed = sum(row["status"] == "PASS" for row in checks)
    failed = len(checks) - passed
    summary = {
        "status": "passed" if failed == 0 else "failed",
        "validation_date": date.today().isoformat(),
        "qiime_release": EXPECTED_QIIME_RELEASE,
        "qiime_version": EXPECTED_QIIME_VERSION,
        "q2cli_version": EXPECTED_Q2CLI_VERSION,
        "q2_cutadapt_version": EXPECTED_Q2_CUTADAPT_VERSION,
        "cutadapt_version": EXPECTED_CUTADAPT_VERSION,
        "qiime_environment_sha256": sha256(root / "env/qiime2.yml"),
        "positive_source_commit": EXPECTED_SOURCE_COMMIT,
        "positive_samples": 4,
        "positive_input_pairs": EXPECTED_INPUT_PAIRS,
        "positive_retained_pairs": EXPECTED_RETAINED_PAIRS,
        "positive_retention_percent": round(
            100 * EXPECTED_RETAINED_PAIRS / EXPECTED_INPUT_PAIRS, 2
        ),
        "positive_r1_primer_matches": EXPECTED_ANCHORED_R1_MATCHES,
        "positive_r2_primer_matches": EXPECTED_ANCHORED_R2_MATCHES,
        "unanchored_retained_pairs": EXPECTED_UNANCHORED_PAIRS,
        "anchor_delta_pairs": EXPECTED_ANCHOR_DELTA,
        "wrong_region_retained_pairs": 0,
        "negative_input_pairs": EXPECTED_ATACAMA_PAIRS,
        "negative_retained_pairs": 0,
        "provenance_actions": 2,
        "checks_total": len(checks),
        "checks_passed": passed,
        "checks_failed": failed,
        "figures": [
            "11-primer-presence-gate",
            "11-trimming-retention-audit",
            "11-primer-read-architecture",
            "11-region-branch-decision",
        ],
    }
    (output_dir / "primer-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    if failed:
        failed_ids = [row["check_id"] for row in checks if row["status"] == "FAIL"]
        raise RuntimeError(f"Article 11 validation failed: {', '.join(failed_ids)}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
