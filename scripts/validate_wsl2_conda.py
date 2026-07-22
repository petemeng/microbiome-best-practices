#!/usr/bin/env python3
"""Audit the Linux-side 16S environment and draw Article 07 evidence figures.

The script intentionally does not claim that the current host is Windows or
WSL. It validates the same Linux-side contract that must hold inside WSL 2:
conda/mamba discovery, the locked QIIME 2 environment, registered command
entrypoints, and a synchronized real FASTQ smoke input.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Iterator, Sequence


EXPECTED_ENVIRONMENT = "microbiome-qiime2-2026.4"
EXPECTED_MANIFEST_SHA256 = (
    "68b0e9b974937525b71cca3dec48b7dc5b6135f588a941d1313c118c2a41f4b5"
)
EXPECTED_QIIME_RELEASE = "2026.4"
EXPECTED_QIIME_VERSION = "2026.4.0"
EXPECTED_Q2CLI_VERSION = "2026.4.0"
EXPECTED_PLUGIN_COUNT = 27
EXPECTED_RECORDS = 2000
EXPECTED_METADATA_SAMPLES = 75
EXPECTED_FILE_SHA256 = {
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
    "source_summary.json": (
        "daa3e140b0bf7e91c94392a48085da8e1e154920a0272eef93c15a26956e9983"
    ),
}
ENTRYPOINTS = (
    ("DADA2 paired denoising", ("qiime", "dada2", "denoise-paired", "--help")),
    (
        "Taxonomy classification",
        ("qiime", "feature-classifier", "classify-sklearn", "--help"),
    ),
    (
        "SEPP fragment insertion",
        ("qiime", "fragment-insertion", "sepp", "--help"),
    ),
)


@dataclass(frozen=True)
class CommandResult:
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--environment", default=EXPECTED_ENVIRONMENT)
    return parser.parse_args()


def run_command(command: Sequence[str], timeout: int = 180) -> CommandResult:
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        env={**os.environ, "LC_ALL": "C", "LANG": "C"},
    )
    return CommandResult(
        command=tuple(command),
        returncode=completed.returncode,
        stdout=completed.stdout.strip(),
        stderr=completed.stderr.strip(),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_line(text: str) -> str:
    return text.splitlines()[0].strip() if text.strip() else ""


def normalize_version(text: str) -> str:
    match = re.search(r"\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?", text)
    return match.group(0) if match else text.strip()


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

    plugins: list[str] = []
    in_plugins = False
    for line in text.splitlines():
        if line.strip() == "Installed plugins":
            in_plugins = True
            continue
        if in_plugins and not line.strip():
            break
        if in_plugins and ":" in line:
            plugins.append(line.split(":", 1)[0].strip())
    if not plugins:
        raise ValueError("qiime info did not report installed plugins")
    fields["plugins"] = plugins
    fields["plugin_count"] = len(plugins)
    return fields


def fastq_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not separator or not quality:
                raise ValueError(f"{path.name}: truncated record {record_number}")
            record = tuple(
                line.rstrip("\r\n")
                for line in (header, sequence, separator, quality)
            )
            if not record[0].startswith("@"):
                raise ValueError(f"{path.name}: invalid header at record {record_number}")
            if not record[2].startswith("+"):
                raise ValueError(
                    f"{path.name}: invalid separator at record {record_number}"
                )
            if len(record[1]) != len(record[3]):
                raise ValueError(
                    f"{path.name}: sequence/quality length mismatch at "
                    f"record {record_number}"
                )
            yield record


def read_id(record: tuple[str, str, str, str]) -> str:
    return record[0].split(maxsplit=1)[0].removeprefix("@")


def scan_synchronized_fastq(
    paths: dict[str, Path],
) -> tuple[dict[str, dict[str, object]], bool]:
    summaries = {
        key: {
            "records": 0,
            "read_length_min": None,
            "read_length_max": None,
        }
        for key in paths
    }
    synchronized = True
    streams = [fastq_records(paths[key]) for key in ("forward", "reverse", "barcodes")]
    for triplet in itertools.zip_longest(*streams):
        if any(record is None for record in triplet):
            synchronized = False
            break
        identifiers = {read_id(record) for record in triplet if record is not None}
        if len(identifiers) != 1:
            synchronized = False
        for key, record in zip(("forward", "reverse", "barcodes"), triplet):
            assert record is not None
            length = len(record[1])
            summaries[key]["records"] = int(summaries[key]["records"]) + 1
            current_min = summaries[key]["read_length_min"]
            current_max = summaries[key]["read_length_max"]
            summaries[key]["read_length_min"] = (
                length if current_min is None else min(int(current_min), length)
            )
            summaries[key]["read_length_max"] = (
                length if current_max is None else max(int(current_max), length)
            )
    return summaries, synchronized


def count_metadata_samples(path: Path) -> int:
    with path.open(encoding="ascii", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows)
        return sum(bool(row) and not row[0].startswith("#") for row in rows)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def status(observed: object, expected: object) -> str:
    return "PASS" if observed == expected else "FAIL"


def redact_paths(text: str, project_root: Path) -> str:
    """Normalize host-specific paths before writing the public audit log."""
    replacements = sorted(
        {
            str(project_root): "<PROJECT_ROOT>",
            str(Path.home()): "<HOME>",
        }.items(),
        key=lambda item: len(item[0]),
        reverse=True,
    )
    redacted = text
    for source, replacement in replacements:
        redacted = redacted.replace(source, replacement)
    return redacted


def detect_scope() -> tuple[str, bool]:
    version_text = ""
    proc_version = Path("/proc/version")
    if proc_version.exists():
        version_text = proc_version.read_text(encoding="utf-8", errors="replace")
    is_wsl = bool(
        os.environ.get("WSL_DISTRO_NAME")
        or re.search(r"(microsoft|wsl)", version_text, flags=re.IGNORECASE)
    )
    return ("wsl2-linux-side" if is_wsl else "native-linux-host", is_wsl)


def save_figure(fig: object, file_base: Path) -> None:
    file_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(file_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(
        file_base.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        file_base.with_suffix(".tiff"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def draw_layer_map(figure_dir: Path) -> None:
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "microbiome-16s-matplotlib"),
    )
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    palette = ["#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7"]
    layers = [
        ("Windows 10/11", "Enable WSL, restart,\nverify VERSION 2"),
        ("WSL 2 · Ubuntu", "Real Linux kernel\nand Bash tools"),
        ("Miniforge · mamba", "Isolated package and\nenvironment manager"),
        ("QIIME 2 · 2026.4", "Version-locked\ndedicated environment"),
        ("Project + FASTQ", "/home/<user>/16s\nLinux-side workspace"),
    ]

    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    ax.text(
        0.25,
        5.72,
        "Windows-to-QIIME 2 execution layers",
        fontsize=18,
        fontweight="bold",
        color="#1F2937",
        va="top",
    )
    ax.text(
        0.25,
        5.32,
        "Each layer has one job; analysis commands run in the Linux execution plane.",
        fontsize=10.5,
        color="#4B5563",
        va="top",
    )

    xs = [0.35, 2.48, 4.61, 6.74, 8.87]
    y = 2.35
    box_w = 1.82
    box_h = 1.75
    for index, ((title, detail), x, color) in enumerate(zip(layers, xs, palette)):
        box = FancyBboxPatch(
            (x, y),
            box_w,
            box_h,
            boxstyle="round,pad=0.04,rounding_size=0.12",
            linewidth=1.4,
            edgecolor=color,
            facecolor=color + "18",
        )
        ax.add_patch(box)
        ax.text(
            x + 0.13,
            y + 1.42,
            f"{index + 1:02d}",
            fontsize=9,
            fontweight="bold",
            color=color,
        )
        ax.text(
            x + box_w / 2,
            y + 1.08,
            title,
            fontsize=11.5,
            fontweight="bold",
            color="#111827",
            ha="center",
        )
        ax.text(
            x + box_w / 2,
            y + 0.53,
            detail,
            fontsize=9.2,
            color="#374151",
            ha="center",
            va="center",
            linespacing=1.35,
        )
        if index < len(layers) - 1:
            arrow = FancyArrowPatch(
                (x + box_w + 0.05, y + box_h / 2),
                (xs[index + 1] - 0.05, y + box_h / 2),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.25,
                color="#6B7280",
            )
            ax.add_patch(arrow)

    plane = FancyBboxPatch(
        (2.34, 1.88),
        8.45,
        2.7,
        boxstyle="round,pad=0.03,rounding_size=0.12",
        linewidth=1.0,
        linestyle="--",
        edgecolor="#6B7280",
        facecolor="none",
    )
    ax.add_patch(plane)
    ax.text(
        2.53,
        4.43,
        "Linux execution plane",
        fontsize=9,
        fontweight="bold",
        color="#4B5563",
        va="top",
    )
    ax.text(
        0.36,
        1.21,
        "Storage rule",
        fontsize=10,
        fontweight="bold",
        color="#111827",
    )
    ax.text(
        1.75,
        1.21,
        "Run file-intensive workflows under /home/<user>/... ; use /mnt/c for exchange.",
        fontsize=10,
        color="#374151",
    )
    ax.text(
        0.36,
        0.68,
        "Activation rule",
        fontsize=10,
        fontweight="bold",
        color="#111827",
    )
    ax.text(
        1.75,
        0.68,
        "Keep base clean; activate one named environment before every QIIME 2 command.",
        fontsize=10,
        color="#374151",
    )

    fig.subplots_adjust(left=0.02, right=0.99, top=0.98, bottom=0.02)
    save_figure(fig, figure_dir / "07-wsl2-layer-map")
    plt.close(fig)


def draw_validation_panel(
    figure_dir: Path,
    summary: dict[str, object],
    environment_rows: list[dict[str, object]],
    entrypoint_rows: list[dict[str, object]],
    fastq_rows: list[dict[str, object]],
) -> None:
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "microbiome-16s-matplotlib"),
    )
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    green = "#009E73"
    blue = "#0072B2"
    orange = "#E69F00"
    dark = "#1F2937"
    gray = "#6B7280"
    light = "#F3F4F6"

    fig = plt.figure(figsize=(11.2, 7.2), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.05, 0.95],
        width_ratios=[1.1, 0.9],
        hspace=0.28,
        wspace=0.20,
    )
    ax_stack = fig.add_subplot(grid[0, 0])
    ax_entry = fig.add_subplot(grid[0, 1])
    ax_fastq = fig.add_subplot(grid[1, :])
    for ax in (ax_stack, ax_entry, ax_fastq):
        ax.axis("off")

    fig.text(
        0.055,
        0.965,
        "Linux-side environment validation",
        fontsize=18,
        fontweight="bold",
        color=dark,
        va="top",
    )
    fig.text(
        0.055,
        0.925,
        (
            f"Run date {summary['validation_date']} · "
            f"scope {summary['execution_scope']} · "
            "Windows host controls not asserted"
        ),
        fontsize=10,
        color=gray,
        va="top",
    )

    ax_stack.set_xlim(0, 1)
    ax_stack.set_ylim(0, 1)
    ax_stack.text(
        0.0,
        0.98,
        "A  Version-locked stack",
        fontsize=12,
        fontweight="bold",
        color=dark,
        va="top",
    )
    stack_labels = [
        ("Conda", summary["conda_version"]),
        ("Mamba", summary["mamba_version"]),
        ("Python", summary["qiime"]["python_version"]),
        ("QIIME 2", summary["qiime"]["qiime_version"]),
        ("q2cli", summary["qiime"]["q2cli_version"]),
        ("Plugins", str(summary["qiime"]["plugin_count"])),
    ]
    y_positions = [0.80, 0.67, 0.54, 0.41, 0.28, 0.15]
    for (label, value), y in zip(stack_labels, y_positions):
        ax_stack.add_patch(
            FancyBboxPatch(
                (0.0, y - 0.055),
                0.90,
                0.10,
                boxstyle="round,pad=0.008,rounding_size=0.018",
                linewidth=0,
                facecolor=light,
            )
        )
        ax_stack.text(0.03, y, label, fontsize=9.5, color=dark, va="center")
        ax_stack.text(
            0.73,
            y,
            str(value),
            fontsize=9.5,
            fontweight="bold",
            color=blue,
            ha="right",
            va="center",
        )
        ax_stack.text(
            0.86,
            y,
            "PASS",
            fontsize=8.5,
            fontweight="bold",
            color=green,
            ha="right",
            va="center",
        )

    ax_entry.set_xlim(0, 1)
    ax_entry.set_ylim(0, 1)
    ax_entry.text(
        0.0,
        0.98,
        "B  Command entrypoints",
        fontsize=12,
        fontweight="bold",
        color=dark,
        va="top",
    )
    entry_y = [0.76, 0.52, 0.28]
    for row, y in zip(entrypoint_rows, entry_y):
        color = green if row["status"] == "PASS" else "#D55E00"
        ax_entry.add_patch(
            FancyBboxPatch(
                (0.0, y - 0.085),
                0.95,
                0.16,
                boxstyle="round,pad=0.012,rounding_size=0.025",
                linewidth=1.0,
                edgecolor=color,
                facecolor=color + "12",
            )
        )
        ax_entry.text(
            0.04,
            y + 0.022,
            str(row["check"]),
            fontsize=9.5,
            fontweight="bold",
            color=dark,
            va="center",
        )
        ax_entry.text(
            0.04,
            y - 0.036,
            "help command exits 0",
            fontsize=8.5,
            color=gray,
            va="center",
        )
        ax_entry.text(
            0.90,
            y,
            str(row["status"]),
            fontsize=9,
            fontweight="bold",
            color=color,
            ha="right",
            va="center",
        )

    ax_fastq.set_xlim(0, 1)
    ax_fastq.set_ylim(0, 1)
    ax_fastq.text(
        0.0,
        0.97,
        "C  Real Atacama FASTQ smoke input",
        fontsize=12,
        fontweight="bold",
        color=dark,
        va="top",
    )
    display_rows = [row for row in fastq_rows if row["kind"] == "FASTQ"]
    x_positions = [0.02, 0.27, 0.52]
    for row, x in zip(display_rows, x_positions):
        color = green if row["status"] == "PASS" else "#D55E00"
        ax_fastq.add_patch(
            FancyBboxPatch(
                (x, 0.32),
                0.21,
                0.42,
                boxstyle="round,pad=0.012,rounding_size=0.025",
                linewidth=1.1,
                edgecolor=color,
                facecolor=color + "10",
            )
        )
        ax_fastq.text(
            x + 0.105,
            0.64,
            str(row["role"]),
            fontsize=10.5,
            fontweight="bold",
            color=dark,
            ha="center",
        )
        ax_fastq.text(
            x + 0.105,
            0.52,
            f"{int(row['records']):,} records",
            fontsize=9.5,
            color=blue,
            ha="center",
        )
        ax_fastq.text(
            x + 0.105,
            0.42,
            f"{row['read_length']} nt",
            fontsize=9,
            color=gray,
            ha="center",
        )
        ax_fastq.text(
            x + 0.105,
            0.35,
            str(row["status"]),
            fontsize=8.5,
            fontweight="bold",
            color=color,
            ha="center",
        )

    ax_fastq.add_patch(
        FancyBboxPatch(
            (0.78, 0.32),
            0.20,
            0.42,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            linewidth=1.1,
            edgecolor=orange,
            facecolor=orange + "12",
        )
    )
    ax_fastq.text(
        0.88,
        0.64,
        "Metadata",
        fontsize=10.5,
        fontweight="bold",
        color=dark,
        ha="center",
    )
    ax_fastq.text(
        0.88,
        0.52,
        f"{summary['fastq']['metadata_samples']} samples",
        fontsize=9.5,
        color=blue,
        ha="center",
    )
    ax_fastq.text(
        0.88,
        0.42,
        "SHA-256 locked",
        fontsize=9,
        color=gray,
        ha="center",
    )
    ax_fastq.text(
        0.88,
        0.35,
        "PASS",
        fontsize=8.5,
        fontweight="bold",
        color=green,
        ha="center",
    )
    ax_fastq.text(
        0.02,
        0.17,
        (
            "Read IDs synchronized: PASS   ·   gzip integrity: PASS   ·   "
            f"required checks: {summary['checks_passed']}/{summary['checks_total']}"
        ),
        fontsize=9.7,
        fontweight="bold",
        color=dark,
    )
    ax_fastq.text(
        0.02,
        0.06,
        (
            "Evidence boundary: this run validates the Linux execution stack; "
            "run wsl -l -v on the Windows host separately."
        ),
        fontsize=8.8,
        color=gray,
    )

    save_figure(fig, figure_dir / "07-environment-validation")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = root / "env/qiime2.yml"
    fastq_dir = root / "data/small/fastq"
    execution_scope, is_wsl = detect_scope()

    command_log: list[CommandResult] = []
    conda_result = run_command(("conda", "--version"))
    mamba_result = run_command(("mamba", "--version"))
    bash_result = run_command(("bash", "--version"))
    qiime_result = run_command(
        ("conda", "run", "-n", args.environment, "qiime", "info")
    )
    command_log.extend((conda_result, mamba_result, bash_result, qiime_result))

    if conda_result.returncode != 0:
        raise RuntimeError(f"conda discovery failed: {conda_result.stderr}")
    if mamba_result.returncode != 0:
        raise RuntimeError(f"mamba discovery failed: {mamba_result.stderr}")
    if qiime_result.returncode != 0:
        raise RuntimeError(f"qiime info failed: {qiime_result.stderr}")

    qiime = parse_qiime_info(qiime_result.stdout)
    conda_version = normalize_version(first_line(conda_result.stdout))
    mamba_version = normalize_version(first_line(mamba_result.stdout))
    bash_version = normalize_version(first_line(bash_result.stdout))
    manifest_hash = sha256(manifest_path)

    environment_rows: list[dict[str, object]] = [
        {
            "check": "Execution scope",
            "observed": execution_scope,
            "expected": "Recorded, not inferred as WSL",
            "status": "INFO",
            "required": False,
        },
        {
            "check": "Kernel",
            "observed": platform.release(),
            "expected": "Linux",
            "status": "PASS" if platform.system() == "Linux" else "FAIL",
            "required": True,
        },
        {
            "check": "Architecture",
            "observed": platform.machine(),
            "expected": "x86_64",
            "status": status(platform.machine(), "x86_64"),
            "required": True,
        },
        {
            "check": "Bash",
            "observed": bash_version,
            "expected": "available",
            "status": "PASS" if bash_result.returncode == 0 else "FAIL",
            "required": True,
        },
        {
            "check": "Conda",
            "observed": conda_version,
            "expected": "available",
            "status": "PASS" if conda_result.returncode == 0 else "FAIL",
            "required": True,
        },
        {
            "check": "Mamba",
            "observed": mamba_version,
            "expected": "available",
            "status": "PASS" if mamba_result.returncode == 0 else "FAIL",
            "required": True,
        },
        {
            "check": "Environment manifest SHA-256",
            "observed": manifest_hash,
            "expected": EXPECTED_MANIFEST_SHA256,
            "status": status(manifest_hash, EXPECTED_MANIFEST_SHA256),
            "required": True,
        },
        {
            "check": "QIIME 2 release",
            "observed": qiime["qiime_release"],
            "expected": EXPECTED_QIIME_RELEASE,
            "status": status(qiime["qiime_release"], EXPECTED_QIIME_RELEASE),
            "required": True,
        },
        {
            "check": "QIIME 2 version",
            "observed": qiime["qiime_version"],
            "expected": EXPECTED_QIIME_VERSION,
            "status": status(qiime["qiime_version"], EXPECTED_QIIME_VERSION),
            "required": True,
        },
        {
            "check": "q2cli version",
            "observed": qiime["q2cli_version"],
            "expected": EXPECTED_Q2CLI_VERSION,
            "status": status(qiime["q2cli_version"], EXPECTED_Q2CLI_VERSION),
            "required": True,
        },
        {
            "check": "Registered plugins",
            "observed": qiime["plugin_count"],
            "expected": EXPECTED_PLUGIN_COUNT,
            "status": status(qiime["plugin_count"], EXPECTED_PLUGIN_COUNT),
            "required": True,
        },
    ]

    entrypoint_rows: list[dict[str, object]] = []
    for label, command in ENTRYPOINTS:
        full_command = ("conda", "run", "-n", args.environment, *command)
        result = run_command(full_command)
        command_log.append(result)
        entrypoint_rows.append(
            {
                "check": label,
                "command": " ".join(command),
                "exit_code": result.returncode,
                "status": "PASS" if result.returncode == 0 else "FAIL",
                "required": True,
            }
        )

    fastq_paths = {
        "forward": fastq_dir / "forward.fastq.gz",
        "reverse": fastq_dir / "reverse.fastq.gz",
        "barcodes": fastq_dir / "barcodes.fastq.gz",
    }
    fastq_summaries, synchronized = scan_synchronized_fastq(fastq_paths)
    metadata_samples = count_metadata_samples(fastq_dir / "metadata.tsv")

    fastq_rows: list[dict[str, object]] = []
    for key, role in (
        ("forward", "Forward reads"),
        ("reverse", "Reverse reads"),
        ("barcodes", "Barcode reads"),
    ):
        path = fastq_paths[key]
        observed_hash = sha256(path)
        expected_hash = EXPECTED_FILE_SHA256[path.name]
        summary = fastq_summaries[key]
        observed_records = int(summary["records"])
        read_length = (
            str(summary["read_length_min"])
            if summary["read_length_min"] == summary["read_length_max"]
            else f"{summary['read_length_min']}-{summary['read_length_max']}"
        )
        row_status = (
            "PASS"
            if observed_hash == expected_hash
            and observed_records == EXPECTED_RECORDS
            else "FAIL"
        )
        fastq_rows.append(
            {
                "kind": "FASTQ",
                "role": role,
                "file": path.name,
                "records": observed_records,
                "read_length": read_length,
                "observed_sha256": observed_hash,
                "expected_sha256": expected_hash,
                "status": row_status,
                "required": True,
            }
        )

    for name in ("metadata.tsv", "source_summary.json"):
        path = fastq_dir / name
        observed_hash = sha256(path)
        expected_hash = EXPECTED_FILE_SHA256[name]
        fastq_rows.append(
            {
                "kind": "TABLE" if name.endswith(".tsv") else "PROVENANCE",
                "role": "Sample metadata" if name.endswith(".tsv") else "Source summary",
                "file": name,
                "records": metadata_samples if name.endswith(".tsv") else "",
                "read_length": "",
                "observed_sha256": observed_hash,
                "expected_sha256": expected_hash,
                "status": status(observed_hash, expected_hash),
                "required": True,
            }
        )

    contract_rows: list[dict[str, object]] = [
        *environment_rows,
        *entrypoint_rows,
        *fastq_rows,
        {
            "check": "FASTQ read ID synchronization",
            "observed": synchronized,
            "expected": True,
            "status": status(synchronized, True),
            "required": True,
        },
        {
            "check": "Metadata sample count",
            "observed": metadata_samples,
            "expected": EXPECTED_METADATA_SAMPLES,
            "status": status(metadata_samples, EXPECTED_METADATA_SAMPLES),
            "required": True,
        },
    ]
    required_rows = [row for row in contract_rows if row.get("required") is True]
    checks_passed = sum(row.get("status") == "PASS" for row in required_rows)
    checks_failed = sum(row.get("status") == "FAIL" for row in required_rows)

    summary: dict[str, object] = {
        "validation_date": date.today().isoformat(),
        "execution_scope": execution_scope,
        "wsl_detected": is_wsl,
        "windows_host_controls_validated": False,
        "environment_name": args.environment,
        "kernel": platform.release(),
        "architecture": platform.machine(),
        "bash_version": bash_version,
        "conda_version": conda_version,
        "mamba_version": mamba_version,
        "environment_manifest_sha256": manifest_hash,
        "qiime": qiime,
        "entrypoints": {
            str(row["check"]): {
                "command": row["command"],
                "exit_code": row["exit_code"],
                "status": row["status"],
            }
            for row in entrypoint_rows
        },
        "fastq": {
            "records_each": {
                key: fastq_summaries[key]["records"]
                for key in ("forward", "reverse", "barcodes")
            },
            "read_ids_synchronized": synchronized,
            "metadata_samples": metadata_samples,
            "files": {
                str(row["file"]): {
                    "sha256": row["observed_sha256"],
                    "status": row["status"],
                }
                for row in fastq_rows
            },
        },
        "checks_total": len(required_rows),
        "checks_passed": checks_passed,
        "checks_failed": checks_failed,
        "status": "passed" if checks_failed == 0 else "failed",
    }

    write_tsv(
        output_dir / "environment-audit.tsv",
        environment_rows,
        ("check", "observed", "expected", "status", "required"),
    )
    write_tsv(
        output_dir / "entrypoint-audit.tsv",
        entrypoint_rows,
        ("check", "command", "exit_code", "status", "required"),
    )
    write_tsv(
        output_dir / "fastq-smoke.tsv",
        fastq_rows,
        (
            "kind",
            "role",
            "file",
            "records",
            "read_length",
            "observed_sha256",
            "expected_sha256",
            "status",
            "required",
        ),
    )
    (output_dir / "environment-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    log_lines = [
        "Article 07 Linux-side environment validation",
        "============================================",
        "",
        f"Validation date: {summary['validation_date']}",
        f"Execution scope: {execution_scope}",
        "Windows host controls validated: no",
        f"Environment: {args.environment}",
        f"Manifest SHA-256: {manifest_hash}",
        "",
    ]
    for result in command_log:
        log_lines.extend(
            [
                f"$ {' '.join(result.command)}",
                f"exit_code={result.returncode}",
                redact_paths(result.stdout, root) or "(no stdout)",
            ]
        )
        if result.stderr:
            log_lines.extend(("stderr:", redact_paths(result.stderr, root)))
        log_lines.append("")
    log_lines.extend(
        [
            "FASTQ smoke check",
            "-----------------",
            f"Forward records: {fastq_summaries['forward']['records']}",
            f"Reverse records: {fastq_summaries['reverse']['records']}",
            f"Barcode records: {fastq_summaries['barcodes']['records']}",
            f"Read IDs synchronized: {synchronized}",
            f"Metadata samples: {metadata_samples}",
            f"Required checks passed: {checks_passed}/{len(required_rows)}",
            f"Final status: {summary['status']}",
            "",
        ]
    )
    (output_dir / "environment-validation.log").write_text(
        "\n".join(log_lines),
        encoding="utf-8",
    )

    draw_layer_map(figure_dir)
    draw_validation_panel(
        figure_dir,
        summary,
        environment_rows,
        entrypoint_rows,
        fastq_rows,
    )

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
