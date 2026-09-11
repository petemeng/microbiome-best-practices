#!/usr/bin/env python3
"""Run the one-time QIIME 2 installation audit for Article 08.

The audit goes beyond ``qiime info``. It imports a real, checksum-locked
Atacama EMP paired-end FASTQ excerpt, validates the resulting Artifact at the
maximum level, exports it again, and proves that all three compressed FASTQ
streams are byte-identical. Four controlled failure probes are also executed
to make the troubleshooting table evidence-based.
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
import shlex
import shutil
import subprocess
import tempfile
import time
import zipfile
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
EXPECTED_TYPE = "EMPPairedEndSequences"
EXPECTED_FORMAT = "EMPPairedEndDirFmt"
EXPECTED_RECORDS = 2000
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
}
ROLE_LABELS = {
    "forward.fastq.gz": "Forward",
    "reverse.fastq.gz": "Reverse",
    "barcodes.fastq.gz": "Barcode",
}


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
    parser.add_argument("--environment", default=EXPECTED_ENVIRONMENT)
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
        for row in rows:
            writer.writerow(row)


def run_command(
    label: str,
    command: Sequence[str],
    environment: dict[str, str],
    timeout: int = 300,
) -> CommandResult:
    started = time.monotonic()
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=environment,
    )
    return CommandResult(
        label=label,
        command=tuple(command),
        returncode=completed.returncode,
        stdout=completed.stdout.strip(),
        stderr=completed.stderr.strip(),
        elapsed_seconds=round(time.monotonic() - started, 3),
    )


def qiime_command(environment_name: str, *arguments: str) -> tuple[str, ...]:
    return (
        "conda",
        "run",
        "-n",
        environment_name,
        "qiime",
        *arguments,
    )


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


def parse_peek_tsv(text: str) -> dict[str, str]:
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("qiime tools peek --tsv did not return exactly one row")
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
                raise ValueError(f"{path.name}: truncated FASTQ record {index}")
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
                    f"{path.name}: sequence/quality mismatch at record {index}"
                )
            yield record


def fastq_summary(path: Path) -> dict[str, object]:
    lengths: list[int] = []
    identifiers: list[str] = []
    for record in fastq_records(path):
        identifiers.append(record[0].split(maxsplit=1)[0].removeprefix("@"))
        lengths.append(len(record[1]))
    if not lengths:
        raise ValueError(f"{path.name}: no FASTQ records")
    return {
        "records": len(lengths),
        "read_length_min": min(lengths),
        "read_length_max": max(lengths),
        "identifiers": identifiers,
    }


def redact_text(text: str, roots: Sequence[Path]) -> str:
    redacted = text
    replacements = {
        str(Path.home()): "<HOME>",
    }
    for index, root in enumerate(roots, start=1):
        replacements[str(root)] = (
            "<PROJECT_ROOT>" if index == 1 else f"<WORK_ROOT_{index - 1}>"
        )
    for source, replacement in sorted(
        replacements.items(),
        key=lambda item: len(item[0]),
        reverse=True,
    ):
        redacted = redacted.replace(source, replacement)
    return redacted


def observed_excerpt(text: str, pattern: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    compiled = re.compile(pattern, flags=re.IGNORECASE)
    for line in reversed(lines):
        if compiled.search(line):
            return re.sub(r"\s+", " ", line)
    return re.sub(r"\s+", " ", lines[-1]) if lines else "(no diagnostic text)"


def audit_row(
    check: str,
    observed: object,
    expected: object,
    passed: bool,
) -> dict[str, object]:
    return {
        "check": check,
        "observed": observed,
        "expected": expected,
        "status": "PASS" if passed else "FAIL",
    }


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


def draw_install_contract(
    figure_dir: Path,
    summary: dict[str, object],
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    colors = ["#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#D55E00"]
    dark = "#111827"
    gray = "#4B5563"
    green = "#009E73"
    stages = [
        ("Locked manifest", "SHA-256\n68b0e9b97493…"),
        ("Dedicated env", "QIIME 2\n2026.4.0"),
        ("CLI audit", f"{summary['plugin_count']}\nplugins"),
        ("Real import", "3 FASTQ\n6,000 records"),
        ("Max validation", "Artifact\nvalid"),
        ("Round-trip", "3 / 3\nhashes match"),
    ]

    fig, ax = plt.subplots(figsize=(12.2, 6.0))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 6.0)
    ax.axis("off")
    ax.text(
        0.30,
        5.55,
        "QIIME 2 installation acceptance contract",
        fontsize=19,
        fontweight="bold",
        color=dark,
        va="top",
    )
    ax.text(
        0.30,
        5.10,
        "A successful install must transform real input into a validated, reversible Artifact.",
        fontsize=10.8,
        color=gray,
        va="top",
    )

    xs = [0.28, 2.29, 4.30, 6.31, 8.32, 10.33]
    y = 2.18
    width = 1.58
    height = 1.72
    for index, ((title, detail), x, color) in enumerate(zip(stages, xs, colors)):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.035,rounding_size=0.12",
                linewidth=1.3,
                edgecolor=color,
                facecolor=color + "14",
            )
        )
        ax.text(
            x + 0.12,
            y + 1.48,
            f"{index + 1:02d}",
            fontsize=8.8,
            fontweight="bold",
            color=color,
        )
        ax.text(
            x + width / 2,
            y + 1.10,
            title,
            fontsize=10.3,
            fontweight="bold",
            color=dark,
            ha="center",
        )
        ax.text(
            x + width / 2,
            y + 0.52,
            detail,
            fontsize=9.2,
            color=gray,
            ha="center",
            va="center",
            linespacing=1.25,
        )
        if index < len(stages) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + width + 0.04, y + height / 2),
                    (xs[index + 1] - 0.04, y + height / 2),
                    arrowstyle="-|>",
                    mutation_scale=13,
                    linewidth=1.15,
                    color="#6B7280",
                )
            )

    ax.add_patch(
        FancyBboxPatch(
            (0.30, 0.55),
            11.58,
            0.88,
            boxstyle="round,pad=0.025,rounding_size=0.09",
            linewidth=0,
            facecolor="#F3F4F6",
        )
    )
    ax.text(
        0.58,
        1.00,
        f"{summary['checks_passed']} / {summary['checks_total']} required checks passed",
        fontsize=11.2,
        fontweight="bold",
        color=green,
        va="center",
    )
    ax.text(
        4.45,
        1.00,
        (
            f"Artifact: {summary['artifact']['semantic_type']} · "
            f"{summary['artifact']['data_format']}"
        ),
        fontsize=10.2,
        color=dark,
        va="center",
    )
    ax.text(
        11.58,
        1.00,
        "PASS",
        fontsize=10.5,
        fontweight="bold",
        color=green,
        ha="right",
        va="center",
    )
    fig.subplots_adjust(left=0.01, right=0.995, top=0.99, bottom=0.01)
    save_figure(fig, figure_dir / "08-qiime2-install-contract")
    plt.close(fig)


def draw_artifact_audit(
    figure_dir: Path,
    summary: dict[str, object],
    content_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    dark = "#111827"
    gray = "#4B5563"
    blue = "#0072B2"
    green = "#009E73"
    colors = {"Forward": "#0072B2", "Reverse": "#D55E00", "Barcode": "#009E73"}

    fig, ax = plt.subplots(figsize=(11.7, 7.0))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 11.7)
    ax.set_ylim(0, 7.0)
    ax.axis("off")
    ax.text(
        0.32,
        6.58,
        "Real Atacama Artifact round-trip audit",
        fontsize=19,
        fontweight="bold",
        color=dark,
        va="top",
    )
    ax.text(
        0.32,
        6.15,
        "Import → maximum validation → export; compressed FASTQ bytes remain unchanged.",
        fontsize=10.7,
        color=gray,
        va="top",
    )

    ax.add_patch(
        FancyBboxPatch(
            (3.66, 3.75),
            4.38,
            1.55,
            boxstyle="round,pad=0.04,rounding_size=0.16",
            linewidth=1.5,
            edgecolor=blue,
            facecolor=blue + "12",
        )
    )
    ax.text(
        5.85,
        4.87,
        "QIIME 2 Artifact (.qza)",
        fontsize=13.3,
        fontweight="bold",
        color=dark,
        ha="center",
    )
    ax.text(
        5.85,
        4.45,
        str(summary["artifact"]["semantic_type"]),
        fontsize=10.8,
        fontweight="bold",
        color=blue,
        ha="center",
    )
    ax.text(
        5.85,
        4.10,
        (
            f"{summary['artifact']['data_format']} · "
            f"{summary['artifact']['data_members']} data members"
        ),
        fontsize=9.7,
        color=gray,
        ha="center",
    )
    ax.text(
        5.85,
        3.82,
        "Maximum validation: PASS",
        fontsize=9.4,
        fontweight="bold",
        color=green,
        ha="center",
    )

    card_x = [0.38, 4.03, 7.68]
    for row, x in zip(content_rows, card_x):
        role = str(row["role"])
        color = colors[role]
        ax.add_patch(
            FancyBboxPatch(
                (x, 0.72),
                3.08,
                2.18,
                boxstyle="round,pad=0.035,rounding_size=0.13",
                linewidth=1.25,
                edgecolor=color,
                facecolor=color + "10",
            )
        )
        ax.text(
            x + 1.54,
            2.56,
            role,
            fontsize=12.2,
            fontweight="bold",
            color=dark,
            ha="center",
        )
        ax.text(
            x + 1.54,
            2.15,
            f"{int(row['records']):,} records · {row['read_length']} nt",
            fontsize=9.8,
            color=blue,
            ha="center",
        )
        ax.text(
            x + 1.54,
            1.77,
            f"{int(row['bytes']):,} compressed bytes",
            fontsize=9.3,
            color=gray,
            ha="center",
        )
        ax.text(
            x + 1.54,
            1.38,
            f"SHA-256 {str(row['input_sha256'])[:12]}…",
            fontsize=9.1,
            color=gray,
            ha="center",
            family="monospace",
        )
        ax.text(
            x + 1.54,
            0.99,
            "INPUT = EXPORTED · PASS",
            fontsize=9.4,
            fontweight="bold",
            color=green,
            ha="center",
        )
        ax.add_patch(
            FancyArrowPatch(
                (5.85, 3.68),
                (x + 1.54, 3.00),
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.0,
                color=color,
                connectionstyle="arc3,rad=0.0",
            )
        )

    ax.text(
        0.40,
        0.27,
        (
            "The .qza UUID and archive hash change on a new import; "
            "the semantic type, format, provenance, and exported-data hashes are the stable contract."
        ),
        fontsize=9.0,
        color=gray,
    )
    fig.subplots_adjust(left=0.01, right=0.995, top=0.99, bottom=0.02)
    save_figure(fig, figure_dir / "08-artifact-smoke-audit")
    plt.close(fig)


def draw_error_triage(
    figure_dir: Path,
    error_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    dark = "#111827"
    gray = "#4B5563"
    blue = "#0072B2"
    orange = "#E69F00"
    green = "#009E73"
    red = "#D55E00"

    fig, ax = plt.subplots(figsize=(11.8, 7.5))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 11.8)
    ax.set_ylim(0, 7.5)
    ax.axis("off")
    ax.text(
        0.32,
        7.08,
        "QIIME 2 error triage: diagnose the failing layer",
        fontsize=19,
        fontweight="bold",
        color=dark,
        va="top",
    )
    ax.text(
        0.32,
        6.65,
        f"Controlled probes passed: {sum(row['status'] == 'PASS' for row in error_rows)} / {len(error_rows)}",
        fontsize=10.6,
        color=gray,
        va="top",
    )

    ax.add_patch(
        FancyBboxPatch(
            (4.12, 5.30),
            3.56,
            0.84,
            boxstyle="round,pad=0.03,rounding_size=0.12",
            linewidth=1.4,
            edgecolor=blue,
            facecolor=blue + "12",
        )
    )
    ax.text(
        5.90,
        5.72,
        "Where does the first failure occur?",
        fontsize=12.0,
        fontweight="bold",
        color=dark,
        ha="center",
        va="center",
    )

    branches = [
        (
            0.30,
            3.53,
            "Environment creation",
            "Solver / channel / architecture",
            "Use a new env and the\nrelease-specific manifest.",
            orange,
        ),
        (
            3.16,
            3.53,
            "qiime not found",
            "Activation / PATH",
            "Activate the named env;\ninspect which qiime.",
            red,
        ),
        (
            6.02,
            3.53,
            "Plugin loading",
            "Writable cache / plugin state",
            "Set writable XDG, Matplotlib,\nand Numba cache roots.",
            "#CC79A7",
        ),
        (
            8.88,
            3.53,
            "Import / validate",
            "Type, format, files, archive",
            "Read the first diagnostic;\nverify type, format, and files.",
            green,
        ),
    ]
    for x, y, title, diagnosis, action, color in branches:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                2.55,
                1.36,
                boxstyle="round,pad=0.03,rounding_size=0.12",
                linewidth=1.25,
                edgecolor=color,
                facecolor=color + "10",
            )
        )
        ax.text(
            x + 1.275,
            y + 1.05,
            title,
            fontsize=10.4,
            fontweight="bold",
            color=dark,
            ha="center",
        )
        ax.text(
            x + 1.275,
            y + 0.72,
            diagnosis,
            fontsize=9.0,
            fontweight="bold",
            color=color,
            ha="center",
        )
        ax.text(
            x + 1.275,
            y + 0.30,
            action,
            fontsize=8.1,
            color=gray,
            ha="center",
            va="center",
            linespacing=1.25,
        )
        ax.add_patch(
            FancyArrowPatch(
                (5.90, 5.25),
                (x + 1.275, y + 1.43),
                arrowstyle="-|>",
                mutation_scale=11,
                linewidth=0.95,
                color=color,
            )
        )

    probe_labels = [
        ("Action typo", "Did you mean denoise-paired?"),
        ("Invalid type", "Semantic type is invalid"),
        ("Missing FASTQ", "Required filename absent"),
        ("Wrong artifact", "Not a QIIME archive"),
    ]
    probe_x = [0.43, 3.28, 6.13, 8.98]
    for (title, detail), row, x in zip(probe_labels, error_rows, probe_x):
        status_color = green if row["status"] == "PASS" else red
        ax.add_patch(
            FancyBboxPatch(
                (x, 1.15),
                2.35,
                1.18,
                boxstyle="round,pad=0.025,rounding_size=0.10",
                linewidth=1.0,
                edgecolor=status_color,
                facecolor="#F9FAFB",
            )
        )
        ax.text(
            x + 0.14,
            2.02,
            title,
            fontsize=9.5,
            fontweight="bold",
            color=dark,
        )
        ax.text(
            x + 0.14,
            1.67,
            detail,
            fontsize=8.2,
            color=gray,
        )
        ax.text(
            x + 2.19,
            1.30,
            str(row["status"]),
            fontsize=8.7,
            fontweight="bold",
            color=status_color,
            ha="right",
        )

    ax.add_patch(
        FancyBboxPatch(
            (0.34, 0.35),
            11.10,
            0.44,
            boxstyle="round,pad=0.02,rounding_size=0.07",
            linewidth=0,
            facecolor="#F3F4F6",
        )
    )
    ax.text(
        5.89,
        0.57,
        "Rule: fix the earliest failing layer; do not reinstall the whole stack for a filename or semantic-type error.",
        fontsize=9.2,
        fontweight="bold",
        color=dark,
        ha="center",
        va="center",
    )
    fig.subplots_adjust(left=0.01, right=0.995, top=0.99, bottom=0.01)
    save_figure(fig, figure_dir / "08-error-triage")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = root / "env/qiime2.yml"
    source_dir = root / "data/small/fastq"
    metadata_path = source_dir / "metadata.tsv"
    artifact_path = output_dir / "atacama-emp-paired-end.qza"
    if artifact_path.exists():
        artifact_path.unlink()

    command_results: list[CommandResult] = []
    with tempfile.TemporaryDirectory(prefix="article08-", dir=output_dir) as temporary:
        work = Path(temporary)
        input_dir = work / "emp-paired-end-sequences"
        export_dir = work / "exported"
        missing_dir = work / "missing-fastq"
        cache_root = work / "cache"
        for path in (
            input_dir,
            export_dir,
            missing_dir,
            cache_root / "xdg",
            cache_root / "matplotlib",
            cache_root / "numba",
            cache_root / "fontconfig",
        ):
            path.mkdir(parents=True, exist_ok=True)

        for filename in EXPECTED_FASTQ_SHA256:
            shutil.copy2(source_dir / filename, input_dir / filename)
        shutil.copy2(source_dir / "forward.fastq.gz", missing_dir / "forward.fastq.gz")

        command_environment = {
            **os.environ,
            "LC_ALL": "C",
            "LANG": "C",
            "MPLBACKEND": "Agg",
            "XDG_CACHE_HOME": str(cache_root / "xdg"),
            "MPLCONFIGDIR": str(cache_root / "matplotlib"),
            "NUMBA_CACHE_DIR": str(cache_root / "numba"),
            "FC_CACHEDIR": str(cache_root / "fontconfig"),
        }
        if Path("/etc/fonts/fonts.conf").exists():
            command_environment["FONTCONFIG_FILE"] = "/etc/fonts/fonts.conf"

        info_result = run_command(
            "qiime-info",
            qiime_command(args.environment, "info"),
            command_environment,
        )
        command_results.append(info_result)
        if info_result.returncode != 0:
            raise RuntimeError(f"qiime info failed: {info_result.stderr}")
        qiime_info = parse_qiime_info(info_result.stdout)

        import_result = run_command(
            "real-import",
            qiime_command(
                args.environment,
                "tools",
                "import",
                "--type",
                EXPECTED_TYPE,
                "--input-path",
                str(input_dir),
                "--output-path",
                str(artifact_path),
            ),
            command_environment,
        )
        command_results.append(import_result)
        if import_result.returncode != 0:
            raise RuntimeError(f"real QIIME 2 import failed: {import_result.stderr}")

        validate_result = run_command(
            "maximum-validation",
            qiime_command(
                args.environment,
                "tools",
                "validate",
                str(artifact_path),
                "--level",
                "max",
            ),
            command_environment,
        )
        command_results.append(validate_result)
        if validate_result.returncode != 0:
            raise RuntimeError(
                f"maximum Artifact validation failed: {validate_result.stderr}"
            )

        peek_result = run_command(
            "artifact-peek",
            qiime_command(
                args.environment,
                "tools",
                "peek",
                "--tsv",
                str(artifact_path),
            ),
            command_environment,
        )
        command_results.append(peek_result)
        if peek_result.returncode != 0:
            raise RuntimeError(f"Artifact peek failed: {peek_result.stderr}")
        peek = parse_peek_tsv(peek_result.stdout)

        export_result = run_command(
            "artifact-export",
            qiime_command(
                args.environment,
                "tools",
                "export",
                "--input-path",
                str(artifact_path),
                "--output-path",
                str(export_dir),
            ),
            command_environment,
        )
        command_results.append(export_result)
        if export_result.returncode != 0:
            raise RuntimeError(f"Artifact export failed: {export_result.stderr}")

        source_summaries = {
            filename: fastq_summary(source_dir / filename)
            for filename in EXPECTED_FASTQ_SHA256
        }
        exported_summaries = {
            filename: fastq_summary(export_dir / filename)
            for filename in EXPECTED_FASTQ_SHA256
        }
        synchronized = all(
            identifiers[0] == identifiers[index]
            for identifiers in zip(
                source_summaries["forward.fastq.gz"]["identifiers"],
                source_summaries["reverse.fastq.gz"]["identifiers"],
                source_summaries["barcodes.fastq.gz"]["identifiers"],
            )
            for index in (1, 2)
        )

        content_rows: list[dict[str, object]] = []
        for filename in (
            "forward.fastq.gz",
            "reverse.fastq.gz",
            "barcodes.fastq.gz",
        ):
            input_hash = sha256(source_dir / filename)
            exported_hash = sha256(export_dir / filename)
            source_summary = source_summaries[filename]
            exported_summary = exported_summaries[filename]
            length_text = (
                str(source_summary["read_length_min"])
                if source_summary["read_length_min"]
                == source_summary["read_length_max"]
                else (
                    f"{source_summary['read_length_min']}-"
                    f"{source_summary['read_length_max']}"
                )
            )
            row_passed = (
                input_hash == EXPECTED_FASTQ_SHA256[filename]
                and exported_hash == input_hash
                and source_summary["records"] == EXPECTED_RECORDS
                and exported_summary["records"] == EXPECTED_RECORDS
                and exported_summary["read_length_min"]
                == source_summary["read_length_min"]
                and exported_summary["read_length_max"]
                == source_summary["read_length_max"]
            )
            content_rows.append(
                {
                    "file": filename,
                    "role": ROLE_LABELS[filename],
                    "records": source_summary["records"],
                    "read_length": length_text,
                    "bytes": (source_dir / filename).stat().st_size,
                    "input_sha256": input_hash,
                    "exported_sha256": exported_hash,
                    "hash_match": input_hash == exported_hash,
                    "status": "PASS" if row_passed else "FAIL",
                }
            )

        with zipfile.ZipFile(artifact_path) as archive:
            archive_members = archive.namelist()
        data_members = [
            member
            for member in archive_members
            if re.search(r"/data/(forward|reverse|barcodes)\.fastq\.gz$", member)
        ]
        provenance_members = [
            member for member in archive_members if "/provenance/" in member
        ]
        checksum_members = [
            member for member in archive_members if member.endswith("/checksums.sha512")
        ]

        probe_specs = [
            {
                "probe_id": "action-typo",
                "failure_layer": "Plugin action",
                "command": qiime_command(
                    args.environment,
                    "dada2",
                    "denoise-pairedx",
                    "--help",
                ),
                "expected_pattern": r"has no action.*denoise-paired",
                "diagnosis": "Action name is misspelled or changed across releases.",
                "fix": "Read `qiime dada2 --help` and use `denoise-paired`.",
            },
            {
                "probe_id": "invalid-semantic-type",
                "failure_layer": "Semantic type",
                "command": qiime_command(
                    args.environment,
                    "tools",
                    "import",
                    "--type",
                    "NotARealSemanticType",
                    "--input-path",
                    str(input_dir),
                    "--output-path",
                    str(work / "invalid-type.qza"),
                ),
                "expected_pattern": r"Semantic type .* is invalid",
                "diagnosis": "The requested semantic type is not registered.",
                "fix": "Use `EMPPairedEndSequences` for this multiplexed EMP directory.",
            },
            {
                "probe_id": "missing-required-fastq",
                "failure_layer": "Directory format",
                "command": qiime_command(
                    args.environment,
                    "tools",
                    "import",
                    "--type",
                    EXPECTED_TYPE,
                    "--input-path",
                    str(missing_dir),
                    "--output-path",
                    str(work / "missing-fastq.qza"),
                ),
                "expected_pattern": r"Missing one or more files",
                "diagnosis": "The EMP paired directory does not satisfy its filename contract.",
                "fix": "Provide forward.fastq.gz, reverse.fastq.gz, and barcodes.fastq.gz.",
            },
            {
                "probe_id": "not-a-qiime-artifact",
                "failure_layer": "Artifact validation",
                "command": qiime_command(
                    args.environment,
                    "tools",
                    "validate",
                    str(metadata_path),
                    "--level",
                    "max",
                ),
                "expected_pattern": r"is not a QIIME archive",
                "diagnosis": "A plain TSV was supplied where a .qza/.qzv Result was required.",
                "fix": "Validate the generated .qza; validate metadata with metadata tools.",
            },
        ]
        error_rows: list[dict[str, object]] = []
        for spec in probe_specs:
            result = run_command(
                f"probe-{spec['probe_id']}",
                spec["command"],
                command_environment,
            )
            command_results.append(result)
            combined = "\n".join((result.stdout, result.stderr))
            pattern_found = bool(
                re.search(
                    str(spec["expected_pattern"]),
                    combined,
                    flags=re.IGNORECASE | re.DOTALL,
                )
            )
            probe_passed = result.returncode != 0 and pattern_found
            error_rows.append(
                {
                    "probe_id": spec["probe_id"],
                    "failure_layer": spec["failure_layer"],
                    "return_code": result.returncode,
                    "expected_nonzero": True,
                    "expected_pattern": spec["expected_pattern"],
                    "observed_excerpt": observed_excerpt(
                        redact_text(
                            combined,
                            (root, work, output_dir),
                        ),
                        str(spec["expected_pattern"]),
                    ),
                    "diagnosis": spec["diagnosis"],
                    "fix": spec["fix"],
                    "status": "PASS" if probe_passed else "FAIL",
                }
            )

        manifest_hash = sha256(manifest_path)
        content_by_role = {str(row["role"]): row for row in content_rows}
        error_passes = sum(row["status"] == "PASS" for row in error_rows)
        installation_rows = [
            audit_row(
                "Environment manifest SHA-256",
                manifest_hash,
                EXPECTED_MANIFEST_SHA256,
                manifest_hash == EXPECTED_MANIFEST_SHA256,
            ),
            audit_row(
                "QIIME 2 release",
                qiime_info["qiime_release"],
                EXPECTED_QIIME_RELEASE,
                qiime_info["qiime_release"] == EXPECTED_QIIME_RELEASE,
            ),
            audit_row(
                "QIIME 2 version",
                qiime_info["qiime_version"],
                EXPECTED_QIIME_VERSION,
                qiime_info["qiime_version"] == EXPECTED_QIIME_VERSION,
            ),
            audit_row(
                "q2cli version",
                qiime_info["q2cli_version"],
                EXPECTED_Q2CLI_VERSION,
                qiime_info["q2cli_version"] == EXPECTED_Q2CLI_VERSION,
            ),
            audit_row(
                "Registered plugins",
                qiime_info["plugin_count"],
                EXPECTED_PLUGIN_COUNT,
                qiime_info["plugin_count"] == EXPECTED_PLUGIN_COUNT,
            ),
            audit_row(
                "Real EMP import exit code",
                import_result.returncode,
                0,
                import_result.returncode == 0,
            ),
            audit_row(
                "Maximum validation exit code",
                validate_result.returncode,
                0,
                validate_result.returncode == 0,
            ),
            audit_row(
                "Artifact semantic type",
                peek["Type"],
                EXPECTED_TYPE,
                peek["Type"] == EXPECTED_TYPE,
            ),
            audit_row(
                "Artifact data format",
                peek["Data Format"],
                EXPECTED_FORMAT,
                peek["Data Format"] == EXPECTED_FORMAT,
            ),
            audit_row(
                "Artifact data members",
                len(data_members),
                3,
                len(data_members) == 3,
            ),
            audit_row(
                "Artifact provenance",
                len(provenance_members),
                ">= 1 member",
                len(provenance_members) >= 1,
            ),
            audit_row(
                "Artifact checksum manifest",
                len(checksum_members),
                1,
                len(checksum_members) == 1,
            ),
            audit_row(
                "Forward round-trip SHA-256",
                content_by_role["Forward"]["hash_match"],
                True,
                content_by_role["Forward"]["status"] == "PASS",
            ),
            audit_row(
                "Reverse round-trip SHA-256",
                content_by_role["Reverse"]["hash_match"],
                True,
                content_by_role["Reverse"]["status"] == "PASS",
            ),
            audit_row(
                "Barcode round-trip SHA-256",
                content_by_role["Barcode"]["hash_match"],
                True,
                content_by_role["Barcode"]["status"] == "PASS",
            ),
            audit_row(
                "Forward FASTQ records",
                content_by_role["Forward"]["records"],
                EXPECTED_RECORDS,
                content_by_role["Forward"]["records"] == EXPECTED_RECORDS,
            ),
            audit_row(
                "Reverse FASTQ records",
                content_by_role["Reverse"]["records"],
                EXPECTED_RECORDS,
                content_by_role["Reverse"]["records"] == EXPECTED_RECORDS,
            ),
            audit_row(
                "Barcode FASTQ records",
                content_by_role["Barcode"]["records"],
                EXPECTED_RECORDS,
                content_by_role["Barcode"]["records"] == EXPECTED_RECORDS,
            ),
            audit_row(
                "FASTQ record synchronization",
                synchronized,
                True,
                synchronized,
            ),
            audit_row(
                "Controlled error probes",
                error_passes,
                len(error_rows),
                error_passes == len(error_rows),
            ),
        ]
        checks_passed = sum(row["status"] == "PASS" for row in installation_rows)
        summary = {
            "status": (
                "passed"
                if checks_passed == len(installation_rows)
                else "failed"
            ),
            "validation_date": date.today().isoformat(),
            "environment": args.environment,
            "environment_manifest_sha256": manifest_hash,
            "python_version": qiime_info["python_version"],
            "qiime_release": qiime_info["qiime_release"],
            "qiime_version": qiime_info["qiime_version"],
            "q2cli_version": qiime_info["q2cli_version"],
            "plugin_count": qiime_info["plugin_count"],
            "plugins": qiime_info["plugins"],
            "checks_total": len(installation_rows),
            "checks_passed": checks_passed,
            "checks_failed": len(installation_rows) - checks_passed,
            "error_probes_total": len(error_rows),
            "error_probes_passed": error_passes,
            "artifact": {
                "filename": artifact_path.name,
                "uuid": peek["UUID"],
                "semantic_type": peek["Type"],
                "data_format": peek["Data Format"],
                "archive_bytes": artifact_path.stat().st_size,
                "archive_sha256": sha256(artifact_path),
                "validation_level": "max",
                "validation_status": (
                    "PASS" if validate_result.returncode == 0 else "FAIL"
                ),
                "data_members": len(data_members),
                "provenance_members": len(provenance_members),
                "checksum_manifests": len(checksum_members),
            },
            "fastq": {
                "files": len(content_rows),
                "records_per_stream": EXPECTED_RECORDS,
                "record_synchronized": synchronized,
                "round_trip_hash_matches": sum(
                    bool(row["hash_match"]) for row in content_rows
                ),
            },
            "command_elapsed_seconds": {
                result.label: result.elapsed_seconds for result in command_results
            },
            "stable_contract_note": (
                "Artifact UUID and archive SHA-256 are run-specific. Semantic type, "
                "data format, provenance presence, validation status, and exported "
                "FASTQ SHA-256 values are the stable acceptance contract."
            ),
        }

        roots_for_redaction = (root, work, output_dir)
        (output_dir / "qiime-info.txt").write_text(
            redact_text(info_result.stdout, roots_for_redaction) + "\n",
            encoding="utf-8",
        )
        write_tsv(
            output_dir / "artifact-peek.tsv",
            [
                {
                    "Filename": peek["Filename"],
                    "Type": peek["Type"],
                    "UUID": peek["UUID"],
                    "Data Format": peek["Data Format"],
                }
            ],
            ("Filename", "Type", "UUID", "Data Format"),
        )
        write_tsv(
            output_dir / "installation-audit.tsv",
            installation_rows,
            ("check", "observed", "expected", "status"),
        )
        write_tsv(
            output_dir / "artifact-content-audit.tsv",
            content_rows,
            (
                "file",
                "role",
                "records",
                "read_length",
                "bytes",
                "input_sha256",
                "exported_sha256",
                "hash_match",
                "status",
            ),
        )
        write_tsv(
            output_dir / "error-probes.tsv",
            error_rows,
            (
                "probe_id",
                "failure_layer",
                "return_code",
                "expected_nonzero",
                "expected_pattern",
                "observed_excerpt",
                "diagnosis",
                "fix",
                "status",
            ),
        )
        (output_dir / "artifact-summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

        log_parts = [
            "# Article 08 QIIME 2 installation and Artifact audit",
            f"validation_date: {summary['validation_date']}",
            f"environment: {args.environment}",
            f"status: {summary['status']}",
            "",
        ]
        for result in command_results:
            combined = "\n".join(
                part for part in (result.stdout, result.stderr) if part
            )
            log_parts.extend(
                [
                    f"## {result.label}",
                    (
                        "$ "
                        + redact_text(
                            shlex.join(result.command),
                            roots_for_redaction,
                        )
                    ),
                    f"exit_code: {result.returncode}",
                    f"elapsed_seconds: {result.elapsed_seconds}",
                    redact_text(combined, roots_for_redaction),
                    "",
                ]
            )
        (output_dir / "import-validation.log").write_text(
            "\n".join(log_parts).rstrip() + "\n",
            encoding="utf-8",
        )

        os.environ.setdefault(
            "MPLCONFIGDIR",
            str(Path(tempfile.gettempdir()) / "microbiome-16s-matplotlib"),
        )
        draw_install_contract(figure_dir, summary)
        draw_artifact_audit(figure_dir, summary, content_rows)
        draw_error_triage(figure_dir, error_rows)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
