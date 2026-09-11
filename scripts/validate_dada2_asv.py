#!/usr/bin/env python3
"""One-time q2-dada2 parameter and ASV audit for Article 12."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import shlex
import shutil
import statistics
import subprocess
import tempfile
import time
import zipfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_QIIME_ENV = "microbiome-qiime2-2026.4"
EXPECTED_QIIME_RELEASE = "2026.4"
EXPECTED_QIIME_VERSION = "2026.4.0"
EXPECTED_Q2CLI_VERSION = "2026.4.0"
EXPECTED_Q2_DADA2_VERSION = "2026.4.0"
EXPECTED_DADA2_VERSION = "1.38.0"
EXPECTED_INPUT_TYPE = "SampleData[PairedEndSequencesWithQuality]"
EXPECTED_INPUT_PAIRS = 8614
EXPECTED_SAMPLES = ("1", "1a", "2", "2a")
EXPECTED_SELECTED_STATS = {
    "1": {
        "input": 1810,
        "filtered": 1539,
        "denoised": 1286,
        "merged": 1170,
        "non-chimeric": 1144,
    },
    "1a": {
        "input": 2219,
        "filtered": 2026,
        "denoised": 1417,
        "merged": 1251,
        "non-chimeric": 1251,
    },
    "2": {
        "input": 2347,
        "filtered": 2033,
        "denoised": 1673,
        "merged": 1521,
        "non-chimeric": 1503,
    },
    "2a": {
        "input": 2238,
        "filtered": 2043,
        "denoised": 1510,
        "merged": 1315,
        "non-chimeric": 1315,
    },
}
EXPECTED_PROFILE_TOTALS = {
    "untruncated": {
        "filtered": 7365,
        "denoised": 5658,
        "merged": 4898,
        "non-chimeric": 4861,
        "asv_count": 340,
    },
    "reverse-tail-kept": {
        "filtered": 7515,
        "denoised": 5799,
        "merged": 5154,
        "non-chimeric": 5117,
        "asv_count": 359,
    },
    "selected": {
        "filtered": 7641,
        "denoised": 5886,
        "merged": 5257,
        "non-chimeric": 5213,
        "asv_count": 366,
    },
    "shortened": {
        "filtered": 8180,
        "denoised": 6456,
        "merged": 5250,
        "non-chimeric": 5223,
        "asv_count": 323,
    },
    "strict-ee": {
        "filtered": 6433,
        "denoised": 4912,
        "merged": 4516,
        "non-chimeric": 4478,
        "asv_count": 358,
    },
}
EXPECTED_ASV_LENGTH_COUNTS = {
    252: 7,
    253: 252,
    254: 91,
    255: 9,
    256: 2,
    285: 1,
    286: 1,
    289: 1,
    296: 1,
    303: 1,
}
EXPECTED_OUTPUT_TYPES = {
    "dada2-table.qza": "FeatureTable[Frequency]",
    "dada2-rep-seqs.qza": "FeatureData[Sequence]",
    "dada2-stats.qza": "SampleData[DADA2Stats]",
    "dada2-base-transitions.qza": "DADA2BaseTransitionStats",
    "dada2-feature-frequencies.qza": "ImmutableMetadata",
    "dada2-sample-frequencies.qza": "ImmutableMetadata",
}


@dataclass(frozen=True)
class Profile:
    profile_id: str
    label: str
    trunc_f: int
    trunc_r: int
    max_ee_f: float
    max_ee_r: float
    decision: str


@dataclass(frozen=True)
class CommandResult:
    label: str
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    elapsed_seconds: float


PROFILES = (
    Profile(
        "untruncated",
        "No 3′ truncation",
        0,
        0,
        2.0,
        2.0,
        "Reject: low-quality terminal bases reduce retention.",
    ),
    Profile(
        "reverse-tail-kept",
        "Keep reverse tail",
        220,
        220,
        2.0,
        2.0,
        "Reject: retaining 20 extra R2 bases loses reads without adding needed overlap.",
    ),
    Profile(
        "selected",
        "Selected",
        220,
        200,
        2.0,
        2.0,
        "Select: high retention, long informative reads, and a large overlap reserve.",
    ),
    Profile(
        "shortened",
        "Shortened",
        180,
        150,
        2.0,
        2.0,
        "Sensitivity only: similar reads retained but 90 informative bases are discarded.",
    ),
    Profile(
        "strict-ee",
        "Strict maxEE",
        220,
        200,
        1.0,
        1.0,
        "Reject: maxEE=1 removes substantially more pairs in this dataset.",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-artifact", type=Path)
    parser.add_argument("--input-summary", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--qiime-env", default=EXPECTED_QIIME_ENV)
    return parser.parse_args()


def run_command(
    label: str,
    command: Sequence[str],
    environment: dict[str, str],
    timeout: int = 1800,
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


def isolated_environment(root: Path, label: str) -> dict[str, str]:
    environment = os.environ.copy()
    safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "-", label)
    cache_root = root / "cache" / safe_label
    xdg = cache_root / "xdg"
    numba = cache_root / "numba"
    mpl = cache_root / "matplotlib"
    for path in (xdg, numba, mpl):
        path.mkdir(parents=True, exist_ok=True)
    environment.update(
        {
            "XDG_CACHE_HOME": str(xdg),
            "NUMBA_CACHE_DIR": str(numba),
            "MPLCONFIGDIR": str(mpl),
            "LC_ALL": "C.UTF-8",
            "LANG": "C.UTF-8",
        }
    )
    return environment


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


def extract_zip_text(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"{path.name}: expected one {suffix}, found {len(matches)}")
        return archive.read(matches[0]).decode("utf-8", errors="replace")


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
    metadata: dict[str, str] = {}
    for key in ("uuid", "type", "format", "data-size"):
        match = re.search(rf"^{re.escape(key)}:\s*(.+)$", text, flags=re.MULTILINE)
        metadata[key] = match.group(1).strip() if match else ""
    return metadata


def parse_stats_artifact(path: Path) -> list[dict[str, object]]:
    text = extract_zip_text(path, "/data/stats.tsv")
    rows: list[dict[str, object]] = []
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        sample_id = str(row["sample-id"])
        if sample_id.startswith("#"):
            continue
        parsed: dict[str, object] = {"sample_id": sample_id}
        for key in ("input", "filtered", "denoised", "merged", "non-chimeric"):
            parsed[key] = int(float(row[key]))
        parsed["filter_percent"] = float(row["percentage of input passed filter"])
        parsed["merged_percent"] = float(row["percentage of input merged"])
        parsed["non_chimeric_percent"] = float(
            row["percentage of input non-chimeric"]
        )
        rows.append(parsed)
    return rows


def parse_rep_sequences(path: Path) -> list[tuple[str, str]]:
    text = extract_zip_text(path, "/data/dna-sequences.fasta")
    records: list[tuple[str, str]] = []
    identifier = ""
    sequence_parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if identifier:
                records.append((identifier, "".join(sequence_parts)))
            identifier = line[1:].split(maxsplit=1)[0]
            sequence_parts = []
        else:
            sequence_parts.append(line.strip())
    if identifier:
        records.append((identifier, "".join(sequence_parts)))
    return records


def parse_quality_summary(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as archive:
        for direction in ("forward", "reverse"):
            suffix = f"/data/{direction}-seven-number-summaries.tsv"
            matches = [name for name in archive.namelist() if name.endswith(suffix)]
            if len(matches) != 1:
                raise ValueError(f"{path.name}: missing {suffix}")
            rows = list(
                csv.reader(
                    io.StringIO(
                        archive.read(matches[0]).decode("utf-8", errors="replace")
                    ),
                    delimiter="\t",
                )
            )
            cycles = [int(value) for value in rows[0][1:]]
            by_stat = {
                row[0]: [float(value) for value in row[1:]]
                for row in rows[1:]
            }
            for index, cycle in enumerate(cycles):
                records.append(
                    {
                        "direction": "R1" if direction == "forward" else "R2",
                        "cycle": cycle,
                        "count": int(by_stat["count"][index]),
                        "q02": by_stat["2%"][index],
                        "q09": by_stat["9%"][index],
                        "q25": by_stat["25%"][index],
                        "q50": by_stat["50%"][index],
                        "q75": by_stat["75%"][index],
                        "q91": by_stat["91%"][index],
                        "q98": by_stat["98%"][index],
                    }
                )
    return records


def parse_read_lengths(path: Path) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if "/data/" in name and name.endswith(".fastq.gz")
        ]
        for member in sorted(members):
            filename = Path(member).name
            direction = "R1" if "_R1_" in filename else "R2"
            sample_id = filename.split("_S", 1)[0]
            compressed = archive.read(member)
            lengths: list[int] = []
            qualities: list[float] = []
            with gzip.open(io.BytesIO(compressed), "rt", encoding="ascii") as handle:
                for index, line in enumerate(handle):
                    if index % 4 == 1:
                        lengths.append(len(line.rstrip("\r\n")))
                    elif index % 4 == 3:
                        qvalues = [ord(character) - 33 for character in line.rstrip()]
                        qualities.append(
                            sum(10 ** (-quality / 10) for quality in qvalues)
                        )
            if not lengths:
                raise ValueError(f"{filename}: no reads")
            ordered_ee = sorted(qualities)
            p95_index = max(0, int(0.95 * len(ordered_ee)) - 1)
            summaries.append(
                {
                    "sample_id": sample_id,
                    "direction": direction,
                    "read_pairs": len(lengths),
                    "length_min": min(lengths),
                    "length_median": float(statistics.median(lengths)),
                    "length_max": max(lengths),
                    "expected_errors_median": round(
                        float(statistics.median(qualities)), 4
                    ),
                    "expected_errors_p95": round(ordered_ee[p95_index], 4),
                }
            )
    return summaries


def provenance_actions(path: Path) -> list[dict[str, str]]:
    actions: list[dict[str, str]] = []
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if "/provenance/" in name and name.endswith("/action/action.yaml")
        ]
        for member in sorted(members):
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
            parameters: dict[str, str] = {}
            for key in (
                "trunc_len_f",
                "trunc_len_r",
                "trim_left_f",
                "trim_left_r",
                "max_ee_f",
                "max_ee_r",
                "min_overlap",
                "pooling_method",
                "chimera_method",
                "n_threads",
                "n_reads_learn",
                "hashed_feature_ids",
            ):
                match = re.search(
                    rf"^\s*-\s+{re.escape(key)}:\s*(.+)$",
                    block,
                    flags=re.MULTILINE,
                )
                if match:
                    parameters[key] = match.group(1).strip()
            actions.append(
                {
                    "member": member,
                    "kind": kind.group(1).strip() if kind else "unknown",
                    "plugin": plugin.group(1).strip() if plugin else "",
                    "action": action.group(1).strip() if action else "import",
                    "parameters": json.dumps(
                        parameters, ensure_ascii=False, sort_keys=True
                    ),
                    "raw_text": text,
                }
            )
    return actions


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


def draw_quality_truncation_gate(
    figure_dir: Path,
    quality_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt

    colors = {"R1": "#0072B2", "R2": "#D55E00"}
    cuts = {"R1": 220, "R2": 200}
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.9), facecolor="white")
    for axis, direction in zip(axes, ("R1", "R2")):
        rows = [row for row in quality_rows if row["direction"] == direction]
        cycles = [int(row["cycle"]) for row in rows]
        median = [float(row["q50"]) for row in rows]
        q25 = [float(row["q25"]) for row in rows]
        q09 = [float(row["q09"]) for row in rows]
        color = colors[direction]
        axis.fill_between(cycles, q09, q25, color=color, alpha=0.13, label="9th–25th percentile")
        axis.plot(cycles, q25, color=color, linewidth=1.25, linestyle="--", label="25th percentile")
        axis.plot(cycles, median, color=color, linewidth=2.0, label="Median")
        axis.axvline(cuts[direction], color="#111827", linewidth=1.4, linestyle=":")
        axis.text(
            cuts[direction] - 3,
            7,
            f"truncate = {cuts[direction]}",
            rotation=90,
            ha="right",
            va="bottom",
            fontsize=9,
            color="#111827",
            fontweight="bold",
        )
        axis.set_xlim(1, max(cycles))
        axis.set_ylim(0, 42)
        axis.set_xlabel("Cycle after primer removal")
        axis.set_ylabel("Phred score")
        axis.set_title(
            f"{'A' if direction == 'R1' else 'B'}  {direction} quality distribution",
            loc="left",
            fontweight="bold",
        )
        axis.grid(color="#E5E7EB", linewidth=0.7)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, loc="lower left", fontsize=8.5)
    fig.suptitle(
        "Truncation is a joint quality-and-overlap decision",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.055,
        0.925,
        "The selected 220/200 cut removes unstable terminal cycles while preserving a large merge reserve.",
        fontsize=10,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.13, wspace=0.24)
    save_figure(fig, figure_dir / "12-quality-truncation-gate")
    plt.close(fig)


def draw_parameter_loss_ledger(
    figure_dir: Path,
    sensitivity_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [str(row["label"]) for row in sensitivity_rows]
    retention = [float(row["retention_percent"]) for row in sensitivity_rows]
    selected_index = next(
        index
        for index, row in enumerate(sensitivity_rows)
        if row["profile_id"] == "selected"
    )
    colors = [
        "#009E73" if index == selected_index else "#9CA3AF"
        for index in range(len(labels))
    ]
    selected = sensitivity_rows[selected_index]
    stages = ["input", "filtered", "denoised", "merged", "non_chimeric"]
    stage_labels = ["Input", "Filtered", "Denoised", "Merged", "Non-chimeric"]
    stage_values = [int(selected[key]) for key in stages]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.3, 6.0),
        facecolor="white",
        gridspec_kw={"width_ratios": [1.12, 1]},
    )
    x = np.arange(len(labels))
    axes[0].bar(x, retention, color=colors, width=0.64)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=18, ha="right")
    axes[0].set_ylim(0, 70)
    axes[0].set_ylabel("Non-chimeric reads / input (%)")
    axes[0].set_title("A  Parameter sensitivity", loc="left", fontweight="bold")
    axes[0].grid(axis="y", color="#E5E7EB", linewidth=0.7)
    axes[0].set_axisbelow(True)
    for index, value in enumerate(retention):
        axes[0].text(
            index,
            value + 1.2,
            f"{value:.2f}%",
            ha="center",
            va="bottom",
            fontsize=8.7,
            fontweight="bold" if index == selected_index else "normal",
        )

    y = np.arange(len(stages))
    axes[1].plot(stage_values, y, color="#0072B2", linewidth=2.2, marker="o")
    axes[1].fill_betweenx(y, 0, stage_values, color="#0072B2", alpha=0.08)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(stage_labels)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 9200)
    axes[1].set_xlabel("Read pairs")
    axes[1].set_title("B  Selected-run loss ledger", loc="left", fontweight="bold")
    axes[1].grid(axis="x", color="#E5E7EB", linewidth=0.7)
    axes[1].set_axisbelow(True)
    for index, value in enumerate(stage_values):
        axes[1].text(
            value + 140,
            index,
            f"{value:,}",
            va="center",
            fontsize=9,
            color="#111827",
            fontweight="bold",
        )
    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle(
        "Judge DADA2 by the full read ledger—not by ASV count alone",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.055,
        0.925,
        "Selected: trunc-len-f 220, trunc-len-r 200, maxEE 2/2, exact overlap, consensus chimera removal.",
        fontsize=10,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.08, right=0.985, top=0.82, bottom=0.19, wspace=0.32)
    save_figure(fig, figure_dir / "12-denoising-loss-ledger")
    plt.close(fig)


def draw_overlap_budget(
    figure_dir: Path,
    modal_length: int,
    maximum_length: int,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    fig, ax = plt.subplots(figsize=(12.2, 6.1), facecolor="white")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.text(
        0.25,
        7.65,
        "Overlap must survive every trimming decision",
        fontsize=18,
        fontweight="bold",
        color="#111827",
        va="top",
    )
    ax.text(
        0.25,
        7.12,
        "Budget = retained R1 + retained R2 − biological insert; DADA2 requires at least 12 nt here.",
        fontsize=10,
        color="#4B5563",
        va="top",
    )

    rows = [
        {
            "y": 5.45,
            "title": "A  Selected 220 + 200",
            "sum": 420,
            "modal_overlap": 420 - modal_length,
            "worst_overlap": 420 - maximum_length,
            "color": "#009E73",
            "status": "PASS",
        },
        {
            "y": 2.38,
            "title": "B  Overtrimmed 130 + 120",
            "sum": 250,
            "modal_overlap": 250 - modal_length,
            "worst_overlap": 250 - maximum_length,
            "color": "#D55E00",
            "status": "FAIL",
        },
    ]
    for row in rows:
        y = float(row["y"])
        color = str(row["color"])
        ax.text(0.35, y + 1.05, row["title"], fontsize=12.5, fontweight="bold")
        ax.add_patch(
            FancyBboxPatch(
                (0.55, y),
                4.3,
                0.72,
                boxstyle="round,pad=0.03",
                facecolor="#0072B218",
                edgecolor="#0072B2",
            )
        )
        ax.text(
            2.7,
            y + 0.36,
            f"Retained read budget = {row['sum']} nt",
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )
        ax.add_patch(
            FancyBboxPatch(
                (5.15, y),
                3.1,
                0.72,
                boxstyle="round,pad=0.03",
                facecolor=color + "18",
                edgecolor=color,
            )
        )
        ax.text(
            6.70,
            y + 0.36,
            f"Modal overlap = {row['modal_overlap']} nt",
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )
        ax.add_patch(
            FancyBboxPatch(
                (8.55, y),
                2.0,
                0.72,
                boxstyle="round,pad=0.03",
                facecolor=color + "18",
                edgecolor=color,
            )
        )
        ax.text(
            9.55,
            y + 0.36,
            str(row["status"]),
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color=color,
        )
        ax.text(
            0.6,
            y - 0.45,
            (
                f"Observed insert range: {modal_length} nt mode to {maximum_length} nt maximum; "
                f"minimum observed overlap estimate = {row['worst_overlap']} nt."
            ),
            fontsize=9.2,
            color="#4B5563",
        )
    ax.text(
        0.55,
        0.52,
        "A quality plot can suggest a cut; only the overlap equation can approve it.",
        fontsize=10.2,
        color="#111827",
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.015, right=0.99, top=0.99, bottom=0.02)
    save_figure(fig, figure_dir / "12-overlap-budget")
    plt.close(fig)


def draw_asv_output_audit(
    figure_dir: Path,
    length_counts: Counter[int],
    sample_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    lengths = sorted(length_counts)
    counts = [length_counts[length] for length in lengths]
    samples = [str(row["sample_id"]) for row in sample_rows]
    retention = [float(row["non_chimeric_percent"]) for row in sample_rows]
    merged = [float(row["merged_percent"]) for row in sample_rows]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.8), facecolor="white")
    colors = ["#0072B2" if length <= 256 else "#E69F00" for length in lengths]
    axes[0].bar([str(length) for length in lengths], counts, color=colors)
    axes[0].set_xlabel("Merged ASV length (nt)")
    axes[0].set_ylabel("ASV count")
    axes[0].set_title("A  Representative-sequence lengths", loc="left", fontweight="bold")
    axes[0].grid(axis="y", color="#E5E7EB", linewidth=0.7)
    axes[0].set_axisbelow(True)
    for index, count in enumerate(counts):
        axes[0].text(index, count + 4, str(count), ha="center", fontsize=8.2)
    axes[0].text(
        0.98,
        0.93,
        "Orange = length-audit candidates",
        transform=axes[0].transAxes,
        ha="right",
        fontsize=8.8,
        color="#4B5563",
    )

    x = np.arange(len(samples))
    width = 0.34
    axes[1].bar(x - width / 2, merged, width, color="#56B4E9", label="Merged")
    axes[1].bar(x + width / 2, retention, width, color="#009E73", label="Non-chimeric")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(samples)
    axes[1].set_xlabel("Sample")
    axes[1].set_ylabel("Read pairs / input (%)")
    axes[1].set_ylim(0, 75)
    axes[1].set_title("B  Sample-level retention", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, loc="lower right")
    axes[1].grid(axis="y", color="#E5E7EB", linewidth=0.7)
    axes[1].set_axisbelow(True)
    for index, value in enumerate(retention):
        axes[1].text(
            index + width / 2,
            value + 1.1,
            f"{value:.1f}%",
            ha="center",
            fontsize=8.4,
            fontweight="bold",
        )
    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle(
        "The ASV table still needs biological and sample-level auditing",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.055,
        0.925,
        "Most V4 ASVs are 252–256 nt; five longer sequences are flagged for later taxonomy and off-target review.",
        fontsize=10,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.13, wspace=0.27)
    save_figure(fig, figure_dir / "12-asv-output-audit")
    plt.close(fig)


def run_profile(
    profile: Profile,
    input_artifact: Path,
    work_root: Path,
    qiime_env: str,
) -> tuple[Profile, CommandResult, Path]:
    profile_root = work_root / "profiles" / profile.profile_id
    profile_root.mkdir(parents=True, exist_ok=True)
    command = qiime_command(
        qiime_env,
        "dada2",
        "denoise-paired",
        "--i-demultiplexed-seqs",
        str(input_artifact),
        "--p-trunc-len-f",
        str(profile.trunc_f),
        "--p-trunc-len-r",
        str(profile.trunc_r),
        "--p-trim-left-f",
        "0",
        "--p-trim-left-r",
        "0",
        "--p-max-ee-f",
        str(profile.max_ee_f),
        "--p-max-ee-r",
        str(profile.max_ee_r),
        "--p-trunc-q",
        "2",
        "--p-min-overlap",
        "12",
        "--p-max-merge-mismatch",
        "0",
        "--p-pooling-method",
        "independent",
        "--p-chimera-method",
        "consensus",
        "--p-n-threads",
        "1",
        "--p-n-reads-learn",
        "1000000",
        "--p-hashed-feature-ids",
        "--p-retain-all-samples",
        "--o-table",
        str(profile_root / "table.qza"),
        "--o-representative-sequences",
        str(profile_root / "rep-seqs.qza"),
        "--o-denoising-stats",
        str(profile_root / "stats.qza"),
        "--o-base-transition-stats",
        str(profile_root / "base-transitions.qza"),
        "--verbose",
    )
    result = run_command(
        f"dada2-{profile.profile_id}",
        command,
        isolated_environment(work_root, profile.profile_id),
    )
    return profile, result, profile_root


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    input_artifact = (
        args.input_artifact.resolve()
        if args.input_artifact
        else project_root / "results/11-primer-trimming/nfcore-v4-trimmed.qza"
    )
    input_summary = (
        args.input_summary.resolve()
        if args.input_summary
        else project_root
        / "results/11-primer-trimming/nfcore-v4-trimmed-summary.qzv"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    checks: list[dict[str, object]] = []
    command_results: list[CommandResult] = []

    with tempfile.TemporaryDirectory(prefix="article12-dada2-") as temporary:
        work_root = Path(temporary)
        validator_mpl = work_root / "cache" / "validator-matplotlib"
        validator_mpl.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(validator_mpl)
        base_environment = isolated_environment(work_root, "base")

        info_result = run_command(
            "qiime-info",
            qiime_command(args.qiime_env, "info"),
            base_environment,
        )
        command_results.append(info_result)
        checks.append(
            audit_row(
                "qiime-info-exit",
                info_result.returncode,
                0,
                info_result.returncode == 0,
            )
        )
        info = parse_qiime_info(info_result.stdout) if info_result.returncode == 0 else {}
        plugins = info.get("plugins", {})
        checks.extend(
            [
                audit_row(
                    "qiime-release",
                    info.get("qiime_release"),
                    EXPECTED_QIIME_RELEASE,
                    info.get("qiime_release") == EXPECTED_QIIME_RELEASE,
                ),
                audit_row(
                    "qiime-version",
                    info.get("qiime_version"),
                    EXPECTED_QIIME_VERSION,
                    info.get("qiime_version") == EXPECTED_QIIME_VERSION,
                ),
                audit_row(
                    "q2cli-version",
                    info.get("q2cli_version"),
                    EXPECTED_Q2CLI_VERSION,
                    info.get("q2cli_version") == EXPECTED_Q2CLI_VERSION,
                ),
                audit_row(
                    "q2-dada2-version",
                    plugins.get("dada2") if isinstance(plugins, dict) else None,
                    EXPECTED_Q2_DADA2_VERSION,
                    isinstance(plugins, dict)
                    and plugins.get("dada2") == EXPECTED_Q2_DADA2_VERSION,
                ),
            ]
        )

        r_version_result = run_command(
            "dada2-r-version-and-randomize",
            conda_command(
                args.qiime_env,
                "Rscript",
                "-e",
                (
                    'cat(as.character(packageVersion("dada2")), "\\n"); '
                    'cat(as.character(formals(dada2::learnErrors)$randomize), "\\n")'
                ),
            ),
            base_environment,
        )
        command_results.append(r_version_result)
        r_lines = [
            line.strip()
            for line in r_version_result.stdout.splitlines()
            if line.strip()
        ]
        dada2_version = r_lines[0] if r_lines else ""
        randomize_default = r_lines[1] if len(r_lines) > 1 else ""
        checks.extend(
            [
                audit_row(
                    "dada2-r-command-exit",
                    r_version_result.returncode,
                    0,
                    r_version_result.returncode == 0,
                ),
                audit_row(
                    "dada2-r-version",
                    dada2_version,
                    EXPECTED_DADA2_VERSION,
                    dada2_version == EXPECTED_DADA2_VERSION,
                ),
                audit_row(
                    "learn-errors-randomize-default",
                    randomize_default,
                    "FALSE",
                    randomize_default == "FALSE",
                ),
            ]
        )

        peek_result = run_command(
            "input-peek",
            qiime_command(
                args.qiime_env,
                "tools",
                "peek",
                "--tsv",
                str(input_artifact),
            ),
            base_environment,
        )
        command_results.append(peek_result)
        input_peek = parse_peek_tsv(peek_result.stdout) if peek_result.returncode == 0 else {}
        checks.extend(
            [
                audit_row(
                    "input-artifact-exists",
                    input_artifact.exists(),
                    True,
                    input_artifact.exists(),
                ),
                audit_row(
                    "input-summary-exists",
                    input_summary.exists(),
                    True,
                    input_summary.exists(),
                ),
                audit_row(
                    "input-peek-exit",
                    peek_result.returncode,
                    0,
                    peek_result.returncode == 0,
                ),
                audit_row(
                    "input-semantic-type",
                    input_peek.get("Type"),
                    EXPECTED_INPUT_TYPE,
                    input_peek.get("Type") == EXPECTED_INPUT_TYPE,
                ),
            ]
        )

        validate_input = run_command(
            "input-maximum-validation",
            qiime_command(
                args.qiime_env,
                "tools",
                "validate",
                "--level",
                "max",
                str(input_artifact),
            ),
            base_environment,
        )
        command_results.append(validate_input)
        checks.append(
            audit_row(
                "input-maximum-validation",
                validate_input.returncode,
                0,
                validate_input.returncode == 0,
            )
        )

        quality_rows = parse_quality_summary(input_summary)
        read_length_rows = parse_read_lengths(input_artifact)
        checks.extend(
            [
                audit_row(
                    "quality-directions",
                    sorted({row["direction"] for row in quality_rows}),
                    ["R1", "R2"],
                    {row["direction"] for row in quality_rows} == {"R1", "R2"},
                ),
                audit_row(
                    "quality-cycle-count",
                    len(quality_rows),
                    463,
                    len(quality_rows) == 463,
                ),
                audit_row(
                    "input-fastq-files",
                    len(read_length_rows),
                    8,
                    len(read_length_rows) == 8,
                ),
                audit_row(
                    "input-pairs-from-fastq",
                    sum(
                        int(row["read_pairs"])
                        for row in read_length_rows
                        if row["direction"] == "R1"
                    ),
                    EXPECTED_INPUT_PAIRS,
                    sum(
                        int(row["read_pairs"])
                        for row in read_length_rows
                        if row["direction"] == "R1"
                    )
                    == EXPECTED_INPUT_PAIRS,
                ),
            ]
        )

        profile_results: dict[str, dict[str, object]] = {}
        profile_roots: dict[str, Path] = {}
        with ThreadPoolExecutor(max_workers=len(PROFILES)) as executor:
            futures = {
                executor.submit(
                    run_profile,
                    profile,
                    input_artifact,
                    work_root,
                    args.qiime_env,
                ): profile
                for profile in PROFILES
            }
            for future in as_completed(futures):
                profile, result, profile_root = future.result()
                command_results.append(result)
                checks.append(
                    audit_row(
                        f"profile-{profile.profile_id}-exit",
                        result.returncode,
                        0,
                        result.returncode == 0,
                    )
                )
                if result.returncode != 0:
                    continue
                sample_rows = parse_stats_artifact(profile_root / "stats.qza")
                sequences = parse_rep_sequences(profile_root / "rep-seqs.qza")
                totals = {
                    key: sum(int(row[key]) for row in sample_rows)
                    for key in ("input", "filtered", "denoised", "merged", "non-chimeric")
                }
                profile_results[profile.profile_id] = {
                    "profile": profile,
                    "sample_rows": sample_rows,
                    "sequences": sequences,
                    **totals,
                    "asv_count": len(sequences),
                }
                profile_roots[profile.profile_id] = profile_root

        for profile in PROFILES:
            observed = profile_results.get(profile.profile_id, {})
            expected = EXPECTED_PROFILE_TOTALS[profile.profile_id]
            checks.append(
                audit_row(
                    f"profile-{profile.profile_id}-input",
                    observed.get("input"),
                    EXPECTED_INPUT_PAIRS,
                    observed.get("input") == EXPECTED_INPUT_PAIRS,
                )
            )
            for metric, expected_value in expected.items():
                checks.append(
                    audit_row(
                        f"profile-{profile.profile_id}-{metric}",
                        observed.get(metric),
                        expected_value,
                        observed.get(metric) == expected_value,
                    )
                )

        selected = profile_results["selected"]
        selected_root = profile_roots["selected"]
        selected_outputs = {
            "dada2-table.qza": selected_root / "table.qza",
            "dada2-rep-seqs.qza": selected_root / "rep-seqs.qza",
            "dada2-stats.qza": selected_root / "stats.qza",
            "dada2-base-transitions.qza": selected_root / "base-transitions.qza",
        }
        for filename, source in selected_outputs.items():
            shutil.copy2(source, output_dir / filename)

        selected_sample_rows = list(selected["sample_rows"])
        selected_by_sample = {
            str(row["sample_id"]): row for row in selected_sample_rows
        }
        for sample_id in EXPECTED_SAMPLES:
            observed = selected_by_sample[sample_id]
            expected = EXPECTED_SELECTED_STATS[sample_id]
            for metric, expected_value in expected.items():
                checks.append(
                    audit_row(
                        f"selected-{sample_id}-{metric}",
                        observed[metric],
                        expected_value,
                        observed[metric] == expected_value,
                    )
                )

        selected_sequences = list(selected["sequences"])
        length_counts = Counter(len(sequence) for _, sequence in selected_sequences)
        checks.extend(
            [
                audit_row(
                    "selected-length-distribution",
                    dict(sorted(length_counts.items())),
                    EXPECTED_ASV_LENGTH_COUNTS,
                    dict(sorted(length_counts.items())) == EXPECTED_ASV_LENGTH_COUNTS,
                ),
                audit_row(
                    "selected-asv-count",
                    len(selected_sequences),
                    366,
                    len(selected_sequences) == 366,
                ),
                audit_row(
                    "selected-modal-length",
                    length_counts.most_common(1)[0][0],
                    253,
                    length_counts.most_common(1)[0][0] == 253,
                ),
                audit_row(
                    "selected-common-v4-lengths",
                    sum(count for length, count in length_counts.items() if 252 <= length <= 256),
                    361,
                    sum(count for length, count in length_counts.items() if 252 <= length <= 256)
                    == 361,
                ),
                audit_row(
                    "selected-length-audit-candidates",
                    sum(count for length, count in length_counts.items() if length > 256),
                    5,
                    sum(count for length, count in length_counts.items() if length > 256)
                    == 5,
                ),
            ]
        )

        modal_length = length_counts.most_common(1)[0][0]
        maximum_length = max(length_counts)
        modal_overlap = 220 + 200 - modal_length
        minimum_overlap = 220 + 200 - maximum_length
        checks.extend(
            [
                audit_row(
                    "selected-modal-overlap",
                    modal_overlap,
                    167,
                    modal_overlap == 167,
                ),
                audit_row(
                    "selected-minimum-observed-overlap",
                    minimum_overlap,
                    117,
                    minimum_overlap == 117,
                ),
                audit_row(
                    "selected-overlap-above-minimum",
                    minimum_overlap >= 12,
                    True,
                    minimum_overlap >= 12,
                ),
            ]
        )

        overtrim_root = work_root / "profiles" / "overtrimmed"
        overtrim_root.mkdir(parents=True, exist_ok=True)
        overtrim_profile = Profile(
            "overtrimmed",
            "Overtrimmed",
            130,
            120,
            2.0,
            2.0,
            "Expected failure: retained lengths cannot provide the required overlap.",
        )
        _, overtrim_result, _ = run_profile(
            overtrim_profile,
            input_artifact,
            work_root,
            args.qiime_env,
        )
        command_results.append(overtrim_result)
        overtrim_text = f"{overtrim_result.stdout}\n{overtrim_result.stderr}"
        checks.extend(
            [
                audit_row(
                    "overtrim-expected-nonzero-exit",
                    overtrim_result.returncode,
                    "non-zero",
                    overtrim_result.returncode != 0,
                ),
                audit_row(
                    "overtrim-no-table",
                    (overtrim_root / "table.qza").exists(),
                    False,
                    not (overtrim_root / "table.qza").exists(),
                ),
                audit_row(
                    "overtrim-error-evidence",
                    "valid sequence table" in overtrim_text,
                    True,
                    "valid sequence table" in overtrim_text,
                ),
            ]
        )

        visualization_outputs = (
            "dada2-stats.qzv",
            "dada2-table-summary.qzv",
            "dada2-feature-frequencies.qza",
            "dada2-sample-frequencies.qza",
            "dada2-rep-seqs.qzv",
            "dada2-base-transitions.qzv",
        )
        for filename in visualization_outputs:
            (output_dir / filename).unlink(missing_ok=True)

        visualization_specs = (
            (
                "dada2-stats-visualization",
                qiime_command(
                    args.qiime_env,
                    "metadata",
                    "tabulate",
                    "--m-input-file",
                    str(output_dir / "dada2-stats.qza"),
                    "--o-visualization",
                    str(output_dir / "dada2-stats.qzv"),
                ),
            ),
            (
                "dada2-table-summary",
                qiime_command(
                    args.qiime_env,
                    "feature-table",
                    "summarize",
                    "--i-table",
                    str(output_dir / "dada2-table.qza"),
                    "--o-feature-frequencies",
                    str(output_dir / "dada2-feature-frequencies.qza"),
                    "--o-sample-frequencies",
                    str(output_dir / "dada2-sample-frequencies.qza"),
                    "--o-summary",
                    str(output_dir / "dada2-table-summary.qzv"),
                ),
            ),
            (
                "dada2-rep-seqs-visualization",
                qiime_command(
                    args.qiime_env,
                    "feature-table",
                    "tabulate-seqs",
                    "--i-data",
                    str(output_dir / "dada2-rep-seqs.qza"),
                    "--o-visualization",
                    str(output_dir / "dada2-rep-seqs.qzv"),
                ),
            ),
            (
                "dada2-base-transitions-visualization",
                qiime_command(
                    args.qiime_env,
                    "dada2",
                    "plot-base-transitions",
                    "--i-base-transition-stats",
                    str(output_dir / "dada2-base-transitions.qza"),
                    "--o-visualization",
                    str(output_dir / "dada2-base-transitions.qzv"),
                ),
            ),
        )
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(
                    run_command,
                    label,
                    command,
                    isolated_environment(work_root, label),
                ): label
                for label, command in visualization_specs
            }
            for future in as_completed(futures):
                result = future.result()
                command_results.append(result)
                checks.append(
                    audit_row(
                        f"{result.label}-exit",
                        result.returncode,
                        0,
                        result.returncode == 0,
                    )
                )

        validation_specs = [
            (
                filename,
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "validate",
                    "--level",
                    "max",
                    str(output_dir / filename),
                ),
            )
            for filename in EXPECTED_OUTPUT_TYPES
        ]
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(
                    run_command,
                    f"validate-{filename}",
                    command,
                    isolated_environment(work_root, f"validate-{filename}"),
                ): filename
                for filename, command in validation_specs
            }
            for future in as_completed(futures):
                result = future.result()
                command_results.append(result)
                filename = futures[future]
                checks.append(
                    audit_row(
                        f"{filename}-maximum-validation",
                        result.returncode,
                        0,
                        result.returncode == 0,
                    )
                )

        artifact_rows: list[dict[str, object]] = []
        for filename, expected_type in EXPECTED_OUTPUT_TYPES.items():
            path = output_dir / filename
            metadata = artifact_metadata(path)
            artifact_rows.append(
                {
                    "artifact": filename,
                    "semantic_type": metadata["type"],
                    "data_format": metadata["format"],
                    "uuid_recorded": bool(metadata["uuid"]),
                    "maximum_validation": "PASS",
                }
            )
            checks.extend(
                [
                    audit_row(
                        f"{filename}-semantic-type",
                        metadata["type"],
                        expected_type,
                        metadata["type"] == expected_type,
                    ),
                    audit_row(
                        f"{filename}-uuid-recorded",
                        bool(metadata["uuid"]),
                        True,
                        bool(metadata["uuid"]),
                    ),
                ]
            )

        actions = provenance_actions(output_dir / "dada2-table.qza")
        action_names = [row["action"] for row in actions]
        dada_action = next(row for row in actions if row["action"] == "denoise_paired")
        dada_parameters = json.loads(dada_action["parameters"])
        expected_parameters = {
            "trunc_len_f": "220",
            "trunc_len_r": "200",
            "trim_left_f": "0",
            "trim_left_r": "0",
            "max_ee_f": "2.0",
            "max_ee_r": "2.0",
            "min_overlap": "12",
            "pooling_method": "independent",
            "chimera_method": "consensus",
            "n_threads": "1",
            "n_reads_learn": "1000000",
            "hashed_feature_ids": "true",
        }
        checks.extend(
            [
                audit_row(
                    "provenance-action-count",
                    len(actions),
                    3,
                    len(actions) == 3,
                ),
                audit_row(
                    "provenance-import",
                    "import" in action_names,
                    True,
                    "import" in action_names,
                ),
                audit_row(
                    "provenance-cutadapt",
                    "trim_paired" in action_names,
                    True,
                    "trim_paired" in action_names,
                ),
                audit_row(
                    "provenance-dada2",
                    "denoise_paired" in action_names,
                    True,
                    "denoise_paired" in action_names,
                ),
                audit_row(
                    "provenance-selected-parameters",
                    dada_parameters,
                    expected_parameters,
                    dada_parameters == expected_parameters,
                ),
            ]
        )

        sensitivity_rows: list[dict[str, object]] = []
        for profile in PROFILES:
            result = profile_results[profile.profile_id]
            sequences = list(result["sequences"])
            lengths = Counter(len(sequence) for _, sequence in sequences)
            sensitivity_rows.append(
                {
                    "profile_id": profile.profile_id,
                    "label": profile.label,
                    "trunc_len_f": profile.trunc_f,
                    "trunc_len_r": profile.trunc_r,
                    "max_ee_f": profile.max_ee_f,
                    "max_ee_r": profile.max_ee_r,
                    "input": result["input"],
                    "filtered": result["filtered"],
                    "denoised": result["denoised"],
                    "merged": result["merged"],
                    "non_chimeric": result["non-chimeric"],
                    "retention_percent": round(
                        100 * int(result["non-chimeric"]) / int(result["input"]), 2
                    ),
                    "asv_count": len(sequences),
                    "modal_asv_length": lengths.most_common(1)[0][0],
                    "maximum_asv_length": max(lengths),
                    "estimated_minimum_overlap": (
                        ""
                        if profile.trunc_f == 0 or profile.trunc_r == 0
                        else profile.trunc_f + profile.trunc_r - max(lengths)
                    ),
                    "decision": profile.decision,
                }
            )

        sample_retention_rows: list[dict[str, object]] = []
        for row in selected_sample_rows:
            input_count = int(row["input"])
            sample_retention_rows.append(
                {
                    **row,
                    "filtered_loss": input_count - int(row["filtered"]),
                    "denoising_loss": int(row["filtered"]) - int(row["denoised"]),
                    "merging_loss": int(row["denoised"]) - int(row["merged"]),
                    "chimera_loss": int(row["merged"]) - int(row["non-chimeric"]),
                }
            )

        asv_length_rows = [
            {
                "asv_length": length,
                "asv_count": count,
                "audit_class": (
                    "V4 common range" if 252 <= length <= 256 else "Review candidate"
                ),
            }
            for length, count in sorted(length_counts.items())
        ]
        overlap_rows = [
            {
                "profile_id": "selected",
                "trunc_len_f": 220,
                "trunc_len_r": 200,
                "retained_sum": 420,
                "insert_length_basis": "observed maximum",
                "insert_length": maximum_length,
                "estimated_overlap": minimum_overlap,
                "required_overlap": 12,
                "margin": minimum_overlap - 12,
                "status": "PASS",
            },
            {
                "profile_id": "selected-modal",
                "trunc_len_f": 220,
                "trunc_len_r": 200,
                "retained_sum": 420,
                "insert_length_basis": "observed mode",
                "insert_length": modal_length,
                "estimated_overlap": modal_overlap,
                "required_overlap": 12,
                "margin": modal_overlap - 12,
                "status": "PASS",
            },
            {
                "profile_id": "overtrimmed-modal",
                "trunc_len_f": 130,
                "trunc_len_r": 120,
                "retained_sum": 250,
                "insert_length_basis": "observed mode",
                "insert_length": modal_length,
                "estimated_overlap": 250 - modal_length,
                "required_overlap": 12,
                "margin": 250 - modal_length - 12,
                "status": "EXPECTED FAIL",
            },
        ]

        write_tsv(
            output_dir / "quality-by-cycle.tsv",
            quality_rows,
            (
                "direction",
                "cycle",
                "count",
                "q02",
                "q09",
                "q25",
                "q50",
                "q75",
                "q91",
                "q98",
            ),
        )
        write_tsv(
            output_dir / "read-length-audit.tsv",
            read_length_rows,
            (
                "sample_id",
                "direction",
                "read_pairs",
                "length_min",
                "length_median",
                "length_max",
                "expected_errors_median",
                "expected_errors_p95",
            ),
        )
        write_tsv(
            output_dir / "parameter-sensitivity.tsv",
            sensitivity_rows,
            (
                "profile_id",
                "label",
                "trunc_len_f",
                "trunc_len_r",
                "max_ee_f",
                "max_ee_r",
                "input",
                "filtered",
                "denoised",
                "merged",
                "non_chimeric",
                "retention_percent",
                "asv_count",
                "modal_asv_length",
                "maximum_asv_length",
                "estimated_minimum_overlap",
                "decision",
            ),
        )
        write_tsv(
            output_dir / "sample-retention.tsv",
            sample_retention_rows,
            (
                "sample_id",
                "input",
                "filtered",
                "filter_percent",
                "denoised",
                "merged",
                "merged_percent",
                "non-chimeric",
                "non_chimeric_percent",
                "filtered_loss",
                "denoising_loss",
                "merging_loss",
                "chimera_loss",
            ),
        )
        write_tsv(
            output_dir / "asv-length-distribution.tsv",
            asv_length_rows,
            ("asv_length", "asv_count", "audit_class"),
        )
        write_tsv(
            output_dir / "overlap-audit.tsv",
            overlap_rows,
            (
                "profile_id",
                "trunc_len_f",
                "trunc_len_r",
                "retained_sum",
                "insert_length_basis",
                "insert_length",
                "estimated_overlap",
                "required_overlap",
                "margin",
                "status",
            ),
        )
        write_tsv(
            output_dir / "artifact-audit.tsv",
            artifact_rows,
            (
                "artifact",
                "semantic_type",
                "data_format",
                "uuid_recorded",
                "maximum_validation",
            ),
        )
        write_tsv(
            output_dir / "provenance-audit.tsv",
            [
                {
                    key: value
                    for key, value in row.items()
                    if key != "raw_text"
                }
                for row in actions
            ],
            ("member", "kind", "plugin", "action", "parameters"),
        )
        write_tsv(
            output_dir / "environment-audit.tsv",
            [
                {
                    "component": "QIIME 2 release",
                    "observed": info.get("qiime_release", ""),
                    "expected": EXPECTED_QIIME_RELEASE,
                    "status": "PASS",
                },
                {
                    "component": "q2-dada2",
                    "observed": plugins.get("dada2", "") if isinstance(plugins, dict) else "",
                    "expected": EXPECTED_Q2_DADA2_VERSION,
                    "status": "PASS",
                },
                {
                    "component": "DADA2 R package",
                    "observed": dada2_version,
                    "expected": EXPECTED_DADA2_VERSION,
                    "status": "PASS",
                },
                {
                    "component": "learnErrors randomize default",
                    "observed": randomize_default,
                    "expected": "FALSE",
                    "status": "PASS",
                },
                {
                    "component": "DADA2 execution threads",
                    "observed": 1,
                    "expected": 1,
                    "status": "PASS",
                },
            ],
            ("component", "observed", "expected", "status"),
        )

        draw_quality_truncation_gate(figure_dir, quality_rows)
        draw_parameter_loss_ledger(figure_dir, sensitivity_rows)
        draw_overlap_budget(figure_dir, modal_length, maximum_length)
        draw_asv_output_audit(figure_dir, length_counts, selected_sample_rows)

        figure_bases = (
            "12-quality-truncation-gate",
            "12-denoising-loss-ledger",
            "12-overlap-budget",
            "12-asv-output-audit",
        )
        for base in figure_bases:
            for suffix in (".pdf", ".png", ".tiff"):
                path = figure_dir / f"{base}{suffix}"
                checks.append(
                    audit_row(
                        f"figure-{base}{suffix}",
                        path.exists() and path.stat().st_size > 0,
                        True,
                        path.exists() and path.stat().st_size > 0,
                    )
                )

        replacements = {
            str(project_root): "<PROJECT_ROOT>",
            str(output_dir): "<OUTPUT_DIR>",
            str(figure_dir): "<FIGURE_DIR>",
            str(work_root): "<WORK_DIR>",
            temporary: "<WORK_DIR>",
        }
        selected_result = next(
            result for result in command_results if result.label == "dada2-selected"
        )
        selected_log = (
            f"$ {shlex.join(selected_result.command)}\n\n"
            f"[exit_code] {selected_result.returncode}\n"
            f"[elapsed_seconds] {selected_result.elapsed_seconds}\n\n"
            f"[stdout]\n{selected_result.stdout}\n\n"
            f"[stderr]\n{selected_result.stderr}\n"
        )
        (output_dir / "qiime-dada2.log").write_text(
            sanitize_text(selected_log, replacements),
            encoding="utf-8",
        )
        overtrim_log = (
            f"$ {shlex.join(overtrim_result.command)}\n\n"
            f"[expected_exit_code] non-zero\n"
            f"[observed_exit_code] {overtrim_result.returncode}\n"
            f"[elapsed_seconds] {overtrim_result.elapsed_seconds}\n\n"
            f"[stdout]\n{overtrim_result.stdout}\n\n"
            f"[stderr]\n{overtrim_result.stderr}\n"
        )
        (output_dir / "overtrim-expected-failure.log").write_text(
            sanitize_text(overtrim_log, replacements),
            encoding="utf-8",
        )
        (output_dir / "input-validation.log").write_text(
            sanitize_text(
                (
                    f"$ {shlex.join(validate_input.command)}\n"
                    f"[exit_code] {validate_input.returncode}\n"
                    f"[stdout]\n{validate_input.stdout}\n"
                    f"[stderr]\n{validate_input.stderr}\n"
                ),
                replacements,
            ),
            encoding="utf-8",
        )

        checks_failed = sum(row["status"] != "PASS" for row in checks)
        checks_passed = len(checks) - checks_failed
        summary = {
            "article": 12,
            "validation_date": date.today().isoformat(),
            "dataset": "nf-core/ampliseq primer-trimmed MiSeq V4 paired-end test data",
            "qiime_environment": args.qiime_env,
            "qiime_release": info.get("qiime_release"),
            "q2_dada2_version": (
                plugins.get("dada2") if isinstance(plugins, dict) else None
            ),
            "dada2_r_version": dada2_version,
            "learn_errors_randomize_default": randomize_default,
            "profiles_compared": len(PROFILES),
            "samples": len(selected_sample_rows),
            "input_pairs": int(selected["input"]),
            "selected_trunc_len_f": 220,
            "selected_trunc_len_r": 200,
            "selected_max_ee_f": 2,
            "selected_max_ee_r": 2,
            "selected_filtered_pairs": int(selected["filtered"]),
            "selected_denoised_pairs": int(selected["denoised"]),
            "selected_merged_pairs": int(selected["merged"]),
            "selected_non_chimeric_pairs": int(selected["non-chimeric"]),
            "selected_retention_percent": round(
                100 * int(selected["non-chimeric"]) / int(selected["input"]), 2
            ),
            "selected_asv_count": len(selected_sequences),
            "modal_asv_length": modal_length,
            "maximum_asv_length": maximum_length,
            "length_audit_candidates": sum(
                count for length, count in length_counts.items() if length > 256
            ),
            "estimated_minimum_overlap": minimum_overlap,
            "provenance_actions": len(actions),
            "overtrim_expected_failure": overtrim_result.returncode != 0,
            "checks_total": len(checks),
            "checks_passed": checks_passed,
            "checks_failed": checks_failed,
            "status": "passed" if checks_failed == 0 else "failed",
        }
        (output_dir / "dada2-summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

        validation_lines = [
            "Article 12 q2-dada2 validation",
            f"date\t{summary['validation_date']}",
            f"status\t{summary['status']}",
            f"checks\t{checks_passed}/{len(checks)}",
            f"selected\t220/200; maxEE 2/2; {summary['selected_non_chimeric_pairs']}/{summary['input_pairs']} non-chimeric",
            "",
        ]
        validation_lines.extend(
            f"{row['status']}\t{row['check_id']}\tobserved={row['observed']}\texpected={row['expected']}"
            for row in checks
        )
        (output_dir / "validation.log").write_text(
            "\n".join(validation_lines) + "\n",
            encoding="utf-8",
        )

    if checks_failed:
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
