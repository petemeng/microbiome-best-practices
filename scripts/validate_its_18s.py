#!/usr/bin/env python3
"""Run and audit the ITS2/18S branch contract for Article 13."""

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
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence


QIIME_ENV = "microbiome-qiime2-2026.4"
QIIME_VERSION = "2026.4.0"
ITSXPRESS_VERSION = "2.1.4"
BIOPYTHON_VERSION = "1.87"
VSEARCH_VERSION = "2.30.4"
HMMER_VERSION = "3.4"
INPUT_SAMPLES = ("sample1", "sample2")
INPUT_PAIRS_PER_SAMPLE = 2_500
INPUT_TOTAL_PAIRS = 5_000
UNITE_EVAL_SHA256 = (
    "f59c439ced1622c372a647464c69c2069ac1d37402559f888c5e4a242c7d20cd"
)
UNITE_CLASSIFIER_SHA256 = (
    "5df2370f89b7b64766c0d969ddc5374d3951fdb8bda598ae54ff15cf86b0fca0"
)
UNITE_CLASSIFIER_TAG = "v10.0-2025-02-19-qiime2-2026.4"
UNITE_CLASSIFIER_COMMIT = "26ca7e07979c230ea9a65565e4781b3538962352"
PR2_SNAPSHOT_SHA256 = (
    "5f8eeee6d8b76bdc951f3620d45c750b0084491cee448553ec0baa944428add4"
)
PR2_SOURCE_SHA256 = (
    "0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--qiime-env", default=QIIME_ENV)
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Reuse existing pipeline artifacts while rebuilding audits and figures.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(
    label: str,
    command: Sequence[str],
    environment: dict[str, str],
    timeout: int = 3600,
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


def require_success(result: dict[str, object]) -> None:
    if int(result["returncode"]) != 0:
        raise RuntimeError(
            f"{result['label']} failed\nSTDOUT:\n{result['stdout']}\n"
            f"STDERR:\n{result['stderr']}"
        )


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
    for key in ("uuid", "type", "format", "data-size"):
        match = re.search(rf"^{key}:\s*(.+)$", text, flags=re.MULTILINE)
        result[key] = match.group(1).strip() if match else ""
    return result


def extract_zip_text(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"{path.name}: expected one {suffix}, found {len(matches)}")
        return archive.read(matches[0]).decode("utf-8", errors="replace")


def normalized_read_id(header: str) -> str:
    token = header.split(maxsplit=1)[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def read_fastq_stream(handle: Iterable[str]) -> dict[str, object]:
    lengths: list[int] = []
    expected_errors: list[float] = []
    read_ids: list[str] = []
    block: list[str] = []
    for raw_line in handle:
        block.append(raw_line.rstrip("\r\n"))
        if len(block) == 4:
            header, sequence, plus, quality = block
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("invalid FASTQ record")
            if len(sequence) != len(quality):
                raise ValueError("FASTQ sequence/quality length mismatch")
            read_ids.append(normalized_read_id(header[1:]))
            lengths.append(len(sequence))
            expected_errors.append(
                sum(10 ** (-(ord(character) - 33) / 10) for character in quality)
            )
            block = []
    if block:
        raise ValueError("truncated FASTQ record")
    if not lengths:
        raise ValueError("empty FASTQ")
    return {
        "records": len(lengths),
        "read_ids": read_ids,
        "lengths": lengths,
        "expected_errors": expected_errors,
    }


def read_fastq_path(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        return read_fastq_stream(handle)


def infer_artifact_sample(member: str) -> tuple[str, str]:
    filename = Path(member).name
    direction = "R1" if re.search(r"(?:_|^)R1(?:_|\\.)", filename) else "R2"
    sample = re.sub(
        r"_\d+_L\d{3}_R[12]_\d+\.fastq\.gz$",
        "",
        filename,
    )
    return sample, direction


def artifact_fastq_records(path: Path) -> dict[tuple[str, str], dict[str, object]]:
    results: dict[tuple[str, str], dict[str, object]] = {}
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if "/data/" in name and name.endswith(".fastq.gz")
        ]
        for member in sorted(members):
            sample, direction = infer_artifact_sample(member)
            with gzip.open(io.BytesIO(archive.read(member)), "rt", encoding="ascii") as handle:
                results[(sample, direction)] = read_fastq_stream(handle)
    return results


def five_number(values: list[int]) -> tuple[float, float, float, float, float]:
    ordered = sorted(values)
    quartiles = statistics.quantiles(ordered, n=4, method="inclusive")
    return (
        float(ordered[0]),
        float(quartiles[0]),
        float(statistics.median(ordered)),
        float(quartiles[2]),
        float(ordered[-1]),
    )


def length_rows(
    stage: str,
    records: dict[tuple[str, str], dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (sample, direction), record in sorted(records.items()):
        minimum, q1, median, q3, maximum = five_number(
            list(record["lengths"])  # type: ignore[arg-type]
        )
        expected_errors = sorted(
            list(record["expected_errors"])  # type: ignore[arg-type]
        )
        rows.append(
            {
                "stage": stage,
                "sample_id": sample,
                "direction": direction,
                "reads": record["records"],
                "length_min": int(minimum),
                "length_q1": round(q1, 2),
                "length_median": round(median, 2),
                "length_q3": round(q3, 2),
                "length_max": int(maximum),
                "expected_errors_median": round(
                    statistics.median(expected_errors), 4
                ),
                "expected_errors_p95": round(
                    expected_errors[max(0, int(0.95 * len(expected_errors)) - 1)],
                    4,
                ),
            }
        )
    return rows


def parse_dada2_stats(path: Path) -> list[dict[str, object]]:
    text = extract_zip_text(path, "/data/stats.tsv")
    rows: list[dict[str, object]] = []
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        sample = str(row["sample-id"])
        if sample.startswith("#"):
            continue
        parsed: dict[str, object] = {"sample_id": sample}
        for key in ("input", "filtered", "denoised", "merged", "non-chimeric"):
            parsed[key] = int(float(row[key]))
        parsed["non_chimeric_percent"] = round(
            float(row["percentage of input non-chimeric"]), 3
        )
        rows.append(parsed)
    return rows


def parse_rep_sequences(path: Path) -> list[tuple[str, str]]:
    text = extract_zip_text(path, "/data/dna-sequences.fasta")
    records: list[tuple[str, str]] = []
    identifier = ""
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if identifier:
                records.append((identifier, "".join(chunks)))
            identifier = line[1:].split(maxsplit=1)[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if identifier:
        records.append((identifier, "".join(chunks)))
    return records


def provenance_actions(path: Path) -> list[dict[str, object]]:
    actions: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as archive:
        members = [
            name
            for name in archive.namelist()
            if "/provenance/" in name and name.endswith("/action/action.yaml")
        ]
        for member in sorted(members):
            text = archive.read(member).decode("utf-8", errors="replace")
            action_match = re.search(
                r"^action:\n(?P<block>.*?)(?=^\S|\Z)",
                text,
                flags=re.MULTILINE | re.DOTALL,
            )
            block = action_match.group("block") if action_match else ""
            kind = re.search(r"^[ \t]+type:\s*(.+)$", block, flags=re.MULTILINE)
            plugin = re.search(r"^[ \t]+plugin:\s*(.+)$", block, flags=re.MULTILINE)
            action = re.search(r"^[ \t]+action:\s*(.+)$", block, flags=re.MULTILINE)
            plugin_value = plugin.group(1).strip() if plugin else ""
            plugin_reference = re.search(
                r"environment:plugins:([^'\"]+)", plugin_value
            )
            if plugin_reference:
                plugin_value = plugin_reference.group(1)
            parameters: dict[str, str] = {}
            for key in (
                "region",
                "taxa",
                "cluster_id",
                "threads",
                "trunc_len_f",
                "trunc_len_r",
                "max_ee_f",
                "max_ee_r",
                "chimera_method",
            ):
                match = re.search(
                    rf"^\s*-\s+{key}:\s*(.+)$", text, flags=re.MULTILINE
                )
                if match:
                    parameters[key] = match.group(1).strip()
            actions.append(
                {
                    "member": member,
                    "kind": kind.group(1).strip() if kind else "unknown",
                    "plugin": plugin_value,
                    "action": action.group(1).strip() if action else "import",
                    "parameters": json.dumps(
                        parameters, ensure_ascii=False, sort_keys=True
                    ),
                }
            )
    return actions


def snapshot_audit(path: Path) -> dict[str, object]:
    records = 0
    depths: Counter[int] = Counter()
    domains: Counter[str] = Counter()
    source_indices: list[int] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            records += 1
            header = line[1:].strip()
            prefix, taxonomy = header.split("|", 1)
            source_indices.append(int(prefix.removeprefix("source_record_")))
            ranks = taxonomy.rstrip(";").split(";")
            depths[len(ranks)] += 1
            domains[ranks[0]] += 1
    return {
        "records": records,
        "depth_counts": dict(sorted(depths.items())),
        "domain_counts": dict(domains.most_common()),
        "first_source_index": min(source_indices),
        "last_source_index": max(source_indices),
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


def draw_marker_branch_map(figure_dir: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    dark, gray = "#111827", "#4B5563"
    green, blue, orange = "#009E73", "#0072B2", "#E69F00"
    fig, ax = plt.subplots(figsize=(12.2, 6.5), facecolor="white")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.text(
        0.35,
        7.65,
        "The marker gene determines the preprocessing branch",
        fontsize=18,
        fontweight="bold",
        color=dark,
    )
    ax.text(
        0.35,
        7.16,
        "Do not route fungal ITS and eukaryotic 18S through one bacterial 16S template.",
        fontsize=10.2,
        color=gray,
    )

    branches = [
        (
            0.45,
            "Fungal ITS1 / ITS2",
            [
                ("Locate variable spacer", "ITSxpress · ITS1/ITS2 · taxa", green),
                ("Denoise variable lengths", "DADA2 · trunc-len 0 if quality permits", green),
                ("Classify fungal ASVs", "UNITE v10.0 · matching scope/date", green),
                ("Audit without forced tree", "Length + chimera + non-target checks", green),
            ],
        ),
        (
            6.25,
            "Eukaryotic 18S SSU",
            [
                ("Remove documented primers", "Cutadapt · no ITSxpress", blue),
                ("Denoise one fixed region", "DADA2 · quality/overlap evidence", blue),
                ("Classify eukaryotic ASVs", "PR2 v5.1.1 · nine ranks", blue),
                ("Filter after inspection", "Plastid + bacteria + off-target audit", orange),
            ],
        ),
    ]
    for start_x, title, nodes in branches:
        ax.text(
            start_x,
            6.55,
            title,
            fontsize=13,
            fontweight="bold",
            color=dark,
        )
        for index, (heading, detail, color) in enumerate(nodes):
            y = 5.55 - index * 1.25
            ax.add_patch(
                FancyBboxPatch(
                    (start_x, y),
                    5.3,
                    0.88,
                    boxstyle="round,pad=0.04,rounding_size=0.12",
                    edgecolor=color,
                    facecolor=color + "12",
                    linewidth=1.25,
                )
            )
            ax.text(
                start_x + 0.2,
                y + 0.57,
                heading,
                fontsize=10.1,
                fontweight="bold",
                color=dark,
            )
            ax.text(start_x + 0.2, y + 0.25, detail, fontsize=8.9, color=gray)
            if index < len(nodes) - 1:
                ax.add_patch(
                    FancyArrowPatch(
                        (start_x + 2.65, y - 0.04),
                        (start_x + 2.65, y - 0.31),
                        arrowstyle="-|>",
                        mutation_scale=12,
                        color=color,
                    )
                )
    fig.subplots_adjust(left=0.02, right=0.99, top=0.98, bottom=0.03)
    save_figure(fig, figure_dir / "13-marker-branch-map")
    plt.close(fig)


def draw_length_audit(
    figure_dir: Path,
    raw_records: dict[tuple[str, str], dict[str, object]],
    trimmed_records: dict[tuple[str, str], dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    colors = {"Raw": "#9CA3AF", "ITSxpress": "#009E73"}
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.9), facecolor="white")
    for axis, direction in zip(axes, ("R1", "R2")):
        for stage, records in (("Raw", raw_records), ("ITSxpress", trimmed_records)):
            lengths: list[int] = []
            for (sample, observed_direction), record in records.items():
                if sample in INPUT_SAMPLES and observed_direction == direction:
                    lengths.extend(record["lengths"])  # type: ignore[arg-type]
            bins = np.arange(0, max(lengths) + 11, 10)
            axis.hist(
                lengths,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=2.0,
                color=colors[stage],
                label=f"{stage} · n={len(lengths):,}",
            )
            axis.axvline(
                statistics.median(lengths),
                color=colors[stage],
                linestyle=":",
                linewidth=1.25,
            )
        axis.set_title(f"{direction} read-length distribution", loc="left", fontweight="bold")
        axis.set_xlabel("Read length (nt)")
        axis.set_ylabel("Density")
        axis.legend(frameon=False)
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.8)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle(
        "ITSxpress converts conserved-flank reads into variable-length ITS2 reads",
        x=0.06,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.06,
        0.92,
        "Vertical dotted lines mark stage-specific medians; all two-sample reads are shown.",
        fontsize=9.8,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.13, wspace=0.25)
    save_figure(fig, figure_dir / "13-itsxpress-length-audit")
    plt.close(fig)


def draw_loss_ledger(
    figure_dir: Path,
    stats_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    stages = ["input", "filtered", "denoised", "merged", "non-chimeric"]
    labels = ["Input", "Filtered", "Denoised", "Merged", "Non-chimeric"]
    colors = ["#9CA3AF", "#56B4E9", "#0072B2", "#E69F00", "#009E73"]
    fig, ax = plt.subplots(figsize=(10.8, 6.0), facecolor="white")
    x = np.arange(len(stages))
    width = 0.34
    maximum = max(int(row[stage]) for row in stats_rows for stage in stages)
    for sample_index, row in enumerate(stats_rows):
        values = [int(row[stage]) for stage in stages]
        offset = (sample_index - 0.5) * width
        ax.bar(
            x + offset,
            values,
            width,
            color=colors,
            alpha=0.72 + sample_index * 0.2,
            edgecolor="#111827",
            linewidth=0.35,
            label=str(row["sample_id"]),
        )
        for index, value in enumerate(values):
            ax.text(
                index + offset,
                value + maximum * 0.025,
                f"{value:,}",
                ha="center",
                va="bottom",
                fontsize=8.1,
            )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Read pairs")
    ax.set_ylim(0, maximum * 1.18)
    ax.set_title(
        "ITS2 denoising keeps a sample-resolved loss ledger",
        loc="left",
        fontsize=17,
        fontweight="bold",
    )
    ax.legend(frameon=False, title="Sample")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.12)
    save_figure(fig, figure_dir / "13-its2-loss-ledger")
    plt.close(fig)


def draw_reference_contract(
    figure_dir: Path,
    pr2_summary: dict[str, object],
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    dark, gray = "#111827", "#4B5563"
    green, blue, orange = "#009E73", "#0072B2", "#E69F00"
    domains = dict(pr2_summary["domain_counts"])
    display = [
        ("Eukaryota", int(domains["Eukaryota"]), blue),
        ("Bacteria", int(domains["Bacteria"]), orange),
        ("Plastid labels", int(domains["Eukaryota:plas"]), green),
        ("Mito labels", int(domains["Eukaryota:mito"]), "#CC79A7"),
        ("Archaea", int(domains["Archaea"]), "#D55E00"),
    ]
    fig = plt.figure(figsize=(12.2, 6.3), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.22)
    ax = fig.add_subplot(grid[0, 0])
    names = [item[0] for item in display]
    values = [item[1] for item in display]
    ax.barh(names, values, color=[item[2] for item in display], height=0.58)
    ax.invert_yaxis()
    ax.set_xlabel("PR2 v5.1.1 SSU records")
    ax.set_title("A  Reference scope is auditable", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    for index, value in enumerate(values):
        ax.text(value + 2500, index, f"{value:,}", va="center", fontsize=9, color=dark)
    ax.set_xlim(0, max(values) * 1.18)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)

    panel = fig.add_subplot(grid[0, 1])
    panel.set_xlim(0, 10)
    panel.set_ylim(0, 10)
    panel.axis("off")
    panel.set_title("B  Match marker, scope, version, and ranks", loc="left", fontweight="bold")
    cards = [
        (
            0.4,
            6.3,
            "Fungal ITS",
            "UNITE v10.0 · 2025-02-19\nRaw release: 8 labels incl. SH\nQ2 classifier removes SH",
            green,
        ),
        (
            5.25,
            6.3,
            "Eukaryotic 18S",
            "PR2 v5.1.1 · 2025-10-27\n240,201 SSU records\n9 taxonomy ranks",
            blue,
        ),
        (
            0.4,
            2.4,
            "Never infer from filename",
            "Verify checksum + semantic type.\nMatch release scope, marker,\nand software compatibility.",
            orange,
        ),
        (
            5.25,
            2.4,
            "Keep off-target evidence",
            "PR2 retains bacterial, plastid,\nmitochondrial, and archaeal\nlabels for off-target review.",
            "#D55E00",
        ),
    ]
    for x, y, title, detail, color in cards:
        panel.add_patch(
            FancyBboxPatch(
                (x, y),
                4.35,
                2.35,
                boxstyle="round,pad=0.05,rounding_size=0.13",
                edgecolor=color,
                facecolor=color + "12",
                linewidth=1.3,
            )
        )
        panel.text(
            x + 0.22,
            y + 1.82,
            title,
            fontsize=10.7,
            fontweight="bold",
            color=dark,
        )
        panel.text(
            x + 0.22,
            y + 1.24,
            detail,
            fontsize=8.9,
            color=gray,
            va="top",
            linespacing=1.35,
        )
    fig.suptitle(
        "UNITE and PR2 answer different marker-gene questions",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color=dark,
    )
    fig.subplots_adjust(left=0.08, right=0.99, top=0.84, bottom=0.12)
    save_figure(fig, figure_dir / "13-reference-database-contract")
    plt.close(fig)


def sanitize_text(text: str, replacements: dict[str, str]) -> str:
    result = text
    for source, replacement in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        if source:
            result = result.replace(source, replacement)
    return result


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    its2_dir = project_root / "data/small/its2"
    reference_dir = project_root / "data/small/its-18s"
    source_summary_path = its2_dir / "source_summary.json"
    pr2_summary_path = reference_dir / "pr2-v5.1.1-source-summary.json"
    pr2_snapshot_path = reference_dir / "pr2-v5.1.1-stratified-snapshot.fasta.gz"
    unite_eval_path = reference_dir / "eval_unite_ver2025-02-19-Q2-2026.4.qzv"
    unite_build_audit_path = (
        reference_dir / "unite-v10-classifier-build-audit.json"
    )
    overlay_path = project_root / "env/qiime2-itsxpress-overlay.yml"

    source_summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
    pr2_summary = json.loads(pr2_summary_path.read_text(encoding="utf-8"))
    unite_build_audit = json.loads(
        unite_build_audit_path.read_text(encoding="utf-8")
    )
    pr2_snapshot = snapshot_audit(pr2_snapshot_path)
    checks: list[dict[str, object]] = []
    commands: list[dict[str, object]] = []

    with tempfile.TemporaryDirectory(prefix="article13-", dir="/tmp") as temporary:
        work_root = Path(temporary)
        environment = isolated_environment(work_root)
        os.environ["MPLCONFIGDIR"] = environment["MPLCONFIGDIR"]

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
        qiime_result = run_command(
            "qiime-info",
            qiime_command(args.qiime_env, "info"),
            environment,
        )
        require_success(qiime_result)
        commands.append(qiime_result)
        vsearch_result = run_command(
            "vsearch-version",
            conda_command(args.qiime_env, "vsearch", "--version"),
            environment,
        )
        require_success(vsearch_result)
        commands.append(vsearch_result)
        hmmer_result = run_command(
            "hmmer-version",
            conda_command(args.qiime_env, "hmmsearch", "-h"),
            environment,
        )
        require_success(hmmer_result)
        commands.append(hmmer_result)

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
                "component": "q2cli",
                "observed_version": packages.get("q2cli", ""),
                "expected_version": QIIME_VERSION,
                "status": "PASS"
                if packages.get("q2cli") == QIIME_VERSION
                else "FAIL",
            },
            {
                "component": "ITSxpress",
                "observed_version": packages.get("itsxpress", ""),
                "expected_version": ITSXPRESS_VERSION,
                "status": "PASS"
                if packages.get("itsxpress") == ITSXPRESS_VERSION
                else "FAIL",
            },
            {
                "component": "Biopython",
                "observed_version": packages.get("biopython", ""),
                "expected_version": BIOPYTHON_VERSION,
                "status": "PASS"
                if packages.get("biopython") == BIOPYTHON_VERSION
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
                "component": "HMMER",
                "observed_version": packages.get("hmmer", ""),
                "expected_version": HMMER_VERSION,
                "status": "PASS"
                if packages.get("hmmer") == HMMER_VERSION
                else "FAIL",
            },
        ]
        write_tsv(
            output_dir / "environment-audit.tsv",
            environment_rows,
            ("component", "observed_version", "expected_version", "status"),
        )
        for row in environment_rows:
            checks.append(
                audit_row(
                    f"environment-{str(row['component']).lower().replace(' ', '-')}",
                    row["observed_version"],
                    row["expected_version"],
                    row["status"] == "PASS",
                )
            )

        input_records: dict[tuple[str, str], dict[str, object]] = {}
        fastq_rows: list[dict[str, object]] = []
        source_by_file = {
            str(row["file"]): row for row in source_summary["files"]
        }
        for sample in INPUT_SAMPLES:
            for direction, suffix in (("R1", "r1"), ("R2", "r2")):
                filename = f"{sample}_{suffix}.fastq.gz"
                path = its2_dir / filename
                record = read_fastq_path(path)
                input_records[(sample, direction)] = record
                observed_sha = sha256_file(path)
                expected_sha = str(source_by_file[filename]["prepared_sha256"])
                fastq_rows.append(
                    {
                        "sample_id": sample,
                        "direction": direction,
                        "file": filename,
                        "records": record["records"],
                        "gzip_bytes": path.stat().st_size,
                        "sha256": observed_sha,
                        "paired_id_mismatches": 0,
                    }
                )
                checks.append(
                    audit_row(
                        f"{sample}-{direction.lower()}-records",
                        record["records"],
                        INPUT_PAIRS_PER_SAMPLE,
                        record["records"] == INPUT_PAIRS_PER_SAMPLE,
                    )
                )
                checks.append(
                    audit_row(
                        f"{sample}-{direction.lower()}-sha256",
                        observed_sha,
                        expected_sha,
                        observed_sha == expected_sha,
                    )
                )
            mismatches = sum(
                left != right
                for left, right in zip(
                    input_records[(sample, "R1")]["read_ids"],
                    input_records[(sample, "R2")]["read_ids"],
                )
            )
            fastq_rows[-2]["paired_id_mismatches"] = mismatches
            fastq_rows[-1]["paired_id_mismatches"] = mismatches
            checks.append(audit_row(f"{sample}-paired-read-ids", mismatches, 0, mismatches == 0))
        write_tsv(
            output_dir / "its2-fastq-audit.tsv",
            fastq_rows,
            (
                "sample_id",
                "direction",
                "file",
                "records",
                "gzip_bytes",
                "sha256",
                "paired_id_mismatches",
            ),
        )

        manifest_path = output_dir / "its2-manifest.tsv"
        manifest_rows: list[dict[str, object]] = []
        for sample in INPUT_SAMPLES:
            manifest_rows.append(
                {
                    "sample-id": sample,
                    "forward-absolute-filepath": str(
                        (its2_dir / f"{sample}_r1.fastq.gz").resolve()
                    ),
                    "reverse-absolute-filepath": str(
                        (its2_dir / f"{sample}_r2.fastq.gz").resolve()
                    ),
                }
            )
        write_tsv(
            manifest_path,
            manifest_rows,
            (
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            ),
        )

        artifacts = {
            "raw": output_dir / "its2-raw.qza",
            "raw_summary": output_dir / "its2-raw-summary.qzv",
            "trimmed": output_dir / "its2-trimmed.qza",
            "trimmed_summary": output_dir / "its2-trimmed-summary.qzv",
            "table": output_dir / "its2-table.qza",
            "rep_seqs": output_dir / "its2-rep-seqs.qza",
            "stats": output_dir / "its2-denoising-stats.qza",
            "base_transitions": output_dir / "its2-base-transitions.qza",
            "table_summary": output_dir / "its2-table-summary.qzv",
            "feature_frequencies": output_dir / "its2-feature-frequencies.qza",
            "sample_frequencies": output_dir / "its2-sample-frequencies.qza",
            "rep_seqs_summary": output_dir / "its2-rep-seqs.qzv",
            "stats_summary": output_dir / "its2-stats.qzv",
        }
        if not args.reuse_existing:
            for path in artifacts.values():
                path.unlink(missing_ok=True)

        pipeline_commands = [
            (
                "import-its2",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "import",
                    "--type",
                    "SampleData[PairedEndSequencesWithQuality]",
                    "--input-path",
                    str(manifest_path),
                    "--input-format",
                    "PairedEndFastqManifestPhred33V2",
                    "--output-path",
                    str(artifacts["raw"]),
                ),
            ),
            (
                "validate-raw",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "validate",
                    str(artifacts["raw"]),
                    "--level",
                    "max",
                ),
            ),
            (
                "summarize-raw",
                qiime_command(
                    args.qiime_env,
                    "demux",
                    "summarize",
                    "--i-data",
                    str(artifacts["raw"]),
                    "--o-visualization",
                    str(artifacts["raw_summary"]),
                ),
            ),
            (
                "itsxpress-trim",
                qiime_command(
                    args.qiime_env,
                    "itsxpress",
                    "trim-pair-output-unmerged",
                    "--i-per-sample-sequences",
                    str(artifacts["raw"]),
                    "--p-region",
                    "ITS2",
                    "--p-taxa",
                    "F",
                    "--p-cluster-id",
                    "1.0",
                    "--p-threads",
                    "1",
                    "--o-trimmed",
                    str(artifacts["trimmed"]),
                    "--verbose",
                ),
            ),
            (
                "validate-trimmed",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "validate",
                    str(artifacts["trimmed"]),
                    "--level",
                    "max",
                ),
            ),
            (
                "summarize-trimmed",
                qiime_command(
                    args.qiime_env,
                    "demux",
                    "summarize",
                    "--i-data",
                    str(artifacts["trimmed"]),
                    "--o-visualization",
                    str(artifacts["trimmed_summary"]),
                ),
            ),
            (
                "dada2-its2",
                qiime_command(
                    args.qiime_env,
                    "dada2",
                    "denoise-paired",
                    "--i-demultiplexed-seqs",
                    str(artifacts["trimmed"]),
                    "--p-trunc-len-f",
                    "0",
                    "--p-trunc-len-r",
                    "0",
                    "--p-trim-left-f",
                    "0",
                    "--p-trim-left-r",
                    "0",
                    "--p-max-ee-f",
                    "2",
                    "--p-max-ee-r",
                    "2",
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
                    str(artifacts["table"]),
                    "--o-representative-sequences",
                    str(artifacts["rep_seqs"]),
                    "--o-denoising-stats",
                    str(artifacts["stats"]),
                    "--o-base-transition-stats",
                    str(artifacts["base_transitions"]),
                    "--verbose",
                ),
            ),
            (
                "summarize-table",
                qiime_command(
                    args.qiime_env,
                    "feature-table",
                    "summarize",
                    "--i-table",
                    str(artifacts["table"]),
                    "--m-metadata-file",
                    str(its2_dir / "metadata.tsv"),
                    "--o-feature-frequencies",
                    str(artifacts["feature_frequencies"]),
                    "--o-sample-frequencies",
                    str(artifacts["sample_frequencies"]),
                    "--o-summary",
                    str(artifacts["table_summary"]),
                ),
            ),
            (
                "tabulate-rep-seqs",
                qiime_command(
                    args.qiime_env,
                    "feature-table",
                    "tabulate-seqs",
                    "--i-data",
                    str(artifacts["rep_seqs"]),
                    "--o-visualization",
                    str(artifacts["rep_seqs_summary"]),
                ),
            ),
            (
                "tabulate-stats",
                qiime_command(
                    args.qiime_env,
                    "metadata",
                    "tabulate",
                    "--m-input-file",
                    str(artifacts["stats"]),
                    "--o-visualization",
                    str(artifacts["stats_summary"]),
                ),
            ),
            (
                "validate-unite-evaluation",
                qiime_command(
                    args.qiime_env,
                    "tools",
                    "validate",
                    str(unite_eval_path),
                    "--level",
                    "max",
                ),
            ),
        ]
        reusable_outputs = {
            "import-its2": artifacts["raw"],
            "summarize-raw": artifacts["raw_summary"],
            "itsxpress-trim": artifacts["trimmed"],
            "summarize-trimmed": artifacts["trimmed_summary"],
            "dada2-its2": artifacts["stats"],
            "summarize-table": artifacts["table_summary"],
            "tabulate-rep-seqs": artifacts["rep_seqs_summary"],
            "tabulate-stats": artifacts["stats_summary"],
        }
        for label, command in pipeline_commands:
            reusable = reusable_outputs.get(label)
            if args.reuse_existing and reusable is not None and reusable.exists():
                result = {
                    "label": label,
                    "command": list(command),
                    "returncode": 0,
                    "stdout": f"Reused existing {reusable.name}",
                    "stderr": "",
                    "elapsed_seconds": 0.0,
                }
            else:
                result = run_command(label, command, environment)
            require_success(result)
            commands.append(result)

        trimmed_records = artifact_fastq_records(artifacts["trimmed"])
        raw_artifact_records = artifact_fastq_records(artifacts["raw"])
        length_audit = length_rows("Raw", raw_artifact_records) + length_rows(
            "ITSxpress", trimmed_records
        )
        write_tsv(
            output_dir / "its2-length-audit.tsv",
            length_audit,
            (
                "stage",
                "sample_id",
                "direction",
                "reads",
                "length_min",
                "length_q1",
                "length_median",
                "length_q3",
                "length_max",
                "expected_errors_median",
                "expected_errors_p95",
            ),
        )

        retention_rows: list[dict[str, object]] = []
        for sample in INPUT_SAMPLES:
            raw_f = int(raw_artifact_records[(sample, "R1")]["records"])
            raw_r = int(raw_artifact_records[(sample, "R2")]["records"])
            trim_f = int(trimmed_records[(sample, "R1")]["records"])
            trim_r = int(trimmed_records[(sample, "R2")]["records"])
            retention_rows.append(
                {
                    "sample_id": sample,
                    "raw_forward": raw_f,
                    "raw_reverse": raw_r,
                    "trimmed_forward": trim_f,
                    "trimmed_reverse": trim_r,
                    "retained_pairs": min(trim_f, trim_r),
                    "retention_percent": round(100 * min(trim_f, trim_r) / raw_f, 3),
                }
            )
        write_tsv(
            output_dir / "its2-retention.tsv",
            retention_rows,
            (
                "sample_id",
                "raw_forward",
                "raw_reverse",
                "trimmed_forward",
                "trimmed_reverse",
                "retained_pairs",
                "retention_percent",
            ),
        )

        stats_rows = parse_dada2_stats(artifacts["stats"])
        write_tsv(
            output_dir / "its2-denoising-ledger.tsv",
            stats_rows,
            (
                "sample_id",
                "input",
                "filtered",
                "denoised",
                "merged",
                "non-chimeric",
                "non_chimeric_percent",
            ),
        )
        rep_sequences = parse_rep_sequences(artifacts["rep_seqs"])
        asv_length_counts = Counter(len(sequence) for _, sequence in rep_sequences)
        asv_length_rows = [
            {"length_nt": length, "asv_count": count}
            for length, count in sorted(asv_length_counts.items())
        ]
        write_tsv(
            output_dir / "its2-asv-length-distribution.tsv",
            asv_length_rows,
            ("length_nt", "asv_count"),
        )

        marker_rows = [
            {
                "marker": "Fungal ITS1/ITS2",
                "target": "Variable internal transcribed spacer",
                "preprocessing": "ITSxpress after primer-aware import",
                "denoising": "DADA2 on unmerged trimmed pairs",
                "reference": "UNITE v10.0 2025-02-19",
                "tree_policy": "Do not force one global alignment/tree",
            },
            {
                "marker": "Eukaryotic 18S",
                "target": "SSU rRNA amplicon",
                "preprocessing": "Documented primer removal; no ITSxpress",
                "denoising": "Standard DADA2 with quality/overlap evidence",
                "reference": "PR2 v5.1.1 2025-10-27",
                "tree_policy": "Region-specific alignment after off-target audit",
            },
        ]
        write_tsv(
            output_dir / "marker-branch-decision.tsv",
            marker_rows,
            (
                "marker",
                "target",
                "preprocessing",
                "denoising",
                "reference",
                "tree_policy",
            ),
        )

        reference_rows = [
            {
                "database": "UNITE",
                "asset": "Raw QIIME fungi release",
                "marker_scope": "Fungal ITS",
                "release": "v10.0 · 2025-02-19",
                "record_summary": "20,295 RefS; 81,842 RepS",
                "taxonomy_contract": "8 labels; raw taxonomy includes sh__",
                "compatibility": "Source taxonomy for classifier training",
                "identifier_or_sha256": "doi:10.15156/BIO/3301241",
            },
            {
                "database": "UNITE",
                "asset": "Pre-trained dynamic fungi classifier",
                "marker_scope": "Fungal ITS",
                "release": "v10.0 · 2025-02-19",
                "record_summary": "Classifier artifact",
                "taxonomy_contract": "7 biological rank labels; sh__ removed before fitting",
                "compatibility": "QIIME 2 2026.4",
                "identifier_or_sha256": UNITE_CLASSIFIER_SHA256,
            },
            {
                "database": "PR2",
                "asset": "Full SSU DADA2 FASTA",
                "marker_scope": "Eukaryotic 18S / SSU",
                "release": "v5.1.1 · 2025-10-27",
                "record_summary": (
                    f"{int(pr2_summary['record_count']):,} records; "
                    f"{int(pr2_summary['unique_taxonomy_strings']):,} "
                    "unique taxonomy strings"
                ),
                "taxonomy_contract": "9 ranks; no empty rank fields",
                "compatibility": "DADA2 FASTA; import/train as required",
                "identifier_or_sha256": str(pr2_summary["source_sha256"]),
            },
        ]
        write_tsv(
            output_dir / "reference-database-audit.tsv",
            reference_rows,
            (
                "database",
                "asset",
                "marker_scope",
                "release",
                "record_summary",
                "taxonomy_contract",
                "compatibility",
                "identifier_or_sha256",
            ),
        )
        pr2_domain_rows = [
            {"domain_label": label, "records": count}
            for label, count in dict(pr2_summary["domain_counts"]).items()
        ]
        write_tsv(
            output_dir / "pr2-domain-audit.tsv",
            pr2_domain_rows,
            ("domain_label", "records"),
        )
        pr2_depth_rows = [
            {"taxonomy_ranks": depth, "records": count}
            for depth, count in dict(pr2_summary["taxonomy_depth_counts"]).items()
        ]
        write_tsv(
            output_dir / "pr2-taxonomy-depth-audit.tsv",
            pr2_depth_rows,
            ("taxonomy_ranks", "records"),
        )

        artifact_expected_types = {
            "raw": "SampleData[PairedEndSequencesWithQuality]",
            "trimmed": "SampleData[PairedEndSequencesWithQuality]",
            "table": "FeatureTable[Frequency]",
            "rep_seqs": "FeatureData[Sequence]",
            "stats": "SampleData[DADA2Stats]",
            "base_transitions": "DADA2BaseTransitionStats",
            "feature_frequencies": "ImmutableMetadata",
            "sample_frequencies": "ImmutableMetadata",
        }
        artifact_rows: list[dict[str, object]] = []
        for key, expected_type in artifact_expected_types.items():
            metadata = artifact_metadata(artifacts[key])
            passed = metadata["type"] == expected_type
            artifact_rows.append(
                {
                    "artifact": artifacts[key].name,
                    "semantic_type": metadata["type"],
                    "expected_type": expected_type,
                    "format": metadata["format"],
                    "status": "PASS" if passed else "FAIL",
                }
            )
            checks.append(
                audit_row(
                    f"artifact-{key}-type",
                    metadata["type"],
                    expected_type,
                    passed,
                )
            )
        write_tsv(
            output_dir / "artifact-audit.tsv",
            artifact_rows,
            ("artifact", "semantic_type", "expected_type", "format", "status"),
        )

        action_rows = provenance_actions(artifacts["table"])
        write_tsv(
            output_dir / "provenance-audit.tsv",
            action_rows,
            ("member", "kind", "plugin", "action", "parameters"),
        )
        action_pairs = {
            (str(row["plugin"]), str(row["action"])) for row in action_rows
        }

        for row in retention_rows:
            sample = str(row["sample_id"])
            checks.append(
                audit_row(
                    f"{sample}-itsxpress-pair-sync",
                    row["trimmed_forward"],
                    row["trimmed_reverse"],
                    row["trimmed_forward"] == row["trimmed_reverse"],
                )
            )
            checks.append(
                audit_row(
                    f"{sample}-itsxpress-retains-reads",
                    row["retained_pairs"],
                    "> 0",
                    int(row["retained_pairs"]) > 0,
                )
            )
        raw_medians = [
            float(row["length_median"])
            for row in length_audit
            if row["stage"] == "Raw"
        ]
        trim_medians = [
            float(row["length_median"])
            for row in length_audit
            if row["stage"] == "ITSxpress"
        ]
        checks.append(
            audit_row(
                "itsxpress-changes-read-geometry",
                round(sum(trim_medians), 2),
                f"!= {round(sum(raw_medians), 2)}",
                trim_medians != raw_medians,
            )
        )
        for row in stats_rows:
            sample = str(row["sample_id"])
            stages = [
                int(row[key])
                for key in ("input", "filtered", "denoised", "merged", "non-chimeric")
            ]
            checks.append(
                audit_row(
                    f"{sample}-dada2-monotonic-ledger",
                    stages,
                    "non-increasing",
                    all(left >= right for left, right in zip(stages, stages[1:])),
                )
            )
            checks.append(
                audit_row(
                    f"{sample}-dada2-nonzero-output",
                    row["non-chimeric"],
                    "> 0",
                    int(row["non-chimeric"]) > 0,
                )
            )
        checks.extend(
            [
                audit_row(
                    "asv-count-positive",
                    len(rep_sequences),
                    "> 0",
                    len(rep_sequences) > 0,
                ),
                audit_row(
                    "asv-lengths-variable",
                    len(asv_length_counts),
                    "> 1",
                    len(asv_length_counts) > 1,
                ),
                audit_row(
                    "provenance-itsxpress",
                    sorted(action_pairs),
                    "itsxpress trim_pair_output_unmerged",
                    any(
                        plugin == "itsxpress"
                        and "trim_pair_output_unmerged" in action
                        for plugin, action in action_pairs
                    ),
                ),
                audit_row(
                    "provenance-dada2",
                    sorted(action_pairs),
                    "dada2 denoise_paired",
                    ("dada2", "denoise_paired") in action_pairs,
                ),
                audit_row(
                    "source-total-pairs",
                    source_summary["observed_total_pairs"],
                    INPUT_TOTAL_PAIRS,
                    int(source_summary["observed_total_pairs"]) == INPUT_TOTAL_PAIRS,
                ),
                audit_row(
                    "source-claim-discrepancy-recorded",
                    source_summary["source_documentation_claim"],
                    "claim retained and actual count audited",
                    "10,000" in str(source_summary["source_documentation_claim"]),
                ),
                audit_row(
                    "overlay-exists",
                    overlay_path.exists(),
                    True,
                    overlay_path.exists(),
                ),
                audit_row(
                    "unite-evaluation-sha256",
                    sha256_file(unite_eval_path),
                    UNITE_EVAL_SHA256,
                    sha256_file(unite_eval_path) == UNITE_EVAL_SHA256,
                ),
                audit_row(
                    "unite-evaluation-type",
                    artifact_metadata(unite_eval_path)["type"],
                    "Visualization",
                    artifact_metadata(unite_eval_path)["type"] == "Visualization",
                ),
                audit_row(
                    "unite-classifier-release-tag",
                    unite_build_audit["classifier_release"]["tag"],
                    UNITE_CLASSIFIER_TAG,
                    unite_build_audit["classifier_release"]["tag"]
                    == UNITE_CLASSIFIER_TAG,
                ),
                audit_row(
                    "unite-classifier-release-commit",
                    unite_build_audit["classifier_release"]["commit"],
                    UNITE_CLASSIFIER_COMMIT,
                    unite_build_audit["classifier_release"]["commit"]
                    == UNITE_CLASSIFIER_COMMIT,
                ),
                audit_row(
                    "unite-raw-taxonomy-has-eight-labels-including-sh",
                    unite_build_audit["raw_qiime_fungi_release"][
                        "taxonomy_labels"
                    ],
                    "8 labels including sh__",
                    len(
                        unite_build_audit["raw_qiime_fungi_release"][
                            "taxonomy_labels"
                        ]
                    )
                    == 8
                    and "sh__"
                    in unite_build_audit["raw_qiime_fungi_release"][
                        "taxonomy_labels"
                    ],
                ),
                audit_row(
                    "unite-classifier-removes-sh-before-fitting",
                    unite_build_audit["classifier_release"]["taxonomy_edit"][
                        "search_strings"
                    ],
                    ";sh__.*",
                    unite_build_audit["classifier_release"]["taxonomy_edit"][
                        "search_strings"
                    ]
                    == ";sh__.*"
                    and unite_build_audit["classifier_release"]["taxonomy_edit"][
                        "replacement_strings"
                    ]
                    == "",
                ),
                audit_row(
                    "unite-classifier-sha256",
                    unite_build_audit["classifier_release"]["asset_sha256"],
                    UNITE_CLASSIFIER_SHA256,
                    unite_build_audit["classifier_release"]["asset_sha256"]
                    == UNITE_CLASSIFIER_SHA256,
                ),
                audit_row(
                    "pr2-source-sha256",
                    pr2_summary["source_sha256"],
                    PR2_SOURCE_SHA256,
                    pr2_summary["source_sha256"] == PR2_SOURCE_SHA256,
                ),
                audit_row(
                    "pr2-record-count",
                    pr2_summary["record_count"],
                    240_201,
                    int(pr2_summary["record_count"]) == 240_201,
                ),
                audit_row(
                    "pr2-nine-rank-schema",
                    pr2_summary["taxonomy_depth_counts"],
                    {"9": 240201},
                    pr2_summary["taxonomy_depth_counts"] == {"9": 240201},
                ),
                audit_row(
                    "pr2-no-empty-ranks",
                    pr2_summary["records_with_empty_rank"],
                    0,
                    int(pr2_summary["records_with_empty_rank"]) == 0,
                ),
                audit_row(
                    "pr2-snapshot-sha256",
                    sha256_file(pr2_snapshot_path),
                    PR2_SNAPSHOT_SHA256,
                    sha256_file(pr2_snapshot_path) == PR2_SNAPSHOT_SHA256,
                ),
                audit_row(
                    "pr2-snapshot-records",
                    pr2_snapshot["records"],
                    1002,
                    int(pr2_snapshot["records"]) == 1002,
                ),
                audit_row(
                    "pr2-snapshot-rank-depth",
                    pr2_snapshot["depth_counts"],
                    {9: 1002},
                    pr2_snapshot["depth_counts"] == {9: 1002},
                ),
                audit_row(
                    "pr2-off-target-labels-audited",
                    sorted(dict(pr2_summary["domain_counts"])),
                    "Bacteria + plastid + mitochondria + Archaea",
                    all(
                        label in dict(pr2_summary["domain_counts"])
                        for label in (
                            "Bacteria",
                            "Eukaryota:plas",
                            "Eukaryota:mito",
                            "Archaea",
                        )
                    ),
                ),
            ]
        )

        draw_marker_branch_map(figure_dir)
        draw_length_audit(figure_dir, raw_artifact_records, trimmed_records)
        draw_loss_ledger(figure_dir, stats_rows)
        draw_reference_contract(figure_dir, pr2_summary)

        total_input = sum(int(row["input"]) for row in stats_rows)
        total_non_chimeric = sum(int(row["non-chimeric"]) for row in stats_rows)
        summary = {
            "status": "passed"
            if all(row["status"] == "PASS" for row in checks)
            else "failed",
            "checks_total": len(checks),
            "checks_passed": sum(row["status"] == "PASS" for row in checks),
            "checks_failed": sum(row["status"] == "FAIL" for row in checks),
            "qiime_version": packages.get("qiime2", ""),
            "itsxpress_version": packages.get("itsxpress", ""),
            "input_samples": len(INPUT_SAMPLES),
            "input_read_pairs": INPUT_TOTAL_PAIRS,
            "itsxpress_retained_pairs": sum(
                int(row["retained_pairs"]) for row in retention_rows
            ),
            "dada2_input_pairs": total_input,
            "dada2_non_chimeric_pairs": total_non_chimeric,
            "dada2_non_chimeric_percent": round(
                100 * total_non_chimeric / total_input, 3
            ),
            "asv_count": len(rep_sequences),
            "asv_unique_lengths": len(asv_length_counts),
            "asv_length_min": min(asv_length_counts),
            "asv_length_median": statistics.median(
                [len(sequence) for _, sequence in rep_sequences]
            ),
            "asv_length_max": max(asv_length_counts),
            "pr2_records": int(pr2_summary["record_count"]),
            "pr2_taxonomy_ranks": 9,
            "pr2_snapshot_records": int(pr2_snapshot["records"]),
            "unite_version": "10.0",
            "unite_release_date": "2025-02-19",
            "unite_dynamic_fungi_qiime2_2026_4_sha256": (
                UNITE_CLASSIFIER_SHA256
            ),
            "unite_raw_taxonomy_labels": 8,
            "unite_classifier_taxonomy_labels": 7,
            "unite_classifier_species_hypothesis_removed": True,
        }
        (output_dir / "its-18s-summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
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
        for result in commands:
            command_text = " ".join(str(token) for token in result["command"])
            log_parts.extend(
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
        (output_dir / "validation.log").write_text(
            sanitize_text("\n".join(log_parts), replacements),
            encoding="utf-8",
        )
        itsxpress_result = next(
            result for result in commands if result["label"] == "itsxpress-trim"
        )
        (output_dir / "qiime-itsxpress.log").write_text(
            sanitize_text(
                f"stdout:\n{itsxpress_result['stdout']}\n\n"
                f"stderr:\n{itsxpress_result['stderr']}\n",
                replacements,
            ),
            encoding="utf-8",
        )

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
