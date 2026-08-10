#!/usr/bin/env python3
"""Validate Articles 46–50, native multi-omics outputs, and release assets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

import yaml
from PIL import Image


CHAPTERS = {
    46: ("46-procrustes-mantel-halla", (
        "46-integration/46-1-procrustes",
        "46-integration/46-2-global-tests",
        "46-integration/46-3-confounding-audit",
        "46-integration/46-4-halla-associations",
    )),
    47: ("47-spls-diablo", (
        "47-spls-diablo/47-1-spls-generalization",
        "47-spls-diablo/47-2-diablo-scores",
        "47-spls-diablo/47-3-test-performance",
        "47-spls-diablo/47-4-feature-stability",
    )),
    48: ("48-mmvec-mofa", (
        "48-mmvec-mofa/48-1-cca-overfit-audit",
        "48-mmvec-mofa/48-2-mmvec-loss",
        "48-mmvec-mofa/48-3-mmvec-top-pairs",
        "48-mmvec-mofa/48-4-mofa-summary",
    )),
    49: ("49-multi-kingdom", (
        "49-multi-kingdom/49-1-three-kingdom-pcoa",
        "49-multi-kingdom/49-2-procrustes-audit",
        "49-multi-kingdom/49-3-shared-niche-audit",
        "49-multi-kingdom/49-4-stable-candidates",
    )),
    50: ("50-source-tracking", (
        "50-source-tracking/50-1-feast-main",
        "50-source-tracking/50-2-method-comparison",
        "50-source-tracking/50-3-loo-calibration",
        "50-source-tracking/50-4-missing-source",
    )),
}

REQUIRED_SECTION_PATTERNS = (
    ("opening", r"\{#sec-(?:target|paper-figure)\}"),
    ("theory", r"\{#sec-theory\}"),
    ("data-environment", r"\{#sec-(?:setup|preparation)\}"),
    ("analysis", r"\{#sec-code\}"),
    ("limitations", r"\{#sec-audit\}"),
    ("presentation", r"\{#sec-(?:publication|polish|beautify|visualization|figure)\}"),
    ("pitfalls", r"常见.*(?:坑|误判)"),
    ("methods", r"\{#sec-methods\}"),
    ("transfer", r"\{#sec-own-data\}"),
    ("references", r"\{#sec-references\}"),
)

EXECUTABLE_TOKENS = {
    46: (
        "data/small/paired-ibd-multiomics", "protest(", "mantel(",
        "from halla import HAllA", "pairwise-group-residual.tsv",
        "46-4-halla-associations",
    ),
    47: (
        "data/small/paired-ibd-multiomics", "spls(",
        "tune.block.splsda(", "block.splsda(", "test-predictions.tsv",
        "47-4-feature-stability",
    ),
    48: (
        "data/small/paired-ibd-multiomics", "from mmvec.multimodal import MMvec",
        "from mofapy2.run.entry_point import entry_point", "ValidationMAE",
        "mmvec-native", "multiomics-native", "48-4-mofa-summary",
    ),
    49: (
        "data/small/multi-kingdom-duran", "fungi-otutab.tsv",
        "oomycete-otutab.tsv", "adonis2(", "protest(",
        "cross-kingdom-associations.tsv", "49-4-stable-candidates",
    ),
    50: (
        "data/small/source-tracking-feast", "library(FEAST)", "FEAST(",
        "from sourcetracker._sourcetracker import _gibbs",
        "feast-missing-source.tsv", "50-4-missing-source",
    ),
}

PACKAGE_VERSIONS = {
    "mixOmics": "6.26.0",
    "FEAST": "0.1.0",
    "XICOR": "0.4.1",
    "eva": "0.2.7",
    "psychTools": "2.6.4",
    "corpcor": "1.6.10",
    "cowplot": "1.1.3",
    "ellipse": "0.5.0",
    "ggthemes": "5.1.0",
    "mnormt": "2.1.1",
    "psych": "2.4.3",
    "rARPACK": "0.11-0",
}

CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--site-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/articles-46-50-validation.json"),
    )
    return parser.parse_args()


def read_frontmatter(text: str) -> dict[str, Any]:
    if not text.startswith("---\n"):
        return {}
    end = text.find("\n---\n", 4)
    if end < 0:
        return {}
    value = yaml.safe_load(text[4:end])
    return value if isinstance(value, dict) else {}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_shape(path: Path) -> tuple[int, int]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        return sum(1 for _ in reader), len(header) - 1


def table_sample_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))[1:]


def close(observed: float, expected: float, tolerance: float = 1e-6) -> bool:
    return math.isclose(observed, expected, rel_tol=tolerance, abs_tol=tolerance)


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    site_root = args.site_root.resolve()
    checks: list[dict[str, Any]] = []

    def check(identifier: str, passed: bool, detail: Any) -> None:
        checks.append({
            "id": identifier,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        })

    tutorial = yaml.safe_load((project / "tutorial.yaml").read_text(encoding="utf-8"))
    manifest_titles = {
        int(row["number"]): row["title"]
        for row in tutorial["series"]["chapters"]
    }

    figure_files = 0
    for number, (slug, stems) in CHAPTERS.items():
        qmd = project / "chapters" / f"{slug}.qmd"
        qmd_text = qmd.read_text(encoding="utf-8") if qmd.is_file() else ""
        metadata = read_frontmatter(qmd_text)
        execute = metadata.get("execute", {})
        check(f"article-{number}-qmd", qmd.is_file(), str(qmd))
        check(
            f"article-{number}-eval",
            isinstance(execute, dict) and execute.get("eval") is True,
            execute,
        )
        check(
            f"article-{number}-freeze",
            isinstance(execute, dict) and execute.get("freeze") == "auto",
            execute,
        )
        check(f"article-{number}-not-draft", metadata.get("draft") is not True, metadata.get("draft"))
        check(f"article-{number}-seed", "set.seed(" in qmd_text, "set.seed present")
        check(
            f"article-{number}-manifest-title",
            manifest_titles[number] in str(metadata.get("title", "")),
            {"manifest": manifest_titles[number], "qmd": metadata.get("title")},
        )
        check(
            f"article-{number}-inline-theme",
            all(token in qmd_text for token in (
                "pal_pub <-", "theme_pub <- function", "save_pub <- function",
            )),
            "pal_pub/theme_pub/save_pub inline",
        )
        check(
            f"article-{number}-no-source-theme",
            re.search(r'(?m)^[ \t]*source\("R/theme_pub\.R"\)', qmd_text) is None,
            "no executable external theme dependency",
        )
        for section, pattern in REQUIRED_SECTION_PATTERNS:
            check(
                f"article-{number}-section-{hashlib.sha1(section.encode()).hexdigest()[:8]}",
                re.search(pattern, qmd_text) is not None,
                pattern,
            )
        for token in EXECUTABLE_TOKENS[number]:
            check(
                f"article-{number}-token-{hashlib.sha1(token.encode()).hexdigest()[:8]}",
                token in qmd_text,
                token,
            )

        html = site_root / "_site" / "chapters" / f"{slug}.html"
        html_ok = html.is_file() and html.stat().st_size > 50_000
        check(f"article-{number}-html", html_ok, str(html))
        if html_ok:
            html_text = html.read_text(encoding="utf-8")
            check(
                f"article-{number}-html-title",
                manifest_titles[number] in html_text,
                manifest_titles[number],
            )

        for stem in stems:
            for extension in ("pdf", "svg", "png", "tiff"):
                figure_files += 1
                path = site_root / "figures" / f"{stem}.{extension}"
                exists = path.is_file() and path.stat().st_size > 1_000
                check(f"figure-{stem}-{extension}", exists, str(path))
                if not exists:
                    continue
                if extension == "pdf":
                    check(f"signature-{stem}-pdf", path.read_bytes()[:4] == b"%PDF", "PDF signature")
                elif extension == "svg":
                    svg = path.read_text(encoding="utf-8")
                    check(f"signature-{stem}-svg", "<svg" in svg, "SVG root")
                    check(f"english-{stem}-svg", CJK_RE.search(svg) is None, "no CJK plot text")
                else:
                    with Image.open(path) as image:
                        dpi = image.info.get("dpi", (0, 0))
                        dpi_x = float(dpi[0]) if isinstance(dpi, tuple) else float(dpi)
                        check(f"raster-{stem}-{extension}-dpi", dpi_x >= 590, dpi)
                        check(
                            f"raster-{stem}-{extension}-pixels",
                            image.width >= 2_000 and image.height >= 1_400,
                            [image.width, image.height],
                        )
                        if extension == "tiff":
                            compression = str(image.info.get("compression", "")).lower()
                            check(f"raster-{stem}-tiff-lzw", "lzw" in compression, compression)

    check("figure-file-count", figure_files == 80, figure_files)

    paired = project / "data" / "small" / "paired-ibd-multiomics"
    paired_summary = json.loads((paired / "source-summary.json").read_text(encoding="utf-8"))
    check("paired-source-commit", paired_summary["commit"] == "89a519d8c832008fbc6e650453e83e2f04858d02", paired_summary["commit"])
    check("paired-otutab-shape", table_shape(paired / "otutab.tsv") == (250, 220), table_shape(paired / "otutab.tsv"))
    check("paired-taxonomy-shape", table_shape(paired / "taxonomy.tsv") == (250, 9), table_shape(paired / "taxonomy.tsv"))
    check("paired-metabolites-shape", table_shape(paired / "metabolites.tsv") == (277, 220), table_shape(paired / "metabolites.tsv"))
    check("paired-metabolite-annotation", table_shape(paired / "metabolite-annotation.tsv") == (277, 11), table_shape(paired / "metabolite-annotation.tsv"))
    paired_metadata = read_tsv(paired / "metadata.tsv")
    paired_groups = Counter(row["StudyGroup"] for row in paired_metadata)
    check("paired-groups", paired_groups == {"CD": 88, "UC": 76, "Control": 56}, paired_groups)
    paired_ids = [row["SampleID"] for row in paired_metadata]
    check("paired-sample-order-microbes", table_sample_ids(paired / "otutab.tsv") == paired_ids, "exact order")
    check("paired-sample-order-metabolites", table_sample_ids(paired / "metabolites.tsv") == paired_ids, "exact order")

    kingdoms = project / "data" / "small" / "multi-kingdom-duran"
    kingdom_summary = json.loads((kingdoms / "source-summary.json").read_text(encoding="utf-8"))
    check("kingdom-source-commit", kingdom_summary["commit"] == "6db5e85cc5d442fd95fcdcb7250b72fa9e2ff900", kingdom_summary["commit"])
    for label, filename, shape in (
        ("bacteria", "otutab.tsv", (772, 36)),
        ("fungi", "fungi-otutab.tsv", (1063, 36)),
        ("oomycete", "oomycete-otutab.tsv", (219, 36)),
    ):
        observed = table_shape(kingdoms / filename)
        check(f"kingdom-{label}-shape", observed == shape, observed)
    kingdom_metadata = read_tsv(kingdoms / "metadata.tsv")
    check("kingdom-soil-balance", Counter(row["Soil"] for row in kingdom_metadata) == {"GE": 12, "PU": 12, "SD": 12}, Counter(row["Soil"] for row in kingdom_metadata))
    check("kingdom-compartment-balance", Counter(row["Compartment"] for row in kingdom_metadata) == {"Rhizosphere": 12, "Root": 12, "Soil": 12}, Counter(row["Compartment"] for row in kingdom_metadata))

    tracking = project / "data" / "small" / "source-tracking-feast"
    tracking_summary = json.loads((tracking / "source-summary.json").read_text(encoding="utf-8"))
    check("tracking-source-commit", tracking_summary["commit"] == "2f8f3df8051e0e08341f597a9f4693bfb76b3bf6", tracking_summary["commit"])
    check("tracking-otutab-shape", table_shape(tracking / "otutab.tsv") == (1839, 10), table_shape(tracking / "otutab.tsv"))
    tracking_metadata = read_tsv(tracking / "metadata.tsv")
    check("tracking-source-sink-count", Counter(row["SourceSink"] for row in tracking_metadata) == {"Source": 9, "Sink": 1}, Counter(row["SourceSink"] for row in tracking_metadata))

    lock = json.loads((project / "env" / "renv.lock").read_text(encoding="utf-8"))
    packages = lock.get("Packages", {})
    check("renv-package-count", len(packages) == 418, len(packages))
    for package, version in PACKAGE_VERSIONS.items():
        observed = packages.get(package, {}).get("Version")
        check(f"renv-{package}", observed == version, observed)
    check(
        "renv-FEAST-source-commit",
        packages.get("FEAST", {}).get("RemoteSha")
        == "2f8f3df8051e0e08341f597a9f4693bfb76b3bf6",
        packages.get("FEAST", {}).get("RemoteSha"),
    )

    mmvec_env_text = (project / "env" / "mmvec-native.yml").read_text(encoding="utf-8")
    modern_env_text = (project / "env" / "multiomics.yml").read_text(encoding="utf-8")
    check("mmvec-env-native-name", "name: mmvec-native" in mmvec_env_text, "mmvec-native")
    for token in ("python=3.7.16", "tensorflow=1.15.0", "mmvec=1.0.5"):
        check(f"mmvec-env-{token}", token in mmvec_env_text, token)
    check("multiomics-env-native-name", "name: multiomics-native" in modern_env_text, "multiomics-native")
    for token in ("HAllA==0.8.40", "mofapy2==0.7.4", "scikit-learn=1.7.2"):
        check(f"multiomics-env-{token}", token in modern_env_text, token)
    native_scope = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (
            project / "chapters" / "48-mmvec-mofa.qmd",
            project / "scripts" / "run_article48_native_mmvec.py",
            project / "scripts" / "run_article48_mmvec_mofa.py",
            project / "env" / "mmvec-native.yml",
            project / "env" / "multiomics.yml",
        )
    )
    check(
        "article-48-native-api-contract",
        "from mmvec.multimodal import MMvec" in native_scope
        and "from mofapy2.run.entry_point import entry_point" in native_scope,
        "native MMvec + mofapy2 entry points",
    )

    result_root = site_root / "results"
    summary46 = json.loads((result_root / "46-integration" / "summary.json").read_text(encoding="utf-8"))
    check("result-46-status", summary46["status"] == "passed", summary46["status"])
    check("result-46-pairs", summary46["tested_pairs"] == 1000 and summary46["group_residual_pairwise_fdr_hits"] == 350, summary46)
    check("result-46-halla", summary46["halla_significant_pairs"] == 350 and summary46["halla_significant_blocks"] == 70, summary46)
    global_tests = {(row["Model"], row["Method"]): float(row["Statistic"]) for row in read_tsv(result_root / "46-integration" / "global-concordance.tsv")}
    expected_global_tests = {
        ("Observed", "Procrustes"): 0.588107476116809,
        ("Observed", "Mantel"): 0.597317089319856,
        ("Group-adjusted", "Procrustes"): 0.528801294555665,
        ("Group-adjusted", "Mantel"): 0.554311306676936,
    }
    check(
        "result-46-global-statistics",
        all(close(global_tests[key], value) for key, value in expected_global_tests.items()),
        {" / ".join(key): value for key, value in global_tests.items()},
    )

    summary47 = json.loads((result_root / "47-spls-diablo" / "summary.json").read_text(encoding="utf-8"))
    check("result-47-split", summary47["train_samples"] == 153 and summary47["test_samples"] == 67, summary47)
    check("result-47-test-performance", close(summary47["test_accuracy"], 0.6567) and close(summary47["test_macro_auc"], 0.7984), summary47)
    check("result-47-stability", summary47["stable_selected_features_rate_ge_0_8"] == 26, summary47)

    summary48 = json.loads((result_root / "48-mmvec-mofa" / "summary.json").read_text(encoding="utf-8"))
    check("result-48-native-implementation", summary48["mmvec_implementation"] == "biocore/mmvec native MMvec class", summary48["mmvec_implementation"])
    check("result-48-split", summary48["mmvec_train_samples"] == 153 and summary48["mmvec_holdout_samples"] == 67, summary48)
    check("result-48-best-epoch", summary48["mmvec_epochs_run"] == 100 and summary48["mmvec_best_epoch_one_based"] == 100, summary48)
    check("result-48-probability-contract", summary48["mmvec_probability_max_row_sum_error"] < 1e-12, summary48["mmvec_probability_max_row_sum_error"])
    check("result-48-mofa", summary48["mofa_factors_retained"] == 4 and summary48["mofa_group_associated_factors_fdr"] == 3, summary48)
    check("result-48-method-overlap", summary48["spearman_mmvec_top50_pair_overlap"] == 0, summary48)
    check("result-48-native-versions", summary48["versions"] == {
        "python": "3.10.19", "mmvec": "1.0.5", "tensorflow": "1.15.0",
        "mofapy2": "0.7.4", "scikit_learn": "1.7.2",
    }, summary48["versions"])
    check("result-48-no-obsolete-wrapper-output", not (result_root / "48-mmvec-mofa" / "mmvec-cooccurrence.tsv").exists(), "obsolete output absent")

    summary49 = json.loads((result_root / "49-multi-kingdom" / "summary.json").read_text(encoding="utf-8"))
    check("result-49-design", summary49["samples"] == 36 and summary49["balanced_design_cells"] == 9 and summary49["replicates_per_cell"] == 4, summary49)
    check("result-49-associations", summary49["tested_cross_kingdom_pairs"] == 1200 and summary49["raw_cross_kingdom_q_lt_0_05"] == 505 and summary49["adjusted_cross_kingdom_q_lt_0_05"] == 2, summary49)
    check("result-49-stable", summary49["sign_changed_after_adjustment"] == 491 and summary49["stable_candidates_q_lt_0_05"] == 2, summary49)

    feast = json.loads((result_root / "50-source-tracking" / "summary-feast.json").read_text(encoding="utf-8"))
    source_tracker = json.loads((result_root / "50-source-tracking" / "summary-source-tracker.json").read_text(encoding="utf-8"))
    check("result-50-largest-source", feast["main_largest_source_class"] == "Infant gut" and source_tracker["main_largest_source_class"] == "Infant gut", {"FEAST": feast["main_largest_source_class"], "SourceTracker2": source_tracker["main_largest_source_class"]})
    check("result-50-unknown", close(feast["main_unknown"], 0.1854, tolerance=1e-4) and close(source_tracker["main_unknown"], 0.0044), {"FEAST": feast["main_unknown"], "SourceTracker2": source_tracker["main_unknown"]})
    check("result-50-loo", feast["leave_one_out_known_class_accuracy"] == 0.5 and source_tracker["leave_one_out_known_class_accuracy"] == 0.5, {"FEAST": feast["leave_one_out_known_class_accuracy"], "SourceTracker2": source_tracker["leave_one_out_known_class_accuracy"]})
    check("result-50-missing-source-nonmonotonic", feast["missing_source_unknown"]["Infant gut"] < feast["missing_source_unknown"]["None"], feast["missing_source_unknown"])

    failed = [item for item in checks if item["status"] == "FAIL"]
    payload = {
        "status": "passed" if not failed else "failed",
        "articles": list(CHAPTERS),
        "article_count": len(CHAPTERS),
        "primary_figure_count": sum(len(stems) for _, stems in CHAPTERS.values()),
        "format_file_count": figure_files,
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "failed_checks": failed,
        "checks": checks,
    }
    output = args.output if args.output.is_absolute() else project / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "checks"}, ensure_ascii=False, indent=2, default=str))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
