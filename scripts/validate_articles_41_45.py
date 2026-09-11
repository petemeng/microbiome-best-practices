#!/usr/bin/env python3
"""Validate rendered Articles 41–45, their data contracts, and release assets."""

from __future__ import annotations

import argparse
import csv
import gzip
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
    41: ("41-picrust2", False, (
        "41-picrust2-workflow", "41-nsti-audit",
        "41-pathway-heatmap", "41-pathway-effect",
    )),
    42: ("42-functional-guilds", True, (
        "42-function-coverage", "42-functional-composition",
        "42-functional-heatmap", "42-functional-sensitivity",
    )),
    43: ("43-random-forest", True, (
        "43-nested-roc", "43-calibration",
        "43-permutation-importance", "43-regression-performance",
    )),
    44: ("44-cross-cohort-validation", True, (
        "44-cohort-pcoa", "44-meta-forest",
        "44-heterogeneity", "44-external-validation",
    )),
    45: ("45-survival-analysis", True, (
        "45-treatment-km", "45-taxon-cox",
        "45-ph-diagnostics", "45-time-dependent-roc",
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
    41: (
        "picrust2_pipeline.py", "--coverage", "--stratified", "maximum NSTI",
        "ggpicrust2", "41-picrust2-workflow", "41-pathway-effect",
    ),
    42: (
        "microeco::trans_func", 'prok_database = "FAPROTAX"',
        "EpsilonSquared", "resolution_coverage", "42-function-coverage",
        "42-functional-sensitivity",
    ),
    43: (
        "ranger::ranger", "nested-random-forest-classification",
        "outer_predictions", "permutation_auc <- numeric(50L)", "pROC::roc",
        "43-nested-roc", "43-regression-performance",
    ),
    44: (
        "metafor::rma.uni", 'method = "REML"',
        "leave-one-study-out-classification", "HeldOutStudy", "pROC::roc",
        "44-meta-forest", "44-external-validation",
    ),
    45: (
        "metadata$SamplingWeek", "metadata$FollowUpWeeks",
        "survival::coxph", "survival::cox.zph", "glmnet::cv.glmnet",
        "timeROC::timeROC", "45-time-dependent-roc",
    ),
}

PACKAGE_VERSIONS = {
    "ggpicrust2": "2.5.17",
    "microeco": "2.0.0",
    "ranger": "0.16.0",
    "pROC": "1.18.5",
    "metafor": "4.6-0",
    "survival": "3.6-4",
    "glmnet": "4.1-8",
    "timeROC": "0.4",
}

CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument(
        "--site-root", type=Path, default=Path("."),
        help="Directory containing _site/, figures/, and results/.",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/articles-41-45-validation.json"),
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def count_table(path: Path) -> tuple[int, int, int]:
    features = 0
    reads = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = len(header) - 1
        for row in reader:
            features += 1
            reads += sum(int(value) for value in row[1:])
    return features, samples, reads


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
    for number, (slug, expected_eval, stems) in CHAPTERS.items():
        qmd = project / "chapters" / f"{slug}.qmd"
        qmd_text = qmd.read_text(encoding="utf-8") if qmd.is_file() else ""
        metadata = read_frontmatter(qmd_text)
        execute = metadata.get("execute", {})
        check(f"article-{number}-qmd", qmd.is_file(), str(qmd))
        check(
            f"article-{number}-eval",
            isinstance(execute, dict) and execute.get("eval") is expected_eval,
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

    picrust_root = project / "data" / "small" / "picrust2-chemerin"
    check(
        "picrust-triad-contract",
        count_table(picrust_root / "otutab.tsv") == (37, 24, 108_718),
        count_table(picrust_root / "otutab.tsv"),
    )
    picrust_summary = json.loads((picrust_root / "source-summary.json").read_text(encoding="utf-8"))
    check(
        "picrust-archive-sha256",
        picrust_summary["upstream_archive"]["sha256"]
        == "f57abdc069b6560f0ddb739cf9a341f5679b12f4947f4db2ad5e2cab893669e9",
        picrust_summary["upstream_archive"]["sha256"],
    )
    with gzip.open(picrust_root / "representative-sequences.fasta.gz", "rt", encoding="utf-8") as handle:
        fasta_records = sum(line.startswith(">") for line in handle)
    check("picrust-fasta-records", fasta_records == 37, fasta_records)
    for filename, record in picrust_summary["files"].items():
        path = picrust_root / filename
        check(f"picrust-file-{filename}-sha", path.is_file() and sha256(path) == record["sha256"], filename)

    crc_expected = {
        "crc_xiang": (1_617, 43, 55_706, {"CRC": 21, "Control": 22}),
        "crc_zhao": (837, 102, 19_810, {"CRC": 46, "Control": 56}),
        "crc_zackular": (60_226, 60, 3_856_096, {"CRC": 30, "Control": 30}),
    }
    crc_summary = json.loads(
        (project / "data" / "small" / "cross-cohort-crc" / "source-summary.json")
        .read_text(encoding="utf-8")
    )
    check("crc-license", crc_summary.get("license") == "CC BY-NC 4.0", crc_summary.get("license"))
    for cohort, (features, samples, reads, group_expected) in crc_expected.items():
        cohort_root = project / "data" / "small" / "cross-cohort-crc" / cohort
        observed = count_table(cohort_root / "otutab.tsv")
        check(f"crc-{cohort}-triad", observed == (features, samples, reads), observed)
        metadata = read_tsv(cohort_root / "metadata.tsv")
        groups = Counter(row["Group"] for row in metadata)
        check(f"crc-{cohort}-groups", dict(groups) == group_expected, groups)
        summary_row = crc_summary["studies"][cohort]
        check(f"crc-{cohort}-summary", summary_row["samples"] == samples and summary_row["reads"] == reads, summary_row)

    survival_root = project / "data" / "small" / "survival-t1d"
    check(
        "survival-triad-contract",
        count_table(survival_root / "otutab.tsv") == (348, 173, 3_073_108),
        count_table(survival_root / "otutab.tsv"),
    )
    survival_summary = json.loads((survival_root / "source-summary.json").read_text(encoding="utf-8"))
    check(
        "survival-source-commit",
        survival_summary["source_commit"] == "692a2ac7079d2dbf581893b9a820c80cde0a7e31",
        survival_summary["source_commit"],
    )
    survival_metadata = read_tsv(survival_root / "metadata.tsv")
    event_counts = Counter(row["Event"] for row in survival_metadata)
    week_counts = Counter(row["SamplingWeek"] for row in survival_metadata)
    check("survival-events", event_counts == {"1": 118, "0": 55}, event_counts)
    check("survival-sampling-weeks", week_counts == {"6": 6, "7": 164, "8": 3}, week_counts)
    followup_ok = all(
        close(float(row["FollowUpWeeks"]), float(row["TimeWeeks"]) - float(row["SamplingWeek"]))
        and float(row["FollowUpWeeks"]) > 0
        for row in survival_metadata
    )
    check("survival-followup-origin", followup_ok, "TimeWeeks - SamplingWeek")

    lock = json.loads((project / "env" / "renv.lock").read_text(encoding="utf-8"))
    packages = lock.get("Packages", {})
    check("renv-package-count", len(packages) == 418, len(packages))
    for package, version in PACKAGE_VERSIONS.items():
        observed = packages.get(package, {}).get("Version")
        check(f"renv-{package}", observed == version, observed)
    picrust_env = yaml.safe_load((project / "env" / "picrust2.yml").read_text(encoding="utf-8"))
    dependencies = set(picrust_env.get("dependencies", []))
    check("picrust-env-version", "picrust2=2.6.3" in dependencies, dependencies)
    check("picrust-env-python", "python=3.12.13" in dependencies, dependencies)

    result_root = site_root / "results"
    run_summary = read_tsv(result_root / "41-picrust2" / "run-summary.tsv")[0]
    check("result-41-versions", run_summary["PICRUSt2"] == "2.6.3" and run_summary["ggpicrust2"] == "2.5.17", run_summary)
    check("result-41-asvs", int(run_summary["ASVs"]) == 37 and int(run_summary["BacterialASVs"]) == 37, run_summary)
    check("result-41-max-nsti", close(float(run_summary["MaxNSTI"]), 0.217549), run_summary)
    check("result-41-weighted-nsti", close(float(run_summary["MedianWeightedNSTI"]), 0.0423077018), run_summary)
    pathways = read_tsv(result_root / "41-picrust2" / "pathway-tests.tsv")
    check("result-41-pathway-tests", len(pathways) == 225, len(pathways))
    adjusted = [float(row["AdjustedP"]) for row in pathways]
    check("result-41-no-fdr-hit", sum(value < 0.05 for value in adjusted) == 0, min(adjusted))
    output_hash_path = result_root / "41-picrust2" / "picrust2-output-sha256.tsv"
    output_hashes = sum(1 for line in output_hash_path.read_text(encoding="utf-8").splitlines() if line.strip())
    check("result-41-output-hashes", output_hashes == 22, output_hashes)
    run_log = (result_root / "41-picrust2" / "picrust2-run.log").read_text(encoding="utf-8")
    check("result-41-run-complete", "Completed PICRUSt2 pipeline in 1213.35 seconds." in run_log, "completion log")

    mapping = read_tsv(result_root / "42-functional-guilds" / "mapping-summary.tsv")[0]
    check("result-42-functions", int(mapping["Functions"]) == 93, mapping)
    check("result-42-mapped-features", int(mapping["MappedFeatures"]) == 3_992, mapping)
    check("result-42-coverage", close(float(mapping["OverallReadCoverage"]), 0.3474504066), mapping)
    resolution = {row["Mapping"]: float(row["ReadCoverage"]) for row in read_tsv(result_root / "42-functional-guilds" / "resolution-coverage.tsv")}
    check("result-42-genus-only-coverage", close(resolution["Genus only"], 0.0047830731), resolution)
    function_tests = read_tsv(result_root / "42-functional-guilds" / "function-tests.tsv")
    check("result-42-test-family", len(function_tests) == 93, len(function_tests))

    rf = read_tsv(result_root / "43-random-forest" / "performance-summary.tsv")[0]
    check("result-43-auc", close(float(rf["ClassificationAUC"]), 0.9947753396), rf)
    check("result-43-brier", close(float(rf["Brier"]), 0.0361179507), rf)
    check("result-43-permutation", close(float(rf["PermutationP"]), 1 / 51), rf)
    check("result-43-regression", close(float(rf["RegressionR2"]), 0.4358175377) and close(float(rf["RegressionRMSE"]), 0.3060877330), rf)
    check("result-43-overfit-gap", float(rf["ApparentRegressionR2"]) > float(rf["RegressionR2"]) + 0.4, rf)
    check("result-43-outer-rows", len(read_tsv(result_root / "43-random-forest" / "classification-outer-predictions.tsv")) == 285, "285")
    check("result-43-sample-rows", len(read_tsv(result_root / "43-random-forest" / "classification-sample-predictions.tsv")) == 95, "95")
    check("result-43-permutation-rows", len(read_tsv(result_root / "43-random-forest" / "label-permutation.tsv")) == 50, "50")

    dimensions = read_tsv(result_root / "44-cross-cohort-validation" / "cohort-dimensions.tsv")
    check("result-44-cohorts", len(dimensions) == 3 and sum(int(row["Samples"]) for row in dimensions) == 205, dimensions)
    meta = read_tsv(result_root / "44-cross-cohort-validation" / "random-effects-meta.tsv")
    check("result-44-meta-family", len(meta) == 58 and all(int(row["Studies"]) == 3 for row in meta), len(meta))
    pepto = next(row for row in meta if row["Genus"] == "Peptostreptococcus")
    check("result-44-pepto-effect", close(float(pepto["PooledEffect"]), 0.6007023920), pepto)
    check("result-44-pepto-fdr", close(float(pepto["AdjustedP"]), 0.000117889996), pepto)
    external = read_tsv(result_root / "44-cross-cohort-validation" / "external-performance.tsv")
    aucs = {row["HeldOutStudy"]: float(row["AUC"]) for row in external}
    check("result-44-external-aucs", all(close(aucs[key], value) for key, value in {
        "Xiang": 0.5876623377, "Zhao": 0.4363354037, "Zackular": 0.5777777778,
    }.items()), aucs)
    check("result-44-macro-auc", close(float(external[0]["MacroMeanAUC"]), 0.5339251731), external)

    split = {row["Set"]: row for row in read_tsv(result_root / "45-survival-analysis" / "prediction-split.tsv")}
    check("result-45-split", int(split["Training"]["Samples"]) == 120 and int(split["Test"]["Samples"]) == 53, split)
    check("result-45-events", int(split["Training"]["Events"]) == 82 and int(split["Test"]["Events"]) == 36, split)
    prediction = read_tsv(result_root / "45-survival-analysis" / "prediction-performance.tsv")
    check("result-45-horizons", [int(row["HorizonWeeks"]) for row in prediction] == [10, 15, 20], prediction)
    check("result-45-cindex", close(float(prediction[0]["TestCIndex"]), 0.5378850958), prediction[0])
    observed_auc = [float(row["TestAUC"]) for row in prediction]
    expected_auc = [0.4888888889, 0.6106060606, 0.6031746032]
    check("result-45-time-auc", all(close(a, b) for a, b in zip(observed_auc, expected_auc)), observed_auc)
    cox = read_tsv(result_root / "45-survival-analysis" / "genus-adjusted-cox.tsv")
    check("result-45-cox-family", len(cox) == 24, len(cox))
    check("result-45-no-fdr-hit", sum(float(row["AdjustedP"]) < 0.05 for row in cox) == 0, min(float(row["AdjustedP"]) for row in cox))
    ruminococcus = next(row for row in cox if row["Genus"] == "Ruminococcus")
    check("result-45-ruminococcus", close(float(ruminococcus["HR"]), 1.3580407711) and close(float(ruminococcus["AdjustedP"]), 0.3795030427), ruminococcus)

    failed = [item for item in checks if item["status"] == "FAIL"]
    payload = {
        "status": "passed" if not failed else "failed",
        "articles": list(CHAPTERS),
        "article_count": len(CHAPTERS),
        "primary_figure_count": sum(len(stems) for _, _, stems in CHAPTERS.values()),
        "format_file_count": figure_files,
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "failed_checks": failed,
        "checks": checks,
    }
    output = args.output
    if not output.is_absolute():
        output = project / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(
        {key: payload[key] for key in payload if key != "checks"},
        ensure_ascii=False,
        indent=2,
        default=str,
    ))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
