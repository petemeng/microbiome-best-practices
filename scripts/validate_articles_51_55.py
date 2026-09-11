#!/usr/bin/env python3
"""Validate Articles 51–55, their fixed public data, results and figures."""

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
    51: ("51-longitudinal-analysis", (
        "51-longitudinal-analysis/51-1-alpha-trajectories",
        "51-longitudinal-analysis/51-2-pcoa-trajectories",
        "51-longitudinal-analysis/51-3-volatility",
        "51-longitudinal-analysis/51-4-model-audit",
    )),
    52: ("52-structural-equation-model", (
        "52-structural-equation-model/52-1-prespecified-dag",
        "52-structural-equation-model/52-2-path-coefficients",
        "52-structural-equation-model/52-3-model-diagnostics",
        "52-structural-equation-model/52-4-model-audit",
    )),
    53: ("53-mediation-analysis", (
        "53-mediation-analysis/53-1-mediation-dag",
        "53-mediation-analysis/53-2-effect-decomposition",
        "53-mediation-analysis/53-3-unmeasured-confounding",
        "53-mediation-analysis/53-4-observed-data",
    )),
    54: ("54-mendelian-randomization", (
        "54-mendelian-randomization/54-1-harmonised-scatter",
        "54-mendelian-randomization/54-2-method-forest",
        "54-mendelian-randomization/54-3-leave-one-out",
        "54-mendelian-randomization/54-4-assumption-audit",
    )),
    55: ("55-causal-evidence", (
        "55-causal-evidence/55-1-evidence-ladder",
        "55-causal-evidence/55-2-threat-matrix",
        "55-causal-evidence/55-3-claim-calibration",
        "55-causal-evidence/55-4-triangulation",
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

ARTICLE_TOKENS = {
    51: (
        "data/small/longitudinal-dietswap", "nlme::lme(", "corAR1(",
        "SubjectID", "vegdist(", "Volatility", "51-4-model-audit",
    ),
    52: (
        "data/small/environment.tsv", "piecewiseSEM::psem(",
        "piecewiseSEM::fisherC(", "plspm::plspm(", "br = 2000",
        "52-4-model-audit",
    ),
    53: (
        "data/small/paired-ibd-multiomics", "mediation::mediate(",
        "mediation::medsens(", "FaecalibacteriumCLR", "sims = simulations",
        "53-4-observed-data",
    ),
    54: (
        "data/small/mr-mibiogen", "TwoSampleMR::harmonise_data(",
        "TwoSampleMR::mr(", "MRPRESSO::mr_presso(", "NbDistribution = 2000",
        "colocalisation", "54-4-assumption-audit",
    ),
    55: (
        "七项真实研究提供了哪些证据", "病例–对照研究",
        "人群随机干预", "转移、救援和机制实验",
        "跨设计三角验证", "55-4-triangulation",
    ),
}

PACKAGE_VERSIONS = {
    "piecewiseSEM": "2.3.0.1",
    "plspm": "0.6.0",
    "mediation": "4.5.1",
    "emmeans": "1.10.2",
    "TwoSampleMR": "0.6.6",
    "MRPRESSO": "1.0",
}

CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--site-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/articles-51-55-validation.json"),
    )
    return parser.parse_args()


def frontmatter(text: str) -> dict[str, Any]:
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


def close(observed: float, expected: float, tolerance: float = 1e-4) -> bool:
    return math.isclose(observed, expected, rel_tol=tolerance, abs_tol=tolerance)


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    site = args.site_root.resolve()
    checks: list[dict[str, Any]] = []

    def check(identifier: str, passed: bool, detail: Any) -> None:
        checks.append({
            "id": identifier,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        })

    tutorial = yaml.safe_load((project / "tutorial.yaml").read_text(encoding="utf-8"))
    titles = {
        int(row["number"]): str(row["title"])
        for row in tutorial["series"]["chapters"]
    }

    figure_files = 0
    for number, (slug, stems) in CHAPTERS.items():
        qmd = project / "chapters" / f"{slug}.qmd"
        text = qmd.read_text(encoding="utf-8") if qmd.is_file() else ""
        metadata = frontmatter(text)
        execute = metadata.get("execute", {})
        expected_eval = number != 55
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
        check(f"article-{number}-not-draft", metadata.get("draft") is not True,
              metadata.get("draft"))
        check(
            f"article-{number}-manifest-title",
            titles[number] in str(metadata.get("title", "")),
            {"manifest": titles[number], "qmd": metadata.get("title")},
        )
        for section, pattern in REQUIRED_SECTION_PATTERNS:
            check(
                f"article-{number}-section-{hashlib.sha1(section.encode()).hexdigest()[:8]}",
                re.search(pattern, text) is not None, pattern,
            )
        if number != 55:
            check(f"article-{number}-copyable-code", "可复制代码" in text,
                  "可复制代码")
            check(f"article-{number}-own-data", "换成你自己的数据怎么做" in text,
                  "换成你自己的数据怎么做")
            check(f"article-{number}-seed", "set.seed(" in text, "set.seed present")
            check(
                f"article-{number}-inline-theme",
                all(token in text for token in (
                    "pal_pub <-", "theme_pub <- function", "save_pub <- function",
                )),
                "pal_pub/theme_pub/save_pub inline",
            )
        else:
            check(
                f"article-{number}-evidence-comparison",
                all(token in text for token in (
                    "七项真实研究提供了哪些证据",
                    "不同研究设计的结论上限",
                    "为什么需要跨设计三角验证",
                )),
                "worked evidence comparison",
            )
            check(
                f"article-{number}-own-project",
                "在自己的研究中怎样选择可辩护措辞" in text,
                "claim calibration for the reader's project",
            )
        check(
            f"article-{number}-no-source-theme",
            re.search(r'(?m)^[ \t]*source\("R/theme_pub\.R"\)', text) is None,
            "no executable external theme dependency",
        )
        check(f"article-{number}-no-omicverse", "omicverse" not in text.lower(),
              "native method stack")
        for token in ARTICLE_TOKENS[number]:
            check(
                f"article-{number}-token-{hashlib.sha1(token.encode()).hexdigest()[:8]}",
                token in text, token,
            )

        html = site / "_site" / "chapters" / f"{slug}.html"
        html_ok = html.is_file() and html.stat().st_size > 30_000
        check(f"article-{number}-html", html_ok, str(html))
        if html_ok:
            html_text = html.read_text(encoding="utf-8")
            check(f"article-{number}-html-title", titles[number] in html_text,
                  titles[number])

        for stem in stems:
            for extension in ("pdf", "svg", "png", "tiff"):
                figure_files += 1
                path = site / "figures" / f"{stem}.{extension}"
                exists = path.is_file() and path.stat().st_size > 1_000
                check(f"figure-{stem}-{extension}", exists, str(path))
                if not exists:
                    continue
                if extension == "pdf":
                    check(f"signature-{stem}-pdf", path.read_bytes()[:4] == b"%PDF",
                          "PDF signature")
                elif extension == "svg":
                    svg = path.read_text(encoding="utf-8")
                    check(f"signature-{stem}-svg", "<svg" in svg, "SVG root")
                    check(f"english-{stem}-svg", CJK_RE.search(svg) is None,
                          "no CJK plot text")
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
                            check(f"raster-{stem}-tiff-lzw", "lzw" in compression,
                                  compression)
    check("figure-file-count", figure_files == 80, figure_files)

    longitudinal = project / "data" / "small" / "longitudinal-dietswap"
    longitudinal_summary = json.loads(
        (longitudinal / "source-summary.json").read_text(encoding="utf-8")
    )
    check("dietswap-source-hash", longitudinal_summary["packaged_source_sha256"]
          == "67cc6d117ce9dfbbc472fba828e1bb0d48c621fb574f4f8d3424f6cc9a77fdab",
          longitudinal_summary["packaged_source_sha256"])
    check("dietswap-otutab-shape", table_shape(longitudinal / "otutab.tsv") == (130, 222),
          table_shape(longitudinal / "otutab.tsv"))
    check("dietswap-taxonomy-shape", table_shape(longitudinal / "taxonomy.tsv") == (130, 8),
          table_shape(longitudinal / "taxonomy.tsv"))
    longitudinal_metadata = read_tsv(longitudinal / "metadata.tsv")
    check("dietswap-metadata", len(longitudinal_metadata) == 222 and
          len({row["SubjectID"] for row in longitudinal_metadata}) == 38,
          {"samples": len(longitudinal_metadata),
           "subjects": len({row["SubjectID"] for row in longitudinal_metadata})})
    check("dietswap-sample-order",
          table_sample_ids(longitudinal / "otutab.tsv")
          == [row["SampleID"] for row in longitudinal_metadata], "exact order")

    wetland = project / "data" / "small"
    check("wetland-otutab-shape", table_shape(wetland / "otutab.tsv") == (13628, 90),
          table_shape(wetland / "otutab.tsv"))
    check("wetland-environment-shape", table_shape(wetland / "environment.tsv") == (200, 11),
          table_shape(wetland / "environment.tsv"))
    environment_rows = read_tsv(wetland / "environment.tsv")
    check("wetland-sem-variables", all(
        token in environment_rows[0]
        for token in ("pH", "Conductivity", "TOC", "Temperature", "Precipitation")
    ), list(environment_rows[0]))

    paired = project / "data" / "small" / "paired-ibd-multiomics"
    paired_summary = json.loads((paired / "source-summary.json").read_text(encoding="utf-8"))
    check("paired-source-commit", paired_summary["commit"]
          == "89a519d8c832008fbc6e650453e83e2f04858d02", paired_summary["commit"])
    check("paired-otutab-shape", table_shape(paired / "otutab.tsv") == (250, 220),
          table_shape(paired / "otutab.tsv"))
    paired_metadata = read_tsv(paired / "metadata.tsv")
    check("paired-groups", Counter(row["StudyGroup"] for row in paired_metadata)
          == {"CD": 88, "UC": 76, "Control": 56},
          Counter(row["StudyGroup"] for row in paired_metadata))

    mr = project / "data" / "small" / "mr-mibiogen"
    mr_summary = json.loads((mr / "source-summary.json").read_text(encoding="utf-8"))
    check("mr-source-hash", mr_summary["exposure_source_sha256"]
          == "37001a83d060596fe0b97b63d6a397f01f43a29add2925d406916b7a50b5883e",
          mr_summary["exposure_source_sha256"])
    exposure_rows = read_tsv(mr / "exposure.tsv")
    outcome_rows = read_tsv(mr / "outcome.tsv")
    check("mr-instrument-count", len(exposure_rows) == len(outcome_rows) == 6,
          [len(exposure_rows), len(outcome_rows)])
    check("mr-cross-chromosome", len({row["chr.exposure"] for row in exposure_rows}) == 6,
          [row["chr.exposure"] for row in exposure_rows])
    check("mr-instrument-strength", min(float(row["FStatistic"]) for row in exposure_rows) > 20,
          min(float(row["FStatistic"]) for row in exposure_rows))

    cases = read_tsv(project / "data" / "small" / "causal-evidence" / "evidence-cases.tsv")
    check("evidence-case-count", len(cases) == 7, len(cases))
    check("evidence-levels", [int(row["EvidenceLevel"]) for row in cases] == list(range(1, 8)),
          [row["EvidenceLevel"] for row in cases])
    check("evidence-primary-dois", all(row["DOI"].startswith("10.") for row in cases),
          [row["DOI"] for row in cases])

    result_root = site / "results"
    summary51 = json.loads((result_root / "51-longitudinal-analysis" / "summary.json").read_text())
    check("result-51-design", summary51["samples"] == 222 and
          summary51["subjects"] == 38 and summary51["complete_consecutive_pairs"] == 180,
          summary51)
    check("result-51-ar1-audit", close(summary51["ar1_phi"], 0.1383) and
          close(summary51["ar1_likelihood_ratio_p"], 0.1525), summary51)

    summary52 = json.loads((result_root / "52-structural-equation-model" / "summary.json").read_text())
    check("result-52-fit", summary52["samples"] == 90 and
          close(summary52["full_dag_p"], 0.105) and summary52["reduced_dag_p"] < 0.001,
          summary52)
    check("result-52-focal-path", close(summary52["shannon_to_carbon_beta"], -0.0123)
          and close(summary52["shannon_to_carbon_p"], 0.8647), summary52)

    summary53 = json.loads((result_root / "53-mediation-analysis" / "summary.json").read_text())
    check("result-53-samples", summary53["samples"] == 108 and
          summary53["steroid_exposed"] == 28 and summary53["steroid_unexposed"] == 80,
          summary53)
    check("result-53-acme", close(summary53["acme"], 0.0102) and
          summary53["acme_lower"] < 0 < summary53["acme_upper"] and
          close(summary53["acme_p"], 0.846), summary53)

    summary54 = json.loads((result_root / "54-mendelian-randomization" / "summary.json").read_text())
    check("result-54-instruments", summary54["candidate_instruments"] == 6 and
          summary54["retained_instruments"] == 5 and
          summary54["dropped_palindromic"] == "rs10841473", summary54)
    check("result-54-ivw", close(summary54["ivw_beta"], -0.0167) and
          close(summary54["ivw_p"], 0.8756) and close(summary54["mr_presso_global_p"], 0.183),
          summary54)

    summary55 = json.loads((result_root / "55-causal-evidence" / "summary.json").read_text())
    check("result-55-framework", summary55["primary_cases"] == 7 and
          summary55["evidence_tiers"] == 7 and summary55["statistical_analysis"] is False,
          summary55)

    lock = json.loads((project / "env" / "renv.lock").read_text(encoding="utf-8"))
    packages = lock.get("Packages", {})
    check("renv-package-count", len(packages) == 418, len(packages))
    for package, version in PACKAGE_VERSIONS.items():
        observed = packages.get(package, {}).get("Version")
        check(f"renv-{package}", observed == version, observed)
    check("renv-TwoSampleMR-sha", packages.get("TwoSampleMR", {}).get("RemoteSha")
          == "78d79889fca1bfd1e29098aa16eaa4d99c61af81",
          packages.get("TwoSampleMR", {}).get("RemoteSha"))
    check("renv-MRPRESSO-sha", packages.get("MRPRESSO", {}).get("RemoteSha")
          == "cece763b47e59763a7916974de43c7cb93843e41",
          packages.get("MRPRESSO", {}).get("RemoteSha"))

    failures = [row for row in checks if row["status"] == "FAIL"]
    payload = {
        "status": "passed" if not failures else "failed",
        "articles": list(CHAPTERS),
        "checks": len(checks),
        "passed": len(checks) - len(failures),
        "failed": len(failures),
        "checks_passed": len(checks) - len(failures),
        "checks_failed": len(failures),
        "primary_figure_count": 20,
        "format_file_count": figure_files,
        "failures": failures,
    }
    output = args.output if args.output.is_absolute() else project / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
