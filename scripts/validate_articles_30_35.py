#!/usr/bin/env python3
"""Validate rendered Articles 30–35 and their release assets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import yaml
from PIL import Image


CHAPTERS = {
    30: ("30-compositional-data", (
        "30-closure-artifact", "30-measurement-scales", "30-reference-frame",
    )),
    31: ("31-da-methods", (
        "31-da-hit-counts", "31-da-jaccard", "31-da-evidence-map",
    )),
    32: ("32-multirank-da", (
        "32-rank-evidence-cascade", "32-multirank-effect-map",
        "32-family-genus-coherence",
    )),
    33: ("33-da-visualization", (
        "33-da-volcano", "33-da-cladogram", "33-da-manhattan", "33-da-forest",
    )),
    34: ("34-absolute-quantification", (
        "34-microbial-load", "34-relative-quantitative-effects", "34-qmp-exemplar",
    )),
    35: ("35-cooccurrence-networks", (
        "35-network-matrix", "35-spiec-network", "35-network-edge-audit",
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

PACKAGE_VERSIONS = {
    "ANCOMBC": "2.6.0",
    "ALDEx2": "1.36.0",
    "Maaslin2": "1.18.0",
    "lefser": "1.14.0",
    "corncob": "0.4.2",
    "reconsi": "1.16.0",
    "SpiecEasi": "1.1.2",
}

CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument(
        "--site-root",
        type=Path,
        default=Path("."),
        help="Directory containing _site/ and figures/.",
    )
    parser.add_argument("--output", type=Path, required=True)
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


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    site_root = args.site_root.resolve()
    checks: list[dict[str, Any]] = []

    def check(identifier: str, passed: bool, detail: Any) -> None:
        checks.append({"id": identifier, "status": "PASS" if passed else "FAIL", "detail": detail})

    figure_files = 0
    for number, (slug, stems) in CHAPTERS.items():
        qmd = project / "chapters" / f"{slug}.qmd"
        qmd_text = qmd.read_text(encoding="utf-8")
        metadata = read_frontmatter(qmd_text)
        execute = metadata.get("execute", {})
        check(f"article-{number}-qmd", qmd.is_file(), str(qmd))
        check(f"article-{number}-eval", isinstance(execute, dict) and execute.get("eval") is True, execute)
        check(f"article-{number}-freeze", isinstance(execute, dict) and execute.get("freeze") == "auto", execute)
        check(f"article-{number}-not-draft", metadata.get("draft") is not True, metadata.get("draft"))
        check(f"article-{number}-seed", "set.seed(" in qmd_text, "set.seed present")
        for section, pattern in REQUIRED_SECTION_PATTERNS:
            check(
                f"article-{number}-section-{hashlib.sha1(section.encode()).hexdigest()[:8]}",
                re.search(pattern, qmd_text) is not None,
                pattern,
            )

        html = site_root / "_site" / "chapters" / f"{slug}.html"
        html_ok = html.is_file() and html.stat().st_size > 50_000
        check(f"article-{number}-html", html_ok, str(html))
        if html_ok:
            html_text = html.read_text(encoding="utf-8")
            check(
                f"article-{number}-html-title",
                str(metadata.get("title", "")) in html_text,
                metadata.get("title"),
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
                            image.width >= 2500 and image.height >= 2000,
                            [image.width, image.height],
                        )
                        if extension == "tiff":
                            compression = str(image.info.get("compression", "")).lower()
                            check(f"raster-{stem}-tiff-lzw", "lzw" in compression, compression)

    check("figure-file-count", figure_files == 76, figure_files)

    wetland = count_table(project / "data" / "small" / "otutab.tsv")
    check("wetland-triad-contract", wetland == (13_628, 90, 1_619_670), wetland)
    absolute = count_table(
        project / "data" / "small" / "absolute-quantification" / "otutab.tsv"
    )
    check("absolute-triad-contract", absolute == (234, 135, 4_080_996), absolute)
    absolute_summary = project / "data" / "small" / "absolute-quantification" / "source-summary.json"
    check(
        "absolute-source-summary-sha256",
        sha256(absolute_summary) == "8ac35a7119d1772d5d1c92f30c051f2fb3bdb2b5d3dca3f83c885528352d9d94",
        sha256(absolute_summary),
    )

    lock = json.loads((project / "env" / "renv.lock").read_text(encoding="utf-8"))
    packages = lock.get("Packages", {})
    check("renv-package-count", len(packages) == 418, len(packages))
    for package, version in PACKAGE_VERSIONS.items():
        observed = packages.get(package, {}).get("Version")
        check(f"renv-{package}", observed == version, observed)
    spiec_sha = packages.get("SpiecEasi", {}).get("RemoteSha")
    check(
        "renv-SpiecEasi-sha",
        spiec_sha == "c463727a51d0df34db0c670d3b170195bb3d4eba",
        spiec_sha,
    )

    result_contracts = {
        "31-da-methods": (
            "ANCOM-BC2     79          44",
            "ALDEx2     79          47",
            "MaAsLin2     79          52",
            "LEfSe + BH     79          51",
            "corncob     79          62",
        ),
        "32-multirank-da": ("67.6%", "37.0%", "0.35%", "FAIL"),
        "34-absolute-quantification": ("Desulfovibrio", "Quantitative only", "0.002948167"),
        "35-cooccurrence-networks": ("264.00000000", "30.00000000", "0.03718974", "0.5909091"),
    }
    for slug, tokens in result_contracts.items():
        html_text = (site_root / "_site" / "chapters" / f"{slug}.html").read_text(encoding="utf-8")
        for token in tokens:
            check(f"result-{slug}-{hashlib.sha1(token.encode()).hexdigest()[:8]}", token in html_text, token)

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
