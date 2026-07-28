#!/usr/bin/env python3
"""Validate rendered Articles 36–40 and their release assets."""

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
    36: ("36-network-robustness", (
        "36-role-cartography", "36-attack-robustness",
        "36-group-rewiring", "36-topology-sensitivity",
    )),
    37: ("37-microbial-wgcna", (
        "37-soft-threshold", "37-module-dendrogram",
        "37-module-trait", "37-hub-sensitivity",
    )),
    38: ("38-community-assembly-bnti", (
        "38-tree-filter", "38-bnti-rcbray",
        "38-process-fractions", "38-phylogenetic-signal",
    )),
    39: ("39-neutral-community-model", (
        "39-neutral-abundance-occupancy", "39-neutral-pool-migration",
        "39-neutral-classification-sensitivity",
        "39-neutral-detection-sensitivity",
    )),
    40: ("40-niche-distance-decay", (
        "40-niche-breadth", "40-rank-abundance-models",
        "40-distance-decay", "40-distance-decay-audit",
    )),
}

REQUIRED_SECTIONS = (
    "这一步对应论文里的哪张图",
    "理论",
    "准备工作",
    "可复制代码",
    "审计与升级",
    "出版级美化",
    "常见坑",
    "这段 Methods 怎么写",
    "换成你自己的数据怎么做",
    "参考",
)

EXECUTABLE_TOKENS = {
    36: (
        "SpiecEasi::sparcc", "SpiecEasi::spiec.easi", "calculate_roles",
        "group_bootstrap_frequency", "36-role-cartography",
        "36-attack-robustness", "36-group-rewiring",
        "36-topology-sensitivity",
    ),
    37: (
        "WGCNA::pickSoftThreshold", "WGCNA::blockwiseModules",
        "WGCNA::signedKME", "module_trait_ledger", "37-soft-threshold",
        "37-module-dendrogram", "37-module-trait", "37-hub-sensitivity",
    ),
    38: (
        "data/small/rooted-tree.nwk.gz", "iCAMP::bNTI.cm", "iCAMP::RC.cm",
        "clusterSetRNGStream", "primary_randomizations <- 999L",
        "38-tree-filter", "38-bnti-rcbray", "38-process-fractions",
        "38-phylogenetic-signal",
    ),
    39: (
        "fit_ncm <- function", "bootstrap_migration", "vegan::rrarefy",
        "39-neutral-abundance-occupancy", "39-neutral-pool-migration",
        "39-neutral-classification-sensitivity",
        "39-neutral-detection-sensitivity",
    ),
    40: (
        "standardized_b", "vegan::radfit", "haversine_km",
        "vegan::mantel", "40-niche-breadth", "40-rank-abundance-models",
        "40-distance-decay", "40-distance-decay-audit",
    ),
}

PACKAGE_VERSIONS = {
    "WGCNA": "1.74",
    "iCAMP": "1.5.12",
    "SpiecEasi": "1.1.2",
    "microeco": "2.0.0",
    "vegan": "2.6-6.1",
}

CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument(
        "--site-root", type=Path, default=Path("."),
        help="Directory containing _site/, figures/ and results/.",
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
        check(
            f"article-{number}-not-draft",
            metadata.get("draft") is not True,
            metadata.get("draft"),
        )
        check(f"article-{number}-seed", "set.seed(" in qmd_text, "set.seed present")
        for section in REQUIRED_SECTIONS:
            check(
                f"article-{number}-section-{hashlib.sha1(section.encode()).hexdigest()[:8]}",
                section in qmd_text,
                section,
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
                    check(
                        f"signature-{stem}-pdf",
                        path.read_bytes()[:4] == b"%PDF",
                        "PDF signature",
                    )
                elif extension == "svg":
                    svg = path.read_text(encoding="utf-8")
                    check(f"signature-{stem}-svg", "<svg" in svg, "SVG root")
                    check(
                        f"english-{stem}-svg",
                        CJK_RE.search(svg) is None,
                        "no CJK plot text",
                    )
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
                            check(
                                f"raster-{stem}-tiff-lzw",
                                "lzw" in compression,
                                compression,
                            )

    check("figure-file-count", figure_files == 80, figure_files)
    check(
        "wetland-triad-contract",
        count_table(project / "data" / "small" / "otutab.tsv")
        == (13_628, 90, 1_619_670),
        count_table(project / "data" / "small" / "otutab.tsv"),
    )

    source_summary_path = project / "data" / "small" / "source_summary.json"
    source_summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
    for key, expected in {
        "tree_tips": 14_096,
        "tree_internal_nodes": 14_095,
        "tree_rooted": True,
        "tree_contains_all_features": True,
        "tree_extra_tips": 468,
    }.items():
        check(f"tree-summary-{key}", source_summary.get(key) == expected, source_summary.get(key))
    tree = project / "data" / "small" / "rooted-tree.nwk.gz"
    tree_sha = sha256(tree)
    check(
        "rooted-tree-sha256",
        tree_sha == "05d64719bfe720714fdf03f5158893f54fc6400ab649de807cbecc5599589f61",
        tree_sha,
    )
    with gzip.open(tree, "rt", encoding="utf-8") as handle:
        newick = handle.read().strip()
    check("rooted-tree-newick", newick.startswith("(") and newick.endswith(";"), len(newick))

    lock = json.loads((project / "env" / "renv.lock").read_text(encoding="utf-8"))
    packages = lock.get("Packages", {})
    check("renv-package-count", len(packages) == 418, len(packages))
    for package, version in PACKAGE_VERSIONS.items():
        observed = packages.get(package, {}).get("Version")
        check(f"renv-{package}", observed == version, observed)
    check(
        "renv-SpiecEasi-sha",
        packages.get("SpiecEasi", {}).get("RemoteSha")
        == "c463727a51d0df34db0c670d3b170195bb3d4eba",
        packages.get("SpiecEasi", {}).get("RemoteSha"),
    )

    result_root = site_root / "results"

    attacks = read_tsv(result_root / "36-network-robustness" / "attack-summary.tsv")
    sparcc_attack = next(row for row in attacks if row["Method"] == "SparCC")
    check("result-36-sparcc-giant", int(sparcc_attack["InitialGiantNodes"]) == 40, sparcc_attack)
    check("result-36-degree-auc", close(float(sparcc_attack["DegreeAttackAUC"]), 0.419375), sparcc_attack)
    roles = read_tsv(result_root / "36-network-robustness" / "sparcc-node-roles.tsv")
    role_counts = Counter(row["Role"] for row in roles)
    check("result-36-role-count", len(roles) == 40, len(roles))
    check("result-36-connectors", role_counts["Connector"] == 3, role_counts)
    jaccard_rows = read_tsv(result_root / "36-network-robustness" / "group-edge-jaccard.tsv")
    mean_jaccard = sum(float(row["Jaccard"]) for row in jaccard_rows) / 3
    check("result-36-mean-jaccard", close(mean_jaccard, 0.052847696), mean_jaccard)

    powers = read_tsv(result_root / "37-microbial-wgcna" / "soft-threshold.tsv")
    power_seven = next(row for row in powers if row["Power"] == "7")
    check("result-37-power-r2", close(float(power_seven["SFT.R.sq"]), 0.876504314), power_seven)
    modules = read_tsv(result_root / "37-microbial-wgcna" / "module-sizes.tsv")
    module_sizes = {row["Module"]: int(row["Genera"]) for row in modules}
    check(
        "result-37-module-sizes",
        module_sizes == {"blue": 72, "brown": 57, "grey": 115, "turquoise": 127, "yellow": 51},
        module_sizes,
    )
    traits = read_tsv(result_root / "37-microbial-wgcna" / "module-trait-ledger.tsv")
    check("result-37-trait-tests", len(traits) == 32, len(traits))
    check("result-37-trait-hits", sum(float(row["QValue"]) < 0.05 for row in traits) == 23, "BH q < 0.05")
    hubs = read_tsv(result_root / "37-microbial-wgcna" / "hub-ledger.tsv")
    check("result-37-hub-candidates", sum(row["CandidateHub"] == "TRUE" for row in hubs) == 19, "candidate hubs")

    tree_contract = read_tsv(result_root / "38-community-assembly-bnti" / "tree-tip-contract.tsv")
    check("result-38-tree-contract", len(tree_contract) >= 1, tree_contract[:2])
    assembly = read_tsv(result_root / "38-community-assembly-bnti" / "pairwise-process-ledger.tsv")
    check("result-38-pair-count", len(assembly) == 4005, len(assembly))
    process_counts = Counter(row["Process"] for row in assembly)
    check("result-38-process-total", sum(process_counts.values()) == 4005, process_counts)
    signal = read_tsv(result_root / "38-community-assembly-bnti" / "phylogenetic-signal.tsv")
    check("result-38-signal-tests", len(signal) == 3, len(signal))
    check("result-38-signal-supported", sum(row["Supported"] == "TRUE" for row in signal) == 2, signal)
    filter_rows = read_tsv(result_root / "38-community-assembly-bnti" / "filter-sensitivity.tsv")
    check("result-38-filter-rows", len(filter_rows) == 10, len(filter_rows))

    pools = read_tsv(result_root / "39-neutral-community-model" / "pool-migration-ledger.tsv")
    pooled = next(row for row in pools if row["Pool"] == "All wetlands")
    check("result-39-migration", close(float(pooled["Migration"]), 0.102903265), pooled)
    check("result-39-r2", close(float(pooled["RSquared"]), 0.711545466), pooled)
    otu_neutral = read_tsv(result_root / "39-neutral-community-model" / "otu-neutral-ledger.tsv")
    class_counts = Counter(row["Classification"] for row in otu_neutral)
    check(
        "result-39-class-counts",
        class_counts == {"Above neutral": 4531, "Below neutral": 1086, "Neutral envelope": 8011},
        class_counts,
    )
    detection = read_tsv(result_root / "39-neutral-community-model" / "detection-sensitivity.tsv")
    check("result-39-detection-grid", len(detection) == 17, len(detection))

    breadth = read_tsv(result_root / "40-niche-distance-decay" / "niche-breadth-ledger.tsv")
    check("result-40-breadth-taxa", len(breadth) == 3208, len(breadth))
    rad = read_tsv(result_root / "40-niche-distance-decay" / "rad-model-ledger.tsv")
    winners = Counter()
    for sample in sorted({row["SampleID"] for row in rad}):
        sample_rows = [row for row in rad if row["SampleID"] == sample]
        winner = min(sample_rows, key=lambda row: float(row["AIC"]))["Model"]
        winners[winner] += 1
    check("result-40-rad-winners", winners == {"Mandelbrot": 75, "Zipf": 15}, winners)
    warnings = read_tsv(result_root / "40-niche-distance-decay" / "rad-warning-ledger.tsv")
    check("result-40-rad-warnings", len(warnings) == 3, warnings)
    slopes = read_tsv(result_root / "40-niche-distance-decay" / "distance-decay-slopes.tsv")
    slope_by_region = {row["Region"]: float(row["Slope"]) for row in slopes}
    check("result-40-negative-slopes", all(value < 0 for value in slope_by_region.values()), slope_by_region)
    mantels = read_tsv(result_root / "40-niche-distance-decay" / "distance-mantel-ledger.tsv")
    tw = next(row for row in mantels if row["Region"] == "TW")
    check("result-40-tw-mantel", float(tw["PermutationP"]) > 0.05, tw)

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
