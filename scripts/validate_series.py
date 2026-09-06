#!/usr/bin/env python3
"""Validate the 55-article contract and the executable formal chapters."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import yaml


PILOT_NUMBERS = {
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    34,
    35,
    36,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
}
PILOT_PRIMARY_INPUT_PATHS = {
    2: "data/small/otutab.tsv",
    3: "data/small/otutab.tsv",
    4: "data/small/decontam/otutab.tsv",
    5: "data/small/decontam/otutab.tsv",
    6: "data/small/fastq",
    7: "env/qiime2.yml",
    8: "env/qiime2.yml",
    9: "data/small/otutab.tsv",
    10: "data/small/fastq",
    11: "data/small/primer-trimming",
    12: "data/raw/nfcore-v4",
    13: "data/raw/its2",
    14: "data/small/taxonomy-databases",
    15: "data/small/taxonomy-databases",
    16: "data/small/taxonomy-databases",
    17: "data/small/phyloseq-import",
    18: "data/small/otutab.tsv",
    19: "data/small/otutab.tsv",
    20: "data/small/otutab.tsv",
    21: "data/small/beta-distances/otutab.tsv",
    22: "data/small/otutab.tsv",
    23: "data/small/permanova-dispersion/otutab.tsv",
    24: "data/small/otutab.tsv",
    25: "data/small/otutab.tsv",
    26: "data/small/otutab.tsv",
    27: "data/small/otutab.tsv",
    28: "data/small/otutab.tsv",
    29: "data/small/community-typing-dmm/otutab.tsv",
    30: "data/small/otutab.tsv",
    31: "data/small/otutab.tsv",
    32: "data/small/otutab.tsv",
    33: "data/small/otutab.tsv",
    34: "data/small/absolute-quantification",
    35: "data/small/otutab.tsv",
    36: "data/small/otutab.tsv",
    37: "data/small/otutab.tsv",
    38: "data/small/otutab.tsv",
    39: "data/small/otutab.tsv",
    40: "data/small/otutab.tsv",
    41: "data/small/picrust2-chemerin",
    42: "data/small/otutab.tsv",
    43: "data/small/absolute-quantification",
    44: "data/small/cross-cohort-crc",
    45: "data/small/survival-t1d",
    46: "data/small/paired-ibd-multiomics",
    47: "data/small/paired-ibd-multiomics",
    48: "data/small/paired-ibd-multiomics",
    49: "data/small/multi-kingdom-duran",
    50: "data/small/source-tracking-feast",
    51: "data/small/longitudinal-dietswap",
    52: "data/small/environment.tsv",
    53: "data/small/paired-ibd-multiomics",
    54: "data/small/mr-mibiogen",
    55: "data/small/causal-evidence",
}
REQUIRED_PILOT_SECTION_PATTERNS = (
    ("opening", r"\{#sec-(?:target|paper-figure)\}"),
    ("theory", r"\{#sec-theory\}"),
    ("data-environment", r"\{#sec-(?:setup|preparation)\}"),
    ("analysis", r"\{#sec-code\}"),
    ("presentation", r"\{#sec-(?:publication|polish|beautify|visualization|figure)\}"),
    ("pitfalls", r"常见.*(?:坑|误判)"),
    ("methods", r"\{#sec-methods\}"),
    ("transfer", r"\{#sec-own-data\}"),
    ("references", r"\{#sec-references\}"),
)
PROHIBITED_PUBLIC_PATTERNS = (
    "vegan::dune",
    "作者代码通常长这样",
    "本篇不依赖前面章节",
    "这正是全系列坚持",
    "Omic" + "Verse",
    "omic" + "verse",
    "审阅草稿",
    "开放审阅",
    "GitHub Draft PR",
    "整仓库",
    "单篇复现",
    "只复制本页",
    "只复制本文",
    "独立运行以上",
    "隐藏决定",
)
PROHIBITED_PUBLIC_REGEXES = (
    r'(?m)^title:\s*"第\s*\d{2}\s*篇\s*·',
    r"(?m)^##\s+这一步对应论文里的哪张图",
    r"(?m)^##\s+理论：",
    r"(?m)^##\s+(?:理论：)?为什么这么做",
    r"不复制[^。\n]{0,24}(?:原图|成图)",
    r"(?m)^\s*number-sections:\s*true\s*$",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def frontmatter(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return {}
    end = text.find("\n---\n", 4)
    if end < 0:
        return {}
    payload = yaml.safe_load(text[4:end])
    return payload if isinstance(payload, dict) else {}


def flatten_book_chapters(items: list[Any]) -> list[str]:
    paths: list[str] = []
    for item in items:
        if isinstance(item, str):
            paths.append(item)
        elif isinstance(item, dict):
            nested = item.get("chapters", [])
            if isinstance(nested, list):
                paths.extend(flatten_book_chapters(nested))
    return paths


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    chapters = manifest.get("series", {}).get("chapters", [])
    errors: list[str] = []
    warnings: list[str] = []

    numbers = [item.get("number") for item in chapters]
    if numbers != list(range(1, 56)):
        errors.append("series.chapters must contain consecutive numbers 1..55")

    files = [item.get("file") for item in chapters]
    if len(files) != len(set(files)):
        errors.append("series.chapters contains duplicate file paths")

    quarto_path = root / "_quarto.yml"
    if not quarto_path.exists():
        errors.append("_quarto.yml is missing")
        quarto_chapters: list[str] = []
    else:
        quarto = yaml.safe_load(quarto_path.read_text(encoding="utf-8"))
        quarto_chapters = flatten_book_chapters(quarto.get("book", {}).get("chapters", []))
        if quarto_chapters != files:
            errors.append("_quarto.yml chapter order does not match tutorial.yaml")

    for item in chapters:
        number = item["number"]
        chapter_path = root / item["file"]
        if not chapter_path.exists():
            errors.append(f"missing chapter file: {item['file']}")
            continue
        metadata = frontmatter(chapter_path)
        title = str(metadata.get("title", ""))
        if item["title"] not in title:
            errors.append(f"title mismatch for article {number:02d}: {item['file']}")
        kind = str(item.get("kind", "computational"))
        if kind not in {"computational", "overview", "evidence-synthesis"}:
            errors.append(f"article {number:02d} has unsupported kind: {kind}")

        if number not in PILOT_NUMBERS:
            continue

        text = chapter_path.read_text(encoding="utf-8")
        if metadata.get("draft") is True:
            errors.append(f"article {number:02d} is still marked draft")
        if kind == "computational":
            for section, pattern in REQUIRED_PILOT_SECTION_PATTERNS:
                if re.search(pattern, text) is None:
                    errors.append(
                        f"article {number:02d} is missing section function "
                        f"{section}: {pattern}"
                    )
        for pattern in PROHIBITED_PUBLIC_PATTERNS:
            if pattern in text:
                errors.append(f"article {number:02d} contains prohibited public text: {pattern}")
        for pattern in PROHIBITED_PUBLIC_REGEXES:
            if re.search(pattern, text):
                errors.append(
                    f"article {number:02d} contains prohibited public pattern: {pattern}"
                )
        if kind == "computational" and "set.seed(" not in text:
            errors.append(f"article {number:02d} does not fix a random seed")
        if kind == "computational":
            expected_input = PILOT_PRIMARY_INPUT_PATHS[number]
            if expected_input not in text:
                errors.append(
                    f"article {number:02d} does not read its declared primary input: "
                    f"{expected_input}"
                )

    intro_path = root / "index.qmd"
    if intro_path.exists():
        intro = intro_path.read_text(encoding="utf-8")
        for token in (
            "01-preview-pcoa",
            "01-preview-cap",
            "01-learning-map",
            "19-hill-profile",
            "26-phylum-stacked",
            "33-da-forest",
            "43-nested-roc",
            "55-1-evidence-ladder",
            "## 参考",
        ):
            if token not in intro:
                errors.append(f"article 01 is missing overview evidence {token}")
        for token in (
            "```",
            "## 准备工作",
            "install.packages(",
            "download.file(",
            "SHA-256",
            "推荐的最短可信主线",
            "## 这段 Methods 怎么写",
        ):
            if token in intro:
                errors.append(f"article 01 contains non-overview material: {token}")
        if len(re.findall(r"!\[[^\]]*\]\([^\)]+\)", intro)) < 7:
            errors.append("article 01 must display at least seven representative figures")

    scope_path = root / "chapters/02-scope-and-limits.qmd"
    if scope_path.exists():
        scope = scope_path.read_text(encoding="utf-8")
        for token in (
            "rank_summary",
            "scope_matrix",
            "02-taxonomy-resolution",
            "02-assay-scope-map",
        ):
            if token not in scope:
                errors.append(f"article 02 is missing executable {token}")

    design_path = root / "chapters/03-study-design.qmd"
    if design_path.exists():
        design = design_path.read_text(encoding="utf-8")
        for token in (
            "design_counts",
            "perfect_group_salinity_confounding",
            "rrarefy(",
            "power.t.test(",
            "03-design-confounding",
            "03-power-sensitivity",
        ):
            if token not in design:
                errors.append(f"article 03 is missing executable {token}")

    contamination_path = root / "chapters/04-contamination-controls.qmd"
    if contamination_path.exists():
        contamination = contamination_path.read_text(encoding="utf-8")
        for token in (
            "isContaminant(",
            "contam_combined",
            "contam_batch",
            "04-control-library-size",
            "04-contaminant-prevalence",
            "04-contaminant-burden",
        ):
            if token not in contamination:
                errors.append(f"article 04 is missing executable {token}")

    batch_path = root / "chapters/05-batch-effects.qmd"
    if batch_path.exists():
        batch = batch_path.read_text(encoding="utf-8")
        for token in (
            "run_permanova(",
            "betadisper(",
            "feature_diagnostics",
            "leave_one_plate_out",
            "plate_contrasts",
            "05-batch-design-audit",
            "05-community-batch-effects",
            "05-feature-batch-effects",
            "05-adjustment-diagnostic",
        ):
            if token not in batch:
                errors.append(f"article 05 is missing executable {token}")

    fastq_path = root / "chapters/06-data-and-fastq.qmd"
    if fastq_path.exists():
        fastq = fastq_path.read_text(encoding="utf-8")
        for token in (
            "read_fastq(",
            "decode_phred33",
            "quality_matrix",
            "ObservedSHA256",
            "read_q2_metadata",
            "06-read-architecture",
            "06-fastq-anatomy",
            "06-quality-profile",
            "06-file-contract",
        ):
            if token not in fastq:
                errors.append(f"article 06 is missing executable {token}")

    wsl_path = root / "chapters/07-wsl2-conda.qmd"
    if wsl_path.exists():
        wsl = wsl_path.read_text(encoding="utf-8")
        for token in (
            "wsl.exe --list --verbose",
            "Miniforge3-26.3.2-2-Linux-x86_64.sh",
            "mamba env create",
            "qiime info",
            "validate_wsl2_conda.py",
            "07-wsl2-layer-map",
            "07-environment-validation",
            "data/small/fastq",
        ):
            if token not in wsl:
                errors.append(f"article 07 is missing executable {token}")
        wsl_metadata = frontmatter(wsl_path)
        execute = wsl_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 07 must retain the one-time upstream eval:false policy")

    qiime_install_path = root / "chapters/08-qiime2-install.qmd"
    if qiime_install_path.exists():
        qiime_install = qiime_install_path.read_text(encoding="utf-8")
        for token in (
            "qiime tools import",
            "qiime tools validate",
            "qiime tools export",
            "EMPPairedEndSequences",
            "EMPPairedEndDirFmt",
            "validate_qiime2_install.py",
            "08-qiime2-install-contract",
            "08-artifact-smoke-audit",
            "08-error-triage",
            "data/small/fastq",
        ):
            if token not in qiime_install:
                errors.append(f"article 08 is missing executable {token}")
        qiime_install_metadata = frontmatter(qiime_install_path)
        execute = qiime_install_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 08 must retain the one-time upstream eval:false policy")

    r_ecosystem_path = root / "chapters/09-r-ecosystem.qmd"
    if r_ecosystem_path.exists():
        r_ecosystem = r_ecosystem_path.read_text(encoding="utf-8")
        for token in (
            "renv::restore",
            'BiocManager::install(',
            "phyloseq::phyloseq(",
            "microeco::microtable$new(",
            "validate_r_ecosystem.R",
            "09-r-ecosystem-map",
            "09-package-version-audit",
            "09-object-parity-audit",
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
        ):
            if token not in r_ecosystem:
                errors.append(f"article 09 is missing executable {token}")
        r_ecosystem_metadata = frontmatter(r_ecosystem_path)
        execute = r_ecosystem_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 09 must retain the one-time upstream eval:false policy")

    import_qc_path = root / "chapters/10-import-provenance-qc.qmd"
    if import_qc_path.exists():
        import_qc = import_qc_path.read_text(encoding="utf-8")
        for token in (
            "qiime tools import",
            "qiime tools validate",
            "qiime demux emp-paired",
            "--p-rev-comp-mapping-barcodes",
            "EMPPairedEndSequences",
            "SampleData[PairedEndSequencesWithQuality]",
            "fastqc",
            "python -m multiqc",
            "qiime tools replay-provenance",
            "validate_import_provenance_qc.py",
            "10-import-provenance-map",
            "10-fastqc-module-audit",
            "10-demultiplexing-qc",
            "data/small/fastq/forward.fastq.gz",
            "data/small/fastq/reverse.fastq.gz",
            "data/small/fastq/barcodes.fastq.gz",
            "data/small/fastq/metadata.tsv",
        ):
            if token not in import_qc:
                errors.append(f"article 10 is missing executable {token}")
        import_qc_metadata = frontmatter(import_qc_path)
        execute = import_qc_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 10 must retain the one-time upstream eval:false policy")

    primer_path = root / "chapters/11-primer-trimming.qmd"
    if primer_path.exists():
        primer = primer_path.read_text(encoding="utf-8")
        for token in (
            "qiime cutadapt trim-paired",
            "--p-front-f '^GTGYCAGCMGCCGCGGTAA'",
            "--p-front-r '^GGACTACNVGGGTWTCTAAT'",
            "--p-discard-untrimmed",
            "--action none",
            '-G "^${REV_PRIMER}"',
            "validate_primer_trimming.py",
            "11-primer-presence-gate",
            "11-trimming-retention-audit",
            "11-primer-read-architecture",
            "11-region-branch-decision",
            "data/small/primer-trimming",
            "data/small/fastq/forward.fastq.gz",
            "data/small/fastq/reverse.fastq.gz",
            "50 / 50 PASS",
            "86.14%",
        ):
            if token not in primer:
                errors.append(f"article 11 is missing executable {token}")
        primer_metadata = frontmatter(primer_path)
        execute = primer_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 11 must retain the one-time upstream eval:false policy")

    dada2_path = root / "chapters/12-dada2-asv.qmd"
    if dada2_path.exists():
        dada2 = dada2_path.read_text(encoding="utf-8")
        for token in (
            "qiime dada2 denoise-paired",
            "--p-trunc-len-f 220",
            "--p-trunc-len-r 200",
            "--p-max-ee-f 2",
            "--p-max-ee-r 2",
            "--p-min-overlap 12",
            "--p-max-merge-mismatch 0",
            "--p-retain-all-samples",
            "--o-base-transition-stats",
            "--o-feature-frequencies",
            "--o-sample-frequencies",
            "--o-summary",
            "qiime dada2 plot-base-transitions",
            "12-quality-truncation-gate",
            "12-denoising-loss-ledger",
            "12-overlap-budget",
            "12-asv-output-audit",
            "data/raw/nfcore-v4",
            "122 / 122 PASS",
            "60.52%",
            "5,213",
            "366",
        ):
            if token not in dada2:
                errors.append(f"article 12 is missing executable {token}")
        dada2_metadata = frontmatter(dada2_path)
        execute = dada2_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 12 must retain the one-time upstream eval:false policy")

    its_18s_path = root / "chapters/13-its-18s.qmd"
    if its_18s_path.exists():
        its_18s = its_18s_path.read_text(encoding="utf-8")
        for token in (
            "qiime itsxpress trim-pair-output-unmerged",
            "--p-region ITS2",
            "--p-taxa F",
            "--p-cluster-id 1.0",
            "qiime dada2 denoise-paired",
            "--p-trunc-len-f 0",
            "--p-trunc-len-r 0",
            "PairedEndFastqManifestPhred33V2",
            "UNITE v10.0",
            "PR2 v5.1.1",
            "916145d3b05e20656fde30c4732dc084de72d876",
            "13-marker-branch-map",
            "13-itsxpress-length-audit",
            "13-its2-loss-ledger",
            "13-reference-database-contract",
            "data/raw/its2",
            "55 / 55 PASS",
            "2,065",
            "1,753",
            "240,201",
        ):
            if token not in its_18s:
                errors.append(f"article 13 is missing executable {token}")
        its_18s_metadata = frontmatter(its_18s_path)
        execute = its_18s_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 13 must retain the one-time upstream eval:false policy")

    taxonomy_database_path = root / "chapters/14-taxonomy-databases.qmd"
    if taxonomy_database_path.exists():
        taxonomy_database = taxonomy_database_path.read_text(encoding="utf-8")
        for token in (
            "qiime feature-classifier classify-consensus-vsearch",
            "--p-maxaccepts 50",
            "--p-maxrejects all",
            "--p-perc-identity 0.8",
            "--p-query-cov 0.8",
            "--p-top-hits-only",
            "--p-min-consensus 0.51",
            "SILVA 138.2",
            "GTDB R11-RS232",
            "Greengenes2 2024.09",
            "14-database-decision-map",
            "14-reference-scope-audit",
            "14-rank-assignment-coverage",
            "14-taxonomy-concordance",
            "data/small/taxonomy-databases",
            "91 / 91 PASS",
            "358 / 366",
            "356 / 366",
            "362 / 366",
        ):
            if token not in taxonomy_database:
                errors.append(f"article 14 is missing executable {token}")
        taxonomy_database_metadata = frontmatter(taxonomy_database_path)
        execute = taxonomy_database_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 14 must retain the one-time upstream eval:false policy")

    classifier_path = root / "chapters/15-train-classifier.qmd"
    if classifier_path.exists():
        classifier = classifier_path.read_text(encoding="utf-8")
        for token in (
            "qiime feature-classifier fit-classifier-naive-bayes",
            "qiime feature-classifier classify-sklearn",
            "qiime feature-classifier classify-consensus-vsearch",
            "--p-confidence 0.50",
            "--p-confidence 0.70",
            "--p-confidence 0.90",
            "--p-read-orientation same",
            "--p-maxaccepts 50",
            "236,516",
            "26,236",
            "262,752",
            "80 / 80 PASS",
            "15-classifier-build-contract",
            "15-holdout-rank-performance",
            "15-query-confidence-sensitivity",
            "15-method-concordance",
            "data/small/taxonomy-databases",
            "validate_region_classifier.py",
        ):
            if token not in classifier:
                errors.append(f"article 15 is missing executable {token}")
        classifier_metadata = frontmatter(classifier_path)
        execute = classifier_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 15 must retain the one-time upstream eval:false policy")

    sepp_path = root / "chapters/16-sepp-phylogeny.qmd"
    if sepp_path.exists():
        sepp = sepp_path.read_text(encoding="utf-8")
        for token in (
            "qiime fragment-insertion sepp",
            "qiime fragment-insertion filter-features",
            "SeppReferenceDatabase",
            "Phylogeny[Rooted]",
            "--p-alignment-subset-size 1000",
            "--p-placement-subset-size 5000",
            "--p-threads 8",
            "sepp-refs-gg-13-8.qza",
            "e252b83d7d5fbf2a9e14e594768e3578b33b557c34e22b6abc83b324689b1360",
            "16-sepp-workflow-contract",
            "16-placement-retention-audit",
            "16-placement-diagnostics",
            "16-query-placement-tree",
            "data/small/taxonomy-databases",
            "validate_sepp_phylogeny.py",
            "74 / 74 PASS",
            "203,818",
            "366 / 5,213",
        ):
            if token not in sepp:
                errors.append(f"article 16 is missing executable {token}")
        sepp_metadata = frontmatter(sepp_path)
        execute = sepp_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 16 must retain the one-time upstream eval:false policy")

    phyloseq_path = root / "chapters/17-phyloseq-import.qmd"
    if phyloseq_path.exists():
        phyloseq_import = phyloseq_path.read_text(encoding="utf-8")
        for token in (
            "data/small/phyloseq-import/otutab.tsv",
            "data/small/phyloseq-import/taxonomy.tsv",
            "data/small/phyloseq-import/metadata.tsv",
            "representative-sequences.fasta.gz",
            "rooted-insertion-tree.nwk.gz",
            "phyloseq::otu_table(",
            "taxa_are_rows = TRUE",
            "ape::keep.tip(",
            "Biostrings::readDNAStringSet(",
            "phyloseq::refseq(",
            "phyloseq::phyloseq(",
            "saveRDS(",
            "readRDS(",
            "validate_phyloseq_import.R",
            "17-phyloseq-object-contract",
            "17-id-reconciliation-audit",
            "17-library-feature-audit",
            "17-taxonomy-coverage-audit",
            "82 / 82 PASS",
            "203,818",
            "318 / 4,217",
        ):
            if token not in phyloseq_import:
                errors.append(f"article 17 is missing executable {token}")
        phyloseq_metadata = frontmatter(phyloseq_path)
        execute = phyloseq_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 17 must retain the downstream eval:true policy")

    publication_graphics_path = root / "chapters/18-publication-graphics.qmd"
    if publication_graphics_path.exists():
        publication_graphics = publication_graphics_path.read_text(encoding="utf-8")
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "colorspace::deutan(",
            "colorspace::protan(",
            "systemfonts::font_info(",
            "grDevices::cairo_pdf",
            "svglite::svglite",
            "ragg::agg_png",
            "ragg::agg_tiff",
            "compression = \"lzw\"",
            "18-encoding-redundancy-audit",
            "18-final-size-community-summary",
            "18-export-decision-map",
            "18-raster-resolution-audit",
            "validate_publication_graphics.R",
            "171 / 171 PASS",
            "89 mm",
            "183 mm",
            "600 ppi",
        ):
            if token not in publication_graphics:
                errors.append(f"article 18 is missing executable {token}")
        publication_graphics_metadata = frontmatter(publication_graphics_path)
        execute = publication_graphics_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 18 must retain the downstream eval:true policy")

    alpha_diversity_path = root / "chapters/19-alpha-diversity.qmd"
    if alpha_diversity_path.exists():
        alpha_diversity = alpha_diversity_path.read_text(encoding="utf-8")
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "phyloseq::phyloseq(",
            "vegan::estimateR(",
            "iNEXT::estimateD(",
            'base = "size"',
            'base = "coverage"',
            "size_target <- 10000L",
            "coverage_target <- 0.90",
            "Hill0",
            "Hill1",
            "Hill2",
            "19-depth-completeness-audit",
            "19-rarefaction-curves",
            "19-standardization-comparison",
            "19-hill-profile",
            "validate_alpha_diversity.R",
            "audit_summary$checks_total == 193L",
            "audit_summary$checks_passed == 193L",
        ):
            if token not in alpha_diversity:
                errors.append(f"article 19 is missing executable {token}")
        alpha_diversity_metadata = frontmatter(alpha_diversity_path)
        execute = alpha_diversity_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 19 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 19 must retain freeze:auto")

    alpha_group_tests_path = root / "chapters/20-alpha-group-tests.qmd"
    if alpha_group_tests_path.exists():
        alpha_group_tests = alpha_group_tests_path.read_text(encoding="utf-8")
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "iNEXT::estimateD(",
            'base = "coverage"',
            "coverage_target <- 0.90",
            "kruskal.test(",
            "p.adjust(",
            'method = "holm"',
            "wilcox.test(",
            "cliffs_delta",
            "bootstrap_replicates <- 5000L",
            "20-design-selection-map",
            "20-alpha-group-distributions",
            "20-standardization-sensitivity",
            "20-pairwise-effect-forest",
            "validate_alpha_group_tests.R",
            "audit_summary$checks_total == 188L",
            "audit_summary$checks_passed == 188L",
            "audit_summary$checks_failed == 0L",
            "lmerTest::lmer(",
        ):
            if token not in alpha_group_tests:
                errors.append(f"article 20 is missing executable {token}")
        alpha_group_tests_metadata = frontmatter(alpha_group_tests_path)
        execute = alpha_group_tests_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 20 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 20 must retain freeze:auto")

    beta_distances_path = root / "chapters/21-beta-distances.qmd"
    if beta_distances_path.exists():
        beta_distances = beta_distances_path.read_text(encoding="utf-8")
        for token in (
            "data/small/beta-distances/otutab.tsv",
            "data/small/beta-distances/taxonomy.tsv",
            "data/small/beta-distances/metadata.tsv",
            "rooted-tree.nwk.gz",
            "vegan::rrarefy(",
            "vegan::vegdist(",
            "phyloseq::UniFrac(",
            "aitchison_distance",
            "run_gemelli_rpca.py",
            "gemelli 0.0.13",
            "Lingoes",
            "21-distance-choice-map",
            "21-distance-correlation",
            "21-pcoa-metric-comparison",
            "21-aitchison-zero-sensitivity",
            "validate_beta_distances.R",
            "audit_summary$checks_total == 239L",
            "audit_summary$checks_passed == 239L",
            "audit_summary$checks_failed == 0L",
        ):
            if token not in beta_distances:
                errors.append(f"article 21 is missing executable {token}")
        beta_distances_metadata = frontmatter(beta_distances_path)
        execute = beta_distances_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 21 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 21 must retain freeze:auto")

    pcoa_path = root / "chapters/22-ordination-unconstrained.qmd"
    if pcoa_path.exists():
        pcoa = pcoa_path.read_text(encoding="utf-8")
        for token in ("adonis2(", "betadisper(", "permutest("):
            if token not in pcoa:
                errors.append(f"article 22 is missing executable {token}")
        if "PNG + TIFF" in pcoa and ".tiff" not in pcoa:
            errors.append("article 22 claims TIFF output but does not create it")

    permanova_dispersion_path = root / "chapters/23-permanova-dispersion.qmd"
    if permanova_dispersion_path.exists():
        permanova_dispersion = permanova_dispersion_path.read_text(
            encoding="utf-8"
        )
        for token in (
            "data/small/permanova-dispersion/otutab.tsv",
            "data/small/permanova-dispersion/taxonomy.tsv",
            "data/small/permanova-dispersion/metadata.tsv",
            "make_splitplot_permutations",
            "permute::how(",
            "permute::numPerms(",
            "vegan::adonis2(",
            "vegan::betadisper(",
            "vegan::permutest(",
            "vegan::anosim(",
            "vegan::mrpp(",
            'method = "holm"',
            "23-exchangeability-design",
            "23-pcoa-centroid-dispersion",
            "23-pairwise-permanova",
            "23-dispersion-by-treatment",
            "validate_permanova_dispersion.R",
            "audit_summary$checks_total == 214L",
            "audit_summary$checks_passed == 214L",
            "audit_summary$checks_failed == 0L",
        ):
            if token not in permanova_dispersion:
                errors.append(f"article 23 is missing executable {token}")
        permanova_dispersion_metadata = frontmatter(permanova_dispersion_path)
        execute = permanova_dispersion_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 23 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 23 must retain freeze:auto")

    cap_path = root / "chapters/24-ordination-constrained-cap.qmd"
    if cap_path.exists():
        cap = cap_path.read_text(encoding="utf-8")
        for token in ('by = "terms"', "vif.cca(", "RsquareAdj("):
            if token not in cap:
                errors.append(f"article 24 is missing executable {token}")
        if not re.search(r"(lab1_mod2|lab1_env|lab1_arrow)", cap):
            errors.append("article 24 does not define model-specific axis labels for the environment model")

    environment_variance_path = root / "chapters/25-environment-variance.qmd"
    if environment_variance_path.exists():
        environment_variance = environment_variance_path.read_text(
            encoding="utf-8"
        )
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "data/small/environment.tsv",
            "permute::shuffleSet(",
            "vegan::decostand(",
            "vegan::envfit(",
            "vegan::mantel(",
            "vegan::mantel.partial(",
            "vegan::varpart(",
            "vegan::RsquareAdj(",
            "vegan::vif.cca(",
            'method = "BH"',
            'method = "holm"',
            "25-environment-correlation",
            "25-envfit-pcoa",
            "25-mantel-distance-blocks",
            "25-variation-partition",
            "validate_environment_variance.R",
            "audit_summary$checks_total == 221L",
            "audit_summary$checks_passed == 221L",
            "audit_summary$checks_failed == 0L",
        ):
            if token not in environment_variance:
                errors.append(f"article 25 is missing executable {token}")
        environment_variance_metadata = frontmatter(environment_variance_path)
        execute = environment_variance_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 25 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 25 must retain freeze:auto")

    community_composition_path = root / "chapters/26-community-composition.qmd"
    if community_composition_path.exists():
        community_composition = community_composition_path.read_text(
            encoding="utf-8"
        )
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "clean_taxon",
            "aggregate_rank(",
            "phylum_display",
            "top_n_sensitivity",
            "denominator_sensitivity",
            "make_flow_polygons",
            "log1p(",
            "26-phylum-stacked",
            "26-genus-bubble",
            "26-group-phylum-alluvial",
            "26-genus-heatmap",
        ):
            if token not in community_composition:
                errors.append(f"article 26 is missing executable {token}")
        community_composition_metadata = frontmatter(community_composition_path)
        execute = community_composition_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 26 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 26 must retain freeze:auto")
        # Verification belongs to maintainer tooling, not the public lesson.
        if not (root / "scripts/validate_community_composition.R").is_file():
            errors.append("article 26 independent analysis validator is missing")
        composition_qa = root / "results/26-community-composition/community-composition-summary.json"
        try:
            summary = json.loads(composition_qa.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            summary = {}
        if any(summary.get(key) != expected for key, expected in (
            ("checks_total", 195), ("checks_passed", 195), ("checks_failed", 0)
        )):
            errors.append("article 26 independent analysis checks must pass")

    multirank_composition_path = root / "chapters/27-multirank-composition.qmd"
    if multirank_composition_path.exists():
        multirank_composition = multirank_composition_path.read_text(
            encoding="utf-8"
        )
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "rank_names <- c(\"Phylum\", \"Class\", \"Order\", \"Family\", \"Genus\")",
            "RankQualifiedTaxon",
            "lineage_key",
            "complete_depth",
            "lineage_gap_audit",
            "label_collision_audit",
            "rank_topn_sensitivity",
            "aggregation_sensitivity",
            "denominator_sensitivity",
            "27-rank-resolution-cascade",
            "27-multirank-bubble",
            "27-lineage-ladder",
            "27-topn-coverage",
            "nrow(validation_checks) == 209L",
            "all(validation_checks$status == \"PASS\")",
        ):
            if token not in multirank_composition:
                errors.append(f"article 27 is missing executable {token}")
        multirank_composition_metadata = frontmatter(multirank_composition_path)
        execute = multirank_composition_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 27 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 27 must retain freeze:auto")

    core_rare_path = root / "chapters/28-core-rare-biosphere.qmd"
    if core_rare_path.exists():
        core_rare = core_rare_path.read_text(encoding="utf-8")
        for token in (
            "data/small/otutab.tsv",
            "data/small/taxonomy.tsv",
            "data/small/metadata.tsv",
            "primary_detection <- 0.0001",
            "primary_prevalence <- 0.80",
            "relative_feature >= primary_detection",
            "global_core",
            "group_core",
            "core_membership_patterns",
            "core_threshold_sensitivity",
            "detection_count_audit",
            "vegan::rrarefy(",
            "rarefied_repeat",
            "feature_state_summary",
            "rare_mass_by_sample",
            "28-occupancy-abundance",
            "28-group-core-membership",
            "28-threshold-depth-sensitivity",
            "28-rare-biosphere-mass",
            "nrow(validation_checks) == 320L",
            "all(validation_checks$status == \"PASS\")",
        ):
            if token not in core_rare:
                errors.append(f"article 28 is missing executable {token}")
        core_rare_metadata = frontmatter(core_rare_path)
        execute = core_rare_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 28 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 28 must retain freeze:auto")

    community_typing_path = root / "chapters/29-community-typing-dmm.qmd"
    if community_typing_path.exists():
        community_typing = community_typing_path.read_text(encoding="utf-8")
        for token in (
            "data/small/community-typing-dmm/otutab.tsv",
            "data/small/community-typing-dmm/taxonomy.tsv",
            "data/small/community-typing-dmm/metadata.tsv",
            "DirichletMultinomial::dmn(",
            "candidate_k <- 1:7",
            "primary_seed <- 20260729L",
            "match_profiles",
            "adjusted_rand_index",
            "MaximumPosterior",
            "NormalizedEntropy",
            "vegan::rrarefy(",
            "participant_repeat_audit",
            "audit-original-vs-upgraded-dmm",
            "29-model-selection",
            "29-posterior-ordination",
            "29-component-profiles",
            "29-stability-audit",
            "nrow(validation_checks) == 222L",
            "all(validation_checks$status == \"PASS\")",
        ):
            if token not in community_typing:
                errors.append(f"article 29 is missing executable {token}")
        community_typing_metadata = frontmatter(community_typing_path)
        execute = community_typing_metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append("article 29 must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 29 must retain freeze:auto")

    downstream_contracts = {
        30: (
            "chapters/30-compositional-data.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "closure_audit",
                "reference_effects",
                "30-closure-artifact",
                "30-measurement-scales",
                "30-reference-frame",
            ),
        ),
        31: (
            "chapters/31-da-methods.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "ANCOMBC::ancombc2(",
                "ALDEx2::aldex.clr(",
                "Maaslin2::Maaslin2(",
                "lefser::lefser(",
                "corncob::differentialTest(",
                "31-da-hit-counts",
                "31-da-jaccard",
                "31-da-evidence-map",
            ),
        ),
        32: (
            "chapters/32-multirank-da.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "aggregate_rank <- function",
                "PrimaryReportingEligible",
                "ReportingGate",
                "32-rank-evidence-cascade",
                "32-multirank-effect-map",
                "32-family-genus-coherence",
            ),
        ),
        33: (
            "chapters/33-da-visualization.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "master_results",
                "build-audited-cladogram",
                "33-da-volcano",
                "33-da-cladogram",
                "33-da-manhattan",
                "33-da-forest",
            ),
        ),
        34: (
            "chapters/34-absolute-quantification.qmd",
            (
                'file.path(data_dir, "taxonomy.tsv")',
                'file.path(data_dir, "metadata.tsv")',
                "cell-load.tsv",
                "calculate-qmp",
                "34-microbial-load",
                "34-relative-quantitative-effects",
                "34-qmp-exemplar",
            ),
        ),
        35: (
            "chapters/35-cooccurrence-networks.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "SpiecEasi::sparccboot(",
                "SpiecEasi::spiec.easi(",
                "group_centered_clr",
                "35-network-matrix",
                "35-spiec-network",
                "35-network-edge-audit",
            ),
        ),
        36: (
            "chapters/36-network-robustness.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "SpiecEasi::sparcc(",
                "SpiecEasi::spiec.easi(",
                "calculate_roles",
                "group_bootstrap_frequency",
                "36-role-cartography",
                "36-attack-robustness",
                "36-group-rewiring",
                "36-topology-sensitivity",
            ),
        ),
        37: (
            "chapters/37-microbial-wgcna.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "data/small/environment.tsv",
                "WGCNA::pickSoftThreshold(",
                "WGCNA::blockwiseModules(",
                "WGCNA::signedKME(",
                "37-soft-threshold",
                "37-module-dendrogram",
                "37-module-trait",
                "37-hub-sensitivity",
            ),
        ),
        38: (
            "chapters/38-community-assembly-bnti.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "data/small/environment.tsv",
                "data/small/rooted-tree.nwk.gz",
                "iCAMP::bNTI.cm",
                "iCAMP::RC.cm",
                "clusterSetRNGStream",
                "38-tree-filter",
                "38-bnti-rcbray",
                "38-process-fractions",
                "38-phylogenetic-signal",
            ),
        ),
        39: (
            "chapters/39-neutral-community-model.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "fit_ncm <- function",
                "bootstrap_migration",
                "vegan::rrarefy(",
                "39-neutral-abundance-occupancy",
                "39-neutral-pool-migration",
                "39-neutral-classification-sensitivity",
                "39-neutral-detection-sensitivity",
            ),
        ),
        40: (
            "chapters/40-niche-distance-decay.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "data/small/environment.tsv",
                "standardized_b",
                "vegan::radfit(",
                "haversine_km",
                "vegan::mantel(",
                "40-niche-breadth",
                "40-rank-abundance-models",
                "40-distance-decay",
                "40-distance-decay-audit",
            ),
        ),
        42: (
            "chapters/42-functional-guilds.qmd",
            (
                "data/small/taxonomy.tsv",
                "data/small/metadata.tsv",
                "microeco::trans_func$new",
                'prok_database = "FAPROTAX"',
                "42-function-coverage",
                "42-functional-composition",
                "42-functional-heatmap",
                "42-functional-sensitivity",
            ),
        ),
        43: (
            "chapters/43-random-forest.qmd",
            (
                "cell-load.tsv",
                "ranger::ranger",
                "nested-random-forest-classification",
                "permutation_auc <- numeric(50L)",
                "43-nested-roc",
                "43-calibration",
                "43-permutation-importance",
                "43-regression-performance",
            ),
        ),
        44: (
            "chapters/44-cross-cohort-validation.qmd",
            (
                "metafor::rma.uni",
                'method = "REML"',
                "leave-one-study-out-classification",
                "HeldOutStudy",
                "44-cohort-pcoa",
                "44-meta-forest",
                "44-heterogeneity",
                "44-external-validation",
            ),
        ),
        45: (
            "chapters/45-survival-analysis.qmd",
            (
                "metadata$SamplingWeek",
                "survival::coxph",
                "survival::cox.zph",
                "glmnet::cv.glmnet",
                "timeROC::timeROC",
                "45-treatment-km",
                "45-taxon-cox",
                "45-ph-diagnostics",
                "45-time-dependent-roc",
            ),
        ),
        46: (
            "chapters/46-procrustes-mantel-halla.qmd",
            (
                "data/small/paired-ibd-multiomics",
                "protest(",
                "mantel(",
                "from halla import HAllA",
                "pairwise-group-residual.tsv",
                "46-1-procrustes",
                "46-4-halla-associations",
            ),
        ),
        47: (
            "chapters/47-spls-diablo.qmd",
            (
                "data/small/paired-ibd-multiomics",
                "spls(",
                "tune.block.splsda(",
                "block.splsda(",
                "test-predictions.tsv",
                "47-1-spls-generalization",
                "47-4-feature-stability",
            ),
        ),
        48: (
            "chapters/48-mmvec-mofa.qmd",
            (
                "data/small/paired-ibd-multiomics",
                "from mmvec.multimodal import MMvec",
                "from mofapy2.run.entry_point import entry_point",
                "mmvec-native",
                "multiomics-native",
                "ValidationMAE",
                "48-1-cca-overfit-audit",
                "48-4-mofa-summary",
            ),
        ),
        49: (
            "chapters/49-multi-kingdom.qmd",
            (
                "data/small/multi-kingdom-duran",
                "fungi-otutab.tsv",
                "oomycete-otutab.tsv",
                "adonis2(",
                "protest(",
                "cross-kingdom-associations.tsv",
                "49-1-three-kingdom-pcoa",
                "49-4-stable-candidates",
            ),
        ),
        50: (
            "chapters/50-source-tracking.qmd",
            (
                "data/small/source-tracking-feast",
                "library(FEAST)",
                "FEAST(",
                "from sourcetracker._sourcetracker import _gibbs",
                "feast-missing-source.tsv",
                "50-1-feast-main",
                "50-4-missing-source",
            ),
        ),
        51: (
            "chapters/51-longitudinal-analysis.qmd",
            (
                "data/small/longitudinal-dietswap",
                "nlme::lme(",
                "corAR1(",
                "SubjectID",
                "vegdist(",
                "Volatility",
                "51-1-alpha-trajectories",
                "51-4-model-audit",
            ),
        ),
        52: (
            "chapters/52-structural-equation-model.qmd",
            (
                "data/small/environment.tsv",
                "piecewiseSEM::psem(",
                "piecewiseSEM::fisherC(",
                "plspm::plspm(",
                "br = 2000",
                "52-1-prespecified-dag",
                "52-4-model-audit",
            ),
        ),
        53: (
            "chapters/53-mediation-analysis.qmd",
            (
                "data/small/paired-ibd-multiomics",
                "mediation::mediate(",
                "mediation::medsens(",
                "FaecalibacteriumCLR",
                "53-1-mediation-dag",
                "53-4-observed-data",
            ),
        ),
        54: (
            "chapters/54-mendelian-randomization.qmd",
            (
                "data/small/mr-mibiogen",
                "TwoSampleMR::harmonise_data(",
                "TwoSampleMR::mr(",
                "MRPRESSO::mr_presso(",
                "NbDistribution = 2000",
                "colocalisation",
                "54-1-harmonised-scatter",
                "54-4-assumption-audit",
            ),
        ),
    }
    for number, (relative_path, tokens) in downstream_contracts.items():
        chapter_path = root / relative_path
        if not chapter_path.exists():
            continue
        chapter_text = chapter_path.read_text(encoding="utf-8")
        for token in tokens:
            if token not in chapter_text:
                errors.append(f"article {number:02d} is missing executable {token}")
        metadata = frontmatter(chapter_path)
        execute = metadata.get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not True:
            errors.append(f"article {number:02d} must retain the downstream eval:true policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append(f"article {number:02d} must retain freeze:auto")

    evidence_path = root / "chapters/55-causal-evidence.qmd"
    if evidence_path.exists():
        evidence = evidence_path.read_text(encoding="utf-8")
        for token in (
            "七项真实研究提供了哪些证据",
            "不同研究设计的结论上限",
            "为什么需要跨设计三角验证",
            "病例–对照研究",
            "人群随机干预",
            "转移、救援和机制实验",
            "三角验证",
            "55-1-evidence-ladder",
            "55-4-triangulation",
        ):
            if token not in evidence:
                errors.append(f"article 55 is missing evidence framework token {token}")
        execute = frontmatter(evidence_path).get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 55 must retain the no-code eval:false policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 55 must retain freeze:auto")

    picrust_path = root / "chapters/41-picrust2.qmd"
    if picrust_path.exists():
        picrust_text = picrust_path.read_text(encoding="utf-8")
        for token in (
            "picrust2_pipeline.py", "--coverage", "--stratified",
            "41-picrust2-workflow", "41-nsti-audit",
            "41-pathway-heatmap", "41-pathway-effect",
        ):
            if token not in picrust_text:
                errors.append(f"article 41 is missing executable {token}")
        execute = frontmatter(picrust_path).get("execute", {})
        if not isinstance(execute, dict) or execute.get("eval") is not False:
            errors.append("article 41 must retain the upstream eval:false policy")
        if not isinstance(execute, dict) or execute.get("freeze") != "auto":
            errors.append("article 41 must retain freeze:auto")

    status = "passed" if not errors else "failed"
    payload = {
        "status": status,
        "chapter_count": len(chapters),
        "expected_chapter_count": 55,
        "pilot_articles": sorted(PILOT_NUMBERS),
        "quarto_chapter_count": len(quarto_chapters),
        "errors": errors,
        "warnings": warnings,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
