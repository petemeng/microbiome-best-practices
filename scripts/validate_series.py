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
}
PILOT_PRIMARY_INPUT_PATHS = {
    1: "data/small/otutab.tsv",
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
}
REQUIRED_PILOT_SECTIONS = (
    "这一步对应论文里的哪张图",
    "理论",
    "准备工作",
    "可复制代码",
    "出版级美化",
    "常见坑",
    "这段 Methods 怎么写",
    "换成你自己的数据怎么做",
    "参考",
)
PROHIBITED_PUBLIC_PATTERNS = (
    "vegan::dune",
    "作者代码通常长这样",
    "本篇不依赖前面章节",
    "这正是全系列坚持",
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

        if number not in PILOT_NUMBERS:
            continue

        text = chapter_path.read_text(encoding="utf-8")
        if metadata.get("draft") is True:
            errors.append(f"article {number:02d} is still marked draft")
        for section in REQUIRED_PILOT_SECTIONS:
            if section not in text:
                errors.append(f"article {number:02d} is missing section: {section}")
        for pattern in PROHIBITED_PUBLIC_PATTERNS:
            if pattern in text:
                errors.append(f"article {number:02d} contains prohibited public text: {pattern}")
        if "set.seed(" not in text:
            errors.append(f"article {number:02d} does not fix a random seed")
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
            "adonis2(",
            "capscale(",
            "01-preview-pcoa",
            "01-preview-cap",
            "01-learning-map",
        ):
            if token not in intro:
                errors.append(f"article 01 is missing executable {token}")

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
