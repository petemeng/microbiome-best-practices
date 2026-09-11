#!/usr/bin/env python3
"""Prepare fixed real-data tables used by Articles 46-50."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "articles-46-50"
SMALL = ROOT / "data" / "small"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path, index_label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=True, index_label=index_label)


def parse_lineage(lineage: str) -> dict[str, str]:
    ranks = {
        "Kingdom": "",
        "Phylum": "",
        "Class": "",
        "Order": "",
        "Family": "",
        "Genus": "",
        "Species": "",
    }
    prefix_map = {
        "d": "Kingdom",
        "k": "Kingdom",
        "p": "Phylum",
        "c": "Class",
        "o": "Order",
        "f": "Family",
        "g": "Genus",
        "s": "Species",
    }
    for token in str(lineage).split(";"):
        token = token.strip()
        match = re.match(r"^([dkpcofgs])__(.*)$", token)
        if match:
            ranks[prefix_map[match.group(1)]] = match.group(2).strip()
    return ranks


def prepare_franzosa() -> dict[str, object]:
    raw_dir = RAW / "franzosa-ibd-2019"
    out_dir = SMALL / "paired-ibd-multiomics"
    out_dir.mkdir(parents=True, exist_ok=True)

    genera = pd.read_csv(raw_dir / "genera.counts.tsv", sep="\t").set_index("Sample")
    metabolites = pd.read_csv(raw_dir / "mtb.tsv", sep="\t").set_index("Sample")
    metadata = pd.read_csv(raw_dir / "metadata.tsv", sep="\t").set_index("Sample")
    metabolite_map = pd.read_csv(raw_dir / "mtb.map.tsv", sep="\t").set_index("Compound")

    samples = metadata.index.intersection(genera.index).intersection(metabolites.index)
    metadata = metadata.loc[samples].copy()
    genera = genera.loc[samples].apply(pd.to_numeric, errors="raise")
    metabolites = metabolites.loc[samples].apply(pd.to_numeric, errors="raise")

    min_samples = int(np.ceil(0.10 * len(samples)))
    genus_prevalence = (genera > 0).sum(axis=0)
    eligible_genera = genus_prevalence[genus_prevalence >= min_samples].index
    genus_rank = genera[eligible_genera].sum(axis=0).sort_values(ascending=False)
    selected_genera = genus_rank.head(250).index.tolist()

    genus_ids = [f"MB_G{i:03d}" for i in range(1, len(selected_genera) + 1)]
    genus_id_map = dict(zip(selected_genera, genus_ids))
    otutab = genera[selected_genera].rename(columns=genus_id_map).T.astype(np.int64)

    taxonomy_rows = []
    for lineage in selected_genera:
        row = parse_lineage(lineage)
        row["Lineage"] = lineage
        row["DisplayName"] = row["Genus"] or row["Family"] or row["Order"] or genus_id_map[lineage]
        taxonomy_rows.append(row)
    taxonomy = pd.DataFrame(taxonomy_rows, index=genus_ids)

    annotation = metabolite_map.reindex(metabolites.columns).copy()
    high_confidence = annotation["High.Confidence.Annotation"].astype(str).str.upper().eq("TRUE")
    annotated = annotation["Compound.Name"].notna() & annotation["Compound.Name"].astype(str).ne("NA")
    eligible_metabolites = annotation.index[high_confidence & annotated]
    log_metabolites = np.log1p(metabolites[eligible_metabolites])
    metabolite_rank = log_metabolites.var(axis=0, ddof=1).sort_values(ascending=False)
    selected_metabolites = metabolite_rank.head(300).index.tolist()

    metabolite_ids = [f"MT_M{i:03d}" for i in range(1, len(selected_metabolites) + 1)]
    metabolite_id_map = dict(zip(selected_metabolites, metabolite_ids))
    metabolite_table = (
        metabolites[selected_metabolites]
        .rename(columns=metabolite_id_map)
        .T.astype(float)
    )
    metabolite_annotation = annotation.loc[selected_metabolites].copy()
    metabolite_annotation.insert(0, "OriginalFeature", selected_metabolites)
    metabolite_annotation.index = metabolite_ids
    metabolite_annotation["DisplayName"] = metabolite_annotation["Compound.Name"].astype(str)

    metadata_out = metadata.rename(
        columns={
            "Study.Group": "StudyGroup",
            "Fecal.Calprotectin": "FecalCalprotectin",
            "Publication.Name": "PublicationName",
        }
    )
    metadata_out.index.name = "SampleID"
    metadata_out["StudyGroup"] = pd.Categorical(
        metadata_out["StudyGroup"], categories=["Control", "CD", "UC"]
    )

    write_tsv(otutab, out_dir / "otutab.tsv", "FeatureID")
    write_tsv(taxonomy, out_dir / "taxonomy.tsv", "FeatureID")
    write_tsv(metadata_out, out_dir / "metadata.tsv", "SampleID")
    write_tsv(metabolite_table, out_dir / "metabolites.tsv", "MetaboliteID")
    write_tsv(metabolite_annotation, out_dir / "metabolite-annotation.tsv", "MetaboliteID")

    summary = {
        "dataset": "FRANZOSA_IBD_2019",
        "citation": "Franzosa et al. Nature Microbiology 2019; doi:10.1038/s41564-018-0306-4",
        "curated_resource": "Muller et al. Scientific Data 2022; doi:10.1038/s41597-022-01644-2",
        "repository": "https://github.com/borenstein-lab/microbiome-metabolome-curated-data",
        "commit": "89a519d8c832008fbc6e650453e83e2f04858d02",
        "raw_sha256": {path.name: sha256(path) for path in sorted(raw_dir.glob("*.tsv"))},
        "samples": len(samples),
        "group_counts": metadata_out["StudyGroup"].value_counts(sort=False).astype(int).to_dict(),
        "original_microbe_features": int(genera.shape[1]),
        "retained_microbe_features": int(otutab.shape[0]),
        "microbe_filter": "prevalence >=10%; top 250 by total count",
        "original_metabolite_features": int(metabolites.shape[1]),
        "retained_metabolite_features": int(metabolite_table.shape[0]),
        "metabolite_filter": "high-confidence named features; up to top 300 by log1p variance",
        "contract": "otutab/taxonomy/metadata plus paired metabolites and metabolite annotation",
    }
    (out_dir / "source-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def read_taxonomy_without_header(path: Path, feature_ids: pd.Index) -> pd.DataFrame:
    parsed: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields:
                continue
            feature_id = fields[0]
            lineage_tokens = [token for token in fields[1:] if re.match(r"^[dkpcofgs]__", token)]
            row = parse_lineage(";".join(lineage_tokens))
            row["Lineage"] = ";".join(lineage_tokens)
            row["DisplayName"] = row["Species"] or row["Genus"] or row["Family"] or feature_id
            parsed[feature_id] = row
    taxonomy = pd.DataFrame.from_dict(parsed, orient="index").reindex(feature_ids)
    taxonomy.index.name = "FeatureID"
    return taxonomy.fillna("")


def load_duran_view(raw_dir: Path, count_name: str, taxonomy_name: str, design_name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(raw_dir / count_name, sep="\t", index_col=0)
    design = pd.read_csv(raw_dir / design_name, sep="\t")
    design["Replicate"] = design["SampleID"].str.extract(r"\.(\d+)$")[0].astype(int)
    design["PairID"] = (
        design["soil"].astype(str)
        + "_"
        + design["compartment"].astype(str)
        + "_R"
        + design["Replicate"].astype(str)
    )
    taxonomy = read_taxonomy_without_header(raw_dir / taxonomy_name, counts.index)
    return counts, taxonomy, design


def prepare_duran() -> dict[str, object]:
    raw_dir = RAW / "duran-interkingdom-2018"
    out_dir = SMALL / "multi-kingdom-duran"
    out_dir.mkdir(parents=True, exist_ok=True)

    view_specs = {
        "bacteria": (
            "otu_tab_filter_001_bac.txt",
            "taxonomy_ref_bac.txt",
            "design_own_bac.txt",
        ),
        "fungi": (
            "otu_tab_filter_001_fun.txt",
            "taxonomy_ref_fun.txt",
            "design_own_fun.txt",
        ),
        "oomycetes": (
            "otu_tab_filter_001_oo.txt",
            "taxonomy_ref_oo.txt",
            "design_own_oomyc.txt",
        ),
    }
    loaded = {name: load_duran_view(raw_dir, *spec) for name, spec in view_specs.items()}
    common_pair_ids = sorted(
        set.intersection(*(set(design["PairID"]) for _, _, design in loaded.values()))
    )
    if len(common_pair_ids) != 36:
        raise RuntimeError(f"Expected 36 three-view pairs, found {len(common_pair_ids)}")

    output_names = {
        "bacteria": ("otutab.tsv", "taxonomy.tsv"),
        "fungi": ("fungi-otutab.tsv", "fungi-taxonomy.tsv"),
        "oomycetes": ("oomycete-otutab.tsv", "oomycete-taxonomy.tsv"),
    }
    paired_designs = {}
    view_shapes = {}
    for name, (counts, taxonomy, design) in loaded.items():
        paired_design = design[design["PairID"].isin(common_pair_ids)].copy()
        if paired_design["PairID"].duplicated().any():
            raise RuntimeError(f"Duplicated paired biological sample in {name}")
        paired_design = paired_design.set_index("PairID").loc[common_pair_ids]
        sample_ids = paired_design["SampleID"].tolist()
        paired_counts = counts.loc[:, sample_ids].copy()
        paired_counts.columns = common_pair_ids
        paired_counts = paired_counts.astype(np.int64)
        count_file, taxonomy_file = output_names[name]
        write_tsv(paired_counts, out_dir / count_file, "FeatureID")
        write_tsv(taxonomy.loc[paired_counts.index], out_dir / taxonomy_file, "FeatureID")
        paired_designs[name] = paired_design
        view_shapes[name] = [int(paired_counts.shape[0]), int(paired_counts.shape[1])]

    bacterial_design = paired_designs["bacteria"]
    metadata = pd.DataFrame(index=common_pair_ids)
    metadata["Soil"] = bacterial_design["soil"]
    compartment_map = {"soil": "Soil", "rhizosphere": "Rhizosphere", "root": "Root"}
    metadata["Compartment"] = bacterial_design["compartment"].map(compartment_map)
    metadata["Replicate"] = bacterial_design["Replicate"].astype(int)
    metadata["BacterialSampleID"] = paired_designs["bacteria"]["SampleID"]
    metadata["FungalSampleID"] = paired_designs["fungi"]["SampleID"]
    metadata["OomyceteSampleID"] = paired_designs["oomycetes"]["SampleID"]
    write_tsv(metadata, out_dir / "metadata.tsv", "SampleID")

    summary = {
        "dataset": "Duran_Thiergart_Arabidopsis_interkingdom_2018",
        "citation": "Duran et al. Cell 2018; doi:10.1016/j.cell.2018.10.020",
        "repository": "https://github.com/ththi/Microbial-Interkingdom-Suppl",
        "commit": "6db5e85cc5d442fd95fcdcb7250b72fa9e2ff900",
        "raw_sha256": {path.name: sha256(path) for path in sorted(raw_dir.glob("*.txt"))},
        "paired_samples": len(common_pair_ids),
        "soil_counts": metadata["Soil"].value_counts().sort_index().astype(int).to_dict(),
        "compartment_counts": metadata["Compartment"].value_counts().sort_index().astype(int).to_dict(),
        "view_shapes_features_by_samples": view_shapes,
        "pairing_key": "soil + compartment + biological replicate",
        "contract": "bacterial three-piece set plus paired fungal and oomycete count/taxonomy tables",
    }
    (out_dir / "source-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def prepare_feast() -> dict[str, object]:
    raw_dir = RAW / "feast-example"
    out_dir = SMALL / "source-tracking-feast"
    out_dir.mkdir(parents=True, exist_ok=True)

    otutab = pd.read_csv(raw_dir / "otu_example.txt", sep="\t", index_col=0)
    metadata = pd.read_csv(raw_dir / "metadata_example.txt", sep="\t").set_index("SampleID")
    otutab = otutab.loc[:, metadata.index].astype(np.int64)

    def source_class(env: str) -> str:
        value = env.lower()
        if value.startswith("infant gut"):
            return "Infant gut"
        if value.startswith("adult gut"):
            return "Adult gut"
        if value.startswith("adult skin"):
            return "Adult skin"
        if value.startswith("soil"):
            return "Soil"
        return "Other"

    metadata["SourceClass"] = metadata["Env"].map(source_class)
    metadata["FeastID"] = 1
    taxonomy = pd.DataFrame(index=otutab.index)
    for rank in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]:
        taxonomy[rank] = ""
    taxonomy["DisplayName"] = taxonomy.index
    taxonomy["TaxonomyStatus"] = "Anonymous taxon in the official FEAST example"

    write_tsv(otutab, out_dir / "otutab.tsv", "FeatureID")
    write_tsv(taxonomy, out_dir / "taxonomy.tsv", "FeatureID")
    write_tsv(metadata, out_dir / "metadata.tsv", "SampleID")

    summary = {
        "dataset": "FEAST_official_source_tracking_example",
        "citation": "Shenhav et al. Nature Methods 2019; doi:10.1038/s41592-019-0431-x",
        "repository": "https://github.com/cozygene/FEAST",
        "commit": "2f8f3df8051e0e08341f597a9f4693bfb76b3bf6",
        "raw_sha256": {path.name: sha256(path) for path in sorted(raw_dir.glob("*.txt"))},
        "samples": int(otutab.shape[1]),
        "features": int(otutab.shape[0]),
        "sink_samples": int((metadata["SourceSink"] == "Sink").sum()),
        "source_samples": int((metadata["SourceSink"] == "Source").sum()),
        "source_class_counts": metadata.loc[metadata["SourceSink"] == "Source", "SourceClass"].value_counts().sort_index().astype(int).to_dict(),
        "taxonomy_scope": "The official demo anonymises taxa; taxonomy rows are retained explicitly as unknown.",
        "contract": "otutab/taxonomy/metadata with source/sink labels and candidate source classes",
    }
    (out_dir / "source-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    summaries = {
        "articles_46_48": prepare_franzosa(),
        "article_49": prepare_duran(),
        "article_50": prepare_feast(),
    }
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
