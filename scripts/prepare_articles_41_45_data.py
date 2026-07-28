#!/usr/bin/env python3
"""Prepare fixed real-data triads for Articles 41, 44 and 45.

The script never downloads data.  It converts already downloaded, checksum-
verified upstream files into the project-wide otutab/taxonomy/metadata contract.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


RANKS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_PREFIX = {
    "k__": "Kingdom",
    "p__": "Phylum",
    "c__": "Class",
    "o__": "Order",
    "f__": "Family",
    "g__": "Genus",
    "s__": "Species",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(frame: pd.DataFrame, path: Path, index_label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=True, index_label=index_label, lineterminator="\n")


def deterministic_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as dst:
            shutil.copyfileobj(src, dst)


def decode_strings(values: np.ndarray) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_biom_hdf5(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        feature_ids = decode_strings(handle["observation/ids"][:])
        sample_ids = decode_strings(handle["sample/ids"][:])
        values = handle["observation/matrix/data"][:]
        indices = handle["observation/matrix/indices"][:]
        indptr = handle["observation/matrix/indptr"][:]
        shape = tuple(int(value) for value in handle.attrs["shape"])
    matrix = sparse.csr_matrix((values, indices, indptr), shape=shape).toarray()
    if not np.allclose(matrix, np.rint(matrix)):
        raise ValueError(f"Non-integer counts found in {path}")
    return pd.DataFrame(np.rint(matrix).astype(np.int64), index=feature_ids, columns=sample_ids)


def prepare_picrust2(source: Path, output_root: Path, archive: Path) -> dict:
    output = output_root / "picrust2-chemerin"
    output.mkdir(parents=True, exist_ok=True)
    biom_path = source / "table.biom"
    fasta_path = source / "seqs.fna"
    metadata_path = source / "metadata.tsv"
    counts = read_biom_hdf5(biom_path)
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str).set_index("SampleID")
    missing = set(counts.columns) - set(metadata.index)
    if missing:
        raise ValueError(f"PICRUSt2 metadata is missing samples: {sorted(missing)}")
    metadata = metadata.loc[counts.columns, ["Facility", "Genotype"]]
    taxonomy = pd.DataFrame("", index=counts.index, columns=RANKS)
    write_tsv(counts, output / "otutab.tsv", "FeatureID")
    write_tsv(taxonomy, output / "taxonomy.tsv", "FeatureID")
    write_tsv(metadata, output / "metadata.tsv", "SampleID")
    shutil.copyfile(fasta_path, output / "representative-sequences.fasta")
    deterministic_gzip(fasta_path, output / "representative-sequences.fasta.gz")
    shutil.copyfile(biom_path, output / "table.biom")
    payload = {
        "dataset": "PICRUSt2 chemerin 16S tutorial data",
        "source_url": "http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/chemerin_16S.zip",
        "tutorial_url": "https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.5.0)",
        "upstream_archive": {"sha256": sha256(archive), "bytes": archive.stat().st_size},
        "dimensions": {
            "features": int(counts.shape[0]),
            "samples": int(counts.shape[1]),
            "reads": int(counts.to_numpy().sum()),
            "genotype_counts": metadata["Genotype"].value_counts().sort_index().to_dict(),
            "facility_counts": metadata["Facility"].value_counts().sort_index().to_dict(),
        },
        "taxonomy_boundary": "The official tutorial archive does not distribute taxonomy; taxonomy.tsv intentionally contains blank ranks because PICRUSt2 uses representative sequences.",
        "files": {},
    }
    for name in (
        "otutab.tsv",
        "taxonomy.tsv",
        "metadata.tsv",
        "representative-sequences.fasta",
        "representative-sequences.fasta.gz",
        "table.biom",
    ):
        path = output / name
        payload["files"][name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    write_json(payload, output / "source-summary.json")
    return payload


def parse_taxonomy_label(label: str, cohort: str, ordinal: int) -> tuple[str, dict[str, str]]:
    parsed = {rank: "" for rank in RANKS}
    denovo = ""
    for field in str(label).split(";"):
        field = field.strip()
        for prefix, rank in RANK_PREFIX.items():
            if field.startswith(prefix):
                parsed[rank] = field[len(prefix):].strip()
                break
        if field.startswith("d__"):
            denovo = field[3:].strip()
    if not denovo:
        denovo = hashlib.sha1(str(label).encode("utf-8")).hexdigest()[:14]
    feature_id = f"{cohort}_{denovo}"
    if not feature_id:
        feature_id = f"{cohort}_feature_{ordinal:06d}"
    parsed["OriginalFeatureLabel"] = str(label)
    return feature_id, parsed


def read_microbiomehd_metadata(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8", encoding_errors="replace")
    first = frame.columns[0]
    frame = frame.rename(columns={first: "SampleID"})
    frame["SampleID"] = frame["SampleID"].astype(str).str.strip()
    return frame.set_index("SampleID")


def prepare_crc_cohort(cohort: str, source: Path, output_root: Path) -> dict:
    table_path = source / "RDP" / f"{cohort}.otu_table.100.denovo.rdp_assigned"
    metadata_path = source / f"{cohort}.metadata.txt"
    raw_counts = pd.read_csv(table_path, sep="\t", index_col=0)
    raw_counts.columns = raw_counts.columns.astype(str).str.strip()
    metadata = read_microbiomehd_metadata(metadata_path)

    if cohort == "crc_zackular":
        group = pd.Series(index=metadata.index, dtype="object")
        group.loc[group.index.str.lower().str.startswith("cancer")] = "CRC"
        group.loc[group.index.str.lower().str.startswith("healthy")] = "Control"
    else:
        disease = metadata["DiseaseState"].astype(str).str.upper()
        group = disease.map({"CRC": "CRC", "H": "Control", "HEALTHY": "Control", "CONTROL": "Control"})
    metadata = metadata.assign(Group=group)
    keep_samples = [sample for sample in raw_counts.columns if sample in metadata.index and pd.notna(metadata.at[sample, "Group"])]
    if not keep_samples:
        raise ValueError(f"No CRC/control samples align for {cohort}")
    raw_counts = raw_counts.loc[:, keep_samples]
    metadata = metadata.loc[keep_samples]

    records: list[dict[str, str]] = []
    feature_ids: list[str] = []
    seen: dict[str, int] = {}
    for ordinal, label in enumerate(raw_counts.index, start=1):
        feature_id, parsed = parse_taxonomy_label(str(label), cohort.replace("crc_", ""), ordinal)
        seen[feature_id] = seen.get(feature_id, 0) + 1
        if seen[feature_id] > 1:
            feature_id = f"{feature_id}_{seen[feature_id]}"
        feature_ids.append(feature_id)
        records.append(parsed)
    raw_counts.index = feature_ids
    taxonomy = pd.DataFrame(records, index=feature_ids)
    raw_counts = raw_counts.loc[raw_counts.sum(axis=1) > 0]
    taxonomy = taxonomy.loc[raw_counts.index]
    metadata_out = pd.DataFrame(index=metadata.index)
    metadata_out["Study"] = cohort
    metadata_out["Group"] = metadata["Group"]
    for column in ("DiseaseState", "Amplicon", "Platform", "PatientGeo"):
        metadata_out[column] = metadata[column] if column in metadata.columns else ""

    destination = output_root / "cross-cohort-crc" / cohort
    write_tsv(raw_counts.astype(np.int64), destination / "otutab.tsv", "FeatureID")
    write_tsv(taxonomy, destination / "taxonomy.tsv", "FeatureID")
    write_tsv(metadata_out, destination / "metadata.tsv", "SampleID")
    return {
        "features": int(raw_counts.shape[0]),
        "samples": int(raw_counts.shape[1]),
        "reads": int(raw_counts.to_numpy().sum()),
        "groups": metadata_out["Group"].value_counts().sort_index().to_dict(),
        "genus_assigned_features": int(taxonomy["Genus"].fillna("").ne("").sum()),
        "files": {
            name: {"sha256": sha256(destination / name), "bytes": (destination / name).stat().st_size}
            for name in ("otutab.tsv", "taxonomy.tsv", "metadata.tsv")
        },
    }


def prepare_microbiomehd(source_root: Path, output_root: Path, archives: dict[str, Path]) -> dict:
    studies = {}
    for cohort in ("crc_xiang", "crc_zhao", "crc_zackular"):
        studies[cohort] = prepare_crc_cohort(cohort, source_root / f"{cohort}_results", output_root)
    payload = {
        "dataset": "MicrobiomeHD standardized CRC 16S cohorts",
        "record_doi": "10.5281/zenodo.1146764",
        "record_url": "https://zenodo.org/records/1146764",
        "database_citation": "Duvallet C, Gibbons SM, Gurry T, et al. Nature Communications. 2017;8:1784.",
        "license": "CC BY-NC 4.0",
        "analysis_boundary": "Only CRC and healthy-control stool samples are retained; adenoma and mock samples are excluded before analysis.",
        "archives": {
            name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for name, path in archives.items()
        },
        "studies": studies,
    }
    write_json(payload, output_root / "cross-cohort-crc" / "source-summary.json")
    return payload


def prepare_survival(source: Path, output_root: Path, commit: str) -> dict:
    counts = pd.read_csv(source / "otu.tab.txt", sep="\t", index_col=0)
    taxonomy = pd.read_csv(source / "tax.tab.txt", sep="\t", index_col=0, dtype=str).fillna("")
    metadata = pd.read_csv(source / "sam.dat.txt", sep="\t", index_col=0)
    counts.columns = counts.columns.astype(str)
    common = [sample for sample in counts.columns if sample in metadata.index]
    counts = counts.loc[:, common]
    metadata = metadata.loc[common]
    taxonomy = taxonomy.loc[counts.index]
    def sampling_week(sample_id: str) -> int:
        for pattern in (r"F([0-9]+)W$", r"FW([0-9]+)$"):
            match = re.search(pattern, sample_id)
            if match:
                return int(match.group(1))
        raise ValueError(f"Cannot parse sampling week from {sample_id}")

    metadata_out = pd.DataFrame(index=metadata.index)
    metadata_out["TimeWeeks"] = pd.to_numeric(metadata["T1Dweek"], errors="raise").astype(float)
    metadata_out["SamplingWeek"] = [sampling_week(sample_id) for sample_id in metadata_out.index]
    metadata_out["FollowUpWeeks"] = metadata_out["TimeWeeks"] - metadata_out["SamplingWeek"]
    if not (metadata_out["FollowUpWeeks"] > 0).all():
        raise ValueError("Every event/censor time must occur after baseline sampling")
    metadata_out["Event"] = pd.to_numeric(metadata["T1D"], errors="raise").astype(int)
    metadata_out["Sex"] = metadata["Sex"].astype(str)
    metadata_out["Treatment"] = metadata["Antibiotics"].astype(str)
    destination = output_root / "survival-t1d"
    write_tsv(counts.astype(np.int64), destination / "otutab.tsv", "FeatureID")
    write_tsv(taxonomy.loc[:, RANKS], destination / "taxonomy.tsv", "FeatureID")
    write_tsv(metadata_out, destination / "metadata.tsv", "SampleID")
    payload = {
        "dataset": "MiSurv processed NOD-mouse gut 16S data with time to T1D onset",
        "source_repository": "https://github.com/wg99526/MiSurvGit",
        "source_commit": commit,
        "qiita_study": "10508",
        "primary_study": "Zhang XS et al. Nature Microbiology. 2018;3:1403-1413. doi:10.1038/s41564-018-0251-7.",
        "processed_data_citation": "Gu W et al. Microbiology Spectrum. 2023;11:e05059-22. doi:10.1128/spectrum.05059-22.",
        "license_boundary": "The source GitHub repository does not state a data license. Cite both studies and verify redistribution terms before republishing the source tables.",
        "outcome_boundary": "The event is T1D onset in NOD mice, not human death or overall survival.",
        "dimensions": {
            "features": int(counts.shape[0]),
            "samples": int(counts.shape[1]),
            "reads": int(counts.to_numpy().sum()),
            "events": int(metadata_out["Event"].sum()),
            "censored": int((metadata_out["Event"] == 0).sum()),
            "sampling_week_counts": metadata_out["SamplingWeek"].value_counts().sort_index().to_dict(),
            "treatment_counts": metadata_out["Treatment"].value_counts().sort_index().to_dict(),
        },
        "source_files": {},
        "files": {},
    }
    for name in ("otu.tab.txt", "tax.tab.txt", "sam.dat.txt", "tree.tre"):
        path = source / name
        payload["source_files"][name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    for name in ("otutab.tsv", "taxonomy.tsv", "metadata.tsv"):
        path = destination / name
        payload["files"][name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    write_json(payload, destination / "source-summary.json")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--picrust-source", type=Path, required=True)
    parser.add_argument("--picrust-archive", type=Path, required=True)
    parser.add_argument("--microbiomehd-source", type=Path, required=True)
    parser.add_argument("--crc-xiang-archive", type=Path, required=True)
    parser.add_argument("--crc-zhao-archive", type=Path, required=True)
    parser.add_argument("--crc-zackular-archive", type=Path, required=True)
    parser.add_argument("--survival-source", type=Path, required=True)
    parser.add_argument("--survival-commit", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    report = {
        "picrust2": prepare_picrust2(args.picrust_source, args.output_root, args.picrust_archive),
        "cross_cohort_crc": prepare_microbiomehd(
            args.microbiomehd_source,
            args.output_root,
            {
                "crc_xiang_results.tar.gz": args.crc_xiang_archive,
                "crc_zhao_results.tar.gz": args.crc_zhao_archive,
                "crc_zackular_results.tar.gz": args.crc_zackular_archive,
            },
        ),
        "survival_t1d": prepare_survival(args.survival_source, args.output_root, args.survival_commit),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
