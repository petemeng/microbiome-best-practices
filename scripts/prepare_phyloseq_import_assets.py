#!/usr/bin/env python3
"""Freeze a self-contained Article 17 phyloseq input bundle from QIIME 2 Artifacts.

Run this script with the Python interpreter from the locked QIIME 2 environment,
where ``biom-format`` and ``DendroPy`` are available.  The generated bundle is
the reader-facing input contract; Article 17 itself does not read an earlier
chapter's RDS or QIIME 2 result directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tempfile
import zipfile
from pathlib import Path
from typing import Iterable

import biom
import dendropy


RANKS = (
    ("d", "Domain"),
    ("p", "Phylum"),
    ("c", "Class"),
    ("o", "Order"),
    ("f", "Family"),
    ("g", "Genus"),
    ("s", "Species"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-artifact", required=True, type=Path)
    parser.add_argument("--taxonomy-artifact", required=True, type=Path)
    parser.add_argument("--sequences-artifact", required=True, type=Path)
    parser.add_argument("--tree-artifact", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def qza_uuid(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        roots = {name.split("/", 1)[0] for name in archive.namelist() if "/" in name}
    if len(roots) != 1:
        raise ValueError(f"Expected one Artifact root in {path}, found {sorted(roots)}")
    return next(iter(roots))


def qza_member(path: Path, suffix: str) -> bytes:
    with zipfile.ZipFile(path) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one member ending with {suffix!r} in {path}, found {matches}"
            )
        return archive.read(matches[0])


def write_tsv(path: Path, header: Iterable[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_deterministic_gzip(path: Path, payload: bytes) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as handle:
            handle.write(payload)


def read_fasta(payload: bytes) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    for raw_line in payload.decode("utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            if current in records:
                raise ValueError(f"Duplicate FASTA identifier: {current}")
            records[current] = []
        elif current is None:
            raise ValueError("FASTA sequence occurred before the first header")
        else:
            records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def format_fasta(records: dict[str, str], order: list[str]) -> bytes:
    lines: list[str] = []
    for name in order:
        sequence = records[name]
        lines.append(f">{name}")
        lines.extend(sequence[i : i + 80] for i in range(0, len(sequence), 80))
    return ("\n".join(lines) + "\n").encode("utf-8")


def parse_taxonomy(value: str) -> list[str]:
    parsed: dict[str, str] = {}
    for token in value.split(";"):
        token = token.strip()
        match = re.match(r"^([dkpcofgs])__(.*)$", token)
        if match:
            rank, label = match.groups()
            parsed[rank] = label.strip() or "Unassigned"
    return [parsed.get(prefix, "Unassigned") for prefix, _ in RANKS]


def main() -> None:
    args = parse_args()
    source_paths = {
        "feature_table": args.table_artifact.resolve(strict=True),
        "taxonomy": args.taxonomy_artifact.resolve(strict=True),
        "representative_sequences": args.sequences_artifact.resolve(strict=True),
        "rooted_tree": args.tree_artifact.resolve(strict=True),
    }
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="article17-biom-") as temp_dir:
        biom_payload = qza_member(source_paths["feature_table"], "/data/feature-table.biom")
        biom_path = Path(temp_dir) / "feature-table.biom"
        biom_path.write_bytes(biom_payload)
        table = biom.load_table(str(biom_path))

    sample_ids = [str(value) for value in table.ids(axis="sample")]
    feature_ids = [str(value) for value in table.ids(axis="observation")]
    counts_by_feature = {
        feature_id: [int(round(value)) for value in table.data(feature_id, axis="observation")]
        for feature_id in feature_ids
    }
    feature_ids = sorted(feature_ids, key=lambda item: (-sum(counts_by_feature[item]), item))

    if table.shape != (366, 4):
        raise ValueError(f"Expected a 366 x 4 table, observed {table.shape}")
    if sum(sum(counts_by_feature[item]) for item in feature_ids) != 5213:
        raise ValueError("Expected exactly 5,213 reads in the frozen table")
    if sample_ids != ["1", "1a", "2", "2a"]:
        raise ValueError(f"Unexpected sample order: {sample_ids}")

    taxonomy_lines = qza_member(
        source_paths["taxonomy"], "/data/taxonomy.tsv"
    ).decode("utf-8").splitlines()
    taxonomy_reader = csv.DictReader(taxonomy_lines, delimiter="\t")
    taxonomy_raw = {
        row["Feature ID"]: {
            "ranks": parse_taxonomy(row["Taxon"]),
            "confidence": float(row["Confidence"]),
            "raw": row["Taxon"],
        }
        for row in taxonomy_reader
    }

    sequences = read_fasta(
        qza_member(source_paths["representative_sequences"], "/data/dna-sequences.fasta")
    )
    tree_payload = qza_member(source_paths["rooted_tree"], "/data/tree.nwk")
    tree_text = tree_payload.decode("utf-8").strip() + "\n"
    tree = dendropy.Tree.get(data=tree_text, schema="newick", preserve_underscores=True)
    tree_tip_ids = {leaf.taxon.label for leaf in tree.leaf_node_iter()}

    feature_set = set(feature_ids)
    if set(taxonomy_raw) != feature_set:
        raise ValueError("Taxonomy feature IDs do not exactly match the count table")
    if set(sequences) != feature_set:
        raise ValueError("Representative-sequence IDs do not exactly match the count table")
    if not feature_set.issubset(tree_tip_ids):
        raise ValueError("The full insertion tree is missing one or more count-table features")
    if len(tree_tip_ids) != 203818:
        raise ValueError(f"Expected 203,818 tree tips, observed {len(tree_tip_ids)}")

    otutab_path = output_dir / "otutab.tsv"
    taxonomy_path = output_dir / "taxonomy.tsv"
    confidence_path = output_dir / "taxonomy-confidence.tsv"
    metadata_path = output_dir / "metadata.tsv"
    sequences_path = output_dir / "representative-sequences.fasta.gz"
    tree_path = output_dir / "rooted-insertion-tree.nwk.gz"

    write_tsv(
        otutab_path,
        ["FeatureID", *sample_ids],
        ([feature_id, *counts_by_feature[feature_id]] for feature_id in feature_ids),
    )
    write_tsv(
        taxonomy_path,
        ["FeatureID", *(name for _, name in RANKS)],
        ([feature_id, *taxonomy_raw[feature_id]["ranks"]] for feature_id in feature_ids),
    )
    write_tsv(
        confidence_path,
        ["FeatureID", "Confidence", "RawTaxonomy"],
        (
            [
                feature_id,
                format(taxonomy_raw[feature_id]["confidence"], ".15g"),
                taxonomy_raw[feature_id]["raw"],
            ]
            for feature_id in feature_ids
        ),
    )
    metadata_rows = [
        ["1", "Library_1", "S103", "base", "16S_V4", "PairedEnd"],
        ["1a", "Library_1a", "S103", "suffix_a", "16S_V4", "PairedEnd"],
        ["2", "Library_2", "S115", "base", "16S_V4", "PairedEnd"],
        ["2a", "Library_2a", "S115", "suffix_a", "16S_V4", "PairedEnd"],
    ]
    write_tsv(
        metadata_path,
        [
            "SampleID",
            "LibraryLabel",
            "FilenameSet",
            "FilenameVariant",
            "AmpliconRegion",
            "ReadLayout",
        ],
        metadata_rows,
    )
    write_deterministic_gzip(sequences_path, format_fasta(sequences, feature_ids))
    write_deterministic_gzip(tree_path, tree_text.encode("utf-8"))

    output_paths = [
        otutab_path,
        taxonomy_path,
        confidence_path,
        metadata_path,
        sequences_path,
        tree_path,
    ]
    source_summary = {
        "schema_version": 1,
        "bundle_id": "nfcore-v4-phyloseq-input-v1",
        "purpose": "Self-contained Article 17 phyloseq import contract",
        "data_boundary": (
            "The four library labels are nf-core test-data file identifiers. "
            "FilenameSet and FilenameVariant describe filenames only and are not biological groups."
        ),
        "counts_and_sequences": {
            "source": "nf-core/test-datasets ampliseq testdata",
            "commit": "f282ad2bafb3ca90190ec87290066ee1985df655",
            "license": "MIT",
            "processing": "QIIME 2/q2-dada2 2026.4.0; 220/200 truncation; maxEE 2/2",
        },
        "taxonomy": {
            "reference": "SILVA 138.2 SSURef NR99 515F/806R V4 LCA reference",
            "license": "CC BY 4.0",
            "method": "q2-feature-classifier classify-sklearn, confidence 0.70",
        },
        "tree": {
            "reference": "QIIME 2 Greengenes 13_8 99% SEPP reference database",
            "method": "q2-fragment-insertion sepp, alignment subset 1000, placement subset 5000",
            "full_tree_tips": len(tree_tip_ids),
            "query_tips": len(feature_set),
            "reference_tips": len(tree_tip_ids - feature_set),
            "note": "The bundle contains the complete insertion tree, not a query-only diagnostic tree.",
        },
        "dimensions": {
            "features": len(feature_ids),
            "samples": len(sample_ids),
            "reads": sum(sum(counts_by_feature[item]) for item in feature_ids),
            "taxonomy_ranks": len(RANKS),
            "representative_sequences": len(sequences),
        },
        "source_artifacts": {
            key: {
                "filename": path.name,
                "sha256": sha256(path),
                "uuid": qza_uuid(path),
            }
            for key, path in source_paths.items()
        },
        "files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in output_paths
        },
    }
    summary_path = output_dir / "source-summary.json"
    summary_path.write_text(
        json.dumps(source_summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps({
        "output_dir": str(output_dir),
        "features": len(feature_ids),
        "samples": len(sample_ids),
        "reads": source_summary["dimensions"]["reads"],
        "tree_tips": len(tree_tip_ids),
        "source_summary_sha256": sha256(summary_path),
    }, indent=2))


if __name__ == "__main__":
    main()
