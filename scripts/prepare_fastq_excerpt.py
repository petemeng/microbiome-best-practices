#!/usr/bin/env python3
"""Create a deterministic, synchronized excerpt of the QIIME 2 Atacama FASTQs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import statistics
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO


SOURCE_FILES = {
    "forward": {
        "url": "https://data.qiime2.org/2024.5/tutorials/atacama-soils/1p/forward.fastq.gz",
        "sha256": "392f178b5f15967fafb9332c0493a8b6ecbf4a7857538023813a630c488c2702",
    },
    "reverse": {
        "url": "https://data.qiime2.org/2024.5/tutorials/atacama-soils/1p/reverse.fastq.gz",
        "sha256": "5407370313a3974e5b70878153841a61df408e49ced133b3ace44b2db286cfae",
    },
    "barcodes": {
        "url": "https://data.qiime2.org/2024.5/tutorials/atacama-soils/1p/barcodes.fastq.gz",
        "sha256": "c07f517e920d9e3924ce1e87102f057232b57041947523b90d560f68a0ff3c0f",
    },
    "metadata": {
        "url": "https://data.qiime2.org/2024.5/tutorials/atacama-soils/sample_metadata.tsv",
        "sha256": "7cff810ad86a621ebc78a16b690255d839ea58e2622ad45a11a838d0f8c5b3bd",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--forward", type=Path, required=True)
    parser.add_argument("--reverse", type=Path, required=True)
    parser.add_argument("--barcodes", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--records", type=int, default=2000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fastq_records(handle: TextIO, source: str) -> Iterator[tuple[str, str, str, str]]:
    record_number = 0
    while True:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        record_number += 1
        if not sequence or not separator or not quality:
            raise ValueError(f"{source}: truncated FASTQ record {record_number}")
        fields = tuple(line.rstrip("\r\n") for line in (header, sequence, separator, quality))
        header_text, sequence_text, separator_text, quality_text = fields
        if not header_text.startswith("@"):
            raise ValueError(f"{source}: record {record_number} header does not start with @")
        if not separator_text.startswith("+"):
            raise ValueError(f"{source}: record {record_number} separator does not start with +")
        if len(sequence_text) != len(quality_text):
            raise ValueError(
                f"{source}: record {record_number} sequence/quality lengths differ "
                f"({len(sequence_text)} != {len(quality_text)})"
            )
        yield fields


def count_records(path: Path, source: str) -> int:
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        return sum(1 for _ in fastq_records(handle, source))


def read_id(record: tuple[str, str, str, str]) -> str:
    return record[0].split(maxsplit=1)[0].removeprefix("@")


def write_record(handle: BinaryIO, record: tuple[str, str, str, str]) -> None:
    handle.write(("\n".join(record) + "\n").encode("ascii"))


def summarize_reads(records: list[tuple[str, str, str, str]]) -> dict[str, float | int]:
    lengths = [len(record[1]) for record in records]
    mean_q = [
        sum(ord(character) - 33 for character in record[3]) / len(record[3])
        for record in records
    ]
    bases = sum(lengths)
    q30_bases = sum(
        sum((ord(character) - 33) >= 30 for character in record[3])
        for record in records
    )
    n_bases = sum(record[1].count("N") for record in records)
    return {
        "records": len(records),
        "read_length_min": min(lengths),
        "read_length_median": statistics.median(lengths),
        "read_length_max": max(lengths),
        "median_read_mean_phred": round(statistics.median(mean_q), 6),
        "q30_base_fraction": round(q30_bases / bases, 8),
        "n_base_fraction": round(n_bases / bases, 8),
    }


def main() -> int:
    args = parse_args()
    if args.records < 2:
        raise ValueError("--records must be at least 2")

    inputs = {
        "forward": args.forward.resolve(),
        "reverse": args.reverse.resolve(),
        "barcodes": args.barcodes.resolve(),
        "metadata": args.metadata.resolve(),
    }
    for key, path in inputs.items():
        observed = sha256(path)
        expected = SOURCE_FILES[key]["sha256"]
        if observed != expected:
            raise ValueError(f"{key} sha256 mismatch: expected {expected}, observed {observed}")

    counts = {
        key: count_records(inputs[key], key)
        for key in ("forward", "reverse", "barcodes")
    }
    if len(set(counts.values())) != 1:
        raise ValueError(f"source FASTQ record counts differ: {counts}")
    total_records = counts["forward"]
    if args.records > total_records:
        raise ValueError(
            f"requested {args.records} records but source contains {total_records}"
        )

    selected_indices = {
        round(index * (total_records - 1) / (args.records - 1))
        for index in range(args.records)
    }
    if len(selected_indices) != args.records:
        raise ValueError("evenly spaced record selection produced duplicate indices")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {
        key: args.output_dir / f"{key}.fastq.gz"
        for key in ("forward", "reverse", "barcodes")
    }
    selected_records: dict[str, list[tuple[str, str, str, str]]] = {
        key: [] for key in ("forward", "reverse", "barcodes")
    }

    with (
        gzip.open(inputs["forward"], "rt", encoding="ascii", newline="") as forward_handle,
        gzip.open(inputs["reverse"], "rt", encoding="ascii", newline="") as reverse_handle,
        gzip.open(inputs["barcodes"], "rt", encoding="ascii", newline="") as barcode_handle,
        output_paths["forward"].open("wb") as forward_raw,
        output_paths["reverse"].open("wb") as reverse_raw,
        output_paths["barcodes"].open("wb") as barcode_raw,
        gzip.GzipFile(fileobj=forward_raw, mode="wb", filename="", mtime=0) as forward_out,
        gzip.GzipFile(fileobj=reverse_raw, mode="wb", filename="", mtime=0) as reverse_out,
        gzip.GzipFile(fileobj=barcode_raw, mode="wb", filename="", mtime=0) as barcode_out,
    ):
        streams = (
            fastq_records(forward_handle, "forward"),
            fastq_records(reverse_handle, "reverse"),
            fastq_records(barcode_handle, "barcodes"),
        )
        for record_index, triplet in enumerate(itertools.zip_longest(*streams)):
            if any(record is None for record in triplet):
                raise ValueError("source FASTQ streams ended at different records")
            forward_record, reverse_record, barcode_record = triplet
            identifiers = {
                read_id(forward_record),
                read_id(reverse_record),
                read_id(barcode_record),
            }
            if len(identifiers) != 1:
                raise ValueError(
                    f"record {record_index + 1} has non-matching read IDs: {identifiers}"
                )
            if record_index not in selected_indices:
                continue
            for key, record, output_handle in (
                ("forward", forward_record, forward_out),
                ("reverse", reverse_record, reverse_out),
                ("barcodes", barcode_record, barcode_out),
            ):
                selected_records[key].append(record)
                write_record(output_handle, record)

    metadata_text = inputs["metadata"].read_text(encoding="ascii")
    metadata_path = args.output_dir / "metadata.tsv"
    metadata_path.write_text(
        "\n".join(metadata_text.splitlines()) + "\n",
        encoding="ascii",
        newline="\n",
    )
    with metadata_path.open(encoding="ascii", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    metadata_samples = sum(
        bool(row) and not row[0].startswith("#")
        for row in rows[1:]
    )

    excerpt_ids = {
        key: [read_id(record) for record in records]
        for key, records in selected_records.items()
    }
    if not (
        excerpt_ids["forward"]
        == excerpt_ids["reverse"]
        == excerpt_ids["barcodes"]
    ):
        raise ValueError("excerpt FASTQ read IDs are not synchronized")

    output_files = {
        path.name: {
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
        }
        for path in (*output_paths.values(), metadata_path)
    }
    summary = {
        "dataset": "QIIME 2 Atacama soils paired-end tutorial, 2024.5 1% subset",
        "study": (
            "Neilson JW et al. Significant Impacts of Increasing Aridity on the "
            "Arid Soil Microbiome. mSystems. 2017;2:e00195-16."
        ),
        "source_files": {
            key: {
                **SOURCE_FILES[key],
                "records": counts.get(key),
            }
            for key in SOURCE_FILES
        },
        "selection": {
            "method": "2000 evenly spaced, order-preserving synchronized records",
            "source_records": total_records,
            "selected_records": args.records,
            "first_source_index_1_based": min(selected_indices) + 1,
            "last_source_index_1_based": max(selected_indices) + 1,
        },
        "metadata_samples": metadata_samples,
        "read_id_synchronization": True,
        "phred_offset": 33,
        "read_summaries": {
            key: summarize_reads(records)
            for key, records in selected_records.items()
        },
        "output_files": output_files,
    }
    summary_path = args.output_dir / "source_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
