#!/usr/bin/env python3
"""Prepare the pinned ITSxpress soil tutorial FASTQs for Article 13.

The upstream repository names its four plain-text FASTQ files ``*.fastq.gz``.
This script verifies the bytes at the pinned commit and writes deterministic,
actually gzip-compressed copies for QIIME 2 import.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


SOURCE_REPOSITORY = "https://github.com/USDA-ARS-GBRU/itsxpress-tutorial"
SOURCE_COMMIT = "916145d3b05e20656fde30c4732dc084de72d876"
SOURCE_DATE = "2024-04-15"
SOURCE_LICENSE = "No separate repository or data license stated at the pinned commit"

FILES = {
    "sample1_r1.fastq.gz": {
        "sample_id": "sample1",
        "direction": "forward",
        "sha256": "5e44d1794ff70790fc2ecf807390e459c60db8cb2d807f5a5cba8532d42ef6c4",
        "bytes": 1_498_149,
        "records": 2_500,
    },
    "sample1_r2.fastq.gz": {
        "sample_id": "sample1",
        "direction": "reverse",
        "sha256": "bffb461483b81e74beeeb65fed3ac3e5321535d1023ed3334ade9a2dd99e3ae7",
        "bytes": 1_533_167,
        "records": 2_500,
    },
    "sample2_r1.fastq.gz": {
        "sample_id": "sample2",
        "direction": "forward",
        "sha256": "00d49fbba7d69f6fad91e128581df2e144faefdd466b844bf0734f67e76a031c",
        "bytes": 1_529_971,
        "records": 2_500,
    },
    "sample2_r2.fastq.gz": {
        "sample_id": "sample2",
        "direction": "reverse",
        "sha256": "c31efec6c90f024122f5f282126ee14b32b4980daab517494cb1686ddba81aca",
        "bytes": 1_564_795,
        "records": 2_500,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def normalized_read_id(header: bytes) -> bytes:
    token = header.split(maxsplit=1)[0]
    if token.endswith(b"/1") or token.endswith(b"/2"):
        token = token[:-2]
    return token


def validate_fastq(payload: bytes, expected_records: int) -> list[bytes]:
    lines = payload.splitlines()
    if len(lines) != expected_records * 4:
        raise ValueError(
            f"FASTQ contains {len(lines) // 4} records; expected {expected_records}"
        )
    read_ids: list[bytes] = []
    for offset in range(0, len(lines), 4):
        header, sequence, plus, quality = lines[offset : offset + 4]
        if not header.startswith(b"@"):
            raise ValueError(f"record {offset // 4 + 1} has an invalid header")
        if not plus.startswith(b"+"):
            raise ValueError(f"record {offset // 4 + 1} has an invalid plus line")
        if len(sequence) != len(quality):
            raise ValueError(
                f"record {offset // 4 + 1} has unequal sequence/quality lengths"
            )
        read_ids.append(normalized_read_id(header))
    return read_ids


def deterministic_gzip(payload: bytes) -> bytes:
    import io

    buffer = io.BytesIO()
    with gzip.GzipFile(
        filename="",
        mode="wb",
        fileobj=buffer,
        compresslevel=9,
        mtime=0,
    ) as handle:
        handle.write(payload)
    return buffer.getvalue()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source_records: dict[str, list[bytes]] = {}
    file_audit: list[dict[str, object]] = []
    for filename, expected in FILES.items():
        source_path = args.source_dir / filename
        payload = source_path.read_bytes()
        if payload.startswith(b"\x1f\x8b"):
            raise ValueError(
                f"{source_path} is gzip-compressed, but the pinned source audit "
                "expects the repository's plain-text bytes"
            )
        observed_sha256 = sha256_bytes(payload)
        if observed_sha256 != expected["sha256"]:
            raise ValueError(
                f"{filename}: SHA-256 {observed_sha256} != {expected['sha256']}"
            )
        if len(payload) != expected["bytes"]:
            raise ValueError(
                f"{filename}: {len(payload)} bytes != {expected['bytes']}"
            )
        read_ids = validate_fastq(payload, int(expected["records"]))
        source_records[filename] = read_ids

        compressed = deterministic_gzip(payload)
        output_path = args.output_dir / filename
        output_path.write_bytes(compressed)
        with gzip.open(output_path, "rb") as handle:
            if handle.read() != payload:
                raise ValueError(f"{filename}: gzip round-trip failed")

        file_audit.append(
            {
                "file": filename,
                "sample_id": expected["sample_id"],
                "direction": expected["direction"],
                "records": expected["records"],
                "source_bytes": len(payload),
                "source_sha256": observed_sha256,
                "source_is_gzip": False,
                "prepared_bytes": len(compressed),
                "prepared_sha256": sha256_bytes(compressed),
                "prepared_is_gzip": True,
            }
        )

    for sample in ("sample1", "sample2"):
        forward = source_records[f"{sample}_r1.fastq.gz"]
        reverse = source_records[f"{sample}_r2.fastq.gz"]
        if forward != reverse:
            mismatches = sum(a != b for a, b in zip(forward, reverse))
            raise ValueError(f"{sample}: {mismatches} paired read IDs are mismatched")

    metadata = (
        "sample-id\tCollection\n"
        "sample1\tAfghanistan collection\n"
        "sample2\tTurkey collection\n"
    )
    (args.output_dir / "metadata.tsv").write_text(metadata, encoding="utf-8")

    prepared_digest = hashlib.sha256()
    for row in sorted(file_audit, key=lambda item: str(item["file"])):
        prepared_digest.update(str(row["prepared_sha256"]).encode("ascii"))

    summary = {
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_commit_date": SOURCE_DATE,
        "source_license": SOURCE_LICENSE,
        "source_documentation_claim": "two samples subsampled to 10,000 read pairs",
        "observed_sample_count": 2,
        "observed_pairs_per_sample": 2_500,
        "observed_total_pairs": 5_000,
        "pair_id_mismatches": 0,
        "preparation": (
            "Verified the pinned plain-text source bytes and wrote deterministic "
            "gzip streams (mtime=0, no embedded filename)."
        ),
        "prepared_collection_sha256": prepared_digest.hexdigest(),
        "files": file_audit,
    }
    (args.output_dir / "source_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        f"Prepared {len(file_audit)} FASTQs, 2 samples, "
        "5,000 synchronized read pairs."
    )


if __name__ == "__main__":
    main()
