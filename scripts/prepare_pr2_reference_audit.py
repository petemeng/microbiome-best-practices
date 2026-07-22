#!/usr/bin/env python3
"""Create a compact, deterministic audit fixture from PR2 v5.1.1 SSU."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import statistics
from collections import Counter
from pathlib import Path


VERSION = "5.1.1"
RELEASE_TAG = "v5.1.1"
RELEASE_DATE = "2025-10-27"
RELEASE_COMMIT = "acd4e3d887d9cb8132b699f886cf7448c947b88e"
RELEASE_DOI = "10.5281/zenodo.17458343"
SOURCE_URL = (
    "https://github.com/pr2database/pr2database/releases/download/v5.1.1/"
    "pr2_version_5.1.1_SSU_dada2.fasta.gz"
)
SOURCE_SHA256 = "0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d"
SOURCE_BYTES = 54_833_479
EXPECTED_RECORDS = 240_201
SNAPSHOT_STRIDE = 240


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def fasta_records(path: Path):
    header: str | None = None
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError("sequence observed before first FASTA header")
                chunks.append(line.upper())
    if header is not None:
        yield header, "".join(chunks)


def deterministic_gzip(payload: bytes) -> bytes:
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
    observed_bytes = args.input.stat().st_size
    observed_sha256 = hashlib.sha256(args.input.read_bytes()).hexdigest()
    if observed_sha256 != SOURCE_SHA256:
        raise ValueError(f"PR2 SHA-256 {observed_sha256} != {SOURCE_SHA256}")
    if observed_bytes != SOURCE_BYTES:
        raise ValueError(f"PR2 size {observed_bytes} != {SOURCE_BYTES}")

    lengths: list[int] = []
    depth_counts: Counter[int] = Counter()
    domain_counts: Counter[str] = Counter()
    supergroup_counts: Counter[str] = Counter()
    fourth_rank_counts: Counter[str] = Counter()
    ambiguous_records = 0
    empty_rank_records = 0
    taxonomy_strings: set[str] = set()
    snapshot_records: list[tuple[int, str, str]] = []
    last_record: tuple[str, str] | None = None

    for index, (header, sequence) in enumerate(fasta_records(args.input), start=1):
        last_record = (header, sequence)
        ranks = header.rstrip(";").split(";")
        lengths.append(len(sequence))
        depth_counts[len(ranks)] += 1
        domain_counts[ranks[0] if ranks else "MISSING"] += 1
        supergroup_counts[ranks[1] if len(ranks) > 1 else "MISSING"] += 1
        fourth_rank_counts[ranks[3] if len(ranks) > 3 else "MISSING"] += 1
        ambiguous_records += int(any(base not in "ACGT" for base in sequence))
        empty_rank_records += int(any(not rank for rank in ranks))
        taxonomy_strings.add(header)
        if index == 1 or (index - 1) % SNAPSHOT_STRIDE == 0:
            snapshot_records.append((index, header, sequence))

    if len(lengths) != EXPECTED_RECORDS:
        raise ValueError(f"PR2 has {len(lengths)} records; expected {EXPECTED_RECORDS}")
    if snapshot_records[-1][0] != EXPECTED_RECORDS:
        if last_record is None:
            raise ValueError("PR2 reference contains no FASTA records")
        snapshot_records.append((EXPECTED_RECORDS, *last_record))

    snapshot_text = "".join(
        f">source_record_{index}|{header}\n{sequence}\n"
        for index, header, sequence in snapshot_records
    ).encode("ascii")
    snapshot_gzip = deterministic_gzip(snapshot_text)
    snapshot_path = args.output_dir / "pr2-v5.1.1-stratified-snapshot.fasta.gz"
    snapshot_path.write_bytes(snapshot_gzip)

    summary = {
        "database": "PR2 SSU DADA2 reference",
        "version": VERSION,
        "release_tag": RELEASE_TAG,
        "release_date": RELEASE_DATE,
        "release_commit": RELEASE_COMMIT,
        "release_doi": RELEASE_DOI,
        "license": "MIT (release metadata)",
        "source_url": SOURCE_URL,
        "source_filename": args.input.name,
        "source_bytes": observed_bytes,
        "source_sha256": observed_sha256,
        "record_count": len(lengths),
        "unique_taxonomy_strings": len(taxonomy_strings),
        "taxonomy_depth_counts": dict(sorted(depth_counts.items())),
        "records_with_empty_rank": empty_rank_records,
        "domain_counts": dict(domain_counts.most_common()),
        "supergroup_counts": dict(supergroup_counts.most_common()),
        "fourth_rank_counts_top20": dict(fourth_rank_counts.most_common(20)),
        "fungi_records": fourth_rank_counts["Fungi"],
        "records_with_ambiguous_bases": ambiguous_records,
        "sequence_length": {
            "minimum": min(lengths),
            "q1": statistics.quantiles(lengths, n=4, method="inclusive")[0],
            "median": statistics.median(lengths),
            "mean": round(statistics.fmean(lengths), 3),
            "q3": statistics.quantiles(lengths, n=4, method="inclusive")[2],
            "maximum": max(lengths),
        },
        "snapshot": {
            "selection": (
                f"record 1, every {SNAPSHOT_STRIDE}th source record thereafter, "
                "and the final record"
            ),
            "record_count": len(snapshot_records),
            "uncompressed_sha256": hashlib.sha256(snapshot_text).hexdigest(),
            "gzip_sha256": hashlib.sha256(snapshot_gzip).hexdigest(),
            "gzip_bytes": len(snapshot_gzip),
        },
        "version_display_note": (
            "The PR2 landing page rendered '5.1.12' on the audit date, while "
            "the Git tag, release title, Zenodo metadata, DOI version, commit, "
            "and repository README all identify v5.1.1."
        ),
    }
    (args.output_dir / "pr2-v5.1.1-source-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"Audited {len(lengths):,} PR2 records and wrote "
        f"{len(snapshot_records):,} deterministic snapshot records."
    )


if __name__ == "__main__":
    main()
