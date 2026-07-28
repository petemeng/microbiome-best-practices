#!/usr/bin/env python3
"""Prepare the frozen MiBioGen -> IBD two-sample MR teaching bundle.

The script reads the official MiBioGen top-hit table already cached in
``data/raw`` and requests six association records from the official NHGRI-EBI
GWAS Catalog summary-statistics API. Raw API responses are cached so later runs
are offline and auditable.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MIBIOGEN = ROOT / "data" / "raw" / "MBG.allHits.p1e4.txt"
RAW_API = ROOT / "data" / "raw" / "mr-mibiogen"
OUT = ROOT / "data" / "small" / "mr-mibiogen"
EXPOSURE_ID = "genus.Bifidobacterium.id.436"
OUTCOME_ID = "GCST004131"
EXPECTED_MIBIOGEN_SHA256 = (
    "37001a83d060596fe0b97b63d6a397f01f43a29add2925d406916b7a50b5883e"
)

# One outcome-available lead per chromosome among MiBioGen P < 5e-6 hits.
# This deliberately does not claim to replace ancestry-matched LD clumping.
SELECTED_RSIDS = (
    "rs182549",
    "rs73797465",
    "rs857444",
    "rs10841473",
    "rs7322849",
    "rs75344046",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def api_url(rsid: str) -> str:
    return (
        "https://www.ebi.ac.uk/gwas/summary-statistics/api/associations/"
        f"{rsid}?study_accession={OUTCOME_ID}&reveal=all"
    )


def request_or_cache(rsid: str) -> tuple[dict, Path]:
    RAW_API.mkdir(parents=True, exist_ok=True)
    target = RAW_API / f"{OUTCOME_ID}-{rsid}.json"
    if not target.exists():
        request = urllib.request.Request(
            api_url(rsid),
            headers={"User-Agent": "microbiome-best-practices/1.0"},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        parsed = json.loads(payload)
        target.write_text(
            json.dumps(parsed, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        time.sleep(0.25)
    return json.loads(target.read_text(encoding="utf-8")), target


def association_rows(payload: dict) -> list[dict]:
    embedded = payload.get("_embedded", {}).get("associations", {})
    if isinstance(embedded, dict):
        return list(embedded.values())
    if isinstance(embedded, list):
        return embedded
    raise RuntimeError("Unexpected GWAS Catalog association payload")


def choose_harmonised_record(rows: list[dict], rsid: str) -> dict:
    rows = [row for row in rows if row.get("study_accession") == OUTCOME_ID]
    if not rows:
        raise RuntimeError(f"No {OUTCOME_ID} association returned for {rsid}")

    # Code 10 = forward strand, alleles correct. Prefer it over a duplicated
    # code-11 (flipped-allele) row. Other valid harmonised records are retained
    # only as a fallback and the full payload remains in data/raw for audit.
    priority = {10: 0, 5: 1, 12: 2, 11: 3, 6: 4, 13: 5}
    usable = [
        row
        for row in rows
        if row.get("hm_effect_allele")
        and row.get("hm_other_allele")
        and row.get("hm_beta") is not None
        and row.get("standard_error") is not None
    ]
    if not usable:
        raise RuntimeError(f"No harmonised effect record returned for {rsid}")
    return sorted(usable, key=lambda row: priority.get(int(row.get("hm_code", 99)), 99))[0]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    observed_sha256 = sha256(MIBIOGEN)
    if observed_sha256 != EXPECTED_MIBIOGEN_SHA256:
        raise RuntimeError(f"Unexpected MiBioGen checksum: {observed_sha256}")

    mbg = pd.read_csv(MIBIOGEN, sep="\t")
    exposure = mbg[
        (mbg["bac"] == EXPOSURE_ID) & mbg["rsID"].isin(SELECTED_RSIDS)
    ].copy()
    exposure = exposure.set_index("rsID").loc[list(SELECTED_RSIDS)].reset_index()
    if exposure.shape[0] != len(SELECTED_RSIDS):
        raise RuntimeError("Not all pre-specified MiBioGen instruments were found")
    if not (exposure["P.weightedSumZ"] < 5e-6).all():
        raise RuntimeError("Instrument outside the pre-specified P < 5e-6 threshold")
    if exposure["chr"].duplicated().any():
        raise RuntimeError("Teaching subset must contain one SNP per chromosome")

    exposure_out = pd.DataFrame(
        {
            "SNP": exposure["rsID"],
            "exposure": "Bifidobacterium genus abundance",
            "id.exposure": EXPOSURE_ID,
            "chr.exposure": exposure["chr"].astype(int),
            "pos.exposure": exposure["bp"].astype(int),
            "effect_allele.exposure": exposure["eff.allele"].str.upper(),
            "other_allele.exposure": exposure["ref.allele"].str.upper(),
            "beta.exposure": exposure["beta"].astype(float),
            "se.exposure": exposure["SE"].astype(float),
            "pval.exposure": exposure["P.weightedSumZ"].astype(float),
            "eaf.exposure": float("nan"),
            "samplesize.exposure": exposure["N"].astype(int),
            "ncohorts.exposure": exposure["Ncohorts"].astype(int),
        }
    )
    exposure_out["FStatistic"] = (
        exposure_out["beta.exposure"] / exposure_out["se.exposure"]
    ) ** 2

    outcome_rows: list[dict] = []
    selection_rows: list[dict] = []
    response_hashes: dict[str, str] = {}
    for rsid in SELECTED_RSIDS:
        payload, cache_path = request_or_cache(rsid)
        rows = association_rows(payload)
        chosen = choose_harmonised_record(rows, rsid)
        response_hashes[cache_path.name] = sha256(cache_path)
        outcome_rows.append(
            {
                "SNP": rsid,
                "outcome": "Inflammatory bowel disease",
                "id.outcome": OUTCOME_ID,
                "chr.outcome": int(chosen["chromosome"]),
                "pos.outcome": int(chosen["base_pair_location"]),
                "effect_allele.outcome": chosen["hm_effect_allele"].upper(),
                "other_allele.outcome": chosen["hm_other_allele"].upper(),
                "beta.outcome": float(chosen["hm_beta"]),
                "se.outcome": float(chosen["standard_error"]),
                "pval.outcome": float(chosen["p_value"]),
                "samplesize.outcome": 59957,
                "ncase.outcome": 25042,
                "ncontrol.outcome": 34915,
                "eaf.outcome": "",
            }
        )
        selection_rows.append(
            {
                "SNP": rsid,
                "returned_records": len(rows),
                "selected_harmonisation_code": int(chosen["hm_code"]),
                "selection_rule": "prefer harmonised alleles-correct records (code 10/5) over flipped duplicates",
                "api_url": api_url(rsid),
                "cache_file": str(cache_path.relative_to(ROOT)),
                "cache_sha256": response_hashes[cache_path.name],
            }
        )

    outcome_out = pd.DataFrame(outcome_rows)
    selection = pd.DataFrame(selection_rows)
    exposure_out.to_csv(OUT / "exposure.tsv", sep="\t", index=False)
    outcome_out.to_csv(OUT / "outcome.tsv", sep="\t", index=False)
    selection.to_csv(OUT / "selection-ledger.tsv", sep="\t", index=False)

    summary = {
        "exposure": "MiBioGen genus.Bifidobacterium.id.436",
        "exposure_citation": "Kurilshikov et al. Nature Genetics 2021; doi:10.1038/s41588-020-00763-1",
        "exposure_source_url": "https://molgenis26.gcc.rug.nl/downloads/MiBioGen/MBG.allHits.p1e4.txt",
        "exposure_source_sha256": observed_sha256,
        "outcome": "Inflammatory bowel disease, GCST004131",
        "outcome_citation": "de Lange et al. Nature Genetics 2017; doi:10.1038/ng.3760; PMID:28067908",
        "outcome_api": "NHGRI-EBI GWAS Catalog summary-statistics API",
        "outcome_api_url_template": "https://www.ebi.ac.uk/gwas/summary-statistics/api/associations/{rsid}?study_accession=GCST004131&reveal=all",
        "outcome_sample_size": {"cases": 25042, "controls": 34915, "total": 59957},
        "instrument_rule": "MiBioGen P < 5e-6; one outcome-available SNP per chromosome",
        "instrument_count": len(SELECTED_RSIDS),
        "instrument_note": "The frozen subset is cross-chromosome and therefore pairwise independent by chromosome, but it is not a substitute for ancestry-matched population LD clumping.",
        "effect_allele_frequency_available": False,
        "steiger_available": False,
        "colocalisation_available": False,
        "api_response_sha256": response_hashes,
        "contract": "TwoSampleMR-format exposure.tsv and outcome.tsv plus an API selection ledger",
    }
    (OUT / "source-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"Prepared {len(SELECTED_RSIDS)} MiBioGen -> IBD instruments; "
        f"minimum F={exposure_out['FStatistic'].min():.2f}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
