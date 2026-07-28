#!/usr/bin/env python3
"""Run a seeded SourceTracker2 audit on the official FEAST example."""

from __future__ import annotations

import json
from importlib import metadata as importlib_metadata
from pathlib import Path

import numpy as np
import pandas as pd
from skbio.stats import subsample_counts

# SourceTracker2 2.0.1 still references the alias removed in NumPy 2.x.
if not hasattr(np, "int"):
    np.int = int  # type: ignore[attr-defined]

from sourcetracker._sourcetracker import _gibbs


SEED = 20260726
COVERAGE = 5000
ROOT = Path.cwd().resolve()
DATA_DIR = ROOT / "data" / "small" / "source-tracking-feast"
OUT_DIR = ROOT / "results" / "50-source-tracking"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_OUTPUTS = {
    "feast-loo-predictions.tsv",
    "feast-main-class.tsv",
    "feast-main-sample.tsv",
    "feast-missing-source.tsv",
    "feast-native",
    "source-tracker-loo-predictions.tsv",
    "source-tracker-main.tsv",
    "summary-feast.json",
    "summary-source-tracker.json",
}
unexpected = {path.name for path in OUT_DIR.iterdir()} - EXPECTED_OUTPUTS
if unexpected:
    raise RuntimeError(
        "Refusing to overwrite unexpected files in result directory: "
        + ", ".join(sorted(unexpected))
    )

feature_by_sample = pd.read_csv(DATA_DIR / "otutab.tsv", sep="\t", index_col=0)
sample_by_feature = feature_by_sample.T.astype(np.int64)
metadata = pd.read_csv(DATA_DIR / "metadata.tsv", sep="\t", index_col=0)
if list(sample_by_feature.index) != list(metadata.index):
    raise ValueError("Count-table samples and metadata rows are not aligned")
if (sample_by_feature.to_numpy() < 0).any():
    raise ValueError("Counts must be non-negative")

source_ids = metadata.index[metadata["SourceSink"] == "Source"].tolist()
sink_ids = metadata.index[metadata["SourceSink"] == "Sink"].tolist()
if len(sink_ids) != 1:
    raise ValueError("The official example should contain exactly one sink")
source_classes = list(dict.fromkeys(metadata.loc[source_ids, "SourceClass"]))


def rarefy_rows(frame: pd.DataFrame, depth: int) -> pd.DataFrame:
    rarefied = np.vstack(
        [
            subsample_counts(row.to_numpy(dtype=np.int64), depth, replace=False)
            for _, row in frame.iterrows()
        ]
    )
    return pd.DataFrame(rarefied, index=frame.index, columns=frame.columns)


def run_sourcetracker(
    training_ids: list[str], sink_id: str, seed: int
) -> tuple[pd.Series, pd.Series]:
    np.random.seed(seed)
    grouped_sources = sample_by_feature.loc[training_ids].groupby(
        metadata.loc[training_ids, "SourceClass"], sort=True
    ).sum()
    if (grouped_sources.sum(axis=1) < COVERAGE).any():
        raise ValueError("At least one source class is shallower than COVERAGE")
    sink_frame = sample_by_feature.loc[[sink_id]]
    if int(sink_frame.sum(axis=1).iloc[0]) < COVERAGE:
        raise ValueError("Sink is shallower than COVERAGE")
    source_rarefied = rarefy_rows(grouped_sources, COVERAGE)
    sink_rarefied = rarefy_rows(sink_frame, COVERAGE)
    mean, standard_deviation = _gibbs(
        source_rarefied,
        sink_rarefied,
        alpha1=0.001,
        alpha2=0.001,
        beta=10,
        restarts=10,
        draws_per_restart=2,
        burnin=100,
        delay=10,
    )
    mean_row = mean.loc[sink_id]
    sd_row = standard_deviation.loc[sink_id]
    if not np.isclose(mean_row.sum(), 1.0, atol=1e-6):
        raise ValueError("SourceTracker2 contributions do not sum to one")
    return mean_row, sd_row


main_mean, main_sd = run_sourcetracker(source_ids, sink_ids[0], SEED)
main_table = pd.DataFrame(
    {
        "SinkSampleID": sink_ids[0],
        "SourceClass": main_mean.index,
        "Contribution": main_mean.to_numpy(),
        "MonteCarloSD": main_sd.reindex(main_mean.index).to_numpy(),
    }
)
main_table.to_csv(OUT_DIR / "source-tracker-main.tsv", sep="\t", index=False)

class_sizes = metadata.loc[source_ids, "SourceClass"].value_counts()
eligible_loo_ids = [
    sample_id
    for sample_id in source_ids
    if class_sizes[metadata.loc[sample_id, "SourceClass"]] >= 2
]
loo_rows: list[dict[str, object]] = []
for index, held_out in enumerate(eligible_loo_ids, start=1):
    training_ids = [sample_id for sample_id in source_ids if sample_id != held_out]
    mean, sd = run_sourcetracker(training_ids, held_out, SEED + index)
    aligned_mean = mean.reindex(source_classes + ["Unknown"], fill_value=0.0)
    aligned_sd = sd.reindex(source_classes + ["Unknown"], fill_value=0.0)
    known = aligned_mean[source_classes]
    true_class = metadata.loc[held_out, "SourceClass"]
    predicted_known = str(known.idxmax())
    predicted_with_unknown = str(aligned_mean.idxmax())
    row: dict[str, object] = {
        "SampleID": held_out,
        "TrueClass": true_class,
        "PredictedKnownClass": predicted_known,
        "PredictedIncludingUnknown": predicted_with_unknown,
        "CorrectKnownClass": predicted_known == true_class,
    }
    for class_name in source_classes + ["Unknown"]:
        safe_name = class_name.replace(" ", "")
        row[safe_name] = float(aligned_mean[class_name])
        row[f"{safe_name}SD"] = float(aligned_sd[class_name])
    loo_rows.append(row)

loo_table = pd.DataFrame(loo_rows)
loo_table.to_csv(
    OUT_DIR / "source-tracker-loo-predictions.tsv", sep="\t", index=False
)

summary = {
    "status": "passed",
    "seed": SEED,
    "samples": int(metadata.shape[0]),
    "sink_samples": len(sink_ids),
    "source_samples": len(source_ids),
    "source_classes": source_classes,
    "taxa": int(feature_by_sample.shape[0]),
    "coverage": COVERAGE,
    "alpha1": 0.001,
    "alpha2": 0.001,
    "beta": 10,
    "restarts": 10,
    "draws_per_restart": 2,
    "burnin": 100,
    "delay": 10,
    "main_largest_source_class": str(main_mean.drop("Unknown").idxmax()),
    "main_unknown": float(main_mean["Unknown"]),
    "leave_one_out_samples": int(loo_table.shape[0]),
    "leave_one_out_known_class_accuracy": float(
        loo_table["CorrectKnownClass"].mean()
    ),
    "versions": {
        "python": ".".join(map(str, __import__("sys").version_info[:3])),
        "sourcetracker": importlib_metadata.version("sourcetracker"),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    },
}
(OUT_DIR / "summary-source-tracker.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
    encoding="utf-8",
)
print(json.dumps(summary, ensure_ascii=False, indent=2))
