#!/usr/bin/env python3
"""Run a checksum-auditable gemelli Robust Aitchison RPCA."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import biom
import gemelli
import numpy as np
import pandas as pd
import scipy
import skbio
import sklearn
import gemelli.optspace as gemelli_optspace
from gemelli.rpca import rpca, rpca_table_processing
from scipy.sparse.linalg import svds as scipy_svds


SEED = 20260722
ROUND_DIGITS = 12


def deterministic_svds(*args, **kwargs):
    """Give gemelli's ARPACK initialization an explicit reproducible RNG."""
    kwargs.setdefault("rng", np.random.default_rng(SEED))
    return scipy_svds(*args, **kwargs)


def canonicalize_results(ordination, distance):
    """Fix arbitrary PC signs and remove sub-machine-precision write drift."""
    sample_scores = ordination.samples.copy().astype(float)
    feature_loadings = ordination.features.copy().astype(float)
    for axis in sample_scores.columns:
        anchor = sample_scores[axis].abs().idxmax()
        if sample_scores.loc[anchor, axis] < 0:
            sample_scores.loc[:, axis] *= -1
            feature_loadings.loc[:, axis] *= -1

    sample_scores = sample_scores.round(ROUND_DIGITS)
    feature_loadings = feature_loadings.round(ROUND_DIGITS)
    eigvals = ordination.eigvals.astype(float).round(ROUND_DIGITS)
    proportion = ordination.proportion_explained.astype(float).round(
        ROUND_DIGITS
    )
    canonical_ordination = skbio.OrdinationResults(
        ordination.short_method_name,
        ordination.long_method_name,
        eigvals,
        samples=sample_scores,
        features=feature_loadings,
        proportion_explained=proportion,
    )
    canonical_distance = skbio.DistanceMatrix(
        np.round(np.asarray(distance.data, dtype=float), ROUND_DIGITS),
        ids=distance.ids,
    )
    return canonical_ordination, canonical_distance


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-components", type=int, default=3)
    parser.add_argument("--min-sample-count", type=int, default=500)
    parser.add_argument("--min-feature-count", type=int, default=10)
    parser.add_argument("--min-feature-frequency", type=float, default=10.0)
    parser.add_argument("--max-iterations", type=int, default=5)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    counts_df = pd.read_csv(input_path, sep="\t", index_col=0)
    if counts_df.index.has_duplicates or counts_df.columns.has_duplicates:
        raise ValueError("Feature and sample identifiers must be unique")
    counts = counts_df.to_numpy(dtype=float)
    if not np.isfinite(counts).all():
        raise ValueError("Count table contains non-finite values")
    if (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
        raise ValueError("Count table must contain non-negative integers")
    if (counts.sum(axis=0) == 0).any() or (counts.sum(axis=1) == 0).any():
        raise ValueError("Count table contains an empty sample or feature")

    table = biom.Table(
        counts,
        observation_ids=counts_df.index.astype(str).tolist(),
        sample_ids=counts_df.columns.astype(str).tolist(),
        table_id="GlobalPatterns feature table",
    )
    processed = rpca_table_processing(
        table,
        min_sample_count=args.min_sample_count,
        min_feature_count=args.min_feature_count,
        min_feature_frequency=args.min_feature_frequency,
    )

    np.random.seed(SEED)
    gemelli_optspace.svds = deterministic_svds
    ordination, distance = rpca(
        table,
        n_components=args.n_components,
        min_sample_count=args.min_sample_count,
        min_feature_count=args.min_feature_count,
        min_feature_frequency=args.min_feature_frequency,
        max_iterations=args.max_iterations,
    )
    ordination, distance = canonicalize_results(ordination, distance)

    distance_path = output_dir / "robust-aitchison-distance.tsv"
    sample_path = output_dir / "robust-aitchison-sample-scores.tsv"
    feature_path = output_dir / "robust-aitchison-feature-loadings.tsv"
    eig_path = output_dir / "robust-aitchison-eigenvalues.tsv"
    ordination_path = output_dir / "robust-aitchison-ordination.txt"
    summary_path = output_dir / "gemelli-summary.json"

    distance.write(distance_path)
    ordination.write(ordination_path)
    ordination.samples.rename_axis("SampleID").reset_index().to_csv(
        sample_path, sep="\t", index=False, float_format="%.12g"
    )
    ordination.features.rename_axis("FeatureID").reset_index().to_csv(
        feature_path, sep="\t", index=False, float_format="%.12g"
    )
    pd.DataFrame(
        {
            "Axis": ordination.eigvals.index.astype(str),
            "Eigenvalue": ordination.eigvals.to_numpy(),
            "ProportionExplained": ordination.proportion_explained.to_numpy(),
        }
    ).to_csv(eig_path, sep="\t", index=False, float_format="%.12g")

    distance_array = np.asarray(distance.data, dtype=float)
    off_diagonal = distance_array[np.triu_indices_from(distance_array, k=1)]
    summary = {
        "seed": SEED,
        "input": {
            "path": str(input_path),
            "sha256": sha256(input_path),
            "features": int(table.shape[0]),
            "samples": int(table.shape[1]),
            "reads": int(counts.sum()),
        },
        "filter": {
            "min_sample_count": args.min_sample_count,
            "min_feature_count": args.min_feature_count,
            "min_feature_frequency_percent": args.min_feature_frequency,
            "features_retained": int(processed.shape[0]),
            "samples_retained": int(processed.shape[1]),
            "reads_retained": int(processed.sum()),
        },
        "model": {
            "n_components": args.n_components,
            "max_iterations": args.max_iterations,
            "arpack_seed": SEED,
            "canonical_round_digits": ROUND_DIGITS,
            "distance_min_off_diagonal": float(off_diagonal.min()),
            "distance_median_off_diagonal": float(np.median(off_diagonal)),
            "distance_max_off_diagonal": float(off_diagonal.max()),
        },
        "versions": {
            "python": platform.python_version(),
            "gemelli": gemelli.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "scikit_bio": skbio.__version__,
            "biom_format": biom.__version__,
        },
        "runtime": {
            "python_executable": sys.executable,
            "python_prefix": sys.prefix,
            "python_hash_seed": os.environ.get("PYTHONHASHSEED", "unset"),
        },
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in (
                distance_path,
                sample_path,
                feature_path,
                eig_path,
                ordination_path,
            )
        },
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "gemelli RPCA complete: "
        f"{processed.shape[0]} features x {processed.shape[1]} samples"
    )


if __name__ == "__main__":
    main()
