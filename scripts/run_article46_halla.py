#!/usr/bin/env python3
"""Run the fixed HAllA analysis for Article 46 on Franzosa 2019."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from halla import HAllA
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "small" / "paired-ibd-multiomics"
OUT = ROOT / "results" / "46-integration" / "halla-v1"
SEED = 20260726


def residualize(matrix: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    design = pd.get_dummies(groups, drop_first=True, dtype=float)
    design.insert(0, "Intercept", 1.0)
    x = design.to_numpy(dtype=float)
    y = matrix.to_numpy(dtype=float)
    fitted = x @ np.linalg.lstsq(x, y, rcond=None)[0]
    return pd.DataFrame(y - fitted, index=matrix.index, columns=matrix.columns)


def pairwise_audit(x: pd.DataFrame, y: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for x_name in x.columns:
        for y_name in y.columns:
            rho, pvalue = spearmanr(x[x_name], y[y_name])
            rows.append((x_name, y_name, float(rho), float(pvalue)))
    result = pd.DataFrame(rows, columns=["Microbe", "Metabolite", "Rho", "PValue"])
    result["QValue"] = multipletests(result["PValue"], method="fdr_bh")[1]
    return result.sort_values(["QValue", "PValue", "Microbe", "Metabolite"])


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing result directory: {OUT}")
    OUT.parent.mkdir(parents=True, exist_ok=True)

    otutab = pd.read_csv(DATA / "otutab.tsv", sep="\t", index_col=0)
    taxonomy = pd.read_csv(DATA / "taxonomy.tsv", sep="\t", index_col=0, keep_default_na=False)
    metadata = pd.read_csv(DATA / "metadata.tsv", sep="\t", index_col=0)
    metabolite = pd.read_csv(DATA / "metabolites.tsv", sep="\t", index_col=0)
    metabolite_annotation = pd.read_csv(
        DATA / "metabolite-annotation.tsv", sep="\t", index_col=0, keep_default_na=False
    )

    samples = metadata.index.tolist()
    if otutab.columns.tolist() != samples or metabolite.columns.tolist() != samples:
        raise RuntimeError("Paired sample order differs across views")

    count_samples = otutab.T.astype(float)
    log_counts = np.log(count_samples + 0.5)
    clr = log_counts.sub(log_counts.mean(axis=1), axis=0)
    log_metabolite = np.log1p(metabolite.T.astype(float))

    microbe_features = clr.var(axis=0).sort_values(ascending=False).head(25).index
    metabolite_features = log_metabolite.var(axis=0).sort_values(ascending=False).head(40).index
    clr = clr.loc[:, microbe_features]
    log_metabolite = log_metabolite.loc[:, metabolite_features]

    microbe_labels = {
        feature: f"{taxonomy.loc[feature, 'DisplayName']} [{feature}]" for feature in microbe_features
    }
    metabolite_labels = {
        feature: f"{metabolite_annotation.loc[feature, 'DisplayName']} [{feature}]"
        for feature in metabolite_features
    }

    clr_labeled = clr.rename(columns=microbe_labels)
    metabolite_labeled = log_metabolite.rename(columns=metabolite_labels)
    groups = metadata.loc[samples, "StudyGroup"]
    clr_residual = residualize(clr_labeled, groups)
    metabolite_residual = residualize(metabolite_labeled, groups)

    raw_audit = pairwise_audit(clr_labeled, metabolite_labeled)
    residual_audit = pairwise_audit(clr_residual, metabolite_residual)
    raw_audit.to_csv(OUT.parent / "pairwise-raw.tsv", sep="\t", index=False)
    residual_audit.to_csv(OUT.parent / "pairwise-group-residual.tsv", sep="\t", index=False)

    x_path = OUT.parent / "halla-microbes-group-residual.tsv"
    y_path = OUT.parent / "halla-metabolites-group-residual.tsv"
    clr_residual.T.to_csv(x_path, sep="\t", index_label="Feature")
    metabolite_residual.T.to_csv(y_path, sep="\t", index_label="Feature")

    model = HAllA(
        pdist_metric="spearman",
        linkage_method="average",
        permute_iters=999,
        permute_func="gpd",
        permute_speedup=True,
        fdr_alpha=0.05,
        fdr_method="fdr_bh",
        fnr_thresh=0.20,
        rank_cluster="best",
        out_dir=str(OUT),
        verbose=True,
        no_progress=True,
        dont_copy=True,
        num_threads=4,
        seed=SEED,
    )
    model.load(str(x_path), str(y_path))
    model.run()

    summary = {
        "status": "passed",
        "seed": SEED,
        "samples": len(samples),
        "microbe_features": int(clr_residual.shape[1]),
        "metabolite_features": int(metabolite_residual.shape[1]),
        "tested_pairs": int(clr_residual.shape[1] * metabolite_residual.shape[1]),
        "raw_pairwise_fdr_hits": int((raw_audit["QValue"] < 0.05).sum()),
        "group_residual_pairwise_fdr_hits": int((residual_audit["QValue"] < 0.05).sum()),
        "halla_significant_pairs": int(model.fdr_reject_table.sum()),
        "halla_significant_blocks": int(len(model.significant_blocks)),
        "halla_metric": "Spearman",
        "permutations": 999,
        "fdr_method": "Benjamini-Hochberg",
        "fnr_threshold": 0.20,
        "group_adjustment": "feature-wise residualization on StudyGroup before HAllA",
    }
    (OUT.parent / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
