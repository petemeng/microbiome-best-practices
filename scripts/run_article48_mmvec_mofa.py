#!/usr/bin/env python3
"""Run native Spearman/CCA/MOFA analyses after native MMvec for Article 48."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from mofapy2 import __version__ as mofapy2_version
from mofapy2.run.entry_point import entry_point
from scipy.stats import kruskal, spearmanr
from sklearn import __version__ as sklearn_version
from sklearn.cross_decomposition import CCA
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "small" / "paired-ibd-multiomics"
OUT = ROOT / "results" / "48-mmvec-mofa"
SEED = 20260726


def zscore(frame: pd.DataFrame) -> pd.DataFrame:
    spread = frame.std(axis=0, ddof=1).replace(0, np.nan)
    return frame.sub(frame.mean(axis=0), axis=1).div(spread, axis=1).fillna(0.0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    native_summary_path = OUT / "summary-mmvec-native.json"
    top_pairs_path = OUT / "mmvec-top-pairs.tsv"
    conditional_path = OUT / "mmvec-conditional-probabilities.tsv"
    if not native_summary_path.is_file() or not top_pairs_path.is_file():
        raise RuntimeError(
            "Run scripts/run_article48_native_mmvec.py in mmvec-native first"
        )

    otutab = pd.read_csv(DATA / "otutab.tsv", sep="\t", index_col=0)
    taxonomy = pd.read_csv(DATA / "taxonomy.tsv", sep="\t", index_col=0, keep_default_na=False)
    metadata = pd.read_csv(DATA / "metadata.tsv", sep="\t", index_col=0)
    metabolites = pd.read_csv(DATA / "metabolites.tsv", sep="\t", index_col=0)
    metabolite_annotation = pd.read_csv(
        DATA / "metabolite-annotation.tsv", sep="\t", index_col=0,
        keep_default_na=False,
    )
    samples = metadata.index.tolist()
    if otutab.columns.tolist() != samples or metabolites.columns.tolist() != samples:
        raise RuntimeError("Paired sample order differs across views")

    microbe_features = otutab.sum(axis=1).sort_values(ascending=False).head(100).index.tolist()
    metabolite_log = np.log1p(metabolites.T.astype(float))
    metabolite_features = (
        metabolite_log.var(axis=0, ddof=1).sort_values(ascending=False).head(120).index.tolist()
    )
    microbe_samples = otutab.loc[microbe_features, samples].T.astype(float)
    metabolite_samples = metabolites.loc[metabolite_features, samples].T.astype(float)

    log_microbe = np.log(microbe_samples + 0.5)
    clr_microbe = log_microbe.sub(log_microbe.mean(axis=1), axis=0)
    log_metabolite = np.log1p(metabolite_samples)

    association_rows: list[dict[str, float | str]] = []
    for microbe in microbe_features:
        for metabolite in metabolite_features:
            rho, pvalue = spearmanr(clr_microbe[microbe], log_metabolite[metabolite])
            association_rows.append({
                "microbe": microbe,
                "metabolite": metabolite,
                "Rho": float(rho),
                "PValue": float(pvalue),
            })
    spearman = pd.DataFrame(association_rows)
    spearman["QValue"] = multipletests(spearman["PValue"], method="fdr_bh")[1]
    spearman["AbsRho"] = spearman["Rho"].abs()
    spearman = spearman.sort_values(
        ["QValue", "AbsRho", "microbe", "metabolite"],
        ascending=[True, False, True, True],
    )
    spearman.to_csv(OUT / "spearman-all-pairs.tsv", sep="\t", index=False)

    x_cca = zscore(clr_microbe)
    y_cca = zscore(log_metabolite)
    cca_model = CCA(n_components=3, scale=False, max_iter=1000, tol=1e-6)
    x_scores, y_scores = cca_model.fit_transform(x_cca, y_cca)
    canonical_correlations = [
        float(np.corrcoef(x_scores[:, index], y_scores[:, index])[0, 1])
        for index in range(3)
    ]
    cca_scores = pd.DataFrame({
        "SampleID": samples,
        "StudyGroup": metadata.loc[samples, "StudyGroup"].to_numpy(),
        **{f"MicrobiomeCC{i + 1}": x_scores[:, i] for i in range(3)},
        **{f"MetabolomeCC{i + 1}": y_scores[:, i] for i in range(3)},
    }).set_index("SampleID")
    cca_scores.to_csv(OUT / "cca-scores.tsv", sep="\t", index=True)

    scaled_microbe = zscore(clr_microbe)
    scaled_metabolite = zscore(log_metabolite)
    mofa = entry_point()
    mofa.set_data_options(
        scale_views=False, scale_groups=False, center_groups=True,
        use_float32=False,
    )
    mofa.set_data_matrix(
        [[scaled_microbe.to_numpy()], [scaled_metabolite.to_numpy()]],
        likelihoods=["gaussian", "gaussian"],
        views_names=["Microbiome", "Metabolome"],
        groups_names=["AllSamples"],
        samples_names=[samples],
        features_names=[microbe_features, metabolite_features],
    )
    mofa.set_model_options(
        factors=5, spikeslab_weights=True, ard_weights=True,
        ard_factors=False,
    )
    mofa.set_train_options(
        iter=1000, convergence_mode="medium", dropR2=0.001,
        verbose=False, quiet=True, seed=SEED, gpu_mode=False,
        outfile=str(OUT / "mofa-model.hdf5"),
    )
    mofa.build()
    mofa.run()
    mofa.save(outfile=str(OUT / "mofa-model.hdf5"), save_data=False)

    z = np.asarray(mofa.model.nodes["Z"].getExpectations()["E"])
    weights = [np.asarray(item["E"]) for item in mofa.model.nodes["W"].getExpectations()]
    n_factors = z.shape[1]
    factor_names = [f"Factor{i + 1}" for i in range(n_factors)]
    factor_scores = pd.DataFrame(z, index=samples, columns=factor_names)
    factor_scores.insert(0, "StudyGroup", metadata.loc[samples, "StudyGroup"])
    factor_scores.to_csv(OUT / "mofa-factor-scores.tsv", sep="\t", index_label="SampleID")
    pd.DataFrame(weights[0], index=microbe_features, columns=factor_names).to_csv(
        OUT / "mofa-loadings-microbiome.tsv", sep="\t", index_label="FeatureID"
    )
    pd.DataFrame(weights[1], index=metabolite_features, columns=factor_names).to_csv(
        OUT / "mofa-loadings-metabolome.tsv", sep="\t", index_label="MetaboliteID"
    )

    r2 = np.asarray(mofa.model.calculate_variance_explained()[0])
    r2_rows = []
    for view_index, view_name in enumerate(["Microbiome", "Metabolome"]):
        for factor_index, factor_name in enumerate(factor_names):
            r2_rows.append({
                "View": view_name,
                "Factor": factor_name,
                "VarianceExplained": float(r2[view_index, factor_index]),
            })
    r2_frame = pd.DataFrame(r2_rows)
    r2_frame.to_csv(OUT / "mofa-variance-explained.tsv", sep="\t", index=False)

    group_tests = []
    for factor_name in factor_names:
        values = [
            factor_scores.loc[
                factor_scores["StudyGroup"] == group, factor_name
            ].to_numpy()
            for group in ["Control", "CD", "UC"]
        ]
        statistic, pvalue = kruskal(*values)
        group_tests.append({
            "Factor": factor_name,
            "KruskalWallisH": statistic,
            "PValue": pvalue,
        })
    group_tests = pd.DataFrame(group_tests)
    group_tests["QValue"] = multipletests(group_tests["PValue"], method="fdr_bh")[1]
    group_tests.to_csv(OUT / "mofa-factor-group-tests.tsv", sep="\t", index=False)

    top_pairs = pd.read_csv(top_pairs_path, sep="\t")
    conditional = pd.read_csv(conditional_path, sep="\t", index_col=0)
    native_summary = json.loads(native_summary_path.read_text(encoding="utf-8"))
    spearman_ranked = spearman.sort_values("AbsRho", ascending=False)
    spearman_top = set(zip(
        spearman_ranked.head(50)["microbe"],
        spearman_ranked.head(50)["metabolite"],
    ))
    mmvec_top = set(zip(top_pairs.head(50)["microbe"], top_pairs.head(50)["metabolite"]))
    overlap = len(spearman_top & mmvec_top)

    summary = {
        "status": "passed",
        "seed": SEED,
        "samples": len(samples),
        "microbe_features": len(microbe_features),
        "metabolite_features": len(metabolite_features),
        "spearman_tested_pairs": int(len(spearman)),
        "spearman_fdr_hits": int((spearman["QValue"] < 0.05).sum()),
        "cca_canonical_correlations": canonical_correlations,
        "mmvec_implementation": native_summary["implementation"],
        "mmvec_train_samples": native_summary["train_samples"],
        "mmvec_holdout_samples": native_summary["test_samples"],
        "mmvec_epochs_run": native_summary["epochs_run"],
        "mmvec_best_epoch_one_based": native_summary["best_epoch_one_based"],
        "mmvec_best_validation_mae": native_summary["best_validation_mae"],
        "mmvec_probability_max_row_sum_error": float(
            np.abs(conditional.sum(axis=1) - 1.0).max()
        ),
        "spearman_mmvec_top50_pair_overlap": overlap,
        "mofa_factors_retained": n_factors,
        "mofa_total_variance_explained_by_view": r2_frame.groupby("View")[
            "VarianceExplained"
        ].sum().to_dict(),
        "mofa_group_associated_factors_fdr": int((group_tests["QValue"] < 0.05).sum()),
        "versions": {
            "python": platform.python_version(),
            "mmvec": native_summary["versions"]["mmvec"],
            "tensorflow": native_summary["versions"]["tensorflow"],
            "mofapy2": mofapy2_version,
            "scikit_learn": sklearn_version,
        },
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
