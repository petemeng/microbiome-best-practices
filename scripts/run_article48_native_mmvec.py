#!/usr/bin/env python
"""Run Article 48 with the official native MMvec TensorFlow 1.x model.

This script is intentionally compatible with Python 3.7 and MMvec 1.0.5.
It uses ``mmvec.multimodal.MMvec`` directly, fixes a stratified holdout, records
one training and validation diagnostic per epoch, and restores the epoch with
the lowest validation error before exporting conditional scores.
"""

import argparse
import json
import platform
import shutil
from pathlib import Path

import mmvec
import numpy as np
import pandas as pd
import tensorflow as tf
from mmvec.multimodal import MMvec
from scipy.sparse import coo_matrix


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "small" / "paired-ibd-multiomics"
SEED = 20260726


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "results" / "48-mmvec-mofa")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--latent-dim", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--learning-rate", type=float, default=1e-5)
    return parser.parse_args()


def softmax_rows(values):
    shifted = values - values.max(axis=1, keepdims=True)
    exponentiated = np.exp(shifted)
    return exponentiated / exponentiated.sum(axis=1, keepdims=True)


def stratified_split(metadata):
    rng = np.random.RandomState(SEED)
    train_ids = []
    for group_name in ("Control", "CD", "UC"):
        ids = np.array(sorted(metadata.index[metadata["StudyGroup"] == group_name]))
        rng.shuffle(ids)
        train_ids.extend(ids[: int(np.floor(0.70 * len(ids)))].tolist())
    train_set = set(train_ids)
    train_ids = [sample for sample in metadata.index if sample in train_set]
    test_ids = [sample for sample in metadata.index if sample not in train_set]
    return train_ids, test_ids


def main():
    args = parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    log_dir = out / "mmvec-native-tensorboard"
    if log_dir.exists():
        shutil.rmtree(str(log_dir))
    log_dir.mkdir(parents=True, exist_ok=True)

    otutab = pd.read_csv(DATA / "otutab.tsv", sep="\t", index_col=0)
    taxonomy = pd.read_csv(DATA / "taxonomy.tsv", sep="\t", index_col=0)
    metadata = pd.read_csv(DATA / "metadata.tsv", sep="\t", index_col=0)
    metabolites = pd.read_csv(DATA / "metabolites.tsv", sep="\t", index_col=0)
    annotation = pd.read_csv(DATA / "metabolite-annotation.tsv", sep="\t", index_col=0)

    samples = metadata.index.tolist()
    if list(otutab.columns) != samples or list(metabolites.columns) != samples:
        raise RuntimeError("Paired sample order differs across views")

    microbe_features = otutab.sum(axis=1).sort_values(ascending=False).head(100).index.tolist()
    metabolite_log = np.log1p(metabolites.T.astype(float))
    metabolite_features = (
        metabolite_log.var(axis=0, ddof=1).sort_values(ascending=False).head(120).index.tolist()
    )
    microbe_samples = otutab.loc[microbe_features, samples].T.astype(float)
    metabolite_samples = metabolites.loc[metabolite_features, samples].T.astype(float)

    train_ids, test_ids = stratified_split(metadata)
    split = pd.DataFrame({
        "SampleID": samples,
        "StudyGroup": metadata.loc[samples, "StudyGroup"].values,
        "MMvecSet": ["Train" if sample in set(train_ids) else "Test" for sample in samples],
    })
    split.to_csv(out / "mmvec-split.tsv", sep="\t", index=False)

    train_microbes = microbe_samples.loc[train_ids]
    test_microbes = microbe_samples.loc[test_ids]
    train_metabolites = metabolite_samples.loc[train_ids]
    test_metabolites = metabolite_samples.loc[test_ids]
    train_x = coo_matrix(train_microbes.values)
    test_x = coo_matrix(test_microbes.values)

    np.random.seed(SEED)
    graph = tf.Graph()
    with graph.as_default():
        tf.set_random_seed(SEED)
        config = tf.ConfigProto(
            intra_op_parallelism_threads=1,
            inter_op_parallelism_threads=1,
            allow_soft_placement=False,
        )
        with tf.Session(config=config) as session:
            model = MMvec(
                latent_dim=args.latent_dim,
                u_scale=1.0,
                v_scale=1.0,
                learning_rate=args.learning_rate,
                beta_1=0.9,
                beta_2=0.95,
                device_name="/cpu:0",
                batch_size=args.batch_size,
                clipnorm=10.0,
                save_path=str(log_dir),
            )
            model(
                session,
                train_x,
                train_metabolites.values,
                test_x,
                test_metabolites.values,
            )
            iterations_per_epoch = max(1, model.nnz // args.batch_size)
            rows = []
            best = None
            for epoch in range(1, args.epochs + 1):
                for _ in range(iterations_per_epoch):
                    session.run(model.train)
                train_loss, validation_mae, u, ub, v, vb = session.run([
                    model.log_loss,
                    model.cv,
                    model.qUmain,
                    model.qUbias,
                    model.qVmain,
                    model.qVbias,
                ])
                rows.append((epoch, float(train_loss), float(validation_mae)))
                if best is None or validation_mae < best["ValidationMAE"]:
                    best = {
                        "Epoch": epoch,
                        "TrainLoss": float(train_loss),
                        "ValidationMAE": float(validation_mae),
                        "U": u.copy(),
                        "Ubias": ub.copy(),
                        "V": v.copy(),
                        "Vbias": vb.copy(),
                    }

            model.U = best["U"]
            model.Ubias = best["Ubias"]
            model.V = best["V"]
            model.Vbias = best["Vbias"]
            ranks = model.ranks()
            model.writer.close()

    loss = pd.DataFrame(rows, columns=["Epoch", "TrainLoss", "ValidationMAE"])
    loss.to_csv(out / "mmvec-loss.tsv", sep="\t", index=False)
    ranks_frame = pd.DataFrame(ranks, index=microbe_features, columns=metabolite_features)
    ranks_frame.index.name = "Microbe"
    ranks_frame.to_csv(out / "mmvec-log-conditionals.tsv", sep="\t")
    probabilities = pd.DataFrame(
        softmax_rows(ranks), index=microbe_features, columns=metabolite_features
    )
    probabilities.index.name = "Microbe"
    probabilities.to_csv(out / "mmvec-conditional-probabilities.tsv", sep="\t")

    pair_rows = []
    for i, microbe in enumerate(microbe_features):
        for j, metabolite in enumerate(metabolite_features):
            pair_rows.append((
                microbe,
                metabolite,
                float(ranks[i, j]),
                float(probabilities.iloc[i, j]),
                taxonomy.loc[microbe, "DisplayName"],
                annotation.loc[metabolite, "DisplayName"],
            ))
    top_pairs = pd.DataFrame(pair_rows, columns=[
        "microbe", "metabolite", "score", "conditional_probability",
        "MicrobeName", "MetaboliteName",
    ]).sort_values(["score", "microbe", "metabolite"], ascending=[False, True, True])
    top_pairs.head(100).to_csv(out / "mmvec-top-pairs.tsv", sep="\t", index=False)

    microbe_embeddings = pd.DataFrame(
        best["U"], index=microbe_features,
        columns=["Axis%d" % (i + 1) for i in range(args.latent_dim)],
    )
    microbe_embeddings["Bias"] = best["Ubias"].ravel()
    metabolite_values = np.hstack((
        np.zeros((args.latent_dim, 1)), best["V"]
    )).T
    metabolite_embeddings = pd.DataFrame(
        metabolite_values,
        index=metabolite_features,
        columns=["Axis%d" % (i + 1) for i in range(args.latent_dim)],
    )
    metabolite_embeddings["Bias"] = np.hstack((0.0, best["Vbias"].ravel()))
    microbe_embeddings.to_csv(out / "mmvec-microbe-embeddings.tsv", sep="\t")
    metabolite_embeddings.to_csv(out / "mmvec-metabolite-embeddings.tsv", sep="\t")

    summary = {
        "status": "passed",
        "implementation": "biocore/mmvec native MMvec class",
        "seed": SEED,
        "samples": len(samples),
        "train_samples": len(train_ids),
        "test_samples": len(test_ids),
        "microbe_features": len(microbe_features),
        "metabolite_features": len(metabolite_features),
        "epochs_run": args.epochs,
        "iterations_per_epoch": iterations_per_epoch,
        "best_epoch_one_based": int(best["Epoch"]),
        "best_validation_mae": best["ValidationMAE"],
        "best_training_loss": best["TrainLoss"],
        "latent_dim": args.latent_dim,
        "batch_size": args.batch_size,
        "learning_rate": args.learning_rate,
        "probability_max_row_sum_error": float(
            np.max(np.abs(probabilities.sum(axis=1).values - 1.0))
        ),
        "versions": {
            "python": platform.python_version(),
            "mmvec": mmvec.__version__,
            "tensorflow": tf.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (out / "summary-mmvec-native.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
