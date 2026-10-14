"""Refit the GRNBoost2 models of the 34 MP3 targets with motif-supported edges: 7 donors x seeds 20261003/20261017 and
3 leave-one-study-out models (seed 20261003), all 1,164 TFs, pooled-network GBM settings; constant targets skipped.
Models run in REGULON_WORKERS processes (default 6); this does not change results.
Inputs: scenic/{expression.npy, genes.tsv, cells.tsv.gz}, regulon/MP3_motif_edges.tsv, cistarget TF list
Outputs: results/communication/regulon_edges/fits/<model>.tsv.gz
"""
import os
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import EXTERNAL, RESULTS, out_dir
import multiprocessing as mp
import numpy as np
import pandas as pd
from arboreto.core import infer_partial_network

DONOR_SEEDS = [20261003, 20261017]
LEAVE_OUT_SEED = 20261003
GBM_KWARGS = {"learning_rate": 0.01, "n_estimators": 5000, "max_features": 0.1, "subsample": 0.9}
EARLY_STOP_WINDOW = 25
WORKERS = int(os.environ.get("REGULON_WORKERS", "6"))
REG = RESULTS / "communication" / "scenic"
FITS = out_dir("communication", "regulon_edges", "fits")


def model_blocks(meta):
    """(name, scope, key, seed) of the 14 donor-seed and 3 leave-one-study-out models."""
    return ([(f"{d}__seed{s}", "donor", d, s) for d in sorted(meta.sample_key.unique()) for s in DONOR_SEEDS]
            + [(f"leave_{s}", "leave_study_out", s, LEAVE_OUT_SEED) for s in sorted(meta.study.unique())])


def model_cells(meta, scope, key):
    return np.flatnonzero(meta.sample_key.eq(key).to_numpy() if scope == "donor" else ~meta.study.eq(key).to_numpy())


def load_inputs():
    X = np.load(REG / "expression.npy", mmap_mode="r")
    genes = pd.read_csv(REG / "genes.tsv", sep="\t").gene.tolist()
    meta = pd.read_csv(REG / "cells.tsv.gz", sep="\t")
    tfs = set((EXTERNAL / "cistarget" / "allTFs_hg38.txt").read_text().split())
    tf_idx = [i for i, g in enumerate(genes) if g in tfs]
    edges = pd.read_csv(RESULTS / "communication/regulon/MP3_motif_edges.tsv", sep="\t")
    targets = sorted(edges.assay_target.unique())
    assert len(edges) == 122 and len(targets) == 34 and len(tf_idx) == 1164
    return X, genes, meta, tf_idx, targets


def fit_block(block, X, genes, meta, tf_idx, targets):
    """Adjacencies of the MP3 targets in one model; constant targets are not fitted."""
    _, scope, key, seed = block
    ix = model_cells(meta, scope, key)
    tf_matrix = np.asarray(X[np.ix_(ix, tf_idx)], order="C")
    tf_names = [genes[i] for i in tf_idx]
    gix = {g: i for i, g in enumerate(genes)}
    out = []
    for gene in targets:
        y = np.asarray(X[ix, gix[gene]])
        if np.std(y) < 1e-12:
            continue
        out.append(infer_partial_network(
            regressor_type="GBM", regressor_kwargs=GBM_KWARGS, tf_matrix=tf_matrix, tf_matrix_gene_names=tf_names,
            target_gene_name=gene, target_gene_expression=y, include_meta=False,
            early_stop_window_length=EARLY_STOP_WINDOW, seed=seed))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["TF", "target", "importance"])


def _run(block):
    fit_block(block, *INPUTS).to_csv(FITS / f"{block[0]}.tsv.gz", sep="\t", index=False)
    return block[0]


def main():
    global INPUTS
    INPUTS = load_inputs()
    blocks = model_blocks(INPUTS[2])
    with mp.get_context("fork").Pool(WORKERS) as pool:
        done = pool.map(_run, blocks, chunksize=1)
    print(f"refitted {len(done)} models x {len(INPUTS[4])} MP3 targets")


if __name__ == "__main__":
    main()
