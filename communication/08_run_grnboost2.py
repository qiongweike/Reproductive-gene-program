"""GRNBoost2 (arboreto) models for all 10,470 targets with the 1,164 measured hg38 TFs as predictors: learning rate
0.01, <=5,000 trees, max_features 0.1, subsample 0.9, early-stop window 25, seed 20261003 for every target.
Targets run in REGULON_WORKERS processes (default 10); models are independent, so this does not change results.
Inputs: scenic/{expression.npy, genes.tsv}, data/external/cistarget/allTFs_hg38.txt
Outputs: results/communication/scenic/full_GRNBoost_adjacencies.tsv.gz
"""
import os
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import EXTERNAL, RESULTS
import multiprocessing as mp
import numpy as np
import pandas as pd
from arboreto.core import infer_partial_network

SEED = 20261003
GBM_KWARGS = {"learning_rate": 0.01, "n_estimators": 5000, "max_features": 0.1, "subsample": 0.9}
EARLY_STOP_WINDOW = 25
WORKERS = int(os.environ.get("REGULON_WORKERS", "10"))
CHUNK = 20
REG = RESULTS / "communication" / "scenic"


def load_inputs():
    X = np.load(REG / "expression.npy", mmap_mode="r")
    genes = pd.read_csv(REG / "genes.tsv", sep="\t").gene.to_numpy()
    tfs = set((EXTERNAL / "cistarget" / "allTFs_hg38.txt").read_text().split())
    tf_idx = [i for i, g in enumerate(genes) if g in tfs]
    return X, genes, np.asarray(X[:, tf_idx], order="C"), genes[tf_idx].tolist()


def fit_targets(targets, X, genes, tf_matrix, tf_names):
    """Adjacencies (TF, target, importance) of the given target columns."""
    return pd.concat([infer_partial_network(
        regressor_type="GBM", regressor_kwargs=GBM_KWARGS, tf_matrix=tf_matrix, tf_matrix_gene_names=tf_names,
        target_gene_name=str(genes[j]), target_gene_expression=X[:, j], include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW, seed=SEED) for j in targets])


def _fit_chunk(targets):
    return fit_targets(targets, *INPUTS)


def main():
    global INPUTS
    INPUTS = load_inputs()
    n = len(INPUTS[1])
    chunks = [list(range(s, min(s + CHUNK, n))) for s in range(0, n, CHUNK)]
    with mp.get_context("fork").Pool(WORKERS) as pool:
        adjacency = pd.concat(pool.map(_fit_chunk, chunks, chunksize=1), ignore_index=True)
    adjacency.to_csv(REG / "full_GRNBoost_adjacencies.tsv.gz", sep="\t", index=False)
    print(f"GRNBoost2: {len(adjacency)} edges, {adjacency.target.nunique()} targets, {len(INPUTS[3])} TFs")


if __name__ == "__main__":
    main()
