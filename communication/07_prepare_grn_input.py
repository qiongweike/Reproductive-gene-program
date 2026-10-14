"""Pooled GRN input: 8,139 discovery spermatogonia (seven NMF donors) x 10,470 genes, log1p CP10k (full-assay depth).
Genes: upper-case symbols unique in and shared by all seven assays, >=5% pooled detection, >=1% detection in >=4 donors.
Inputs: results/programs/nmf/donor_manifest.tsv and the per-donor count, gene and cell files it lists
Outputs: results/communication/scenic/{expression.npy, genes.tsv, cells.tsv.gz}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import ROOT, RESULTS, out_dir
import numpy as np
import pandas as pd
from scipy import io

N_CELLS = 8139
MIN_POOLED_DETECTION = 0.05
MIN_DONOR_DETECTION = 0.01
MIN_DONORS = 4
OUT = out_dir("communication", "scenic")


def unique_upper_index(gene_file):
    genes = pd.read_csv(gene_file, sep="\t").gene.astype(str).str.upper()
    counts = genes.value_counts()
    return {g: i for i, g in enumerate(genes) if counts[g] == 1}


def main():
    manifest = pd.read_csv(RESULTS / "programs/nmf/donor_manifest.tsv", sep="\t").query("eligible")
    assert len(manifest) == 7 and manifest.n_spg.sum() == N_CELLS
    maps = [unique_upper_index(ROOT / row.gene_file) for row in manifest.itertuples()]
    common = sorted(set.intersection(*(set(m) for m in maps)))
    mats, cells, detection = [], [], []
    for row, gmap in zip(manifest.itertuples(), maps):
        x = io.mmread(ROOT / row.matrix_file).T.tocsr()
        meta = pd.read_csv(ROOT / row.cell_file, sep="\t")
        depth = np.asarray(x.sum(1)).ravel()
        assert np.allclose(depth, meta.full_assay_counts_recomputed)
        z = x[:, [gmap[g] for g in common]].toarray().astype(np.float32)
        detection.append(np.mean(z > 0, 0))
        mats.append(np.log1p(z * (1e4 / depth[:, None])).astype(np.float32))
        meta["sample_key"] = row.sample_key
        meta["study"] = row.study
        cells.append(meta)
    X = np.concatenate(mats)
    detection = np.asarray(detection)
    keep = (np.mean(X > 0, 0) >= MIN_POOLED_DETECTION) & ((detection >= MIN_DONOR_DETECTION).sum(0) >= MIN_DONORS)
    X = X[:, keep]
    np.save(OUT / "expression.npy", X)
    pd.DataFrame({"gene": np.array(common)[keep], "overall_detection": np.mean(X > 0, 0),
                  "n_donors_detection01": (detection[:, keep] >= MIN_DONOR_DETECTION).sum(0)}).to_csv(
        OUT / "genes.tsv", sep="\t", index=False)
    pd.concat(cells).to_csv(OUT / "cells.tsv.gz", sep="\t", index=False)
    print(f"GRN input: {X.shape[0]} cells x {X.shape[1]} genes")


if __name__ == "__main__":
    main()
