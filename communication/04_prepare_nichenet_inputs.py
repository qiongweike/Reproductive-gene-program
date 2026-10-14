"""NicheNet inputs: detection and mean log1p CP10k (full-assay depth) of every gene with a unique upper-case symbol,
per donor x communication group eligible for CellChat (>=30 cells), and the assay symbol/Ensembl identity table.
Inputs: data/processed/GSE202647_annotated.h5ad, cellchat_inputs/all_cell_eligibility.tsv.gz (step 02)
Outputs: results/communication/nichenet_inputs/donor_group_expression.tsv.gz,
  results/communication/nichenet/assay_gene_identity.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, out_dir
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

OUT = out_dir("communication", "nichenet_inputs")
IDENTITY = out_dir("communication", "nichenet") / "assay_gene_identity.tsv"


def main():
    meta = pd.read_csv(RESULTS / "communication/cellchat_inputs/all_cell_eligibility.tsv.gz",
                       sep="\t", low_memory=False)
    rows = []
    with h5py.File(PROCESSED / "GSE202647_annotated.h5ad", "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        symbols = var.index.astype(str)
        upper = pd.Series(symbols).str.upper()
        pd.DataFrame({"assay_symbol": symbols, "ensembl": var.gene_ids.astype(str).to_numpy(),
                      "symbol_upper": upper.to_numpy()}).to_csv(IDENTITY, sep="\t", index=False)
        layer = sparse_dataset(f["layers/counts"])
        gi = np.flatnonzero(~upper.duplicated(keep=False))
        for (donor, group), m in meta[meta.communication_eligible].groupby(["donor", "communication_group"]):
            idx = obs.index.get_indexer(m.source_cell_id)
            assert (idx >= 0).all()
            a = layer[idx, :].tocsr()
            depth = np.asarray(a.sum(1)).ravel()
            assert np.allclose(depth, m.full_assay_total_counts)
            x = a[:, gi].tocsr()
            detection = np.asarray((x > 0).mean(0)).ravel()
            norm = x.multiply(1e4 / depth[:, None]).tocsr()
            norm.data = np.log1p(norm.data)
            rows.append(pd.DataFrame(dict(donor=donor, diagnosis=m.diagnosis.iloc[0], group=group, n_cells=len(m),
                                          gene=upper.iloc[gi].to_numpy(), detection=detection,
                                          mean_logCP10k=np.asarray(norm.mean(0)).ravel())))
    expr = pd.concat(rows)
    expr.to_csv(OUT / "donor_group_expression.tsv.gz", sep="\t", index=False)
    print(f"{len(rows)} donor groups x {len(gi)} genes")


if __name__ == "__main__":
    main()
