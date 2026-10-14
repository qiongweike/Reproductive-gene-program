"""Export per-donor spermatogonial count matrices: discovery donors for NMF and GSE326372 NOA cases for program scoring.

Inputs:  data/processed/normal_spg_<study>.h5ad (GSE124263, GSE153947, GSE254315), data/processed/disease_spg.h5ad,
         data/processed/disease_gene_presence_by_study.tsv
Outputs: results/programs/nmf/input/<sample>.{counts.mtx.gz,genes.tsv,cells.tsv} (genes x cells, full-assay counts),
         results/programs/nmf/donor_manifest.tsv (all reference donors; >=200 SPG are discovery donors),
         results/programs/nmf/NOA_assignment_manifest.tsv
"""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

from utils.paths import PROCESSED, out_dir

DISCOVERY_STUDIES = ["GSE124263", "GSE153947", "GSE254315"]
DISCOVERY_DIAGNOSES = ["Normal", "OA"]
MIN_SPG_PER_DONOR = 200
SYMBOL_FIXES = {"GSE124263": {"C12orf44": "ATG101"}}  # source-verified rename, applied before upper-casing
NOA_STUDY = "GSE326372"

OUT = out_dir("programs", "nmf")
INPUT = out_dir("programs", "nmf", "input")


def merge_symbols(x, symbols):
    """Sum count columns (cells x genes) that share an upper-case symbol; genes come back sorted."""
    genes, idx = np.unique(np.asarray(symbols), return_inverse=True)
    mapper = sparse.csr_matrix((np.ones(len(idx)), (np.arange(len(idx)), idx)), shape=(len(idx), len(genes)))
    return x @ mapper, genes


def export(sample_key, x, genes, cells):
    """Write one donor (cells x genes counts) as genes x cells Matrix Market plus gene and cell tables."""
    base = INPUT / sample_key
    paths = dict(matrix_file=base.with_suffix(".counts.mtx.gz"), gene_file=base.with_suffix(".genes.tsv"),
                 cell_file=base.with_suffix(".cells.tsv"))
    with gzip.open(paths["matrix_file"], "wb") as f:
        io.mmwrite(f, x.T.tocoo(), field="real", precision=17)
    pd.DataFrame({"gene": genes}).to_csv(paths["gene_file"], sep="\t", index=False)
    cells.to_csv(paths["cell_file"], sep="\t", index=False)
    return {k: str(v) for k, v in paths.items()}


def discovery_donors():
    rows = []
    for study in DISCOVERY_STUDIES:
        a = ad.read_h5ad(PROCESSED / f"normal_spg_{study}.h5ad")
        assert a.obs.diagnosis.astype(str).isin(DISCOVERY_DIAGNOSES).all()
        fixes = SYMBOL_FIXES.get(study, {})
        x, genes = merge_symbols(a.layers["counts"].tocsr().astype(np.float64),
                                 [fixes.get(g, g).upper() for g in a.var_names.astype(str)])
        obs = a.obs.copy()
        obs["cell_id"] = a.obs_names.astype(str)
        obs["sample_key"] = study + "__" + obs.donor.astype(str)
        obs["full_assay_counts_recomputed"] = np.asarray(x.sum(axis=1)).ravel()
        assert (obs.full_assay_counts_recomputed > 0).all()
        for sample_key, n in obs.groupby("sample_key", observed=True).size().items():
            ids = np.flatnonzero(obs.sample_key.to_numpy() == sample_key)
            cells = obs.iloc[ids]
            row = dict(sample_key=sample_key, study=study, donor=str(cells.donor.iloc[0]),
                       diagnosis=str(cells.diagnosis.iloc[0]), n_spg=int(n), eligible=bool(n >= MIN_SPG_PER_DONOR),
                       n_features=len(genes), matrix_file="", gene_file="", cell_file="")
            if row["eligible"]:
                row.update(export(sample_key, x[ids], genes, cells))
            rows.append(row)
    return pd.DataFrame(rows)


def noa_cases():
    """GSE326372 NOA spermatogonia restricted to genes measured in that study; depth stays the full-assay total."""
    a = ad.read_h5ad(PROCESSED / "disease_spg.h5ad")
    a = a[(a.obs.study.astype(str) == NOA_STUDY) & (a.obs.diagnosis.astype(str) == "NOA")].copy()
    full = np.asarray(a.layers["counts"].sum(axis=1)).ravel().astype(float)
    assert np.allclose(full, a.obs.total_counts.to_numpy(float), rtol=1e-9, atol=1e-6)
    presence = pd.read_csv(PROCESSED / "disease_gene_presence_by_study.tsv", sep="\t")
    a = a[:, a.var_names.isin(set(presence.loc[presence[NOA_STUDY].astype(int) == 1, "gene_symbol"]))].copy()
    x, genes = merge_symbols(a.layers["counts"].tocsr().astype(np.float64), [g.upper() for g in a.var_names])
    rows = []
    for donor in sorted(a.obs.donor.astype(str).unique()):
        ids = np.flatnonzero(a.obs.donor.astype(str).to_numpy() == donor)
        cells = a.obs.iloc[ids].copy()
        cells["cell_id"] = a.obs_names[ids].astype(str)
        cells["sample_key"] = donor
        row = dict(sample_key=donor, study=NOA_STUDY, donor=donor, diagnosis="NOA", n_spg=len(ids), eligible=True,
                   n_features=len(genes))
        row.update(export(donor, x[ids], genes, cells))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    donors = discovery_donors()
    donors.to_csv(OUT / "donor_manifest.tsv", sep="\t", index=False)
    noa = noa_cases()
    noa.to_csv(OUT / "NOA_assignment_manifest.tsv", sep="\t", index=False)
    print(f"{donors.eligible.sum()} discovery donors ({donors.n_spg[donors.eligible].sum()} SPG); "
          f"{len(noa)} NOA cases ({noa.n_spg.sum()} SPG)")


if __name__ == "__main__":
    main()
