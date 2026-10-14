"""Cell eligibility and per-donor CellChat inputs (log1p CP10k on full-assay depth, measured CellChat genes).

Eligible: QC pass and confident atlas identity; spermatogonia also need the strict primary-SPG flag and an SPG call and
enter as U/T/D groups; groups >=30 cells. Unassigned-state spermatogonia (SPG_unknown) are kept out of the input.
Inputs: GSE202647_annotated.h5ad, results/atlas/cell_annotations.tsv.gz, database_genes.tsv (step 01)
Outputs: results/communication/cellchat_inputs/ (eligibility, coverage, manifest, <donor>/ inputs)
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, out_dir
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy.io import mmwrite

DONORS = ["GSE202647_OA1", "GSE202647_OA2", "GSE202647_NOA1", "GSE202647_NOA2"]
BASE_SEED = 20261003  # donor j (order above) uses BASE_SEED + j
MIN_CELLS = 30
STATE_GROUPS = {"undifferentiated": "SPG_U", "transition": "SPG_T", "differentiating": "SPG_D",
                "unknown": "SPG_unknown"}
UNASSIGNED = "SPG_unknown"
OUT = out_dir("communication", "cellchat_inputs")


def eligibility(meta):
    meta = meta[meta.donor.isin(DONORS)].copy()
    assert meta.cell_id.is_unique and meta.source_cell_id.is_unique
    spg = meta.primary_spg
    meta["communication_group"] = meta.atlas_cell_type
    meta.loc[spg, "communication_group"] = meta.loc[spg, "primary_spg_state"].map(STATE_GROUPS)
    assert meta.loc[spg, "communication_group"].notna().all()
    meta["identity_eligible"] = ((spg & meta.annotation_confident & meta.atlas_cell_type.eq("SPG"))
                                 | (~spg & meta.annotation_confident & ~meta.atlas_cell_type.isin(["SPG", "Uncertain"])))
    meta["pre_group_eligible"] = meta.qc_analysis_pass & meta.identity_eligible
    counts = (meta[meta.pre_group_eligible].groupby(["donor", "diagnosis", "communication_group"]).size()
              .rename("n_cells").reset_index())
    counts["eligible_ge30"] = counts.n_cells.ge(MIN_CELLS)
    meta = meta.merge(counts[["donor", "communication_group", "n_cells", "eligible_ge30"]],
                      on=["donor", "communication_group"], how="left")
    large = meta.eligible_ge30.eq(True)  # groups absent from the counts (not pre-eligible) are not large
    meta["communication_eligible"] = meta.pre_group_eligible & large
    meta["exclusion_reason"] = np.select([~meta.qc_analysis_pass, ~meta.identity_eligible, ~large],
                                         ["atlas_QC_failed", "identity_uncertain_or_nonprimary_SPG", "group_under30"],
                                         default="")
    return meta, counts


def main():
    db = pd.read_csv(OUT / "database_genes.tsv", sep="\t").gene.astype(str)
    meta, counts = eligibility(pd.read_csv(RESULTS / "atlas/cell_annotations.tsv.gz", sep="\t", low_memory=False))
    meta.to_csv(OUT / "all_cell_eligibility.tsv.gz", sep="\t", index=False)
    counts.to_csv(OUT / "donor_group_coverage.tsv", sep="\t", index=False)
    manifest = []
    with h5py.File(PROCESSED / "GSE202647_annotated.h5ad", "r") as f:
        obs = read_elem(f["obs"])
        genes = read_elem(f["var"]).index.astype(str).to_numpy()
        assert len(set(genes)) == len(genes)
        layer = sparse_dataset(f["layers/counts"])
        gi = np.flatnonzero(np.isin(genes, db))
        pd.DataFrame({"gene": genes[gi]}).to_csv(OUT / "input_genes.tsv", sep="\t", index=False)
        pd.DataFrame({"gene": db, "assay_present": db.isin(genes)}).to_csv(
            OUT / "database_assay_coverage.tsv", sep="\t", index=False)
        for j, donor in enumerate(DONORS):
            m = meta[meta.donor.eq(donor) & meta.communication_eligible & meta.communication_group.ne(UNASSIGNED)]
            idx = obs.index.get_indexer(m.source_cell_id)
            assert (idx >= 0).all()
            a = layer[idx, :].tocsr()
            depth = np.asarray(a.sum(1)).ravel()
            assert np.allclose(depth, m.full_assay_total_counts, rtol=1e-10, atol=1e-6)
            norm = a[:, gi].tocsr().multiply(1e4 / depth[:, None]).tocsr()
            norm.data = np.log1p(norm.data)
            d = OUT / donor
            d.mkdir(exist_ok=True)
            mmwrite(d / "normalized.mtx", norm.T)
            m.to_csv(d / "cells.tsv", sep="\t", index=False)
            manifest.append(dict(donor=donor, diagnosis=m.diagnosis.iloc[0], seed=BASE_SEED + j, n_cells=len(m),
                                 n_groups=m.communication_group.nunique()))
    manifest = pd.DataFrame(manifest)
    manifest.to_csv(OUT / "donor_manifest.tsv", sep="\t", index=False)
    print(manifest.to_string(index=False), f"\nsignalling genes measured: {len(gi)}/{len(db)}")


if __name__ == "__main__":
    main()
