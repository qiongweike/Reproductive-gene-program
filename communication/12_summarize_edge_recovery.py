"""Recovery of the 122 motif-supported TF -> MP3 edges in the 17 refitted models (Fig 5g). An edge is evaluable when TF
and target are each detected in >=5% of the model's cells and the target was fitted (unselected = importance 0);
a donor recovers an edge only when both seeds select it.
Inputs: regulon_edges/fits/ (step 11), scenic/ GRN input,
  regulon/{MP3_motif_edges.tsv, structural_TF_selection.tsv}
Outputs: results/communication/regulon_edges/ (MP3_edges_by_model.tsv and three summaries)
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir
import numpy as np
import pandas as pd

DONOR_SEEDS = [20261003, 20261017]
LEAVE_OUT_SEED = 20261003
MIN_DETECTION = 0.05
REG = RESULTS / "communication" / "scenic"
OUT = out_dir("communication", "regulon_edges")


def model_blocks(meta):
    """(name, scope, key, seed) of the 14 donor-seed and 3 leave-one-study-out models (as in step 11)."""
    return ([(f"{d}__seed{s}", "donor", d, s) for d in sorted(meta.sample_key.unique()) for s in DONOR_SEEDS]
            + [(f"leave_{s}", "leave_study_out", s, LEAVE_OUT_SEED) for s in sorted(meta.study.unique())])


def model_cells(meta, scope, key):
    return np.flatnonzero(meta.sample_key.eq(key).to_numpy() if scope == "donor" else ~meta.study.eq(key).to_numpy())


def edge_rows(block, fit, edges, X, gix, meta):
    name, scope, key, seed = block
    ix = model_cells(meta, scope, key)
    needed = sorted(set(edges.TF) | set(edges.target))
    det = dict(zip(needed, (X[np.ix_(ix, [gix[g] for g in needed])] > 0).mean(0)))
    assert not fit.duplicated(["TF", "target"]).any()
    fit["importance_share"] = fit.importance / fit.groupby("target").importance.transform("sum")
    fit["positive_edge_rank"] = fit.groupby("target").importance.rank(ascending=False, method="min")
    lookup = fit.set_index(["TF", "target"])
    fitted = set(fit.target)
    rows = []
    for e in edges.itertuples():
        ok = det[e.TF] >= MIN_DETECTION and det[e.target] >= MIN_DETECTION and e.target in fitted
        found = (e.TF, e.target) in lookup.index
        hit = lookup.loc[(e.TF, e.target)] if ok and found else None
        rows.append(dict(block=name, scope=scope, key=key, seed=seed, n_cells=len(ix), TF=e.TF, target=e.target,
                         TF_detection=det[e.TF], target_detection=det[e.target], eligible=ok,
                         edge_selected=bool(found) if ok else np.nan,
                         importance=float(hit.importance) if hit is not None else (0.0 if ok else np.nan),
                         importance_share=float(hit.importance_share) if hit is not None else (0.0 if ok else np.nan),
                         positive_edge_rank=float(hit.positive_edge_rank) if hit is not None else np.nan,
                         reason="PASS" if ok else "low_detection_or_no_target_model"))
    return rows


def main():
    X = np.load(REG / "expression.npy", mmap_mode="r")
    gix = {g: i for i, g in enumerate(pd.read_csv(REG / "genes.tsv", sep="\t").gene)}
    meta = pd.read_csv(REG / "cells.tsv.gz", sep="\t")
    edges = (pd.read_csv(RESULTS / "communication/regulon/MP3_motif_edges.tsv", sep="\t")
             [["TF", "assay_target"]].rename(columns={"assay_target": "target"}))
    rows = []
    for block in model_blocks(meta):
        fit = pd.read_csv(OUT / "fits" / f"{block[0]}.tsv.gz", sep="\t")
        rows += edge_rows(block, fit, edges, X, gix, meta)
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "MP3_edges_by_model.tsv", sep="\t", index=False, na_rep="NA")
    q = tab[tab.scope.eq("donor") & tab.eligible].copy()
    q["edge_selected"] = q.edge_selected.astype(bool)
    donor = q.groupby(["TF", "target", "key"]).agg(n_seeds=("seed", "nunique"), both_seed_selected=("edge_selected", "all"),
                                                   mean_importance_share=("importance_share", "mean")).reset_index()
    assert donor.n_seeds.eq(len(DONOR_SEEDS)).all()
    donor.to_csv(OUT / "edge_donor_seed_concordance.tsv", sep="\t", index=False)
    summary = donor.groupby(["TF", "target"]).agg(n_eligible_donors=("key", "nunique"),
                                                  n_donors_both_seeds=("both_seed_selected", "sum")).reset_index()
    left_out = (tab[tab.scope.eq("leave_study_out") & tab.eligible].groupby(["TF", "target"])
                .agg(n_eligible_leave_study=("block", "nunique"), n_leave_study_selected=("edge_selected", "sum"))
                .reset_index())
    summary = summary.merge(left_out, on=["TF", "target"], how="outer")
    summary.to_csv(OUT / "edge_recovery_summary.tsv", sep="\t", index=False, na_rep="NA")
    selection = pd.read_csv(RESULTS / "communication/regulon/structural_TF_selection.tsv", sep="\t")
    shown = summary[summary.TF.isin(selection.query("selected").TF)]
    shown.to_csv(OUT / "structural_edge_recovery_summary.tsv", sep="\t", index=False)
    print(f"{tab.block.nunique()} models x {len(edges)} edges; evaluable rows {int(tab.eligible.sum())}")


if __name__ == "__main__":
    main()
