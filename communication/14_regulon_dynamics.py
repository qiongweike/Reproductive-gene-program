"""Regulon AUCell along the primary Slingshot lineages (Fig 5h, Supplementary Fig 2c): bin means of every regulon,
non-MP curves of the 12 structural TFs in all qualified lineages, and in each donor's selected lineage the regulon
peak bin minus the MP3 score peak bin (median of tied maxima; range spans all ties).
Inputs: results/trajectory/ (04: lineage_bin_cells, trajectory_bin_scores, selected_primary_lineages), regulon/ (step 13)
Outputs: results/trajectory/complete_regulon_bin_dynamics.tsv.gz, results/communication/{regulon_trajectory_curves.tsv,
  regulon_peak_offsets.tsv}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir
import numpy as np
import pandas as pd

SCENARIO = "primary"
N_BINS = 10
NON_MP = "__exclude_all_MP"
TRAJ = RESULTS / "trajectory"
OUT = out_dir("communication")


def is_max(values):
    return np.isclose(values, values.max(), rtol=1e-9, atol=1e-12)


def bin_dynamics(auc):
    cells = pd.read_csv(TRAJ / "lineage_bin_cells.tsv.gz", sep="\t")
    cells = cells[cells.scenario.eq(SCENARIO)]
    bins = pd.read_csv(TRAJ / "trajectory_bin_scores.tsv", sep="\t")
    bins = bins[bins.scenario.eq(SCENARIO) & bins.metric.eq("MP_2")]  # one row per lineage bin
    rows = []
    for b in bins.itertuples():
        ids = cells.loc[cells.sample_key.eq(b.sample_key) & cells.lineage.eq(b.lineage) & cells.bin.eq(b.bin), "cell_id"]
        assert len(ids) == b.n_cells
        for name, value in auc.loc[ids].mean(0).items():
            rows.append(dict(sample_key=b.sample_key, lineage=b.lineage, bin=b.bin, time_scaled=b.time_scaled,
                             n_cells=b.n_cells, trajectory_eligible=b.eligible, bin_eligible=b.bin_eligible,
                             regulon=name, mean_AUCell=value, exclude_all_MP=NON_MP in name))
    return pd.DataFrame(rows)


def peak_offsets(curves, selected):
    ref = pd.read_csv(TRAJ / "trajectory_bin_scores.tsv", sep="\t")
    ref = ref[ref.scenario.eq(SCENARIO) & ref.eligible & ref.bin_eligible & ref.metric.eq("MP_3")]
    lineages = selected[["sample_key", "lineage", "n_cells"]].rename(columns={"n_cells": "lineage_n_cells"})
    main = curves.merge(lineages, on=["sample_key", "lineage"], validate="many_to_one")
    rows = []
    for (sample, lineage, tf), x in main.groupby(["sample_key", "lineage", "TF"]):
        r = ref[ref.sample_key.eq(sample) & ref.lineage.eq(lineage)]
        assert len(x) == N_BINS and len(r) == N_BINS
        ok = x.mean_AUCell.max() - x.mean_AUCell.min() > 1e-12 and r.mean_score.max() - r.mean_score.min() > 1e-12
        tp = x.loc[is_max(x.mean_AUCell), "bin"].to_numpy()
        rp = r.loc[is_max(r.mean_score), "bin"].to_numpy()
        rows.append(dict(sample_key=sample, donor=x.donor.iloc[0], selected_lineage=lineage,
                         lineage_n_cells=int(x.lineage_n_cells.iloc[0]), TF=tf,
                         peak_offset=float(np.median(tp) - np.median(rp)) if ok else np.nan,
                         offset_min=int(tp.min() - rp.max()) if ok else np.nan,
                         offset_max=int(tp.max() - rp.min()) if ok else np.nan,
                         regulon_peak_bins=";".join(map(str, tp)), MP3_peak_bins=";".join(map(str, rp)), n_bins=len(x),
                         status="observed_bin_peak" if ok else "constant_not_estimable"))
    return pd.DataFrame(rows)


def main():
    auc = pd.read_csv(RESULTS / "communication/regulon/normal_full_and_nonMP_AUCell.tsv.gz", sep="\t", index_col=0)
    dynamics = bin_dynamics(auc)
    dynamics.to_csv(TRAJ / "complete_regulon_bin_dynamics.tsv.gz", sep="\t", index=False, na_rep="NA")
    tfs = (pd.read_csv(RESULTS / "communication/regulon/structural_TF_selection.tsv", sep="\t")
           .query("selected").TF.tolist())
    a = dynamics.copy()
    a["TF"] = a.regulon.str.split("(").str[0]
    a = a[a.exclude_all_MP & a.trajectory_eligible & a.bin_eligible & a.TF.isin(tfs)].copy()
    a["donor"] = a.sample_key.str.split("__").str[-1]
    a.to_csv(OUT / "regulon_trajectory_curves.tsv", sep="\t", index=False)
    selected = pd.read_csv(TRAJ / "selected_primary_lineages.tsv", sep="\t")
    offsets = peak_offsets(a, selected)
    assert len(offsets) == len(selected) * len(tfs)
    offsets.to_csv(OUT / "regulon_peak_offsets.tsv", sep="\t", index=False, na_rep="NA")
    print(f"{dynamics.regulon.nunique()} regulons x {len(dynamics) // dynamics.regulon.nunique()} lineage bins; "
          f"{len(offsets)} peak offsets")


if __name__ == "__main__":
    main()
