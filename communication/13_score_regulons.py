"""Regulon AUCell without program members (top 5% of seeded-tie rankings), its within-state association with the
programs (Figs 5f, 5i) and its projection onto GSE202647 spermatogonia (Fig 4g). Association: per donor x state with
>=50 cells, partial rank correlation of AUCell with raw program means given ranked depth, detected genes, G1/S and G2/M.
Inputs: scenic/ (steps 07, 09), regulon/structural_TF_selection.tsv, results/programs/,
  resources/candidate_143_genes.txt, all_cell_eligibility.tsv.gz, GSE202647_annotated.h5ad
Outputs: results/communication/regulon/ (AUCell, association, bubble and projection tables; projection_U_AUCell.tsv)
"""
import os
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"  # single-threaded BLAS keeps the least-squares residuals bit-reproducible
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESOURCES, RESULTS, ROOT, out_dir
import pickle
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import io, stats
np.__dict__.setdefault("object", object)  # pySCENIC 0.12.1 still refers to np.object (removed in NumPy 1.24)
from ctxcore.genesig import Regulon
from pyscenic.aucell import aucell4r, create_rankings

SEED = 20261003
AUC_THRESHOLD = 0.05
MIN_TARGETS = 10
MIN_STATE_CELLS = 50
MIN_MEASURED_FRACTION = 0.8
STATES = ["undifferentiated", "transition", "differentiating"]
PROGRAMS = ["MP_1", "MP_2", "MP_3"]
CYCLE = ["G1/S", "G2/M"]
NON_MP = "__exclude_all_MP"
EXTENDED = "__exclude_all_MP_and_fixed143"
CASE_STATES = {"SPG_U": "undifferentiated", "SPG_T": "transition", "SPG_D": "differentiating"}
U_CELLS_PER_DONOR = {"GSE202647_NOA1": 160, "GSE202647_NOA2": 28, "GSE202647_OA1": 137, "GSE202647_OA2": 107}
REG = RESULTS / "communication" / "scenic"
OUT = out_dir("communication", "regulon")


def save(d, name, index=False):
    d.to_csv(OUT / name, sep="\t", index=index, na_rep="NA")


def canonical_symbol():
    aliases = pd.read_csv(RESULTS / "programs/identity/verified_aliases.tsv", sep="\t")
    amap = dict(zip(aliases.old_symbol.str.upper(), aliases.approved_symbol))
    return lambda g: amap.get(str(g).upper(), str(g).upper())


def regulon_copy(reg, weights, suffix=""):
    return Regulon(name=reg.name + suffix, gene2weight=weights, gene2occurrence={},
                   transcription_factor=reg.transcription_factor, context=reg.context)


def program_scores(mp, allmp, canon):
    """Per-cell raw program means and non-MP cycle means (mean log1p CP10k, full-assay depth), depth and genes."""
    manifest = pd.read_csv(RESULTS / "programs/donor_manifest.tsv", sep="\t").query("eligible")
    marks = pd.read_csv(RESULTS / "programs/identity/state_markers.tsv", sep="\t")
    sets = {p: mp.loc[mp.MP.eq(p), "gene"].map(canon).tolist() for p in PROGRAMS}
    sets.update({s: sorted(set(marks.loc[marks.state.eq(s), "gene"].map(canon)) - allmp) for s in CYCLE})
    rows, cover = [], []
    for m in manifest.itertuples():
        x = io.mmread(ROOT / m.matrix_file).T.tocsr()
        g = pd.read_csv(ROOT / m.gene_file, sep="\t").gene.map(canon)
        assert g.is_unique
        gi = {v: i for i, v in enumerate(g)}
        c = pd.read_csv(ROOT / m.cell_file, sep="\t")
        depth = np.asarray(x.sum(1)).ravel()
        assert np.allclose(depth, c.full_assay_counts_recomputed)
        d = pd.DataFrame(dict(cell_id=c.cell_id, sample_key=m.sample_key, study=m.study, state=c.state,
                              logUMI=np.log1p(depth), log_detected=np.log1p(np.asarray((x > 0).sum(1)).ravel())))
        for label, members in sets.items():
            use = [z for z in members if z in gi]
            v = x[:, [gi[z] for z in use]].toarray()
            d[label] = np.log1p(v / depth[:, None] * 1e4).mean(1)
            cover.append(dict(sample_key=m.sample_key, set=label, fixed_members=len(members), measured=len(use),
                              missing=";".join(sorted(set(members) - set(use)))))
        rows.append(d)
    out = pd.concat(rows, ignore_index=True).set_index("cell_id")
    save(out, "fixed_program_scores_by_cell.tsv.gz", True)
    save(pd.DataFrame(cover), "fixed_program_score_coverage.tsv")
    return out


def association(auc, meta, prefix, tfmap):
    """Within-state partial rank correlations, donor means over states and the summary over donors."""
    rows = []
    for (sample, state), ids in meta.groupby(["sample_key", "state"]).groups.items():
        if state not in STATES:
            continue
        z = meta.loc[ids]
        n = len(z)
        a = auc.loc[ids]
        valid = a.notna().all(0).to_numpy()
        cor = np.full((a.shape[1], len(PROGRAMS)), np.nan)
        if n >= MIN_STATE_CELLS:
            cv = np.column_stack([np.ones(n), stats.rankdata(z[["logUMI", "log_detected"] + CYCLE], axis=0) / n])
            y = stats.rankdata(z[PROGRAMS], axis=0) / n
            yr = y - cv @ np.linalg.lstsq(cv, y, rcond=1e-10)[0]
            yn = np.linalg.norm(yr, axis=0)
            yr /= np.maximum(yn, 1e-12)
            av = stats.rankdata(a.loc[:, valid], axis=0) / n
            ar = av - cv @ np.linalg.lstsq(cv, av, rcond=1e-10)[0]
            an = np.linalg.norm(ar, axis=0)
            ar /= np.maximum(an, 1e-12)
            rr = ar.T @ yr
            rr[an <= 1e-10, :] = np.nan
            rr[:, yn <= 1e-10] = np.nan
            cor[valid] = rr
        for j, name in enumerate(a.columns):
            for k, p in enumerate(PROGRAMS):
                rows.append(dict(sample_key=sample, study=z.study.iloc[0], state=state, n_cells=n, regulon=name,
                                 TF=tfmap[name], program=p, rho=cor[j, k], eligible_n50=n >= MIN_STATE_CELLS))
    z = pd.DataFrame(rows)
    save(z, prefix + "_within_state_associations.tsv")
    d = z.groupby(["sample_key", "study", "regulon", "TF", "program"]).agg(
        rho=("rho", "mean"), eligible_states=("rho", "count")).reset_index()
    save(d, prefix + "_donor_associations.tsv")
    s = d.groupby(["regulon", "TF", "program"]).agg(
        mean_rho=("rho", "mean"), median_rho=("rho", "median"), min_rho=("rho", "min"), max_rho=("rho", "max"),
        n_donors=("rho", "count"), n_positive=("rho", lambda x: int(x.gt(0).sum())),
        n_negative=("rho", lambda x: int(x.lt(0).sum()))).reset_index()
    save(s, prefix + "_association_summary.tsv")
    return s


def state_bubbles(cells, x, normal, regulons, genes):
    """Mean non-MP AUCell and mean target detection for every donor x state, including states without cells."""
    gix = {g: i for i, g in enumerate(genes)}
    rows = []
    for sample, m in cells.groupby("sample_key"):
        for st in STATES:
            ii = m.index[m.state.eq(st)].to_numpy()
            n = len(ii)
            det = (x[ii] > 0).mean(0) if n else None
            for r in regulons:
                av = normal.iloc[ii][r.name] if n else None
                ji = [gix[g] for g in r.genes]
                rows.append(dict(sample_key=sample, study=m.study.iloc[0], state=st, n_cells=n, regulon=r.name,
                                 TF=r.transcription_factor, n_targets=len(ji),
                                 mean_AUCell=float(av.mean()) if n else np.nan,
                                 median_AUCell=float(av.median()) if n else np.nan,
                                 mean_detected_target_fraction=float(det[ji].mean()) if n else np.nan,
                                 eligible_n50=n >= MIN_STATE_CELLS))
    return pd.DataFrame(rows)


def case_expression(genes, canon):
    """log1p CP10k (full-assay depth) of the GRN genes measured in pre-group-eligible GSE202647 U/T/D spermatogonia."""
    cm = pd.read_csv(RESULTS / "communication/cellchat_inputs/all_cell_eligibility.tsv.gz",
                     sep="\t", low_memory=False)
    cm = cm[cm.pre_group_eligible & cm.communication_group.isin(list(CASE_STATES))].copy()
    with h5py.File(PROCESSED / "GSE202647_annotated.h5ad", "r") as f:
        obs = read_elem(f["obs"])
        names = read_elem(f["var"]).index.map(canon)
        assert names.is_unique
        col = {g: i for i, g in enumerate(names)}
        common = [g for g in genes if g in col]
        ii = obs.index.get_indexer(cm.source_cell_id)
        assert np.all(ii >= 0)
        raw = sparse_dataset(f["layers/counts"])[ii].tocsr()
        depth = np.asarray(raw.sum(1)).ravel()
        assert np.allclose(depth, cm.full_assay_total_counts)
        z = raw[:, [col[g] for g in common]].toarray().astype(np.float32)
        case = np.log1p(z * (1e4 / depth[:, None]))
    cm["state"] = cm.communication_group.map(CASE_STATES)
    return cm, common, case


def projection(full, allmp, ranks, cells, x, genes, canon, selected):
    """Fixed regulons scored on normal and GSE202647 cells over the shared genes, in the same tie order."""
    cm, common, case = case_expression(genes, canon)
    cset = set(common)
    gix = {g: i for i, g in enumerate(genes)}
    ids = ["normal::" + s for s in cells.cell_id] + ["case::" + s for s in cm.source_cell_id]
    px = pd.DataFrame(np.concatenate([x[:, [gix[g] for g in common]], case]), index=ids, columns=common)
    order = [g for g in ranks.columns if g in cset]
    pr = px[order].rank(axis=1, ascending=False, method="first", na_option="bottom").astype("uint32") - 1
    scored, cover = [], []
    for r in full:
        for drop_mp in [False, True]:
            dd = {g: w for g, w in r.gene2weight.items() if not drop_mp or g not in allmp}
            measured = {g: w for g, w in dd.items() if g in cset}
            frac = len(measured) / len(dd) if dd else 0
            ok = frac >= MIN_MEASURED_FRACTION and len(measured) >= MIN_TARGETS
            suffix = NON_MP if drop_mp else ""
            cover.append(dict(regulon=r.name + suffix, TF=r.transcription_factor, exclude_all_MP=drop_mp,
                              original_targets=len(dd), measured_targets=len(measured), measured_fraction=frac,
                              eligible=ok, missing_targets=";".join(sorted(set(dd) - cset))))
            if ok:
                scored.append(regulon_copy(r, measured, suffix))
    cover = pd.DataFrame(cover)
    save(cover, "projection_regulon_coverage.tsv")
    pa = aucell4r(pr, scored, auc_threshold=AUC_THRESHOLD, num_workers=1).reindex(index=ids, columns=cover.regulon)
    save(pa, "projection_AUCell.tsv.gz", True)
    nm = pd.DataFrame(dict(projection_cell_id=ids[:len(cells)], sample_key=cells.sample_key, state=cells.state,
                           study=cells.study, diagnosis="Normal", primary_eligible=True))
    dm = pd.DataFrame(dict(projection_cell_id=ids[len(cells):], sample_key=cm.donor.to_numpy(), state=cm.state.to_numpy(),
                           study=cm.study.to_numpy(), diagnosis=cm.diagnosis.to_numpy(),
                           primary_eligible=cm.communication_eligible.to_numpy()))
    pm = pd.concat([nm, dm], ignore_index=True).set_index("projection_cell_id")
    save(pm, "projection_cell_metadata.tsv.gz", True)
    rows = []
    groups = pm.groupby(["sample_key", "state", "diagnosis", "primary_eligible"]).groups
    for (sample, st, diagnosis, eligible), ii in groups.items():
        vals = pa.loc[ii]
        for name in pa:
            rows.append(dict(sample_key=sample, state=st, diagnosis=diagnosis, primary_eligible=eligible, n_cells=len(ii),
                             regulon=name, mean_AUCell=vals[name].mean(), median_AUCell=vals[name].median(),
                             n_scored=int(vals[name].notna().sum())))
    out = pd.DataFrame(rows)
    save(out, "projection_donor_state.tsv")
    f4 = out.merge(cover, on="regulon").query('state == "undifferentiated" and exclude_all_MP and TF in @selected')
    f4 = f4[f4.sample_key.str.startswith("GSE202647_")]
    assert f4.groupby("sample_key").n_cells.first().to_dict() == U_CELLS_PER_DONOR
    save(f4, "projection_U_AUCell.tsv")
    return len(common), len(cm)


def main():
    canon = canonical_symbol()
    genes = list(map(canon, pd.read_csv(REG / "genes.tsv", sep="\t").gene))
    assert len(set(genes)) == len(genes) == 10470
    cells = pd.read_csv(REG / "cells.tsv.gz", sep="\t")
    x = np.load(REG / "expression.npy")
    mp = pd.read_csv(RESULTS / "programs/MP_genes.tsv", sep="\t")
    allmp = set(mp.gene.map(canon))
    cand = {canon(g) for g in (RESOURCES / "candidate_143_genes.txt").read_text().split()}
    with open(REG / "full_SCENIC_regulons.pkl", "rb") as f:
        learned = pickle.load(f)
    assert all(len(r.genes) == len(set(map(canon, r.genes))) for r in learned)
    full = [regulon_copy(r, {canon(g): w for g, w in r.gene2weight.items()}) for r in learned]
    non_mp, extended, coverage, tf_of = [], [], [], {}
    for r in full:
        for label, suffix, excluded in [("full", "", set()), ("nonMP", NON_MP, allmp), ("extended", EXTENDED, allmp | cand)]:
            dd = {g: w for g, w in r.gene2weight.items() if g not in excluded}
            eligible = label == "full" or len(dd) >= MIN_TARGETS
            tf_of[r.name + suffix] = r.transcription_factor
            coverage.append(dict(regulon=r.name + suffix, base_regulon=r.name, TF=r.transcription_factor, version=label,
                                 n_targets=len(dd), eligible=eligible))
            if eligible and label != "full":
                (non_mp if label == "nonMP" else extended).append(regulon_copy(r, dd, suffix))
    save(pd.DataFrame(coverage), "all_regulon_scoring_qualification.tsv")
    ranks = create_rankings(pd.DataFrame(x, index=cells.cell_id, columns=genes), seed=SEED)
    auc = aucell4r(ranks, full + non_mp + extended, auc_threshold=AUC_THRESHOLD, num_workers=1)
    normal = auc.loc[:, [r.name for r in full + non_mp]]
    save(normal, "normal_full_and_nonMP_AUCell.tsv.gz", True)
    ext = auc.reindex(columns=[r.name + EXTENDED for r in full])
    ext.columns = [r.name for r in full]
    save(ext, "extended_AUCell.tsv.gz", True)

    meta = program_scores(mp, allmp, canon).loc[normal.index]
    nonmp_summary = association(normal[[r.name for r in non_mp]], meta, "nonMP", tf_of)
    ext_summary = association(ext, meta, "extended", {r.name: r.transcription_factor for r in full})
    selected = pd.read_csv(OUT / "structural_TF_selection.tsv", sep="\t").query("selected").TF.tolist()
    nonmp_summary["regulon"] = nonmp_summary.regulon.str.replace(NON_MP, "", regex=False)
    cmp = ext_summary.merge(nonmp_summary, on=["regulon", "TF", "program"], how="left", suffixes=("_extended", "_nonMP"))
    cmp["structurally_selected_TF"] = cmp.TF.isin(selected)
    save(cmp, "nonMP_vs_extended_associations.tsv")
    save(state_bubbles(cells, x, normal, non_mp, genes), "nonMP_AUCell_by_state.tsv")
    n_common, n_case = projection(full, allmp, ranks, cells, x, genes, canon, selected)
    print(f"AUCell: {len(full)} regulons, {len(non_mp)} non-MP, {len(extended)} extended; "
          f"projection on {n_common} shared genes, {n_case} GSE202647 cells")


if __name__ == "__main__":
    main()
