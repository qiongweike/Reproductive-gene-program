"""Donor-level Hotspot pair summaries (Fig. 2d,e) and their comparison with within-state residual correlations.
Pair C is averaged over a donor's available contexts per scope (whole donor or states), then over donors; pairs need >=4
donors and >=2 studies. Residual reference: within-state member Spearman correlation after regressing ranks on log depth,
log detected genes and non-MP G1/S and G2/M scores (members detected in >=5% of cells), on pairs detected in >=5% of cells.
Inputs: results/trajectory/{donor_manifest.tsv, hotspot/*}, results/programs/{MP_genes.tsv, identity/state_markers.tsv}, count matrices.
Outputs: results/trajectory/hotspot/{donor_residual_comparison.tsv, MP3_pair_consensus.tsv}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir
import itertools

import numpy as np
import pandas as pd
from scipy import io, stats

MIN_DETECTION = 0.05
MIN_CYCLE_GENES = 5
MIN_DONORS, MIN_STUDIES = 4, 2
CYCLE_STATES = ['G1/S', 'G2/M']
PROGRAMS = RESULTS / 'programs'
TRAJ = RESULTS / 'trajectory'
OUT = out_dir('trajectory', 'hotspot')


def donor_files(sample):
    return (PROGRAMS / 'nmf/input' / f'{sample}.counts.mtx.gz', PROGRAMS / 'input' / f'{sample}.genes.tsv',
            PROGRAMS / 'nmf/input' / f'{sample}.cells.tsv')


def mean_available(a):
    a = np.asarray(a)
    n = np.isfinite(a).sum(0)
    return np.divide(np.nansum(a, axis=0), n, out=np.full(a.shape[1:], np.nan), where=n > 0)


def comparison(x, y):
    x, y = np.asarray(x), np.asarray(y)
    ok = np.isfinite(x) & np.isfinite(y)
    a, b = x[ok], y[ok]
    return dict(n_common_pairs=int(ok.sum()),
                pair_rank_spearman=float(stats.spearmanr(a, b).statistic) if len(a) > 2 and np.ptp(a) > 0 and np.ptp(b) > 0 else np.nan,
                sign_agreement=float(np.mean(np.sign(a) == np.sign(b))) if len(a) else np.nan,
                Hotspot_mean_C=float(a.mean()) if len(a) else np.nan, residual_mean_rho=float(b.mean()) if len(b) else np.nan,
                Hotspot_positive_fraction=float(np.mean(a > 0)) if len(a) else np.nan,
                residual_positive_fraction=float(np.mean(b > 0)) if len(a) else np.nan)


def residual_correlations(sample, contexts, target, cycle):
    """Within-state residual Spearman correlations of the members, as vectors over the 1,176 pairs."""
    counts_file, gene_file, cell_file = donor_files(sample)
    raw = io.mmread(counts_file).T.tocsr()
    genes = pd.read_csv(gene_file, sep='\t').gene.astype(str).str.upper()
    cells = pd.read_csv(cell_file, sep='\t')
    dup = genes.duplicated(False)
    gi = {g: i for i, g in enumerate(genes) if not dup.iloc[i]}
    lib = np.asarray(raw.sum(1)).ravel()
    ng = np.asarray((raw > 0).sum(1)).ravel()
    cov_genes = {k: [g for g in v if g in gi] for k, v in cycle.items()}
    assert min(map(len, cov_genes.values())) >= MIN_CYCLE_GENES
    sel = sorted(set(target).union(*map(set, cov_genes.values())) & set(gi))
    mat = np.log1p(raw[:, [gi[g] for g in sel]].multiply(1e4 / lib[:, None]).toarray())
    di = {g: i for i, g in enumerate(sel)}
    pa, pb = map(np.array, zip(*itertools.combinations(range(len(target)), 2)))
    out = {}
    for ctx in contexts.itertuples():
        corr = np.full((len(target), len(target)), np.nan)
        if ctx.scope == 'within_state' and ctx.state_eligible:
            ix = np.flatnonzero(cells.state.eq(ctx.state))
            n = len(ix)
            cov = np.column_stack([np.log1p(lib[ix]), np.log1p(ng[ix])] + [mat[ix][:, [di[g] for g in gs]].mean(1) for gs in cov_genes.values()])
            cov = np.column_stack([np.ones(n), stats.rankdata(cov, axis=0) / n])
            y = np.column_stack([mat[ix, di[g]] if g in di else np.full(n, np.nan) for g in target])
            ids = np.flatnonzero(np.array([g in di for g in target]) & (np.mean(y > 0, axis=0) >= MIN_DETECTION))
            rk = stats.rankdata(y[:, ids], axis=0) / n
            rr = rk - cov @ np.linalg.lstsq(cov, rk, rcond=1e-10)[0]
            norm = np.linalg.norm(rr, axis=0)
            ok = norm > 1e-8
            ids, rr = ids[ok], rr[:, ok] / norm[ok]
            corr[np.ix_(ids, ids)] = rr.T @ rr
        out[ctx.context] = corr[pa, pb]
    return out


def main():
    manifest = pd.read_csv(TRAJ / 'donor_manifest.tsv', sep='\t')
    mp = pd.read_csv(PROGRAMS / 'MP_genes.tsv', sep='\t')
    target = mp[mp.MP.eq('MP_3')].sort_values('rank').gene.str.upper().tolist()
    all_mp = set(mp.gene.str.upper())
    markers = pd.read_csv(PROGRAMS / 'identity/state_markers.tsv', sep='\t')
    cycle = {s: sorted(set(markers.loc[markers.state.eq(s), 'gene'].str.upper()) - all_mp) for s in CYCLE_STATES}
    names = pd.DataFrame(list(itertools.combinations(target, 2)), columns=['gene1', 'gene2'])
    contexts = pd.read_csv(OUT / 'context_eligibility.tsv', sep='\t')
    pairs = pd.read_csv(OUT / 'MP3_pair_local_correlations.tsv.gz', sep='\t')
    split = {k: p.set_index(['gene1', 'gene2']).reindex(pd.MultiIndex.from_frame(names))
             for k, p in pairs.groupby(['context', 'model'], sort=False)}
    assert all(len(p) == len(names) for p in split.values())

    donor_rows, comparisons = [], []
    for m in manifest.itertuples():
        ctx = contexts[contexts.sample_key.eq(m.sample_key)]
        residual = residual_correlations(m.sample_key, ctx, target, cycle)
        for model in ['bernoulli', 'danb']:
            for scope in ['whole_donor', 'within_state']:
                pv, av, bv = [], [], []
                for c in ctx[ctx.scope.eq(scope)].context:
                    p = split[(c, model)]
                    v, r = p.C.to_numpy(), residual[c]
                    ok = np.isfinite(v) & np.isfinite(r) & p.both_detection_ge05.to_numpy()
                    pv.append(v); av.append(np.where(ok, v, np.nan)); bv.append(np.where(ok, r, np.nan))
                v, a, b = mean_available(pv), mean_available(av), mean_available(bv)
                donor_rows.append(names.assign(mean_Hotspot_C=v, mean_Hotspot_C_matched_residual_contexts=a,
                                               mean_residual_rho_matched_contexts=b, study=m.study, scope=scope, model=model))
                finite = np.isfinite(v)
                comparisons.append(dict(sample_key=m.sample_key, study=m.study, donor=m.donor, scope=scope, model=model,
                                        n_pairs_available=int(finite.sum()),
                                        mean_Hotspot_C_all_available=float(np.nanmean(v)) if finite.any() else np.nan,
                                        positive_fraction_all_available=float(np.mean(v[finite] > 0)) if finite.any() else np.nan,
                                        **comparison(a, b)))
    pd.DataFrame(comparisons).to_csv(OUT / 'donor_residual_comparison.tsv', sep='\t', index=False, na_rep='NA')

    donor = pd.concat(donor_rows, ignore_index=True)
    consensus = []
    for (scope, model), d in donor.groupby(['scope', 'model'], sort=False):
        for (g1, g2), p in d.groupby(['gene1', 'gene2'], sort=False):
            ok = p.mean_Hotspot_C.notna()
            n_studies = p.loc[ok, 'study'].nunique()
            consensus.append(dict(scope=scope, model=model, gene1=g1, gene2=g2, n_donors=int(ok.sum()), n_studies=n_studies,
                                  mean_C_equal_donor=p.mean_Hotspot_C.mean(), n_positive_donors=int((p.mean_Hotspot_C > 0).sum()),
                                  minimum_donor_C=p.mean_Hotspot_C.min(), maximum_donor_C=p.mean_Hotspot_C.max(),
                                  n_joint_donors=int(p.mean_Hotspot_C_matched_residual_contexts.notna().sum()),
                                  mean_Hotspot_C_joint_equal_donor=p.mean_Hotspot_C_matched_residual_contexts.mean(),
                                  mean_residual_rho_joint_equal_donor=p.mean_residual_rho_matched_contexts.mean(),
                                  covered_ge4_donors_ge2_studies=bool(ok.sum() >= MIN_DONORS and n_studies >= MIN_STUDIES)))
    consensus = pd.DataFrame(consensus)
    consensus.to_csv(OUT / 'MP3_pair_consensus.tsv', sep='\t', index=False, na_rep='NA')
    c = pd.DataFrame(comparisons).groupby(['model', 'scope']).mean_Hotspot_C_all_available.mean()
    print('Donor-mean local pair C:', '; '.join(f'{k[0]} {k[1]} {v:.3f}' for k, v in c.items()))


if __name__ == '__main__':
    main()
