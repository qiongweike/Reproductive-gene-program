"""Residual MP3 member correlations within donor x state, before and after removing published state markers,
and the K = 4 average-linkage partitions (B1-B4) of the donor-equal consensus, used as enrichment gene sets.
Ranks of member log1p(CP10k) are residualized against ranked depth, detected genes and G1/S and G2/M scores
in states with >=100 cells; members need >=5% detection and nonconstant residuals.
Inputs: results/programs/{donor_manifest,MP_genes}.tsv, results/programs/identity/state_markers.tsv
Outputs: results/coordination/marker_removal/{known_marker_sensitivity,state_eligibility,MP3_branches}.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import ROOT, RESULTS, out_dir

import numpy as np
import pandas as pd
from scipy import io, stats
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

STATES = ['undifferentiated', 'transition', 'differentiating']
CYCLE_SETS = ['G1/S', 'G2/M']
MIN_CELLS = 100
MIN_DETECTION = 0.05
MIN_CYCLE_GENES = 5
MIN_PAIR_DONORS = 4
MIN_PAIR_STUDIES = 2
MIN_CLUSTERED_GENES = 10
N_PARTITIONS = 4


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def save(rows, path):
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False, na_rep='NA')


def unique_map(genes):
    """Column index of each upper-case symbol; symbols occurring twice are left out."""
    s = pd.Series([str(x).upper() for x in genes])
    dup = s.duplicated(False)
    return {g: i for i, g in enumerate(s) if not dup.iloc[i]}


def mean_available(a):
    n = np.isfinite(a).sum(axis=0)
    return np.divide(np.nansum(a, axis=0), n, out=np.full(a.shape[1:], np.nan), where=n > 0)


def donor_correlations(r, target, cycle, markers, strata, rows):
    """Member x member residual correlations per eligible state, averaged over states (missing pairs stay NaN)."""
    raw = io.mmread(resolve(r.matrix_file)).T.tocsr()
    gi = unique_map(pd.read_csv(resolve(r.gene_file), sep='\t').gene.astype(str))
    meta = pd.read_csv(resolve(r.cell_file), sep='\t')
    lib = np.asarray(raw.sum(1)).ravel()
    ng = np.asarray((raw > 0).sum(1)).ravel()
    assert np.allclose(lib, meta.full_assay_counts_recomputed)
    covg = {k: [g for g in v if g in gi] for k, v in cycle.items()}
    assert min(map(len, covg.values())) >= MIN_CYCLE_GENES
    sel = sorted(set(target).union(*map(set, covg.values())) & set(gi))
    mat = np.log1p(raw[:, [gi[g] for g in sel]].multiply(1e4 / lib[:, None]).toarray())
    di = {g: i for i, g in enumerate(sel)}
    per_state = []
    for state in STATES:
        ix = np.flatnonzero(meta.state.eq(state))
        n = len(ix)
        strata.append(dict(sample_key=r.sample_key, study=r.study, state=state, n_cells=n, eligible=n >= MIN_CELLS))
        if n < MIN_CELLS:
            continue
        cov = np.column_stack([np.log1p(lib[ix]), np.log1p(ng[ix])] +
                              [mat[ix][:, [di[g] for g in gs]].mean(1) for gs in covg.values()])
        cov = np.column_stack([np.ones(n), stats.rankdata(cov, axis=0) / n])
        y = np.column_stack([mat[ix, di[g]] if g in di else np.full(n, np.nan) for g in target])
        ids = np.flatnonzero(np.array([g in di for g in target]) & (np.mean(y > 0, axis=0) >= MIN_DETECTION))
        rk = stats.rankdata(y[:, ids], axis=0) / n
        rr = rk - cov @ np.linalg.lstsq(cov, rk, rcond=1e-10)[0]
        norm = np.linalg.norm(rr, axis=0)
        ok = norm > 1e-8
        ids, rr = ids[ok], rr[:, ok] / norm[ok]
        corr = np.full((len(target), len(target)), np.nan)
        corr[np.ix_(ids, ids)] = rr.T @ rr
        per_state.append(corr)
        for label, sub in [('fixed_MP3', list(ids)),
                           ('remove_published_state_markers', [i for i in ids if target[i] not in markers])]:
            z = corr[np.ix_(sub, sub)]
            rows.append(dict(sample_key=r.sample_key, study=r.study, state=state, sensitivity=label, n_genes=len(sub),
                             mean_rho=np.mean(z[np.triu_indices(len(sub), 1)]) if len(sub) > 1 else np.nan))
    return mean_available(np.array(per_state))


def partitions(arr, manifest, target, markers):
    """Average-linkage K = 4 partition of the consensus over member pairs covered by >=4 donors and >=2 studies."""
    donors = np.isfinite(arr).sum(0)
    studies = np.zeros(donors.shape, int)
    for ii in manifest.groupby('study').indices.values():
        studies += np.any(np.isfinite(arr[list(ii)]), axis=0)
    valid = (donors >= MIN_PAIR_DONORS) & (studies >= MIN_PAIR_STUDIES)
    ids = list(range(len(target)))
    while not np.all(valid[np.ix_(ids, ids)]):
        ids.pop(int(np.argmin(valid[np.ix_(ids, ids)].sum(0))))
        if len(ids) < MIN_CLUSTERED_GENES:
            raise RuntimeError('Fewer than ten members have adequate pair coverage')
    a = np.clip(mean_available(arr)[np.ix_(ids, ids)], -1, 1)
    np.fill_diagonal(a, 1)
    tree = linkage(squareform(np.maximum(0, 1 - a), checks=False), method='average')
    cluster, leaves = fcluster(tree, N_PARTITIONS, criterion='maxclust'), leaves_list(tree)
    rows = []
    for i, gene in enumerate(target):
        j = ids.index(i) if i in ids else None
        rows.append(dict(gene=gene, original_rank=i + 1,
                         branch=f'B{cluster[j]}' if j is not None else 'unassigned_coverage',
                         published_state_marker=gene in markers,
                         leaf_order=int(np.where(leaves == j)[0][0]) if j is not None else np.nan))
    return rows


def main():
    out = out_dir('coordination', 'marker_removal')
    manifest = pd.read_csv(RESULTS / 'programs' / 'donor_manifest.tsv', sep='\t').query('eligible')
    mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t')
    target = mp.loc[mp.MP.eq('MP_3')].sort_values('rank').gene.str.upper().tolist()
    prior = pd.read_csv(RESULTS / 'programs' / 'identity' / 'state_markers.tsv', sep='\t')
    cycle = {s: sorted(set(prior.loc[prior.state.eq(s), 'gene'].str.upper()) - set(mp.gene.str.upper()))
             for s in CYCLE_SETS}
    markers = set(prior.loc[~prior.state.isin(CYCLE_SETS), 'gene'].str.upper())
    strata, rows = [], []
    arr = np.array([donor_correlations(r, target, cycle, markers, strata, rows) for r in manifest.itertuples()])
    branches = partitions(arr, manifest, target, markers)
    save(rows, out / 'known_marker_sensitivity.tsv')
    save(strata, out / 'state_eligibility.tsv')
    save(branches, out / 'MP3_branches.tsv')
    print(len(rows) // 2, 'donor-state strata;', pd.DataFrame(branches).branch.value_counts().sort_index().to_dict())


if __name__ == '__main__':
    main()
