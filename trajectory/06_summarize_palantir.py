"""Palantir qualification and MP3-minus-MP2 median member peak bins for Palantir and the Slingshot lineages (Fig. 3f).
Whole-donor Palantir pseudotime is qualified with the Slingshot branch rules (no curve weight); each Slingshot lineage is
summarized on its cells (curve weight >=0.5) that also have a Palantir pseudotime. Ten quantile bins; peaks as in step 04.
Inputs: results/trajectory/{donor_manifest.tsv, fixed_cell_scores.tsv.gz, lineage_cell_pseudotime.tsv.gz,
        palantir/donor_start_eligibility.tsv, palantir/<donor>.cell_results.tsv.gz}, results/programs/MP_genes.tsv, count matrices.
Outputs: results/trajectory/palantir/{palantir_qualification.tsv, MP2_MP3_member_peaks.tsv, MP2_MP3_peak_order.tsv}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir

import numpy as np
import pandas as pd
from scipy import io, stats

MIN_CURVE_WEIGHT = 0.5
MIN_CELLS = 100
MIN_STATE_CELLS = 30
N_BINS = 10
MIN_BIN_CELLS = 10
CONSTANT_RANGE = 1e-12
PROGRAMS = RESULTS / 'programs'
TRAJ = RESULTS / 'trajectory'
OUT = out_dir('trajectory', 'palantir')


def donor_files(sample):
    return (PROGRAMS / 'nmf/input' / f'{sample}.counts.mtx.gz', PROGRAMS / 'input' / f'{sample}.genes.tsv',
            PROGRAMS / 'nmf/input' / f'{sample}.cells.tsv')


def save(rows, name):
    pd.DataFrame(rows).to_csv(OUT / name, sep='\t', index=False, na_rep='NA')


def rho(x, y):
    return float(stats.spearmanr(x, y).statistic) if len(x) > 2 and np.nanstd(x) > 0 and np.nanstd(y) > 0 else np.nan


def qualify(z, root_ok):
    n_u = int(z.state.eq('undifferentiated').sum())
    n_d = int(z.state.eq('differentiating').sum())
    e, d = rho(z.pseudotime, z.early_marker), rho(z.pseudotime, z.late_marker)
    u = z.loc[z.state.eq('undifferentiated'), 'pseudotime'].median()
    l = z.loc[z.state.eq('differentiating'), 'pseudotime'].median()
    reason = []
    if len(z) < MIN_CELLS: reason.append('fewer_than_100_cells')
    if n_u < MIN_STATE_CELLS: reason.append('fewer_than_30_U')
    if n_d < MIN_STATE_CELLS: reason.append('fewer_than_30_D')
    if not root_ok: reason.append('root_marker_inconsistent')
    if not l > u: reason.append('D_not_after_U')
    if not (e < 0 and d > 0): reason.append('marker_gradient_inconsistent')
    return dict(n_cells=len(z), n_U=n_u, n_D=n_d, rho_early_marker=e, rho_late_marker=d, median_U=u, median_D=l,
                eligible=not reason, reason=';'.join(reason) if reason else 'PASS')


def member_expression(sample, genes):
    """log1p(CP10k) per cell for uniquely mapped members, with a cell_id -> row index map."""
    counts_file, gene_file, cell_file = donor_files(sample)
    raw = io.mmread(counts_file).T.tocsr().astype(float)
    upper = pd.read_csv(gene_file, sep='\t').gene.astype(str).str.upper()
    dup = upper.duplicated(False).to_numpy()
    gi = {g: i for i, g in enumerate(upper) if not dup[i]}
    cells = pd.read_csv(cell_file, sep='\t')
    lib = np.asarray(raw.sum(1)).ravel()
    assert np.allclose(lib, cells.full_assay_counts_recomputed, rtol=1e-10)
    expr = {g: np.log1p(raw[:, gi[g]].toarray().ravel() * 1e4 / lib) for g in genes if g in gi}
    return {c: i for i, c in enumerate(cells.cell_id)}, expr


def context_peaks(z, meta, mp, ci, expr):
    z = z.copy()
    z['bin'] = pd.qcut(z.pseudotime.rank(method='first'), N_BINS, labels=False) + 1
    z['time_scaled'] = (z.pseudotime - z.pseudotime.min()) / (z.pseudotime.max() - z.pseudotime.min())
    ix = np.array([ci[c] for c in z.cell_id.to_numpy()])
    groups = [(int(b), np.flatnonzero(z.bin.to_numpy() == b), float(a.time_scaled.mean())) for b, a in z.groupby('bin', sort=True)]
    peaks = []
    for row in mp.itertuples():
        v = expr.get(row.gene)
        means = [float(v[ix[bix]].mean()) if v is not None else np.nan for _, bix, _ in groups]
        if all(np.isfinite(x) and len(bix) >= MIN_BIN_CELLS for x, (_, bix, _) in zip(means, groups)):
            top, low = max(means), min(means)
            k = means.index(top)
            peaks.append(dict(**meta, MP=row.MP, gene=row.gene, peak_bin=groups[k][0], peak_time=groups[k][2],
                              n_tied_peak=means.count(top), dynamic_range=top - low, nonconstant=top - low > CONSTANT_RANGE,
                              reason='PASS' if top - low > CONSTANT_RANGE else 'CONSTANT_MEMBER'))
        else:
            peaks.append(dict(**meta, MP=row.MP, gene=row.gene, peak_bin=np.nan, peak_time=np.nan, n_tied_peak=0,
                              dynamic_range=np.nan, nonconstant=False, reason='MISSING_OR_INSUFFICIENT_BIN'))
    return peaks


def main():
    manifest = pd.read_csv(TRAJ / 'donor_manifest.tsv', sep='\t')
    starts = pd.read_csv(OUT / 'donor_start_eligibility.tsv', sep='\t').set_index('sample_key')
    scores = pd.read_csv(TRAJ / 'fixed_cell_scores.tsv.gz', sep='\t').set_index('cell_id')
    sl = pd.read_csv(TRAJ / 'lineage_cell_pseudotime.tsv.gz', sep='\t')
    mp = pd.read_csv(PROGRAMS / 'MP_genes.tsv', sep='\t')
    mp = mp[mp.MP.isin(['MP_2', 'MP_3'])].copy()
    mp.gene = mp.gene.str.upper()

    quals, peaks = [], []
    for m in manifest.itertuples():
        st = starts.loc[m.sample_key]
        if not st.run_eligible:
            quals.append(dict(sample_key=m.sample_key, status='ROOT_INELIGIBLE', eligible=False, reason=st.reason, n_input_cells=m.n_spg))
            continue
        z = pd.read_csv(OUT / f'{m.sample_key}.cell_results.tsv.gz', sep='\t').set_index('cell_id', drop=False)
        zz = z[np.isfinite(z.pseudotime)].copy()
        q = qualify(zz, st.root_marker_ok)
        quals.append(dict(sample_key=m.sample_key, status='FIT_RETURNED', n_input_cells=m.n_spg, **q))
        contexts = [('global_all_cells', 'Palantir', 'global', zz, q)]
        for lineage, t in sl[sl.sample_key.eq(m.sample_key)].groupby('lineage', sort=False):
            a = t[t.pseudotime.notna() & t.curve_weight.ge(MIN_CURVE_WEIGHT)]
            a = a[a.cell_id.isin(zz.index)]
            ss = scores.loc[a.cell_id, ['state', 'early_marker', 'late_marker']].copy()
            ss['cell_id'] = a.cell_id.to_numpy()
            ss['pseudotime'] = a.pseudotime.to_numpy()
            contexts.append(('paired_primary_shared_cells', 'Slingshot', lineage, ss, qualify(ss, st.root_marker_ok)))
        ci, expr = member_expression(m.sample_key, mp.gene)
        for scope, method, lineage, cz, cq in contexts:
            meta = dict(sample_key=m.sample_key, scope=scope, method=method, lineage=lineage,
                        trajectory_eligible=cq['eligible'], trajectory_reason=cq['reason'])
            peaks += context_peaks(cz, meta, mp, ci, expr)
    save(quals, 'palantir_qualification.tsv')
    pk = pd.DataFrame(peaks)
    save(pk, 'MP2_MP3_member_peaks.tsv')

    order = []
    keys = ['sample_key', 'scope', 'method', 'lineage']
    for k, a in pk.groupby(keys, sort=False):
        p2, p3 = (float(a.loc[a.MP.eq(x) & a.nonconstant, 'peak_bin'].median()) for x in ['MP_2', 'MP_3'])
        order.append(dict(zip(keys, k), MP2_median_peak=p2, MP3_median_peak=p3, MP3_minus_MP2_peak=p3 - p2,
                          MP2_precedes_MP3=p3 > p2, trajectory_eligible=bool(a.trajectory_eligible.iloc[0]),
                          trajectory_reason=a.trajectory_reason.iloc[0]))
    order = pd.DataFrame(order)
    save(order, 'MP2_MP3_peak_order.tsv')
    g = order[order.method.eq('Palantir') & order.trajectory_eligible]
    print('Palantir MP3-MP2 median member peak:', ', '.join(f'{s} {v:g}' for s, v in zip(g.sample_key, g.MP3_minus_MP2_peak)))


if __name__ == '__main__':
    main()
