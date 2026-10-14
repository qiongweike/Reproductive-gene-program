"""Residual program coordination with added late germ-cell (spermatid) and somatic marker covariates and with published
state markers removed; base model as in 03 (states with >=100 cells, averaged within donor); donor sign-flip tests.
Inputs: results/programs/{donor_manifest,MP_genes}.tsv, results/programs/identity/state_markers.tsv,
        results/atlas/annotation_prespecified_rules.json
Outputs: results/coordination/specificity/background_{coexpression_by_donor_state,coexpression_by_donor,
         donor_differences,donor_effects}.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import ROOT, RESULTS, out_dir
import itertools
import json

import numpy as np
import pandas as pd
from scipy import io, stats
from statsmodels.stats.multitest import multipletests

STATES = ['undifferentiated', 'transition', 'differentiating']
CYCLE_SETS = ['G1/S', 'G2/M']
MIN_CELLS = 100
MIN_DETECTION = 0.05
LATE_PANEL = 'Spermatid'
SOMATIC_PANELS = ['Sertoli', 'Leydig', 'Peritubular_smooth_muscle', 'Endothelial', 'Immune', 'Erythroid']
MODELS = ['base', 'plus_late', 'plus_somatic', 'plus_late_somatic',
          'remove_published_state_markers', 'remove_published_plus_late_somatic']
RESIDUAL_MODEL = {'remove_published_state_markers': 'base', 'remove_published_plus_late_somatic': 'plus_late_somatic'}


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def save(rows, name):
    pd.DataFrame(rows).to_csv(out_dir('coordination', 'specificity') / name, sep='\t', index=False, na_rep='NA')


def two_sided(d):
    """Exact two-sided donor sign-flip test of the mean difference, with a descriptive t interval."""
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    n = len(d)
    if n < 2:
        return dict(n_donors=n, mean_difference=float(d.mean()) if n else np.nan,
                    median_difference=float(np.median(d)) if n else np.nan, ci95_low=np.nan, ci95_high=np.nan,
                    p_two_sided=np.nan, status='insufficient_donors')
    null = np.asarray(list(itertools.product([-1., 1.], repeat=n))) @ d / n
    p = float(np.mean(np.abs(null) >= abs(d.mean()) - 1e-12))
    half = stats.t.ppf(.975, n - 1) * stats.sem(d)
    return dict(n_donors=n, mean_difference=float(d.mean()), median_difference=float(np.median(d)),
                ci95_low=float(d.mean() - half), ci95_high=float(d.mean() + half), p_two_sided=p,
                status='estimated_descriptive', n_exact_sign_flips=2 ** n)


def bh(rows, group):
    d = pd.DataFrame(rows)
    d['q_BH'] = np.nan
    for _, ix in d.groupby(group, observed=True, dropna=False).groups.items():
        ids = d.loc[ix].index[d.loc[ix, 'p_two_sided'].notna()]
        d.loc[ids, 'q_BH'] = multipletests(d.loc[ids, 'p_two_sided'], method='fdr_bh')[1]
    return d


def donor_coexpression(row, programs, panels, cycle, known):
    raw = io.mmread(resolve(row.matrix_file)).T.tocsr()
    genes = pd.read_csv(resolve(row.gene_file), sep='\t').gene.astype(str).str.upper()
    assert genes.is_unique
    meta = pd.read_csv(resolve(row.cell_file), sep='\t')
    depth = np.asarray(raw.sum(1)).ravel()
    ng = np.asarray((raw > 0).sum(1)).ravel()
    assert np.allclose(depth, meta.full_assay_counts_recomputed)
    gi = {g: i for i, g in enumerate(genes)}
    wanted = sorted(set().union(*programs.values(), *panels.values(), *cycle.values()) & set(genes))
    mi = {g: i for i, g in enumerate(wanted)}
    mat = np.log1p(raw[:, [gi[g] for g in wanted]].multiply(1e4 / depth[:, None]).toarray())
    rows = []
    for state in STATES:
        ix = np.flatnonzero(meta.state.eq(state))
        n = len(ix)
        base = dict(sample_key=row.sample_key, study=row.study, donor=row.donor, state=state)
        if n < MIN_CELLS:
            rows += [dict(base, MP=mp, model=m, n_cells=n, n_genes=0, mean_rho=np.nan, status='fewer_than_100_cells')
                     for mp in programs for m in MODELS]
            continue
        cov = np.column_stack([np.log1p(depth[ix]), np.log1p(ng[ix])] +
                              [mat[ix][:, [mi[g] for g in gs if g in mi]].mean(1) for gs in cycle.values()])
        covs = {'base': np.column_stack([np.ones(n), stats.rankdata(cov, axis=0) / n])}
        score = {k: stats.rankdata(mat[ix][:, [mi[g] for g in gs if g in mi]].mean(1)) / n for k, gs in panels.items()}
        late = score[LATE_PANEL][:, None]
        somatic = np.column_stack([score[k] for k in SOMATIC_PANELS])
        covs.update(plus_late=np.column_stack([covs['base'], late]),
                    plus_somatic=np.column_stack([covs['base'], somatic]),
                    plus_late_somatic=np.column_stack([covs['base'], late, somatic]))
        for mp, gs in programs.items():
            ug = [g for g in gs if g in mi and np.mean(mat[ix, mi[g]] > 0) >= MIN_DETECTION]
            y = stats.rankdata(mat[ix][:, [mi[g] for g in ug]], axis=0) / n
            res = {k: y - x @ np.linalg.lstsq(x, y, rcond=1e-10)[0] for k, x in covs.items()}
            ok = np.logical_and.reduce([np.linalg.norm(v, axis=0) > 1e-8 for v in res.values()])
            for model in MODELS:
                src = RESIDUAL_MODEL.get(model, model)
                select = ok.copy()
                if model.startswith('remove_published'):
                    select &= np.array([g not in known for g in ug])
                rr = res[src][:, select]
                norm = np.linalg.norm(rr, axis=0)
                rr = rr / norm if len(norm) else rr
                k = rr.shape[1]
                rho = float((rr.T @ rr)[np.triu_indices(k, 1)].mean()) if k > 1 else np.nan
                rows.append(dict(base, MP=mp, model=model, n_cells=n, n_genes=k, mean_rho=rho,
                                 status='estimated_descriptive' if k > 1 else 'fewer_than_2_residual_genes',
                                 covariate_rank=int(np.linalg.matrix_rank(covs[src])),
                                 n_covariate_columns=covs[src].shape[1]))
    return rows


def main():
    manifest = pd.read_csv(RESULTS / 'programs' / 'donor_manifest.tsv', sep='\t').query('eligible')
    mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t')
    programs = {k: g.sort_values('rank').gene.str.upper().tolist() for k, g in mp.groupby('MP', observed=True)}
    union = set().union(*map(set, programs.values()))
    rules = json.loads((RESULTS / 'atlas' / 'annotation_prespecified_rules.json').read_text())
    panels = {t: sorted({g.upper() for g in rules['assignment_markers'][t] + rules['heldout_markers'][t]} - union)
              for t in [LATE_PANEL] + SOMATIC_PANELS}
    prior = pd.read_csv(RESULTS / 'programs' / 'identity' / 'state_markers.tsv', sep='\t')
    cycle = {s: sorted(set(prior.loc[prior.state.eq(s), 'gene'].str.upper()) - union) for s in CYCLE_SETS}
    known = set(prior.loc[~prior.state.isin(CYCLE_SETS), 'gene'].str.upper())
    coh = pd.DataFrame([r for row in manifest.itertuples()
                        for r in donor_coexpression(row, programs, panels, cycle, known)])
    don = coh[coh.mean_rho.notna()].groupby(['sample_key', 'study', 'donor', 'MP', 'model'], observed=True).agg(
        mean_rho=('mean_rho', 'mean'), n_states=('state', 'size'), states=('state', lambda x: ';'.join(sorted(x))),
        min_genes=('n_genes', 'min')).reset_index()
    diffs, tests = [], []
    for p in programs:
        base = don[(don.MP == p) & don.model.eq('base')].set_index('sample_key')
        for model in MODELS[1:]:
            alt = don[(don.MP == p) & don.model.eq(model)].set_index('sample_key')
            shared = base.index.intersection(alt.index)
            assert base.loc[shared, 'states'].eq(alt.loc[shared, 'states']).all()
            dd = alt.loc[shared, 'mean_rho'] - base.loc[shared, 'mean_rho']
            tests.append(dict(MP=p, model=model, contrast='sensitivity minus base',
                              BH_family='all_15_background_changes', **two_sided(dd)))
            diffs += [dict(sample_key=s, MP=p, model=model, base_rho=base.loc[s, 'mean_rho'],
                           sensitivity_rho=alt.loc[s, 'mean_rho'], difference=dd[s], n_states=base.loc[s, 'n_states'])
                      for s in shared]
    save(coh, 'background_coexpression_by_donor_state.tsv')
    save(don, 'background_coexpression_by_donor.tsv')
    save(diffs, 'background_donor_differences.tsv')
    save(bh(tests, ['BH_family']), 'background_donor_effects.tsv')
    full = don[don.MP.eq('MP_3') & don.model.eq('remove_published_plus_late_somatic')].mean_rho
    print(f'MP3 without state markers, plus late and somatic covariates: donor-equal mean {full.mean():.4f}, '
          f'{int((full > 0).sum())}/{len(full)} donors > 0')


if __name__ == '__main__':
    main()
