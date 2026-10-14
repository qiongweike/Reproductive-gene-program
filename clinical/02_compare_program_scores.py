"""Fixed-background program scores in donor x state groups and within-study case-control contrasts
(Fig 4b, 4d; Supplementary Table 3).
Inputs: normal, cryptozoospermia and GSE202647/GSE326372 spermatogonia (data/processed), results/programs/MP_genes.tsv
and donor_manifest.tsv, the shared measured-gene universe (assay_genes.py), results/atlas/cell_annotations.tsv.gz.
Outputs: results/clinical/coverage.tsv, donor_state_scores.tsv, all_effects.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, RESOURCES, out_dir
import json

import anndata as ad
import numpy as np
import pandas as pd

from assay_genes import canonical, log_cp10k, measured_genes, symbol_mappings
from donor_contrasts import bh, exact_welch

SEED = 20261003
N_BINS = 20                # expression bins over equally weighted discovery-donor means
N_CONTROLS = 30            # background genes drawn per member, without replacement, then pooled
MIN_SHARED_MEMBERS = 25
STATES = ['undifferentiated', 'transition', 'differentiating']
MIN_CELLS = [30, 10]       # primary and sensitivity donor-state thresholds
CASE_FILES = ['crypto_spg.h5ad', 'disease_spg.h5ad']
CONTRASTS = [('GSE202647', 'NOA', 'OA'), ('GSE153947', 'Crypto', 'OA')]
PRIMARY = 'primary'                   # every scored spermatogonium
CONCORDANT = 'atlas_SPG_agreement'   # GSE202647 cells also classified as SPG in the atlas
CONCORDANT_STUDY = 'GSE202647'
FAMILIES = ['study', 'scenario', 'min_cells']   # BH families: study x cell selection x threshold
UNITS = 'log1p_CP10k_full_assay_denominator'
KEYS = ['study', 'donor', 'diagnosis', 'state', 'score_units']


def load(fn, mappings):
    """Cell metadata, canonical genes and ln(1 + CP10k) expression of one spermatogonial assay."""
    a = ad.read_h5ad(PROCESSED / fn)
    studies = a.obs.study.astype(str).unique()
    alias = mappings.get(studies[0], {}) if len(studies) == 1 else {}
    names = [canonical(alias.get(g, g)) for g in map(str, a.var_names)]
    genes, x = log_cp10k(a.layers['counts'], names, pd.to_numeric(a.obs.total_counts).to_numpy(float))
    return a.obs.copy(), genes, x


def gene_columns(genes, wanted):
    j = genes.get_indexer(wanted)
    assert (j >= 0).all(), 'scoring gene missing from assay'
    return j


mappings = symbol_mappings()
mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t').sort_values(['MP', 'rank'])
programs = {k: list(dict.fromkeys(map(canonical, v.gene))) for k, v in mp.groupby('MP')}
common = sorted(set.intersection(*measured_genes().values()))
common_set, common_arr = set(common), np.array(common)
candidates = (RESOURCES / 'candidate_143_genes.txt').read_text().split()
modules = json.loads((RESOURCES / 'excluded_module_genes.json').read_text())
# Program members, the 143 candidate genes and the predefined modules never enter a background pool.
excluded = {canonical(g) for g in [*sum(programs.values(), []), *candidates, *sum(modules.values(), [])]}

# Background pools: expression bins and detection from the seven discovery donors, donors weighted equally.
manifest = pd.read_csv(RESULTS / 'programs' / 'donor_manifest.tsv', sep='\t').query('eligible')
means, detected = [], []
for study in manifest.study.unique():
    obs, genes, x = load(f'normal_spg_{study}.h5ad', mappings)
    j = gene_columns(genes, common)
    for r in manifest[manifest.study.eq(study)].itertuples():
        use = obs.donor.astype(str).eq(r.donor).to_numpy()
        assert use.sum() == r.n_spg, r.sample_key
        z = x[use][:, j]
        means.append(np.asarray(z.mean(0)).ravel())
        detected.append(np.asarray((z > 0).mean(0)).ravel())
bins = pd.qcut(pd.Series(np.mean(means, axis=0)).rank(method='first'), N_BINS, labels=False).to_numpy()
poolable = (np.mean(detected, axis=0) > 0) & ~np.isin(common_arr, list(excluded))

rng = np.random.default_rng(SEED)
members, background, coverage = {}, {}, []
for k, genes_k in programs.items():
    members[k] = [g for g in genes_k if g in common_set]
    assert len(members[k]) >= MIN_SHARED_MEMBERS, k
    drawn = set()
    for g in members[k]:
        pool = np.flatnonzero((bins == bins[common.index(g)]) & poolable)
        assert len(pool) >= N_CONTROLS, g
        drawn.update(common_arr[rng.choice(pool, N_CONTROLS, replace=False)])
    background[k] = sorted(drawn)
    coverage.append(dict(program=k, full_members=len(genes_k), shared_measured_members=len(members[k]),
                         omitted=';'.join(sorted(set(genes_k) - set(members[k]))), background_members=len(drawn)))
out = out_dir('clinical')
pd.DataFrame(coverage).to_csv(out / 'coverage.tsv', sep='\t', index=False)

# Cell scores: member mean minus background mean.
parts = []
for fn in [f'normal_spg_{s}.h5ad' for s in manifest.study.unique()] + CASE_FILES:
    obs, genes, x = load(fn, mappings)
    d = obs[['study', 'donor', 'diagnosis', 'state']].astype(str)
    d['cell_id'] = obs.index.astype(str)
    d['score_units'] = UNITS
    for k in members:
        raw = np.asarray(x[:, gene_columns(genes, members[k])].mean(1)).ravel()
        d[k + '_raw'] = raw
        d[k + '_adjusted'] = raw - np.asarray(x[:, gene_columns(genes, background[k])].mean(1)).ravel()
    parts.append(d.reset_index(drop=True))
cells = pd.concat(parts, ignore_index=True)
assert cells.cell_id.is_unique

ann = pd.read_csv(RESULTS / 'atlas' / 'cell_annotations.tsv.gz', sep='\t',
                  usecols=['study', 'source_cell_id', 'atlas_cell_type'])
atlas_spg = set(ann.loc[ann.study.eq(CONCORDANT_STUDY) & ann.atlas_cell_type.eq('SPG'), 'source_cell_id'])
selections = {PRIMARY: cells, CONCORDANT: cells[cells.study.eq(CONCORDANT_STUDY) & cells.cell_id.isin(atlas_spg)]}
metrics = [k + s for k in members for s in ['_raw', '_adjusted']]
scores, tests = [], []
for scenario, cc in selections.items():
    g = cc[cc.state.isin(STATES)].groupby(KEYS, observed=True)
    ds = g[metrics].mean().reset_index().merge(g.size().reset_index(name='n_cells'), on=KEYS)
    ds['scenario'] = scenario
    scores.append(ds)
    for study, case, control in CONTRASTS:
        if scenario == CONCORDANT and study != CONCORDANT_STUDY:
            continue
        for state in STATES:
            for min_cells in MIN_CELLS:
                z = ds[ds.study.eq(study) & ds.state.eq(state) & ds.n_cells.ge(min_cells)]
                assert not z.donor.duplicated().any()
                for metric in metrics:
                    row = dict(scenario=scenario, study=study, state=state, min_cells=min_cells, metric=metric,
                               case_label=case, control_label=control)
                    row.update(exact_welch(z.loc[z.diagnosis.eq(case), metric], z.loc[z.diagnosis.eq(control), metric]))
                    tests.append(row)
pd.concat(scores).to_csv(out / 'donor_state_scores.tsv', sep='\t', index=False)
effects = bh(pd.DataFrame(tests), FAMILIES)
effects.to_csv(out / 'all_effects.tsv', sep='\t', index=False)
print('universe', len(common), 'genes;', len(cells), 'cells;', len(effects), 'tests')
