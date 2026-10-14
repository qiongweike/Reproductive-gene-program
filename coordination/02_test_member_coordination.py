"""Within-state coordination of fixed program members against 1,000 expression-matched gene sets, per donor x state
(>=50 cells), for the discovery donors and the external NOA cases. Member ranks are residualized against nothing, depth
and detected genes, or these plus G1/S and G2/M scores; one-sided competitive P and BH within model x level.
Inputs: results/programs/{donor_manifest,MP_genes}.tsv, results/programs/identity/state_markers.tsv,
        resources/candidate_143_genes.txt, results/coordination/external_inputs/donor_manifest.tsv
Outputs: results/coordination/coherence/{coherence_results,donor_state_eligibility}.tsv (and in external_NOA/)
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import ROOT, RESULTS, RESOURCES, out_dir
import hashlib

import numpy as np
import pandas as pd
from scipy import io, stats
from statsmodels.stats.multitest import multipletests

SEED = 20261003
N_RANDOM = 1000
N_NEAREST = 100                 # matching candidates per member
MIN_CLASS_POOL = 20
PRIMARY_MIN_CELLS = 100
SENSITIVITY_MIN_CELLS = 50
MIN_DETECTION = 0.05
MIN_MODULE_GENES = 20
MIN_MODULE_FRACTION = 0.60
MIN_CYCLE_GENES = 5
STATES = ['undifferentiated', 'transition', 'differentiating']
CYCLE_SETS = ['G1/S', 'G2/M']
MODELS = ['raw', 'technical', 'technical_cycle']
COHORTS = {'discovery': (RESULTS / 'programs' / 'donor_manifest.tsv', ()),
           'external_NOA': (RESULTS / 'coordination' / 'external_inputs' / 'donor_manifest.tsv', ('external_NOA',))}


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def eligible_donors(path):
    m = pd.read_csv(path, sep='\t')
    return m[m.eligible]


def assay_genes(r):
    return set(pd.read_csv(resolve(r.gene_file), sep='\t').gene.astype(str))


def fixed_gene_sets():
    mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t')
    sets = {k: list(v.gene.str.upper()) for k, v in mp.groupby('MP')}
    sets['candidate_143'] = [g.upper() for g in (RESOURCES / 'candidate_143_genes.txt').read_text().split()]
    return sets


def cycle_genes(exclude, common):
    """Published G1/S and G2/M markers outside all tested sets and inside the common assay."""
    markers = pd.read_csv(RESULTS / 'programs' / 'identity' / 'state_markers.tsv', sep='\t')
    cycle = {k: sorted((set(markers.loc[markers.state.eq(k), 'gene'].str.upper()) - exclude) & set(common))
             for k in CYCLE_SETS}
    assert min(map(len, cycle.values())) >= MIN_CYCLE_GENES, cycle
    return cycle


def eligibility(manifest):
    rows = []
    for r in manifest.itertuples():
        state = pd.read_csv(resolve(r.cell_file), sep='\t').state
        for s in STATES + ['unknown']:
            n = int(state.eq(s).sum())
            level = 'not_estimated' if s == 'unknown' or n < SENSITIVITY_MIN_CELLS else \
                'primary' if n >= PRIMARY_MIN_CELLS else 'sensitivity_50_99'
            reason = 'unknown_state' if s == 'unknown' else 'fewer_than_50_cells' if n < SENSITIVITY_MIN_CELLS else ''
            rows.append(dict(sample_key=r.sample_key, study=r.study, donor=r.donor, state=s, n_cells=n,
                             level=level, reason=reason))
    return pd.DataFrame(rows)


def standardized_residual(y, c):
    q, _, rank, _ = np.linalg.lstsq(c, y, rcond=1e-10)
    r = y - c @ q
    norm = np.linalg.norm(r, axis=0)
    keep = norm > 1e-8
    r[:, keep] /= norm[keep]
    r[:, ~keep] = 0
    return r, keep, int(rank)


def mean_corr(z, ix):
    """Mean off-diagonal correlation of unit-norm columns: (|sum of columns|^2 - G) / (G(G-1))."""
    n = len(ix)
    return float((np.square(z[:, ix].sum(1)).sum() - n) / (n * (n - 1)))


def gene_class(g):
    return 'mitochondrial' if g.startswith('MT-') else 'ribosomal' if g.startswith(('RPL', 'RPS')) else 'other'


def stratum_seed(sample_key, state, module):
    return SEED + int(hashlib.sha256(f'{sample_key}|{state}|{module}'.encode()).hexdigest()[:7], 16)


def matched_sets(selected, pool, classes, feats, seed):
    """Random sets drawn member by member from the 100 nearest same-class genes, without replacement."""
    rng = np.random.default_rng(seed)
    candidates = []
    for j in selected:
        pp = pool[classes[pool] == classes[j]]
        assert len(pp) >= MIN_CLASS_POOL
        d = ((feats[pp] - feats[j]) ** 2).sum(1)
        candidates.append(pp[np.argsort(d, kind='stable')[:N_NEAREST]])
    sets = np.empty((N_RANDOM, len(selected)), dtype=int)
    for b in range(N_RANDOM):
        used = set()
        for j, pp in enumerate(candidates):
            v = int(rng.choice(np.array([x for x in pp if x not in used])))
            sets[b, j] = v
            used.add(v)
    return sets


def test_donor(r, sets, cycle, common, exclude):
    counts = io.mmread(resolve(r.matrix_file)).tocsr().T.tocsr()
    genes = pd.read_csv(resolve(r.gene_file), sep='\t').gene.astype(str).to_numpy()
    meta = pd.read_csv(resolve(r.cell_file), sep='\t')
    assert counts.shape == (len(meta), len(genes)) and len(set(genes)) == len(genes)
    depth = np.asarray(counts.sum(1)).ravel()
    ng = np.asarray((counts > 0).sum(1)).ravel()
    assert np.allclose(depth, meta.full_assay_counts_recomputed, rtol=1e-10, atol=1e-8)
    gi = {g: i for i, g in enumerate(genes)}
    common = [g for g in common if g in gi]
    rows = []
    for state in STATES:
        ids = np.flatnonzero(meta.state.eq(state).to_numpy())
        n = len(ids)
        if n < SENSITIVITY_MIN_CELLS:
            continue
        level = 'primary' if n >= PRIMARY_MIN_CELLS else 'sensitivity_50_99'
        s = counts[ids][:, [gi[g] for g in common]].tocsr()
        s = s.multiply(1e4 / depth[ids, None]).tocsr()
        det = np.asarray((s > 0).mean(0)).ravel()
        means = np.asarray(s.mean(0)).ravel()
        expr = s.toarray()
        np.log1p(expr, out=expr)
        cov = [np.log1p(depth[ids]), np.log1p(ng[ids])] + \
              [expr[:, [common.index(g) for g in gs]].mean(1) for gs in cycle.values()]
        c = np.column_stack([stats.rankdata(x) / n for x in cov])
        c -= c.mean(0)
        sd = c.std(0)
        active = sd > 1e-12
        c[:, active] /= sd[active]
        c[:, ~active] = 0
        eligible = (det >= MIN_DETECTION) & (means > 0)
        allg = np.asarray(common)[eligible]
        mmean, mdet = means[eligible], det[eligible]
        y = stats.rankdata(expr[:, eligible], axis=0) / n
        del expr, s
        raw = y - y.mean(0)
        raw /= np.maximum(np.linalg.norm(raw, axis=0), 1e-30)
        tech, kt, _ = standardized_residual(y, np.column_stack([np.ones(n), c[:, :2]]))
        full, kf, rf = standardized_residual(y, np.column_stack([np.ones(n), c]))
        del y
        usable = kt & kf
        ixg = {g: i for i, g in enumerate(allg)}
        feats = np.column_stack([stats.rankdata(mmean) / len(mmean), stats.rankdata(mdet) / len(mdet)])
        classes = np.array([gene_class(g) for g in allg])
        pool = np.array([i for i, g in enumerate(allg) if usable[i] and g not in exclude])
        for module, gs in sets.items():
            selected = [ixg[g] for g in gs if g in ixg and usable[ixg[g]]]
            base = dict(sample_key=r.sample_key, study=r.study, donor=r.donor, state=state, level=level, n_cells=n,
                        module=module, n_fixed_genes=len(gs), n_retained_genes=len(selected),
                        retained_fraction=len(selected) / len(gs), n_pairs=len(selected) * (len(selected) - 1) // 2,
                        residual_df=n - rf)
            if len(selected) < MIN_MODULE_GENES or len(selected) / len(gs) < MIN_MODULE_FRACTION:
                rows += [dict(base, model=m, status='insufficient_retained_members') for m in MODELS]
                continue
            random_sets = matched_sets(selected, pool, classes, feats, stratum_seed(r.sample_key, state, module))
            for model, z in zip(MODELS, [raw, tech, full]):
                observed = mean_corr(z, selected)
                null = np.array([mean_corr(z, rs) for rs in random_sets])
                rows.append(dict(base, model=model, status='estimated_descriptive_competitive_test',
                                 mean_pairwise_rho=observed, null_mean=float(null.mean()),
                                 null_sd=float(null.std(ddof=1)), null_q025=float(np.quantile(null, .025)),
                                 null_q975=float(np.quantile(null, .975)), excess_rho=observed - float(null.mean()),
                                 p_competitive=float((1 + np.sum(null >= observed)) / (len(null) + 1)),
                                 n_random=len(null)))
        del raw, tech, full
    return rows


def add_bh(res):
    res['q_BH'] = np.nan
    res['BH_family'] = ''
    for (model, level), g in res[res.p_competitive.notna()].groupby(['model', 'level']):
        res.loc[g.index, 'q_BH'] = multipletests(g.p_competitive, method='fdr_bh')[1]
        res.loc[g.index, 'BH_family'] = f'{model}__{level}'
    return res


def main():
    sets = fixed_gene_sets()
    union = set().union(*map(set, sets.values()))
    discovery = eligible_donors(COHORTS['discovery'][0])
    common = sorted(set.intersection(*[assay_genes(r) for r in discovery.itertuples()]))
    cycle = cycle_genes(union, common)
    exclude = union | set().union(*map(set, cycle.values()))
    for name in sys.argv[1:] or COHORTS:
        manifest = eligible_donors(COHORTS[name][0])
        genes = sorted(set(common).intersection(*[assay_genes(r) for r in manifest.itertuples()]))
        assert all(set(v) <= set(genes) for v in cycle.values())
        out = out_dir('coordination', 'coherence', *COHORTS[name][1])
        eligibility(manifest).to_csv(out / 'donor_state_eligibility.tsv', sep='\t', index=False)
        rows = [row for r in manifest.itertuples() for row in test_donor(r, sets, cycle, genes, exclude)]
        res = add_bh(pd.DataFrame(rows))
        res.to_csv(out / 'coherence_results.tsv', sep='\t', index=False)
        primary = res[res.model.eq('technical_cycle') & res.level.eq('primary') & res.q_BH.notna()]
        print(name, len(genes), 'common genes;', len(primary), 'primary technical_cycle tests,',
              int((primary.q_BH < .05).sum()), 'with q < 0.05')


if __name__ == '__main__':
    main()
