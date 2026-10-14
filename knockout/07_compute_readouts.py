"""Knockout readouts and target-minus-control comparisons (five matched controls) for all five configurations. A set's
response is the rank-biserial effect of a tie-corrected two-sided rank-sum test on unsigned non-self distances (noise tied
at zero): MP1-3 and the 143-gene set in every knockout (one BH family); for focal targets and their controls also GO/Reactome
(one family), nine fixed modules (separate family) and gene percentiles, all with the 49 MP3 genes removed.
Inputs: results/knockout/{matching,networks,perturbations}/<config>/, fixed_gene_sets.tsv, results/coordination/pathways/membership.tsv.gz
Outputs: results/knockout/readouts/<config>/*.tsv[.gz]
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS
import numpy as np
import pandas as pd
from scipy import sparse, stats

from ko_common import (FOCAL_TARGETS, KO, N_CONTROLS, NETWORK_CONFIGS, PAIRED_CONFIGS, PAIRED_KEYS, PAIRED_TARGETS,
                       SAMPLE_KEYS, all_targets, fixed_gene_sets, hgnc_resolver, read_tsv, write_tsv)

PROGRAMS = ['MP_1', 'MP_2', 'MP_3', 'fixed_143']
MIN_SET = 10            # programs and ontology sets
MIN_MODULE = 3          # fixed small modules
MIN_BACKGROUND = 10
MIN_RESOLVED = 10       # noise-resolved genes in the universe
MIN_RESOLVED_SET = 5    # noise-resolved set members (3 for modules)


def bh(p):
    p = np.asarray(p, float)
    ix = np.argsort(p)
    q = np.empty(len(p))
    q[ix] = np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1].clip(0, 1)
    return q


def rank_tests(genes, distance, resolved, gene_sets, min_sizes):
    """Competitive rank-sum tests of each gene set within one knockout; also returns gene response percentiles."""
    N = len(genes)
    lookup = {g: i for i, g in enumerate(genes)}
    rows, cols = [], []
    for k, members in enumerate(gene_sets):
        ids = [lookup[g] for g in members if g in lookup]
        rows.extend([k] * len(ids))
        cols.extend(ids)
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(gene_sets), N))
    rank = stats.rankdata(distance, method='average')
    n = np.asarray(A.sum(1)).ravel()
    n_resolved = np.asarray(A @ resolved.astype(float)).ravel()
    u = np.asarray(A @ rank).ravel() - n * (n + 1) / 2
    den = n * (N - n)
    with np.errstate(divide='ignore', invalid='ignore'):
        effect = 2 * u / den - 1
        _, tie_counts = np.unique(distance, return_counts=True)
        tie = np.sum(tie_counts.astype(float) ** 3 - tie_counts)
        var = den / 12 * (N + 1 - tie / (N * (N - 1)))
        p = 2 * stats.norm.sf(np.abs(u - den / 2) / np.sqrt(var))
    min_sizes = np.asarray(min_sizes)
    min_resolved = np.where(min_sizes == MIN_MODULE, MIN_MODULE, MIN_RESOLVED_SET)
    checks = [n < min_sizes, N - n < MIN_BACKGROUND, np.full(len(n), resolved.sum() < MIN_RESOLVED),
              n_resolved < min_resolved, ~np.isfinite(p)]
    ok = ~np.logical_or.reduce(checks)
    reason = np.select(checks, ['SMALL_SET_COVERAGE', 'SMALL_BACKGROUND', 'NOISE_LT10_NONSELF',
                                'RESOLVED_MEMBERS_TOO_FEW', 'DEGENERATE_RANK'], default='PASS')
    table = pd.DataFrame(dict(n_universe=N, n_set=n.astype(int), n_resolved_universe=int(resolved.sum()),
                              n_resolved_set=n_resolved.astype(int), raw_rank_biserial=effect,
                              rank_biserial=np.where(ok, effect, np.nan), conditional_gene_rank_P=np.where(ok, p, np.nan),
                              eligible=ok, reason=reason))
    return table, (rank - .5) / N


def function_library(sets):
    """GO/Reactome sets (canonical symbols) followed by the nine fixed modules."""
    resolve = hgnc_resolver()
    ontology = pd.read_csv(RESULTS / 'coordination' / 'pathways' / 'membership.tsv.gz', sep='\t')
    ontology['gene'] = ontology.gene.map(lambda g: resolve(g)[0] if isinstance(g, str) else g)
    ontology = ontology.drop_duplicates(['pathway', 'gene'])
    meta = ontology[['pathway', 'database', 'name']].drop_duplicates().reset_index(drop=True)
    members = {k: set(v.gene) for k, v in ontology.groupby('pathway', sort=False)}
    modules = [k for k in sets if k.startswith('module_')]
    library = pd.concat([meta, pd.DataFrame({'pathway': modules, 'database': 'module', 'name': modules})], ignore_index=True)
    return library, [members[k] for k in meta.pathway] + [set(sets[k]) for k in modules], [MIN_SET] * len(meta) + [MIN_MODULE] * len(modules)


def knockout_readouts(cfg, keys, focal, sets, library, matches, jobs):
    """Program tests for every knockout; function tests and gene percentiles for focal targets and their controls."""
    library_meta, library_sets, library_min = library
    mp3 = set(sets['MP_3'])
    program_meta = pd.DataFrame({'pathway': PROGRAMS, 'database': 'program', 'name': PROGRAMS})
    program_tabs, function_tabs, percentiles = [], [], {}
    for key in keys:
        nodes = read_tsv(KO / 'networks' / cfg / key / 'WT_node_features.tsv')
        assert not nodes.gene.duplicated().any() and not nodes.mapping_ambiguous.any()
        fm = matches[(matches.sample_key == key) & matches.target.isin(focal) & matches.match_pass]
        with_functions = set(focal) | set(fm.control)
        percentiles[key] = {}
        for job in jobs[jobs.sample_key == key].itertuples():
            g = job.perturbed_gene
            dr = read_tsv(KO / 'perturbations' / cfg / key / f'{g}_diffRegulation.tsv')
            for scope in ['all_nonself'] + (['outside_MP3'] if g in with_functions else []):
                x = dr[~dr.gene.isin({g} | (mp3 if scope == 'outside_MP3' else set()))]
                genes = x.gene.to_numpy()
                resolved = x.above_noise.to_numpy(bool)
                distance = np.where(resolved, x.distance, 0)
                if scope == 'all_nonself':
                    tab, _ = rank_tests(genes, distance, resolved, [set(sets[p]) - {g} for p in PROGRAMS], [MIN_SET] * 4)
                    tab = pd.concat([program_meta, tab], axis=1)
                else:
                    tab, pct = rank_tests(genes, distance, resolved, library_sets, library_min)
                    tab = pd.concat([library_meta, tab], axis=1)
                    percentiles[key][g] = pd.Series(pct if resolved.sum() >= MIN_RESOLVED else np.full(len(pct), np.nan), index=genes)
                tab['sample_key'], tab['perturbed_gene'], tab['scope'], tab['role'] = key, g, scope, job.role
                (program_tabs if scope == 'all_nonself' else function_tabs).append(tab)
    programs = pd.concat(program_tabs, ignore_index=True)
    functions = pd.concat(function_tabs, ignore_index=True)
    programs['q_conditional_family'] = np.nan
    functions['q_conditional_family'] = np.nan
    ok = programs.eligible
    programs.loc[ok, 'q_conditional_family'] = bh(programs.loc[ok, 'conditional_gene_rank_P'])
    for is_module in [False, True]:
        ok = functions.eligible & functions.database.eq('module').eq(is_module)
        functions.loc[ok, 'q_conditional_family'] = bh(functions.loc[ok, 'conditional_gene_rank_P'])
    by_ko = functions.groupby(['sample_key', 'perturbed_gene'])
    functions['effect_rank_within_KO'] = by_ko.rank_biserial.rank(ascending=False, method='min')
    functions['q_rank_within_KO'] = by_ko.q_conditional_family.rank(method='min')
    return programs, functions, percentiles


def compare_with_controls(tests, key, target, controls, status, pathways, scope):
    """Target statistic minus the median (and max/min) of its five control knockouts, per pathway."""
    ss = tests[tests.sample_key.eq(key)]
    t = ss[ss.perturbed_gene.eq(target)].set_index('pathway').reindex(pathways)
    c = (ss[ss.perturbed_gene.isin(controls)].pivot(index='pathway', columns='perturbed_gene', values='rank_biserial')
         .reindex(index=pathways, columns=controls))
    n, med, mn, mx = c.notna().sum(1), c.median(1), c.min(1), c.max(1)
    matched = bool(status.control_matching_eligible)
    te = t.rank_biserial.notna()
    valid = te & (n == N_CONTROLS) & matched
    result = pd.DataFrame(dict(
        pathway=pathways, sample_key=key, target=target, scope=scope, target_rank_biserial=t.rank_biserial.to_numpy(),
        target_eligible=te.to_numpy(), target_reason=t.reason.fillna(status.KO_reason).to_numpy(),
        n_controls_matching=len(controls), n_controls_estimable=n.to_numpy(), control_median=med.to_numpy(),
        control_min=mn.to_numpy(), control_max=mx.to_numpy(), controls=';'.join(controls),
        matching_reason=status.matching_reason, comparison_eligible=valid.to_numpy(),
        delta_vs_control_median=(t.rank_biserial - med).where(valid).to_numpy(),
        delta_vs_control_max=(t.rank_biserial - mx).where(valid).to_numpy(),
        delta_vs_control_min=(t.rank_biserial - mn).where(valid).to_numpy(),
        exceeds_all_controls=(t.rank_biserial > mx).where(valid).to_numpy()))
    result['comparison_reason'] = np.select(
        [np.full(len(pathways), not matched), ~te.to_numpy(), n.to_numpy() < N_CONTROLS],
        [status.matching_reason, 'TARGET_READOUT_INELIGIBLE', 'FEWER_THAN_FIVE_ELIGIBLE_CONTROL_READOUTS'], default='PASS')
    return result


def gene_responses(percentiles, key, target, controls, matched):
    """Per-gene response percentile of the target minus the median of its controls (program-external genes)."""
    pr = percentiles[key][target]
    cc = pd.concat([percentiles[key][g].rename(g) for g in controls], axis=1) if controls else pd.DataFrame(index=pr.index)
    cc = cc.reindex(pr.index)
    n = cc.notna().sum(1)
    z = pd.DataFrame(dict(gene=pr.index, target_rank_percentile=pr.to_numpy(),
                          control_median_rank_percentile=cc.median(1).to_numpy(), n_controls_with_gene=n.to_numpy(),
                          delta_rank_percentile=(pr - cc.median(1)).where((n == N_CONTROLS) & matched).to_numpy()))
    z['sample_key'], z['target'], z['matching_eligible'] = key, target, matched
    return z


def run_config(cfg, sets, library):
    paired = cfg in PAIRED_CONFIGS
    keys = PAIRED_KEYS if paired else SAMPLE_KEYS
    focal = PAIRED_TARGETS if paired else FOCAL_TARGETS
    targets = all_targets(sets) if cfg == 'primary' else focal
    out = KO / 'readouts' / cfg
    out.mkdir(parents=True, exist_ok=True)
    matches = read_tsv(KO / 'matching' / cfg / 'matched_controls.tsv')
    quality = read_tsv(KO / 'matching' / cfg / 'target_eligibility.tsv')
    jobs = read_tsv(KO / 'matching' / cfg / 'perturbation_jobs.tsv')
    programs, functions, percentiles = knockout_readouts(cfg, keys, focal, sets, library, matches, jobs)
    write_tsv(programs, out / 'all_perturbations_four_programs.tsv.gz')
    write_tsv(functions, out / 'complete_GO_Reactome_and_modules.tsv.gz')

    library_meta = library[0]
    program_cmp, function_cmp, points, genes = [], [], [], []
    for key in keys:
        for target in targets:
            status = quality[(quality.sample_key == key) & quality.target.eq(target)].iloc[0]
            matched = bool(status.control_matching_eligible)
            m = matches[(matches.sample_key == key) & matches.target.eq(target)]
            controls = m.control.tolist() if matched else []
            program_cmp.append(compare_with_controls(programs, key, target, controls, status, PROGRAMS, 'all_nonself'))
            if target not in focal:
                continue
            fc = compare_with_controls(functions, key, target, controls, status, library_meta.pathway.tolist(), 'outside_MP3')
            function_cmp.append(fc.merge(library_meta, on='pathway', how='left', validate='one_to_one'))
            for rank, g in enumerate([target] + controls):
                z = programs[(programs.sample_key == key) & programs.perturbed_gene.eq(g) & programs.pathway.eq('MP_3')]
                points.append(dict(sample_key=key, target=target, probe=g,
                                   probe_role='candidate_reference' if rank == 0 else 'matched_control', match_rank=rank,
                                   rank_biserial=z.rank_biserial.iloc[0] if len(z) else np.nan,
                                   readout_eligible=bool(z.eligible.iloc[0]) if len(z) else False,
                                   matching_eligible=matched))
            if target in percentiles[key]:
                genes.append(gene_responses(percentiles, key, target, controls, matched))
    pc = pd.concat(program_cmp, ignore_index=True)
    write_tsv(pc, out / 'all_target_program_control_comparisons.tsv.gz')
    write_tsv(pc[pc.pathway.eq('MP_3') & pc.target.isin(sets['MP_3'])], out / 'MP3_member_readouts.tsv')
    write_tsv(pd.concat(function_cmp, ignore_index=True), out / 'complete_functional_control_comparisons.tsv.gz')
    write_tsv(pd.DataFrame(points), out / 'focal_MP3_control_points.tsv')
    write_tsv(pd.concat(genes, ignore_index=True), out / 'focal_gene_relative_controls.tsv.gz')
    print(cfg, programs.perturbed_gene.groupby(programs.sample_key).nunique().sum(), 'knockouts;',
          int((pc.pathway.eq('MP_3') & pc.comparison_eligible).sum()), 'eligible MP3 comparisons')


sets = fixed_gene_sets()
library = function_library(sets)
for cfg in NETWORK_CONFIGS + PAIRED_CONFIGS:
    run_config(cfg, sets, library)
