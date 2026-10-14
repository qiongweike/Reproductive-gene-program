"""Figure 7 / Supplementary Table 4-5 / Supplementary Fig. 3 tables, the within-donor member summary and the
response-Z profiles (five focal targets and their matched controls, per donor) used by the response-profile PCA.
The Fig. 7d, e and g display sets were fixed before this scan and are listed below.
Inputs: results/knockout/readouts/<config>/, matching/primary/, perturbations/primary/, fixed_gene_sets.tsv
Outputs: results/knockout/panels/{Figure 7 and Supplementary Fig. 3 tables, within_donor_member_summary.tsv, <key>_response_Z.tsv.gz,
         PCA_profile_eligibility.tsv}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd

from ko_common import (FOCAL_TARGETS, KO, NETWORK_CONFIGS, PAIRED_CONFIGS, PAIRED_TARGETS, SAMPLE_KEYS, N_MP3,
                       fixed_gene_sets, read_tsv, write_tsv)

DISPLAY_TARGETS = PAIRED_TARGETS  # FGFR3 and COMMD3
FIG7D_FUNCTIONS = ['GO:1902850', 'GO:0034502', 'R-HSA-141424', 'R-HSA-2467813', 'R-HSA-69620', 'GO:0000776', 'GO:0051304',
                   'GO:0033047', 'GO:0007088', 'GO:0050673', 'GO:0007091', 'GO:0051303', 'GO:0000070', 'R-HSA-69206']
FIG7E_GENES = ['GAN', 'IGFBP6', 'MFHAS1', 'ZWINT', 'PAGE1', 'B2M', 'KLF6', 'MTCH1', 'CDC45', 'GINS2', 'HOXC9', 'TMED9']
FIG7G_FUNCTIONS = ['GO:0008283', 'GO:0042127', 'GO:0006260', 'GO:0008285', 'GO:0000075', 'GO:0000086']
MIN_RESOLVED = 10

P = KO / 'panels'
P.mkdir(parents=True, exist_ok=True)
primary = KO / 'readouts' / 'primary'


def with_donor(z):
    z = z.copy()
    z['donor'] = z.sample_key.str.split('__').str[-1]
    return z


def member_scan(mp3):
    """Fig. 7a and Supplementary Table 4: MP3 excess of all 49 members in four donors, ordering and missingness."""
    a = with_donor(read_tsv(primary / 'MP3_member_readouts.tsv'))
    assert len(a) == 4 * N_MP3 and a.target.nunique() == N_MP3 and set(a.target) == set(mp3)
    order = a.groupby('target', sort=False).agg(n_estimable=('comparison_eligible', 'sum'),
                                                median_delta=('delta_vs_control_median', 'median')).reset_index()
    order['all_four_donors'] = order.n_estimable.eq(4)
    order = order.sort_values(['all_four_donors', 'median_delta', 'target'], ascending=[False, False, True])
    write_tsv(a, P / 'MP3_member_responses.tsv')
    write_tsv(order, P / 'MP3_member_order.tsv')
    quality = read_tsv(KO / 'matching' / 'primary' / 'target_eligibility.tsv')
    write_tsv(a.merge(quality, on=['sample_key', 'target'], how='left', suffixes=('', '_network'), validate='one_to_one'),
              P / 'member_response_status.tsv')
    return a


def within_donor_summary(a):
    """Supplementary Methods: members above the control median / all five controls, median excess and FGFR3 rank."""
    rows = []
    for donor in ['A1', 'A2', 'N2', 'N3']:
        eligible = a[a.donor.eq(donor) & a.comparison_eligible].copy()
        eligible['member_rank'] = eligible.delta_vs_control_median.rank(ascending=False, method='min')
        fg = eligible[eligible.target.eq('FGFR3')].iloc[0]
        rows.append({'donor': donor, 'n_eligible_members': len(eligible),
                     'above_control_median': int(eligible.delta_vs_control_median.gt(0).sum()),
                     'above_all_five': int(eligible.exceeds_all_controls.sum()),
                     'median_member_excess': eligible.delta_vs_control_median.median(),
                     'FGFR3_excess': fg.delta_vs_control_median, 'FGFR3_member_rank': int(fg.member_rank)})
    pd.DataFrame(rows).to_csv(P / 'within_donor_member_summary.tsv', sep='\t', index=False)


def focal_panels(mp3):
    """Fig. 7b-f: focal knockouts with their controls, program and function excesses, gene-level response excesses."""
    write_tsv(with_donor(read_tsv(primary / 'focal_MP3_control_points.tsv')), P / 'focal_and_control_responses.tsv')
    programs = read_tsv(primary / 'all_target_program_control_comparisons.tsv.gz')
    write_tsv(with_donor(programs[programs.target.isin(FOCAL_TARGETS)]), P / 'focal_program_responses.tsv')
    func = read_tsv(primary / 'complete_functional_control_comparisons.tsv.gz')
    d = with_donor(func[func.pathway.isin(FIG7D_FUNCTIONS) & func.target.isin(DISPLAY_TARGETS)])
    assert len(d) == len(FIG7D_FUNCTIONS) * len(DISPLAY_TARGETS) * len(SAMPLE_KEYS)
    write_tsv(d, P / 'function_responses.tsv')
    names = func[['pathway', 'name']].drop_duplicates().set_index('pathway').name
    write_tsv(pd.DataFrame({'pathway': FIG7D_FUNCTIONS, 'name': names.loc[FIG7D_FUNCTIONS].to_numpy(),
                            'display_order': range(len(FIG7D_FUNCTIONS))}), P / 'function_definitions.tsv')
    assert not set(FIG7E_GENES) & set(mp3)
    write_tsv(pd.DataFrame({'gene': FIG7E_GENES}), P / 'member_selection.tsv')
    g = with_donor(read_tsv(primary / 'focal_gene_relative_controls.tsv.gz'))
    ix = pd.MultiIndex.from_product([SAMPLE_KEYS, DISPLAY_TARGETS, FIG7E_GENES], names=['sample_key', 'target', 'gene'])
    e = g[g.target.isin(DISPLAY_TARGETS) & g.gene.isin(FIG7E_GENES)].set_index(['sample_key', 'target', 'gene']).reindex(ix).reset_index()
    e['donor'] = e.sample_key.str.split('__').str[-1]
    write_tsv(e, P / 'member_donor_responses.tsv')
    f = g[g.target.isin(DISPLAY_TARGETS)].pivot(index=['sample_key', 'donor', 'gene'], columns='target',
                                                 values='delta_rank_percentile').reset_index()
    write_tsv(f, P / 'donor_gene_responses.tsv')


def paired_functions():
    """Fig. 7g and Supplementary Table 5: Full-SPG versus U-only excesses; NA unless both conditions are comparable."""
    parts = []
    for cfg in PAIRED_CONFIGS:
        z = read_tsv(KO / 'readouts' / cfg / 'complete_functional_control_comparisons.tsv.gz')
        parts.append(z[z.pathway.isin(FIG7G_FUNCTIONS)].assign(config=cfg))
    pair = with_donor(pd.concat(parts, ignore_index=True))
    by_position = pair.groupby(['sample_key', 'target', 'pathway'])
    assert by_position.controls.nunique().eq(1).all()  # identical controls in both networks
    pair['both_conditions_pass'] = by_position.comparison_eligible.transform('sum').eq(2)
    pair['raw_delta_vs_control_median'] = pair.delta_vs_control_median
    for col in ['delta_vs_control_median', 'delta_vs_control_max', 'delta_vs_control_min']:
        pair.loc[~pair.both_conditions_pass, col] = np.nan
    write_tsv(pair, P / 'paired_function_responses.tsv')


def sampling_configurations():
    """Supplementary Fig. 3: MP3 excess of the five focal targets in the three network configurations."""
    parts = []
    for cfg in NETWORK_CONFIGS:
        z = read_tsv(KO / 'readouts' / cfg / 'all_target_program_control_comparisons.tsv.gz')
        parts.append(z[z.pathway.eq('MP_3') & z.target.isin(FOCAL_TARGETS)].assign(config=cfg))
    write_tsv(with_donor(pd.concat(parts, ignore_index=True)), P / 'configuration_sensitivity.tsv')


def response_profiles():
    """Fig. 7h inputs: Z responses of focal and matched-control knockouts over genes common to all, targets removed."""
    matches = read_tsv(KO / 'matching' / 'primary' / 'matched_controls.tsv')
    matches = matches[matches.target.isin(FOCAL_TARGETS) & matches.match_pass]
    jobs = read_tsv(KO / 'matching' / 'primary' / 'perturbation_jobs.tsv')
    meta = []
    for key in SAMPLE_KEYS:
        mm = matches[matches.sample_key.eq(key)]
        requested = sorted(set(FOCAL_TARGETS) | set(mm.control))
        run = jobs[jobs.sample_key.eq(key) & jobs.perturbed_gene.isin(requested)].perturbed_gene
        drs = {g: read_tsv(KO / 'perturbations' / 'primary' / key / f'{g}_diffRegulation.tsv').set_index('gene') for g in run}
        common = set.intersection(*(set(z.index) for z in drs.values()))
        genes = sorted(common - set(requested))
        Z = pd.DataFrame({g: z.loc[genes, 'Z'] for g, z in drs.items()}, index=genes)
        cols = [g for g in FOCAL_TARGETS if g in Z] + [g for g in sorted(Z) if g not in FOCAL_TARGETS]
        Z = Z[cols]
        assert np.isfinite(Z.to_numpy()).all()
        for g in cols:
            v = Z[g].to_numpy()
            above = int(drs[g].loc[genes, 'above_noise'].sum())
            meta.append(dict(sample_key=key, donor=key.split('__')[1], target=g,
                             role=('Reference' if g == 'TCF3' else 'Candidate') if g in FOCAL_TARGETS else 'Control',
                             associated_candidates=';'.join(mm.loc[mm.control.eq(g), 'target'].sort_values().unique()),
                             n_common=len(genes), max_abs_Z_before=float(drs[g].Z.abs().max()),
                             max_abs_Z_after=float(abs(v).max()), n_above_noise=above,
                             noise_guard_pass=above >= MIN_RESOLVED, eligible=bool(abs(v).max() > 0 and np.std(v) > 0)))
        eligible = [x['target'] for x in meta if x['sample_key'] == key and x['eligible']]
        Z[eligible].rename_axis('gene').to_csv(P / f'{key}_response_Z.tsv.gz', sep='\t', compression='gzip')
        print(key, len(eligible), 'response profiles over', len(genes), 'genes')
    write_tsv(pd.DataFrame(meta), P / 'PCA_profile_eligibility.tsv')


mp3 = fixed_gene_sets()['MP_3']
within_donor_summary(member_scan(mp3))
focal_panels(mp3)
paired_functions()
sampling_configurations()
response_profiles()
