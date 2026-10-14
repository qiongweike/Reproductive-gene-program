"""Hotspot (hotspotsc 1.1.3) autocorrelation and pairwise local correlation of the 49 MP3 members in each donor and state.
Contexts with >=100 cells; exact unweighted 30-neighbour graphs in the saved 20-PC coordinates (rows subset within states,
PCA not refitted); full-library depth. Bernoulli models in all donors, DANB only for all-integer counts (no rounding).
Members enter a fit when uniquely mapped and non-constant; other members and pairs are kept as NA with a reason.
Inputs: results/trajectory/{donor_manifest.tsv, <donor>.primary.pca.tsv.gz}, results/programs/MP_genes.tsv, count matrices.
Outputs: results/trajectory/hotspot/{context_eligibility.tsv, fit_status.tsv, MP3_member_autocorrelations.tsv.gz, MP3_pair_local_correlations.tsv.gz}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir
import itertools

import anndata as ad
import hotspot
import numpy as np
import pandas as pd
from scipy import io, sparse
from statsmodels.stats.multitest import multipletests

SEED = 20261003
PCS = [f'PC{i}' for i in range(1, 21)]
N_NEIGHBORS = 30
MIN_STATE_CELLS = 100
MIN_DETECTION = 0.05
INTEGER_TOLERANCE = 1e-10
STATES = ['undifferentiated', 'transition', 'differentiating']
MODELS = ['bernoulli', 'danb']
META = ['context', 'sample_key', 'study', 'donor', 'scope', 'state', 'n_cells']
PROGRAMS = RESULTS / 'programs'
TRAJ = RESULTS / 'trajectory'
OUT = out_dir('trajectory', 'hotspot')


def donor_files(sample):
    return (PROGRAMS / 'nmf/input' / f'{sample}.counts.mtx.gz', PROGRAMS / 'input' / f'{sample}.genes.tsv',
            PROGRAMS / 'nmf/input' / f'{sample}.cells.tsv')


def load_donor(sample, target):
    counts_file, gene_file, cell_file = donor_files(sample)
    raw = io.mmread(counts_file).T.tocsr().astype(np.float64)
    upper = pd.read_csv(gene_file, sep='\t').gene.astype(str).str.upper()
    dup = upper.duplicated(False).to_numpy()
    gi = {g: i for i, g in enumerate(upper) if not dup[i]}
    cells = pd.read_csv(cell_file, sep='\t')
    pca = pd.read_csv(TRAJ / f'{sample}.primary.pca.tsv.gz', sep='\t').set_index('cell_id').loc[cells.cell_id]
    lib = np.asarray(raw.sum(1)).ravel()
    assert np.allclose(lib, cells.full_assay_counts_recomputed, rtol=1e-10) and np.all(lib > 0)
    integer = not (np.abs(raw.data - np.rint(raw.data)) > INTEGER_TOLERANCE).any()
    counts = np.column_stack([raw[:, gi[g]].toarray().ravel() if g in gi else np.full(len(cells), np.nan) for g in target])
    return cells, counts, lib, pca[PCS].to_numpy(), integer


def run_hotspot(counts, lib, pc, cells, names, model):
    a = ad.AnnData(X=sparse.csc_matrix(counts), obs=pd.DataFrame({'full_depth': lib}, index=cells),
                   var=pd.DataFrame(index=names), obsm={'independent_PCA': pc})
    hs = hotspot.Hotspot(a, model=model, latent_obsm_key='independent_PCA', umi_counts_obs_key='full_depth')
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS, approx_neighbors=False)
    auto = hs.compute_autocorrelations(jobs=1)
    hs.compute_local_correlations(list(names), jobs=1)
    return auto, hs.local_correlation_c, hs.local_correlation_z


def main():
    manifest = pd.read_csv(TRAJ / 'donor_manifest.tsv', sep='\t')
    mp = pd.read_csv(PROGRAMS / 'MP_genes.tsv', sep='\t')
    target = mp[mp.MP.eq('MP_3')].sort_values('rank').gene.str.upper().tolist()
    assert len(target) == len(set(target)) == 49
    pairs = list(itertools.combinations(range(len(target)), 2))
    np.random.seed(SEED)
    contexts, status, member_rows, pair_rows = [], [], [], []
    for m in manifest.itertuples():
        cells, all_counts, all_lib, all_pc, integer = load_donor(m.sample_key, target)
        mapped = np.isfinite(all_counts).all(0)
        for scope, st in [('whole_donor', 'ALL')] + [('within_state', x) for x in STATES]:
            ix = np.arange(len(cells)) if st == 'ALL' else np.flatnonzero(cells.state.eq(st))
            n = len(ix)
            ctx = dict(context=f'{m.sample_key}__{st}', sample_key=m.sample_key, study=m.study, donor=m.donor, scope=scope,
                       state=st, n_cells=n, state_eligible=n >= MIN_STATE_CELLS, DANB_integer_eligible=integer,
                       bernoulli_run_eligible=n >= MIN_STATE_CELLS, danb_run_eligible=n >= MIN_STATE_CELLS and integer,
                       reason='PASS' if n >= MIN_STATE_CELLS else 'FEWER_THAN_100_CELLS', PCA_reestimated_within_state=False)
            contexts.append(ctx)
            counts, lib, pc, ids = all_counts[ix], all_lib[ix], all_pc[ix], cells.cell_id.to_numpy(dtype=str)[ix]
            for model in MODELS:
                meta = {k: ctx[k] for k in META}
                meta['model'] = model
                run_ok = ctx[model + '_run_eligible']
                reasons = [] if ctx['state_eligible'] else ['FEWER_THAN_100_CELLS']
                if model == 'danb' and not integer:
                    reasons.append('NONINTEGER_RAW_COUNTS')
                vals = (counts > 0).astype(float) if model == 'bernoulli' else counts
                varying = np.ptp(np.where(np.isfinite(vals), vals, 0), axis=0) > 0 if n else np.zeros(len(target), dtype=bool)
                eligible = mapped & varying & run_ok
                det = np.mean(counts > 0, axis=0) if n else np.full(len(target), np.nan)
                det[~mapped] = np.nan
                gene_reason = [';'.join(reasons) if reasons else
                               ('MISSING_OR_AMBIGUOUS' if not mapped[j] else
                                ('CONSTANT_DETECTION' if model == 'bernoulli' else 'CONSTANT_COUNTS') if not varying[j] else 'PASS')
                               for j in range(len(target))]
                ac = pd.DataFrame(index=target, columns=['C', 'Z', 'Pval', 'FDR'], dtype=float)
                pair_c = pd.DataFrame(np.nan, index=target, columns=target)
                pair_z = pair_c.copy()
                rec = dict(**meta, status='SKIPPED_INELIGIBLE', reason=';'.join(reasons), n_computational_members=int(eligible.sum()))
                if run_ok:
                    assert eligible.sum() >= 2
                    names = np.array(target)[eligible]
                    auto, c, z = run_hotspot(counts[:, eligible], lib, pc, ids, names, model)
                    ac.loc[names] = auto.loc[names]
                    pair_c.loc[names, names] = c.loc[names, names]
                    pair_z.loc[names, names] = z.loc[names, names]
                    n_auto = int(np.isfinite(ac.C).sum())
                    n_pairs = int(np.isfinite(pair_c.to_numpy()[np.triu_indices(len(target), 1)]).sum())
                    complete = n_auto == len(names) and n_pairs == len(names) * (len(names) - 1) // 2
                    rec.update(status='FIT_RETURNED' if complete else 'RETURNED_WITH_NONFINITE',
                               reason='PASS' if complete else 'NUMERICAL_NONFINITE_VALUES_RETAINED_NA',
                               n_finite_auto=n_auto, n_finite_pairs=n_pairs)
                status.append(rec)
                p = ac.Pval.to_numpy()
                q = multipletests(np.where(np.isfinite(p), p, 1), method='fdr_bh')[1]
                q[~np.isfinite(p)] = np.nan
                for j, g in enumerate(target):
                    estimated = bool(np.isfinite(ac.loc[g, ['C', 'Z', 'Pval']]).all())
                    reason = 'PASS' if estimated else ('NUMERICAL_NONFINITE' if eligible[j] else gene_reason[j])
                    member_rows.append(dict(**meta, gene=g, rank=j + 1, mapped=bool(mapped[j]), model_outcome_nonconstant=bool(varying[j]),
                                            detection=det[j], comparison_detection_eligible=bool(mapped[j] and det[j] >= MIN_DETECTION),
                                            estimable=estimated, reason=reason, C=ac.at[g, 'C'], Z=ac.at[g, 'Z'],
                                            Pval_model_conditional=ac.at[g, 'Pval'], FDR_author_tested=ac.at[g, 'FDR'], q_BH_fixed49=q[j]))
                for i, j in pairs:
                    cij, zij = pair_c.iat[i, j], pair_z.iat[i, j]
                    estimated = bool(np.isfinite(cij) and np.isfinite(zij))
                    reason = 'PASS' if estimated else (
                        ';'.join(dict.fromkeys(gene_reason[k] for k in (i, j) if gene_reason[k] != 'PASS')) or 'NUMERICAL_NONFINITE')
                    pair_rows.append(dict(**meta, gene1=target[i], gene2=target[j], estimable=estimated, reason=reason,
                                          both_detection_ge05=bool(mapped[i] and mapped[j] and det[i] >= MIN_DETECTION and det[j] >= MIN_DETECTION),
                                          C=cij, Z=zij))
    pd.DataFrame(contexts).to_csv(OUT / 'context_eligibility.tsv', sep='\t', index=False)
    status = pd.DataFrame(status)
    status.to_csv(OUT / 'fit_status.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(member_rows).to_csv(OUT / 'MP3_member_autocorrelations.tsv.gz', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(pair_rows).to_csv(OUT / 'MP3_pair_local_correlations.tsv.gz', sep='\t', index=False, na_rep='NA')
    fitted = status[status.status.eq('FIT_RETURNED')].model.value_counts()
    print(f"Hotspot fits: {fitted.get('bernoulli', 0)} Bernoulli, {fitted.get('danb', 0)} DANB; "
          f"{status.status.eq('SKIPPED_INELIGIBLE').sum()} skipped")


if __name__ == '__main__':
    main()
