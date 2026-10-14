"""SPG-trained program scores (Gavish/scalop, saved matched controls, centred on the donor's SPG) in spermatocytes and
spermatogonia of the same donors; strict atlas and author-labelled spermatocytes kept apart; donor sign-flip tests.
Inputs: results/programs/{donor_manifest,MP_genes}.tsv, results/programs/scores/, results/atlas/cell_annotations.tsv.gz,
        results/coordination/specificity/matched_control_genes.tsv, full-library counts (data/raw, data/processed)
Outputs: results/coordination/specificity/{MP_paired_donor_values,MP_paired_effects,donor_state_scores,
         state_paired_donor_differences,state_paired_effects}.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import ROOT, RAW, METADATA, PROCESSED, RESULTS, out_dir
import gzip
import itertools

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse, stats
from statsmodels.stats.multitest import multipletests

STATES = ['undifferentiated', 'transition', 'differentiating']
MIN_PAIR_CELLS = 30             # cells per class for a paired donor
MIN_STATE_CELLS = 30
MIN_HELDOUT_SPC_MARKERS = 2     # withheld spermatocyte markers detected (strict labels)
AUTHOR_SPC_LABELS = ['Leptotene', 'Zygotene', 'Pachytene', 'Diplotene', 'Meiotic divisions']
STRICT_STRATUM, AUTHOR_STRATUM = 'strict_label_primary', 'author_label_auxiliary'   # identity_stratum labels
# HGNC-verified previous symbols, applied before summing full-library rows that name the same gene.
ALIASES = {'C12ORF44': 'ATG101', 'FAM26F': 'CALHM6', 'HIST1H1D': 'H1-3', 'FAM207A': 'SLX9', 'C4ORF3': 'ARLN',
           'RP11-325I22.2': 'LINC01511'}


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def save(rows, name):
    pd.DataFrame(rows).to_csv(out_dir('coordination', 'specificity') / name, sep='\t', index=False, na_rep='NA')


def canonical(x):
    s = str(x).upper()
    return ALIASES.get(s, s)


def collapse(x, genes):
    """Sum rows (genes) that map to the same canonical symbol."""
    names = pd.Index([canonical(g) for g in genes])
    unique = pd.Index(pd.unique(names))
    if len(unique) != len(names):
        tr = sparse.csr_matrix((np.ones(len(names)), (np.arange(len(names)), unique.get_indexer(names))),
                               shape=(len(names), len(unique)))
        x = (x @ tr).tocsr()
    x = sparse.csr_matrix(x, dtype=np.float64)
    x.eliminate_zeros()
    assert np.isfinite(x.data).all() and (x.data >= 0).all()
    return x, unique


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


def author_counts(path, barcodes):
    """Selected barcode columns of a dense gene x barcode table, as cells x genes."""
    inds, vals, indptr, genes = [], [], [0], []
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        look = {v: i for i, v in enumerate(header)}
        positions = np.array([look[b.rsplit('_', 1)[0]] for b in barcodes])
        for line in f:
            if not line.strip():
                continue
            gene, values = line.rstrip('\r\n').split('\t', 1)
            v = np.fromstring(values, sep='\t')
            assert len(v) == len(header)
            v = v[positions]
            nz = np.flatnonzero(v)
            inds.extend(nz.tolist())
            vals.extend(v[nz].tolist())
            indptr.append(len(inds))
            genes.append(gene)
    return sparse.csr_matrix((vals, inds, indptr), shape=(len(genes), len(barcodes))).T.tocsr(), genes


def strict_spermatocytes(row, spg, atlas):
    """Confident atlas spermatocytes (>=2 withheld markers) that are not program-input SPG, with full-library counts."""
    a = atlas[(atlas.study == row.study) & (atlas.donor == row.donor)]
    sp = a[(a.atlas_cell_type == 'Spermatocyte') & a.annotation_confident &
           (a.heldout_detected_Spermatocyte >= MIN_HELDOUT_SPC_MARKERS) & a.qc_analysis_pass].copy()
    sp = sp[~(sp.primary_spg | sp.source_cell_id.isin(spg.cell_id) | sp.cell_id.isin(spg.cell_id))]
    sp['cell_id'] = sp.source_cell_id
    if row.study == 'GSE124263':
        obj = ad.read_h5ad(PROCESSED / 'normal_allcells_GSE124263.h5ad')
        obs = obj.obs[obj.obs.donor.astype(str).eq(row.donor)]
        x, genes, ids = obj[obs.index].layers['counts'], obj.var_names.copy(), obs.index.astype(str)
    else:
        path = next((RAW / row.study).glob(str(a.library.iloc[0]) + '_*_matrix.mtx.gz'))
        stem = path.name.replace('_matrix.mtx.gz', '')
        genes = pd.read_csv(path.parent / f'{stem}_features.tsv.gz', sep='\t', header=None)[1]
        barcodes = pd.read_csv(path.parent / f'{stem}_barcodes.tsv.gz', sep='\t', header=None)[0]
        x, ids = io.mmread(path).T.tocsr(), row.study + ':' + row.donor + ':' + barcodes
    return sp, x, genes, pd.Index(ids), STRICT_STRATUM


def author_spermatocytes(row, spg):
    """Author-labelled meiotic cells passing whole-library QC; 'Diff.spermatogonia/Preleptotene' is not used."""
    author = pd.read_csv(METADATA / 'GSE153947_Cell_metadata.tsv.gz', sep='\t', index_col=0)
    author = author[author.Sample_ID.eq(row.donor)].copy()
    author['raw_barcode'] = author.index.str.rsplit('_', n=1).str[0]
    qc = pd.read_csv(PROCESSED / f'GSE153947_{row.donor}_whole_library_doublet_QC.tsv.gz', sep='\t')
    author = author.join(qc.set_index('raw_barcode'), on='raw_barcode', validate='one_to_one')
    assert author.qc_basic_pass.notna().all()
    keep = author.qc_basic_pass & author.qc_paper_pass & ~author.predicted_doublet & author.doublet_screen_succeeded
    sp = author[keep & author.Cluster_identity.isin(AUTHOR_SPC_LABELS)].copy()
    sp['cell_id'] = row.study + ':' + row.donor + ':' + sp.index
    sp['barcode'] = sp.index
    sp['library'] = str(spg.library.iloc[0])
    sp = sp[~sp.cell_id.isin(spg.cell_id)]
    path = next(iter(list((RAW / 'GSE153947').glob(f'*_{row.donor}_counts.tsv.gz')) +
                     list(RAW.glob(f'*_{row.donor}_counts.tsv.gz'))))
    selected = pd.concat([spg[['cell_id', 'barcode']], sp[['cell_id', 'barcode']]], ignore_index=True)
    assert selected.cell_id.is_unique
    x, genes = author_counts(path, selected.barcode)
    return sp, x, genes, pd.Index(selected.cell_id), AUTHOR_STRATUM


def spg_and_spc(row, spg, atlas):
    """Full-library counts (canonical genes) and metadata of the donor's SPG followed by its SPC."""
    if row.study == 'GSE153947':
        sp, x, genes, ids, stratum = author_spermatocytes(row, spg)
    else:
        sp, x, genes, ids, stratum = strict_spermatocytes(row, spg, atlas)
    sp = sp.assign(cell_class='Spermatocyte', state='not_SPG', identity_stratum=stratum)
    meta = pd.concat([spg.assign(cell_class='SPG', identity_stratum=stratum), sp], ignore_index=True)
    assert meta.cell_id.is_unique
    loc = ids.get_indexer(meta.cell_id)
    assert min(loc) >= 0
    x, genes = collapse(x[loc], genes)
    depth = np.asarray(x.sum(1)).ravel()
    assert (depth > 0).all()
    assert np.allclose(depth[:len(spg)], spg.full_assay_counts_recomputed, atol=1e-5)
    assert np.allclose(depth[len(spg):], sp.total_counts, atol=1e-5)
    meta['full_library_counts'] = depth
    return x, genes, meta


def score_donor(row, programs, control, atlas):
    """Per-cell SPG-centred scores of SPG and SPC, checked against the program-scoring output for SPG."""
    raw = io.mmread(resolve(row.matrix_file)).T.tocsr()
    genes = pd.Index(pd.read_csv(resolve(row.gene_file), sep='\t').gene.map(canonical))
    assert genes.is_unique
    spg = pd.read_csv(resolve(row.cell_file), sep='\t')
    assert np.allclose(np.asarray(raw.sum(1)).ravel(), spg.full_assay_counts_recomputed)
    full, fg, meta = spg_and_spc(row, spg, atlas)
    fi = {g: i for i, g in enumerate(fg)}
    assert set(genes) == set(fg)
    delta = full[:len(spg), [fi[g] for g in genes]] - raw
    assert delta.nnz == 0 or np.max(np.abs(delta.data)) < 1e-6
    saved = pd.read_csv(RESULTS / 'programs' / 'scores' / f'{row.sample_key}_scores.tsv.gz', sep='\t')
    saved = saved.set_index('cell_id').loc[spg.cell_id]
    logx = full.multiply(1e5 / meta.full_library_counts.to_numpy()[:, None]).tocsr()
    logx.data = np.log2(1 + logx.data)
    cells, strata = [], []
    for mp, gs in programs.items():
        target = [fi[g] for g in gs if g in fi]
        ct = control[(control.sample_key == row.sample_key) & (control.MP == mp)]
        assert ct.gene.isin(fg).all() and ct.multiplicity.sum() == ct.n_control_draws.iloc[0]
        target_mean = np.asarray(logx[:, target].mean(1)).ravel()
        weights = ct.multiplicity.to_numpy() / ct.multiplicity.sum()
        rawscore = target_mean - np.asarray(logx[:, [fi[g] for g in ct.gene]] @ weights).ravel()
        score = rawscore - rawscore[:len(spg)].mean()
        assert np.max(np.abs(score[:len(spg)] - saved[mp].to_numpy())) < 1e-9, (row.sample_key, mp)
        cells.append(meta[['cell_id', 'identity_stratum', 'cell_class']].assign(
            sample_key=row.sample_key, study=row.study, donor=row.donor, MP=mp, score=score,
            un_centered_matched_score=rawscore, target_mean_log2CP100K=target_mean))
        for state in STATES:
            ii = np.flatnonzero(spg.state.eq(state))
            strata.append(dict(sample_key=row.sample_key, study=row.study, donor=row.donor, state=state, MP=mp,
                               n_cells=len(ii), eligible=len(ii) >= MIN_STATE_CELLS,
                               mean_score=float(score[ii].mean()) if len(ii) else np.nan))
    return cells, strata


def paired_values(cells):
    keys = ['identity_stratum', 'sample_key', 'study', 'donor', 'MP']
    means = cells.groupby(keys + ['cell_class'], observed=True).agg(
        n_cells=('cell_id', 'size'), mean_score=('score', 'mean'), mean_raw_score=('un_centered_matched_score', 'mean'),
        mean_target=('target_mean_log2CP100K', 'mean')).reset_index()
    pairs = []
    for key, g in means.groupby(keys, observed=True):
        s = g.set_index('cell_class').reindex(['SPG', 'Spermatocyte'])
        s['n_cells'] = s.n_cells.fillna(0)
        pairs.append(dict(zip(keys, key), eligible=bool(min(s.n_cells) >= MIN_PAIR_CELLS),
                          n_SPG=s.loc['SPG', 'n_cells'], n_Spermatocyte=s.loc['Spermatocyte', 'n_cells'],
                          SPG_mean=s.loc['SPG', 'mean_score'], SPC_mean=s.loc['Spermatocyte', 'mean_score'],
                          SPG_minus_SPC=s.loc['SPG', 'mean_score'] - s.loc['Spermatocyte', 'mean_score'],
                          raw_difference=s.loc['SPG', 'mean_raw_score'] - s.loc['Spermatocyte', 'mean_raw_score'],
                          target_only_difference=s.loc['SPG', 'mean_target'] - s.loc['Spermatocyte', 'mean_target']))
    return pd.DataFrame(pairs)


def state_contrasts(st, programs):
    diffs, tests = [], []
    for mp in programs:
        for a, b in itertools.combinations(STATES, 2):
            aa = st[(st.MP == mp) & st.state.eq(a) & st.eligible].set_index('sample_key')
            bb = st[(st.MP == mp) & st.state.eq(b) & st.eligible].set_index('sample_key')
            shared = aa.index.intersection(bb.index)
            v = bb.loc[shared, 'mean_score'] - aa.loc[shared, 'mean_score']
            tests.append(dict(MP=mp, contrast=b + ' minus ' + a, BH_family='all_9_state_contrasts', **two_sided(v)))
            diffs += [dict(MP=mp, sample_key=s, contrast=b + ' minus ' + a, difference=v[s],
                           n_first=aa.loc[s, 'n_cells'], n_second=bb.loc[s, 'n_cells']) for s in shared]
    return diffs, tests


def main():
    manifest = pd.read_csv(RESULTS / 'programs' / 'donor_manifest.tsv', sep='\t').query('eligible')
    mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t')
    programs = {k: g.sort_values('rank').gene.map(canonical).tolist() for k, g in mp.groupby('MP', observed=True)}
    assert {k: len(set(v)) for k, v in programs.items()} == {'MP_1': 50, 'MP_2': 50, 'MP_3': 49}
    control = pd.read_csv(RESULTS / 'coordination' / 'specificity' / 'matched_control_genes.tsv', sep='\t')
    atlas = pd.read_csv(RESULTS / 'atlas' / 'cell_annotations.tsv.gz', sep='\t')
    cells, strata = [], []
    for row in manifest.itertuples():
        c, s = score_donor(row, programs, control, atlas)
        cells += c
        strata += s
    cells = pd.concat(cells, ignore_index=True)
    pairs = paired_values(cells)
    tests = [dict(identity_stratum=stratum, MP=mp, contrast='SPG minus Spermatocyte',
                  studies=';'.join(sorted(g.study.unique())), **two_sided(g.SPG_minus_SPC))
             for (stratum, mp), g in pairs[pairs.eligible].groupby(['identity_stratum', 'MP'], observed=True)]
    st = pd.DataFrame(strata)
    diffs, state_tests = state_contrasts(st, programs)
    save(pairs, 'MP_paired_donor_values.tsv')
    save(bh(tests, ['identity_stratum']), 'MP_paired_effects.tsv')
    save(st, 'donor_state_scores.tsv')
    save(diffs, 'state_paired_donor_differences.tsv')
    save(bh(state_tests, ['BH_family']), 'state_paired_effects.tsv')
    print('paired donors:', pairs[pairs.eligible].groupby('identity_stratum').sample_key.nunique().to_dict())


if __name__ == '__main__':
    main()
