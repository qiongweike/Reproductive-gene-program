"""Adult testis atlas: QC-passing cells of GSE124263, GSE202647, GSE326372 and GSE254315 on commonly measured genes.

Sources keep their full-library QC, doublet calls and depth; a fixed multi-marker classifier assigns compartments
(>= 2 detected markers, standardized score >= 0.5, margin >= 0.25, donor-equal Normal/OA moments, z in [-3, 5]).
Outputs: data/processed/atlas/atlas_allcells_annotated.h5ad; results/atlas/cell_annotations.tsv.gz,
QC tables, identity-marker manifest and donor summaries, withheld-marker agreement.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from utils.paths import ROOT, RAW, PROCESSED, out_dir
from common import ATLAS_ASSIGN, ATLAS_HELDOUT, read_mtx

PRIMARY_SPG = ['normal_spg_GSE124263.h5ad', 'normal_spg_GSE254315.h5ad', 'disease_spg.h5ad']
MIN_MARKERS, MIN_SCORE, MIN_MARGIN, Z_CLIP = 2, 0.5, 0.25, (-3, 5)
OBS_COLUMNS = ['study', 'donor', 'donor_key', 'diagnosis', 'age', 'library', 'barcode', 'enrichment', 'author_cell_type',
               'source_cell_id', 'qc_source_file', 'count_source_file', 'count_type', 'qc_rule_kind',
               'total_counts', 'n_genes', 'pct_counts_mt', 'qc_basic_pass', 'qc_paper_available', 'qc_paper_pass',
               'doublet_score', 'predicted_doublet', 'doublet_screen_succeeded', 'qc_analysis_pass',
               'prior_cell_type', 'prior_state', 'key', 'primary_spg', 'primary_spg_state', 'primary_spg_source']
RES = out_dir('atlas')


def canonical(names):
    return pd.Index(['ATG101' if x == 'C12orf44' else x for x in names], name='gene_symbol')


def cell_key(o):
    return o.study.astype(str) + '|' + o.library.astype(str) + '|' + o.barcode.astype(str)


def primary_spg():
    """Membership and state of the primary spermatogonia, keyed by study|library|barcode."""
    frames = []
    for name in PRIMARY_SPG:
        a = ad.read_h5ad(PROCESSED / name, backed='r')
        o = a.obs.copy()
        a.file.close()
        o['primary_spg_source'] = name
        o['key'] = cell_key(o)
        frames.append(o)
    primary = pd.concat(frames)
    assert primary.key.is_unique
    return primary.set_index('key')


def qc_piece(a, source, kind, primary, measured=None, count_source=None):
    """Inherited QC flags for all barcodes of one source; integer counts of its QC-passing cells."""
    a.var_names = canonical(a.var_names)
    assert a.var_names.is_unique
    o = a.obs.copy()
    o['source_cell_id'] = o.index.astype(str)
    o['donor_key'] = o.study.astype(str) + ':' + o.donor.astype(str)
    o['qc_source_file'] = source
    o['count_source_file'] = count_source or source
    o['count_type'] = 'integer_UMI'
    o['qc_rule_kind'] = kind
    if kind == 'disease_shared_primary':
        o['qc_basic_pass'] = o.qc_pass_basic.astype(bool)
        o['qc_paper_available'] = o.qc_author_thresholds_available.astype(bool)
        o['qc_paper_pass'] = o.qc_author_thresholds.astype(bool)
        o['doublet_screen_succeeded'] = o.doublet_tested.astype(bool)
        o['qc_analysis_pass'] = o.qc_pass.astype(bool)
        o['prior_cell_type'] = o.cell_type.astype(str)
        o['prior_state'] = o.state.astype(str)
    else:
        o['qc_paper_available'] = kind == 'age_paper_and_basic'
        o['qc_analysis_pass'] = (o.qc_basic_pass.astype(bool) & o.qc_paper_pass.astype(bool)
                                 & ~o.predicted_doublet.astype(bool) & o.doublet_screen_succeeded.astype(bool))
        o['prior_cell_type'] = o.get('reference_cell_type', o.author_cell_type).astype(str)
        o['prior_state'] = o.get('state', pd.Series('unknown', index=o.index)).astype(str)
    o['key'] = cell_key(o)
    assert o.key.is_unique
    o['primary_spg'] = o.key.isin(primary.index)
    o['primary_spg_state'] = o.key.map(primary.state.astype(str)).fillna('not_primary_SPG')
    o['primary_spg_source'] = o.key.map(primary.primary_spg_source).fillna('not_primary_SPG')
    raw = sparse.csr_matrix(a.layers['counts'] if 'counts' in a.layers else a.X)
    raw.eliminate_zeros()
    assert np.allclose(np.asarray(raw.sum(1)).ravel(), o.total_counts.to_numpy(), atol=1e-5), source
    assert np.isfinite(raw.data).all() and (raw.data >= 0).all() and np.equal(raw.data, np.rint(raw.data)).all()
    assert ((~o.predicted_doublet.astype(bool)) | (~o.qc_analysis_pass)).all()
    assert o.loc[o.primary_spg, 'qc_analysis_pass'].all(), source
    measured = {'ATG101' if g == 'C12orf44' else g for g in (a.var_names if measured is None else measured)}
    assert measured <= set(a.var_names)
    for c in OBS_COLUMNS:
        if c not in o:
            o[c] = np.nan if c in ['age', 'doublet_score'] else 'not_available'
    keep = o.qc_analysis_pass.to_numpy()
    obs = o.loc[keep, OBS_COLUMNS].copy()
    obs.index = pd.Index(obs.key, name='cell_id')
    return ad.AnnData(raw[keep].astype(np.int32), obs=obs, var=a.var.copy()), o, measured


def gse254315_pieces():
    """Per-donor GSE254315 QC tables joined to the raw integer counts."""
    for table in sorted((PROCESSED / 'normal_age_samples').glob('*_allcells_qc_annotation.tsv.gz')):
        o = pd.read_csv(table, sep='\t', index_col=0, low_memory=False)
        matrix = next((RAW / 'GSE254315').glob(table.name.split('_')[0] + '_*_matrix.mtx.gz'))
        stem = matrix.name.replace('_matrix.mtx.gz', '')
        features = pd.read_csv(matrix.parent / (stem + '_features.tsv.gz'), sep='\t', header=None)
        barcodes = pd.read_csv(matrix.parent / (stem + '_barcodes.tsv.gz'), sep='\t', header=None)[0].astype(str)
        assert barcodes.tolist() == o.barcode.astype(str).tolist()
        x = read_mtx(matrix).astype(np.int32)
        names = features[1].astype(str).replace({'C12orf44': 'ATG101'})
        unique = pd.Index(pd.unique(names), name='gene_symbol')
        mapper = sparse.csr_matrix((np.ones(len(names), dtype=np.int32), (np.arange(len(names)), unique.get_indexer(names))),
                                   shape=(len(names), len(unique)))
        yield ad.AnnData((x @ mapper).tocsr(), obs=o, var=pd.DataFrame(index=unique)), table, matrix


def ingest():
    primary = primary_spg()
    pieces, allobs, gene_sets = [], [], {}

    def add(*args, **kwargs):
        piece, o, measured = qc_piece(*args, primary=primary, **kwargs)
        for study in o.study.astype(str).unique():
            gene_sets[study] = measured if study not in gene_sets else gene_sets[study] & measured
        pieces.append(piece)
        allobs.append(o)

    add(ad.read_h5ad(PROCESSED / 'normal_allcells_GSE124263.h5ad'), 'data/processed/normal_allcells_GSE124263.h5ad',
        'author_retained_and_basic')
    disease = ad.read_h5ad(PROCESSED / 'disease_allcells.h5ad')
    presence = pd.read_csv(PROCESSED / 'disease_gene_presence_by_study.tsv', sep='\t', index_col=0)
    for study in ['GSE202647', 'GSE326372']:
        add(disease[disease.obs.study.astype(str).eq(study), :].copy(), 'data/processed/disease_allcells.h5ad',
            'disease_shared_primary', measured=presence.index[presence[study].eq(1)])
    del disease
    for b, table, matrix in gse254315_pieces():
        add(b, str(table.relative_to(ROOT)), 'age_paper_and_basic', count_source=str(matrix.relative_to(ROOT)))
    common = sorted(set.intersection(*gene_sets.values()))
    pd.DataFrame({'gene_symbol': common}).to_csv(RES / 'common_measured_genes.tsv', sep='\t', index=False)
    A = ad.concat([a[:, common].copy() for a in pieces], join='inner', merge='same')
    pieces.clear()
    assert A.obs_names.is_unique and int(A.obs.primary_spg.sum()) == len(primary)
    A.layers['counts'] = A.X.copy()
    norm = (sparse.diags(1e4 / A.obs.total_counts.to_numpy()) @ A.X).astype(np.float32).tocsr()
    norm.data = np.log1p(norm.data)
    A.X = norm
    # Categorical obs with natural-sorted categories (as stored in h5ad) fix the donor order of the moments.
    A.strings_to_categoricals()
    Q = pd.concat(allobs, ignore_index=True)
    Q.to_csv(RES / 'all_barcodes_qc_inheritance.tsv.gz', sep='\t', index=False)
    Q.groupby(['study', 'donor_key', 'diagnosis'], observed=True).agg(
        raw_cells=('key', 'size'), basic=('qc_basic_pass', 'sum'), singlets_retained=('qc_analysis_pass', 'sum'),
        predicted_doublets=('predicted_doublet', 'sum'), primary_spg=('primary_spg', 'sum')).reset_index().to_csv(
        RES / 'donor_qc_summary.tsv', sep='\t', index=False)
    return A


def annotate(A):
    """Fixed marker classifier; writes the identity-marker tables and returns A with annotations in obs."""
    assert not (set(sum(ATLAS_ASSIGN.values(), [])) & set(sum(ATLAS_HELDOUT.values(), [])))
    rules = {'assignment_markers': ATLAS_ASSIGN, 'heldout_markers': ATLAS_HELDOUT, 'min_detected_assignment_genes': MIN_MARKERS,
             'min_top_standardized_score': MIN_SCORE, 'min_top_minus_second_score': MIN_MARGIN, 'gene_z_clip': list(Z_CLIP),
             'normalizer': 'donor-equal normal/OA full-library logCP10k means and variances',
             'assignment_unit': 'individual cell; ambiguous cells retained as Uncertain',
             'heldout_role': 'diagnostic only, never used to choose label'}
    json.dump(rules, open(RES / 'annotation_prespecified_rules.json', 'w'), indent=2)
    o = A.obs
    labels = list(ATLAS_ASSIGN)
    genes = list(dict.fromkeys(sum(ATLAS_ASSIGN.values(), []) + sum(ATLAS_HELDOUT.values(), [])))
    manifest = pd.DataFrame([{'cell_type': ct, 'role': role, 'gene': g, 'measured': g in A.var_names}
                             for role, panel in [('assignment', ATLAS_ASSIGN), ('heldout', ATLAS_HELDOUT)]
                             for ct, gs in panel.items() for g in gs])
    manifest.to_csv(RES / 'identity_marker_manifest.tsv', sep='\t', index=False)
    genes = [g for g in genes if g in A.var_names]
    X = A[:, genes].X.toarray().astype(np.float64)
    det = A[:, genes].layers['counts'].toarray() > 0
    # Reference moments: unenriched Normal/OA libraries, each donor weighted equally.
    normal = o.diagnosis.astype(str).isin(['Normal', 'OA']) & ~o.enrichment.astype(str).str.contains('enriched', case=False)
    moments, seconds = [], []
    for _, ix in o.loc[normal].groupby('donor_key', observed=True).indices.items():
        ids = np.flatnonzero(normal.to_numpy())[ix]
        moments.append(X[ids].mean(0))
        seconds.append((X[ids] ** 2).mean(0))
    assert len(moments) >= 2
    mu = np.mean(moments, 0)
    sd = np.sqrt(np.maximum(np.mean(seconds, 0) - mu ** 2, 1e-4))
    Z = np.clip((X - mu) / sd, *Z_CLIP)
    scores, ns = [], []
    for ct in labels:
        ix = [genes.index(g) for g in ATLAS_ASSIGN[ct] if g in genes]
        assert len(ix) >= 2
        scores.append(Z[:, ix].mean(1))
        ns.append(det[:, ix].sum(1))
    S, N = np.column_stack(scores), np.column_stack(ns)
    order = np.argsort(S, axis=1)
    top, second, r = order[:, -1], order[:, -2], np.arange(A.n_obs)
    margin = S[r, top] - S[r, second]
    ok = (S[r, top] >= MIN_SCORE) & (margin >= MIN_MARGIN) & (N[r, top] >= MIN_MARKERS)
    ann = np.array(labels, dtype=object)[top]
    ann[~ok] = 'Uncertain'
    o['atlas_cell_type'] = ann
    o['annotation_top_type'] = np.array(labels)[top]
    o['annotation_score'] = S[r, top]
    o['annotation_margin'] = margin
    o['assignment_markers_detected'] = N[r, top]
    o['annotation_source'] = 'prespecified_individual_multimarker'
    o['annotation_confident'] = ok
    for j, ct in enumerate(labels):
        o['assignment_score_' + ct] = S[:, j]
        o['assignment_detected_' + ct] = N[:, j]
    V = []
    for ct in labels:
        ix = [genes.index(g) for g in ATLAS_HELDOUT[ct] if g in genes]
        V.append(Z[:, ix].mean(1))
        o['heldout_score_' + ct] = V[-1]
        o['heldout_detected_' + ct] = det[:, ix].sum(1)
    o['heldout_top_type'] = np.array(labels)[np.column_stack(V).argmax(1)]
    o['heldout_type_agrees'] = o.heldout_top_type.eq(o.atlas_cell_type)
    o['state'] = o.primary_spg_state.astype(str)
    o['full_assay_total_counts'] = o.total_counts.astype(float)
    o['primary_spg_annotation_conflict'] = o.primary_spg & (o.atlas_cell_type != 'SPG')
    # Donor x compartment x gene summaries (the donor is the summary unit), then donor-equal means.
    rows = []
    for (study, donor, diagnosis, ct), ix in o.groupby(['study', 'donor_key', 'diagnosis', 'atlas_cell_type'], observed=True).indices.items():
        mean_log, mean_cp, freq = X[ix].mean(0), np.expm1(X[ix]).mean(0), det[ix].mean(0)
        for j, g in enumerate(genes):
            rows.append({'study': study, 'donor_key': donor, 'diagnosis': diagnosis, 'cell_type': ct, 'gene': g, 'n_cells': len(ix),
                         'mean_log1pCP10k': mean_log[j], 'mean_CP10k': mean_cp[j], 'detection_fraction': freq[j]})
    M = pd.DataFrame(rows).merge(manifest[['cell_type', 'role', 'gene']].rename(columns={'cell_type': 'marker_for'}), on='gene', how='left')
    M.to_csv(RES / 'identity_marker_donor_summary.tsv.gz', sep='\t', index=False)
    M[M.diagnosis.isin(['Normal', 'OA'])].groupby(['cell_type', 'gene', 'role', 'marker_for'], observed=True).agg(
        mean_CP10k=('mean_CP10k', 'mean'), mean_log1pCP10k=('mean_log1pCP10k', 'mean'), detection_fraction=('detection_fraction', 'mean'),
        n_donors=('donor_key', 'nunique'), n_cells=('n_cells', 'sum')).reset_index().to_csv(
        RES / 'identity_marker_donor_equal_summary.tsv', sep='\t', index=False)
    o.groupby(['study', 'donor_key', 'diagnosis', 'atlas_cell_type'], observed=True).agg(
        n_cells=('study', 'size'), heldout_type_agreement=('heldout_type_agrees', 'mean'),
        median_margin=('annotation_margin', 'median')).reset_index().to_csv(RES / 'annotation_heldout_diagnostics.tsv', sep='\t', index=False)
    return A


if __name__ == '__main__':
    A = annotate(ingest())
    (PROCESSED / 'atlas').mkdir(parents=True, exist_ok=True)
    A.write_h5ad(PROCESSED / 'atlas' / 'atlas_allcells_annotated.h5ad', compression='gzip')
    A.obs.to_csv(RES / 'cell_annotations.tsv.gz', sep='\t', index_label='cell_id')
    print(f'atlas: {A.n_obs} cells x {A.n_vars} genes, {int((A.obs.atlas_cell_type != "Uncertain").sum())} assigned', flush=True)
