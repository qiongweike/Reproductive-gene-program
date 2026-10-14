"""QC, Scrublet, Leiden clusters and marker scores for GSE202647 (OA/NOA) and GSE326372 (NOA) libraries.

Cells need >= 200 detected genes and mitochondrial fraction < 20% and must not be Scrublet doublets
(scanpy, rate 0.05, 30 PCs, per library). GSE235321 author log2TPM spermatogonia are only stored:
the study is not analysed, but its measured genes enter the shared gene universe used for clinical scoring.
Outputs: data/processed/GSE202647_qc.h5ad, GSE326372_qc.h5ad, gse235321_spg.h5ad.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gzip
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '4')

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse
from threadpoolctl import threadpool_limits

from utils.paths import RAW, PROCESSED
from common import SEED, MIN_GENES, MAX_MT_PCT, DOUBLET_RATE, DOUBLET_PCS, DISEASE_PANELS, DISEASE_STATE_MARKERS

sc.settings.n_jobs = 4
LIBRARIES = {'GSE202647': [('NOA1', 'GSM6127766'), ('NOA2', 'GSM6127767'), ('NOA3', 'GSM6127768'),
                           ('OA1', 'GSM6127769'), ('OA2', 'GSM6127770')],
             'GSE326372': [(f'testis_{i + 1}', f'GSM{9630067 + i}') for i in range(5)]}
MT_THRESHOLDS = (10, 20, 30)            # 20 is the main rule; 10 and 30 are stored for sensitivity analyses
LATE_GENES = ('PRM1', 'PRM2', 'TNP1', 'TNP2', 'AKAP4', 'ODF1')
N_HVG, N_PCS, N_NEIGHBORS, LEIDEN_RESOLUTION = 3000, 40, 20, 0.6


def collapse_features(matrix, features):
    """Sum features sharing a symbol (float32); join their gene IDs."""
    symbols = features.iloc[:, 1].astype(str).to_numpy()
    indices, unique = pd.factorize(symbols, sort=False)
    mapping = sparse.csr_matrix((np.ones(len(indices)), (np.arange(len(indices)), indices)), shape=(len(indices), len(unique)))
    collapsed = (matrix @ mapping).astype(np.float32).tocsr()
    gene_ids = features.iloc[:, 0].astype(str).groupby(symbols, sort=False).agg(';'.join)
    return collapsed, pd.DataFrame({'gene_ids': gene_ids.reindex(unique).values}, index=pd.Index(unique, name='gene_symbol'))


def read_library(study, sample, gsm):
    base = RAW / study
    prefix = gsm + '_' + sample + ('.' if study == 'GSE202647' else '_')
    features = pd.read_csv(base / (prefix + ('genes.tsv.gz' if study == 'GSE202647' else 'features.tsv.gz')), sep='\t', header=None)
    barcodes = pd.read_csv(base / (prefix + 'barcodes.tsv.gz'), sep='\t', header=None)[0].astype(str)
    with gzip.open(base / (prefix + 'matrix.mtx.gz'), 'rb') as stream:
        matrix = io.mmread(stream).tocsr().T.tocsr()
    assert matrix.shape == (len(barcodes), len(features))
    assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all() and np.equal(matrix.data, np.rint(matrix.data)).all()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix, var = collapse_features(matrix, features)
    donor = study + '_' + sample
    obs = pd.DataFrame({'barcode': barcodes.values, 'study': study, 'sample': sample, 'donor': donor, 'library': gsm,
                        'diagnosis': 'OA' if study == 'GSE202647' and sample.startswith('OA') else 'NOA',
                        'pathology': 'not_individually_resolved', 'author_cell_type': 'not_available',
                        'author_state': 'not_available', 'expression_type': 'integer_UMI_counts', 'enrichment': 'unselected'},
                       index=pd.Index([donor + ':' + b for b in barcodes], name='cell_id'))
    return ad.AnnData(matrix, obs=obs, var=var)


def library_qc(adata):
    """Per-library QC flags, Scrublet on QC-passing cells, late-germ-cell RNA fraction, then ln(1 + CP10k)."""
    counts = adata.X
    total = np.asarray(counts.sum(axis=1)).ravel()
    genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    mt = np.asarray(counts[:, adata.var_names.str.upper().str.startswith('MT-')].sum(axis=1)).ravel()
    adata.obs['total_counts'] = total
    adata.obs['n_genes'] = genes
    adata.obs['pct_counts_mt'] = 100 * mt / np.maximum(total, 1)
    for threshold in MT_THRESHOLDS:
        adata.obs['qc_mt%s' % threshold] = (genes >= MIN_GENES) & (adata.obs.pct_counts_mt < threshold)
    adata.obs['qc_pass_basic'] = adata.obs['qc_mt%s' % MAX_MT_PCT]
    lc = np.log1p(total)
    adata.obs['high_count_flag'] = lc > np.median(lc) + 5 * max(np.median(np.abs(lc - np.median(lc))), 0.01)
    adata.obs['predicted_doublet'] = False
    adata.obs['doublet_score'] = np.nan
    passed = adata.obs.qc_pass_basic.to_numpy()
    sub = adata[passed].copy()
    sc.pp.scrublet(sub, expected_doublet_rate=DOUBLET_RATE, random_state=SEED, n_prin_comps=DOUBLET_PCS, verbose=False)
    adata.obs.loc[sub.obs_names, 'doublet_score'] = sub.obs.doublet_score
    adata.obs.loc[sub.obs_names, 'predicted_doublet'] = sub.obs.predicted_doublet.astype(bool)
    adata.obs['doublet_tested'] = passed
    adata.obs['qc_pass'] = passed & ~adata.obs.predicted_doublet.to_numpy()
    late = counts[:, adata.var_names.get_indexer([g for g in LATE_GENES if g in adata.var_names])]
    adata.obs['late_rna_fraction'] = np.asarray(late.sum(axis=1)).ravel() / np.maximum(total, 1)
    adata.obs['late_genes_detected'] = np.asarray((late > 0).sum(axis=1)).ravel()
    adata.obs['PRM1_detected'] = np.asarray(counts[:, adata.var_names.get_loc('PRM1')].todense()).ravel() > 0 if 'PRM1' in adata.var_names else False
    adata.layers['counts'] = counts.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def cluster_study(adata):
    """Leiden clusters of QC-passing cells; marker scores of all cells; cluster label = best compartment panel."""
    valid = adata.obs.qc_pass.to_numpy()
    sub = adata[valid].copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, flavor='seurat', batch_key='library')
    sc.tl.pca(sub, n_comps=N_PCS, mask_var='highly_variable', svd_solver='arpack', random_state=SEED)
    sc.pp.neighbors(sub, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)
    sc.tl.leiden(sub, resolution=LEIDEN_RESOLUTION, random_state=SEED, key_added='leiden', flavor='igraph', n_iterations=2)
    adata.obs['leiden'] = 'QC_excluded'
    adata.obs.loc[sub.obs_names, 'leiden'] = sub.obs['leiden'].astype(str)
    for label, panel in dict(DISEASE_PANELS, **DISEASE_STATE_MARKERS).items():
        genes = [gene for gene in panel if gene in adata.var_names]
        sc.tl.score_genes(adata, genes, score_name='marker_' + label, random_state=SEED, use_raw=False)
        raw = adata.layers['counts'][:, adata.var_names.get_indexer(genes)]
        adata.obs['detected_' + label] = np.asarray((raw > 0).sum(axis=1)).ravel()
    means = adata.obs.loc[valid].groupby('leiden', observed=True)[['marker_' + name for name in DISEASE_PANELS]].mean()
    best = means.idxmax(axis=1).str.replace('marker_', '', regex=False)
    adata.obs['cluster_cell_type'] = adata.obs.leiden.map(best).fillna('QC_excluded').astype(str)
    return adata


def process_counts(study):
    libraries = [library_qc(read_library(study, sample, gsm)) for sample, gsm in LIBRARIES[study]]
    combined = ad.concat(libraries, join='outer', merge='first', fill_value=0)
    combined.varm['feature_available_by_library'] = pd.DataFrame(
        {a.obs.library.iloc[0]: combined.var_names.isin(a.var_names) for a in libraries}, index=combined.var_names)
    del libraries
    combined = cluster_study(combined)
    combined.write_h5ad(PROCESSED / (study + '_qc.h5ad'), compression='gzip')
    print(f'{study}: {combined.n_obs} cells, {int(combined.obs.qc_pass.sum())} pass QC', flush=True)


def process_tpm():
    """GSE235321 author log2TPM cells; author spermatogonia mapped to U/D (outlier SPG kept as unknown)."""
    objects = []
    for group in ('NOA1', 'NOA2', 'Normal'):
        meta = pd.read_csv(RAW / 'GSE235321' / f'GSE235321_Info.{group}.tsv.gz', sep='\t')
        meta['Cell'] = meta.Cell.str.strip()
        meta = meta.set_index('Cell')
        exp = pd.read_csv(RAW / 'GSE235321' / f'GSE235321_exp.log2TPM.{group}.tsv.gz', sep='\t')
        exp.columns = exp.columns.str.strip()
        assert exp.index.is_unique and exp.columns.is_unique and meta.index.is_unique and set(exp.columns) == set(meta.index)
        exp = exp.loc[:, meta.index]
        assert np.isfinite(exp.values).all() and np.min(exp.values) >= 0
        obs = meta.copy()
        obs['author_cell_type'] = np.where(obs.Stage.str.contains('SPG', regex=False), 'SPG', obs.Stage)
        obs['author_state'] = obs.Stage
        obs['study'] = 'GSE235321'
        obs['donor'] = 'GSE235321_' + obs.Batch.astype(str)
        obs['diagnosis'] = 'Normal' if group == 'Normal' else 'NOA'
        obs['state'] = np.where(obs.Stage.str.startswith('Undiff.SPG'), 'undifferentiated',
                                np.where(obs.Stage == 'Diff.ing SPG', 'differentiating', 'unknown'))
        obs['state_including_author_outlier_sensitivity'] = obs['state'].copy()
        obs['author_outlier'] = obs.Stage.str.startswith('Undiff.SPG-out')
        obs.loc[obs.author_outlier, 'state'] = 'unknown'
        obs['state_annotation_note'] = np.where(obs.author_outlier, 'author_outlier_SPG_retained_but_not_assigned_to_main_state',
                                                'mapped_from_author_state')
        obs['library'] = obs.Batch.astype(str)
        obs['expression_type'] = 'author_log2TPM'
        obs['total_counts'] = np.nan
        obs['author_total_UMI'] = obs.Total_UMI
        obs['n_genes'] = (exp.to_numpy() > 0).sum(axis=0)
        obs['pct_counts_mt'] = np.nan
        obs['qc_pass'] = True
        obs['qc_status'] = 'author_processed_cells; count_based_QC_not_reapplied_to_log2TPM'
        obs['predicted_doublet'] = False
        obs['doublet_tested'] = False
        obs['doublet_score'] = np.nan
        obs.index = pd.Index('GSE235321:' + obs.index, name='cell_id')
        a = ad.AnnData(sparse.csr_matrix(exp.to_numpy(dtype=np.float32).T), obs=obs,
                       var=pd.DataFrame(index=pd.Index(exp.index, name='gene_symbol')))
        a.layers['author_log2TPM'] = a.X.copy()
        objects.append(a)
    allcells = ad.concat(objects, join='inner', merge='first', uns_merge='same')
    spg = allcells[allcells.obs.author_cell_type == 'SPG'].copy()
    spg.write_h5ad(PROCESSED / 'gse235321_spg.h5ad', compression='gzip')
    print(f'GSE235321: {spg.n_obs} author spermatogonia', flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=4):
        for study in LIBRARIES:
            process_counts(study)
        process_tpm()
