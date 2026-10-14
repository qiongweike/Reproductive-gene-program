"""Unintegrated UMAPs of the adult testis reference (unenriched Normal/OA atlas cells) and of the 9,080 normal spermatogonia.

Batch-aware Seurat-flavour selection of 3,000 variable genes (study batches) without withheld markers or
mitochondrial genes, scaling capped at 10, 50 principal components, 30 neighbours, UMAP min_dist 0.3 with
200 epochs, seed 20261001. Inputs: data/processed/atlas/atlas_allcells_annotated.h5ad,
data/processed/normal_reference_geneID_aligned.h5ad. Outputs: results/atlas/{adult_reference,normal_spg}_umap.tsv.gz.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import scanpy as sc
from threadpoolctl import threadpool_limits

from utils.paths import PROCESSED, out_dir
from common import ATLAS_HELDOUT

SEED = 20261001
N_HVG, N_PCS, N_NEIGHBORS, MIN_DIST, N_EPOCHS, SCALE_MAX = 3000, 50, 30, 0.3, 200, 10
sc.settings.n_jobs = 4


def embed(a, stem):
    sc.pp.highly_variable_genes(a, flavor='seurat', n_top_genes=N_HVG, batch_key='study')
    withheld = set(sum(ATLAS_HELDOUT.values(), []))
    hv = a.var.highly_variable.to_numpy() & ~a.var_names.isin(withheld) & ~a.var_names.str.startswith('MT-')
    h = a[:, hv].copy()
    h.layers.clear()
    sc.pp.scale(h, max_value=SCALE_MAX)
    with threadpool_limits(4):
        sc.tl.pca(h, n_comps=N_PCS, svd_solver='arpack', random_state=SEED)
    a.obsm['X_pca'] = h.obsm['X_pca'].copy()
    del h
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep='X_pca', random_state=SEED)
    sc.tl.umap(a, min_dist=MIN_DIST, random_state=SEED, maxiter=N_EPOCHS)
    o = a.obs.copy()
    o['UMAP_1'], o['UMAP_2'] = a.obsm['X_umap'][:, 0], a.obsm['X_umap'][:, 1]
    o.index.name = 'cell_id'
    o.to_csv(out_dir('atlas') / (stem + '_umap.tsv.gz'), sep='\t')
    print(f'{stem}: {a.n_obs} cells, {int(hv.sum())} variable genes', flush=True)


if __name__ == '__main__':
    A = ad.read_h5ad(PROCESSED / 'atlas' / 'atlas_allcells_annotated.h5ad')
    keep = A.obs.diagnosis.astype(str).isin(['Normal', 'OA']) & ~A.obs.enrichment.astype(str).str.contains('enriched', case=False)
    R = A[keep].copy()
    del A
    R.obs['cell_type'] = R.obs.atlas_cell_type.astype(str)
    embed(R, 'adult_reference')
    del R
    S = ad.read_h5ad(PROCESSED / 'normal_reference_geneID_aligned.h5ad')
    S.obs['donor_key'] = S.obs.study.astype(str) + ':' + S.obs.donor.astype(str)
    S.obs['full_assay_total_counts'] = S.obs.total_counts.astype(float)
    embed(S, 'normal_spg')
