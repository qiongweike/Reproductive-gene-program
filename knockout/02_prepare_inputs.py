"""Per-donor scTenifold inputs: full-assay counts of the network nodes (2,000 Seurat HVGs plus eligible fixed priors).
Eligible genes: detection > 5%, positive log-CP10k variance, no MT-/RPL/RPS name, unambiguous identity.
GSE124263 donors also get same-node inputs (nodes eligible in U cells too) for all spermatogonia and U cells only.
Inputs: data/processed/normal_spg_<study>.h5ad, results/knockout/<study>_resolved_identity.tsv, fixed_gene_sets.tsv
Outputs: results/knockout/inputs/<key>/{counts.mtx.gz,genes.tsv,cells.tsv,all_gene_features.tsv.gz},
         results/knockout/inputs_paired/<config>/<key>/{counts.mtx.gz,genes.tsv,cells.tsv}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED
import gzip

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse

from ko_common import KO, PAIRED_CONFIGS, PAIRED_KEYS, STUDY_DONORS, fixed_gene_sets, prior_genes

N_HVG = 2000
MIN_DETECTION = 0.05
EXCLUDED_NAMES = r'^(MT-|RP[SL][0-9])'
U_STATE = 'undifferentiated'


def expression_stats(X, depth):
    """Detection fraction, mean CP10k, variance of log1p CP10k (cells x genes counts, full-assay depth) and log1p CP10k."""
    cp = X.multiply((1e4 / depth)[:, None]).tocsr()
    lg = cp.copy()
    lg.data = np.log1p(lg.data)
    mean_log = np.asarray(lg.mean(0)).ravel()
    variance = np.maximum(np.asarray(lg.power(2).mean(0)).ravel() - mean_log ** 2, 0)
    return np.asarray((X > 0).mean(0)).ravel(), np.asarray(cp.mean(0)).ravel(), variance, lg


def save_input(folder, counts, genes, cells):
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename=str(folder / 'counts.mtx.gz'), mode='wb', mtime=0) as f:
        io.mmwrite(f, counts.T, precision=17)
    genes.to_csv(folder / 'genes.tsv', sep='\t', index=False)
    cells.to_csv(folder / 'cells.tsv', sep='\t', index=False)


def network_nodes(identity, X, depth, prior):
    g = identity.copy()
    det, mean, variance, logcp = expression_stats(X, depth)
    g['gene'] = g.canonical
    g['detection_fraction'], g['mean_CP10k'], g['variance_logCP10k'] = det, mean, variance
    g['eligible'] = (det > MIN_DETECTION) & (variance > 0) & ~g.gene.str.match(EXCLUDED_NAMES) & ~g.mapping_ambiguous
    assert not g[g.eligible].gene.duplicated().any()
    hvg = ad.AnnData(logcp[:, g.eligible.to_numpy()], var=pd.DataFrame(index=g[g.eligible].gene))
    sc.pp.highly_variable_genes(hvg, n_top_genes=N_HVG, flavor='seurat')
    g['is_hvg'] = g.gene.isin(set(hvg.var_names[hvg.var.highly_variable])) & g.eligible
    g['prior_forced'] = g.gene.isin(prior)
    g['in_network'] = g.eligible & (g.is_hvg | g.prior_forced)
    return g


def save_paired_inputs(key, g, X, depth, cells, is_u):
    """Same nodes and cell denominators for the Full-SPG and U-only networks."""
    ui = np.flatnonzero(is_u)
    u_det, _, u_var, _ = expression_stats(X[ui, :], depth[ui])
    same = np.flatnonzero(g.in_network.to_numpy() & (u_det > MIN_DETECTION) & (u_var > 0))
    for cfg, cell_index in zip(PAIRED_CONFIGS, [np.arange(X.shape[0]), ui]):
        xx = X[cell_index, :][:, same]
        nodes = g.iloc[same].copy()
        nodes['detection_fraction'], nodes['mean_CP10k'], nodes['variance_logCP10k'], _ = expression_stats(xx, depth[cell_index])
        save_input(KO / 'inputs_paired' / cfg / key, xx, nodes, cells.iloc[cell_index])
        print(cfg, key, len(cell_index), 'cells', len(same), 'nodes')


prior = prior_genes(fixed_gene_sets())
for study, donors in STUDY_DONORS.items():
    a = ad.read_h5ad(PROCESSED / f'normal_spg_{study}.h5ad')
    identity = pd.read_csv(KO / f'{study}_resolved_identity.tsv', sep='\t').fillna('')
    assert np.array_equal(a.var_names.astype(str), identity.original_gene)
    for donor in donors:
        key = f'{study}__{donor}'
        b = a[a.obs.donor.astype(str).eq(donor)].copy()
        X = sparse.csr_matrix(b.layers['counts'], dtype=float)
        depth = np.asarray(X.sum(1)).ravel()
        assert np.allclose(depth, b.obs.total_counts, rtol=1e-10, atol=1e-6)
        g = network_nodes(identity, X, depth, prior)
        nodes = np.flatnonzero(g.in_network)
        cells = b.obs.reset_index()
        cells.rename(columns={cells.columns[0]: 'cell_id'}, inplace=True)
        cells['full_assay_total_counts'] = depth
        save_input(KO / 'inputs' / key, X[:, nodes], g.iloc[nodes].copy(), cells)
        g.to_csv(KO / 'inputs' / key / 'all_gene_features.tsv.gz', sep='\t', index=False)
        print(key, len(cells), 'cells', len(nodes), 'nodes')
        if key in PAIRED_KEYS:
            save_paired_inputs(key, g, X, depth, cells, b.obs.state.astype(str).eq(U_STATE).to_numpy())
