"""Spatial autocorrelation C of every measured meta-program member (Hotspot DANB on raw counts with full bead UMI,
unweighted exact 30-NN graph) in whole arrays and in HC-SPG beads on the induced subgraph of the same whole-array graph
(no new edges; isolated beads kept), and the cross-array agreement of MP3 member C (Fig 6f). Needs hotspotsc 1.1.3.
Inputs: slideseq h5ad, results/programs/nmf/MP_genes.tsv, results/spatial/neighbourhood/<array>_bead_annotation.tsv.gz.
Outputs: results/spatial/domains/<array>_hotspot_autocorrelations.tsv,
results/spatial/neighbourhood/<array>_HC_SPG_Hotspot.tsv and MP3_spatial_member_consistency.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import hotspot
import numpy as np
import pandas as pd
from hotspot.knn import make_weights_non_redundant
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

from spatial_common import ARRAYS, DOMAINS, NEIGHBOURHOOD, map_genes, programs, read_array

K = 30


def danb_hotspot(counts, umi, barcodes, genes, xy):
    aa = ad.AnnData(X=counts.tocsc(), obs=pd.DataFrame({'UMI': umi}, index=barcodes), var=pd.DataFrame(index=genes),
                    obsm={'spatial': xy})
    return hotspot.Hotspot(aa, model='danb', latent_obsm_key='spatial', umi_counts_obs_key='UMI')


def whole_array(a, x, xy, umi):
    ix, names = map_genes(a, list(dict.fromkeys(sum(programs().values(), []))))
    xx = x[:, ix]
    detected = np.asarray(xx.max(0).toarray()).ravel() > 0
    hs = danb_hotspot(xx[:, detected], umi, a.obs_names, np.array(names)[detected].tolist(), xy)
    hs.create_knn_graph(weighted_graph=False, n_neighbors=K, approx_neighbors=False)
    result = hs.compute_autocorrelations(jobs=1)
    result.index.name = 'gene'
    return result


def hc_spg(a, x, xy, umi, keep):
    """Hotspot on the HC-SPG induced subgraph; unused neighbour slots point to the bead itself with zero weight."""
    nn = NearestNeighbors(n_neighbors=K, algorithm='ball_tree').fit(xy).kneighbors(return_distance=False)
    W = sparse.csr_matrix((np.ones(nn.size), (np.repeat(np.arange(len(xy)), K), nn.ravel())), shape=(len(xy), len(xy)))
    inds = np.flatnonzero(keep)
    S = W[inds][:, inds].tocsr()
    n = len(inds)
    nb = np.tile(np.arange(n)[:, None], (1, K))
    wt = np.zeros((n, K))
    for i in range(n):
        js = S.indices[S.indptr[i]:S.indptr[i + 1]]
        nb[i, :len(js)] = js
        wt[i, :len(js)] = 1.
    wn = make_weights_non_redundant(nb, wt)
    assert np.all(wn[nb == np.arange(n)[:, None]] == 0)
    ix, names = map_genes(a, list(dict.fromkeys(sum(programs(sort=True).values(), []))))
    xx = x[inds][:, ix]
    variable = np.asarray(xx.max(0).toarray()).ravel() > np.asarray(xx.min(0).toarray()).ravel()
    hs = danb_hotspot(xx[:, variable], umi[inds], a.obs_names[inds], np.array(names)[variable].tolist(), xy[inds])
    hs.neighbors = pd.DataFrame(nb, index=a.obs_names[inds])
    hs.weights = pd.DataFrame(wn, index=a.obs_names[inds])
    result = hs.compute_autocorrelations(jobs=1).reindex(names)
    result['estimable'] = variable
    result.index.name = 'gene'
    return result


for array in ARRAYS:
    a, x, xy, umi = read_array(array)
    meta = pd.read_csv(NEIGHBOURHOOD / f'{array}_bead_annotation.tsv.gz', sep='\t').set_index('barcode').loc[a.obs_names]
    whole = whole_array(a, x, xy, umi)
    whole.to_csv(DOMAINS / f'{array}_hotspot_autocorrelations.tsv', sep='\t')
    hc = hc_spg(a, x, xy, umi, meta.HC_SPG.to_numpy())
    hc.to_csv(NEIGHBOURHOOD / f'{array}_HC_SPG_Hotspot.tsv', sep='\t')
    print(array, len(whole), 'members tested;', int(hc.estimable.sum()), 'estimable in', int(meta.HC_SPG.sum()), 'HC-SPG beads')

# Cross-array agreement of MP3 member C, assembled from the per-array tables written above. Population labels
# as read by the Fig 6f figure code: whole array, HC-SPG beads on the induced original graph.
TABLES = {'whole_array': DOMAINS / '{}_hotspot_autocorrelations.tsv',
          'HC_SPG': NEIGHBOURHOOD / '{}_HC_SPG_Hotspot.tsv'}
mp3 = programs(sort=True)['MP_3']
rows = []
for population, table in TABLES.items():
    out = pd.concat([pd.read_csv(str(table).format(p), sep='\t').set_index('gene').reindex(mp3)[['C', 'FDR']]
                     .rename(columns=lambda c: f'{p}_{c}') for p in ARRAYS], axis=1)
    out.index.name = 'gene'
    out['population'] = population
    ok = out.Puck5_C.notna() & out.Puck6_C.notna()
    rho = spearmanr(out.Puck5_C[ok], out.Puck6_C[ok]).statistic
    print(population, f'MP3 member C, cross-array Spearman rho {rho:.6f} over {int(ok.sum())} members')
    rows.append(out.reset_index())
pd.concat(rows).to_csv(NEIGHBOURHOOD / 'MP3_spatial_member_consistency.tsv', sep='\t', index=False)
