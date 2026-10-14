"""Constants and helpers shared by the spatial Python scripts (imported, not run).
Slide-seq arrays Puck5/Puck6 ("array 1/2"); expression is ln(1 + count / full bead UMI * 1e4).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from utils.paths import PROCESSED, RESULTS

ARRAYS = ['Puck5', 'Puck6']
TYPES = ['Endothelial', 'Immune', 'Leydig', 'Peritubular_smooth_muscle', 'Sertoli', 'Spermatid', 'Spermatocyte', 'SPG']
HC_CLASSES = ['singlet', 'doublet_certain']   # high-confidence RCTD spot classes
# Atlas held-out (second) marker panel: reference QC, withheld from the RCTD fit, identity-marker families of Fig 6i.
HELDOUT_MARKERS = {
    'SPG': ['PIWIL4', 'DMRT1', 'UCHL1', 'NANOS3'],
    'Spermatocyte': ['DMC1', 'REC8', 'STAG3', 'PIWIL1', 'HORMAD1', 'MEI1'],
    'Spermatid': ['TNP2', 'IZUMO1', 'SPATA3', 'TSSK6', 'AKAP4', 'ODF1'],
    'Sertoli': ['GATA4', 'INHBB', 'CTSL', 'SHBG', 'FSHR', 'AR'],
    'Leydig': ['CYP11A1', 'DLK1', 'HSD17B3', 'NR5A1', 'LHCGR'],
    'Peritubular_smooth_muscle': ['ACTG2', 'CNN1', 'MYL9', 'WFDC1'],
    'Endothelial': ['EGFL7', 'RAMP2', 'KDR', 'EMCN', 'ESAM'],
    'Immune': ['C1QB', 'AIF1', 'CD68', 'FCER1G', 'LST1'],
    'Erythroid': ['ALAS2', 'AHSP', 'GYPA', 'SLC4A1', 'HBM']}

EXPORTS = PROCESSED / 'atlas' / 'spatial'          # count matrices and coordinates read by the R steps
RCTD_ORIGINAL = RESULTS / 'spatial' / 'rctd_original'
SCORES = RESULTS / 'spatial' / 'program_scores'
DOMAINS = RESULTS / 'spatial' / 'domains'
NEIGHBOURHOOD = RESULTS / 'spatial' / 'neighbourhood'
SUMMARIES = RESULTS / 'spatial' / 'summaries'
SPATIAL = RESULTS / 'spatial'


def read_array(array):
    """AnnData, integer counts (beads x genes, CSR), native coordinates and full bead UMI of one array."""
    a = ad.read_h5ad(PROCESSED / f'slideseq_{array}.h5ad')
    x = sparse.csr_matrix(a.layers['counts'])
    xy = np.asarray(a.obsm['spatial'])
    umi = np.asarray(x.sum(1)).ravel()
    assert np.array_equal(umi, a.obs.total_counts.to_numpy())
    return a, x, xy, umi


def map_genes(a, genes):
    """Column indices and requested names of the measured genes (unique case-insensitive symbol match)."""
    by = {}
    for j, g in enumerate(a.var_names):
        by.setdefault(str(g).upper(), []).append(j)
    ix, found = [], []
    for g in genes:
        hits = by.get(g.upper(), [])
        if len(hits) > 1:
            raise ValueError(f'ambiguous gene {g}: {hits}')
        if hits:
            ix.append(hits[0])
            found.append(g)
    return np.array(ix, dtype=int), found


def programs(sort=False):
    """Meta-program members {MP_1, MP_2, MP_3} in rank order (or alphabetical with sort=True)."""
    mp = pd.read_csv(RESULTS / 'programs' / 'nmf' / 'MP_genes.tsv', sep='\t')
    return {k: sorted(z.gene) if sort else list(z.gene) for k, z in mp.groupby('MP')}


def log_cp10k(counts, umi):
    return np.log1p(counts / umi[:, None] * 1e4)


def median_nn(xy):
    """Median distance to the nearest other bead (the array's distance unit)."""
    return float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))


def save(rows, path):
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


for d in (EXPORTS, RCTD_ORIGINAL, SCORES, DOMAINS, NEIGHBOURHOOD, SUMMARIES):
    d.mkdir(parents=True, exist_ok=True)
