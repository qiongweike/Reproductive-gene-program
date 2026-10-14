"""Export the Slide-seq count matrices and coordinates, and build the strict adult RCTD reference.
Reference cells: confidently annotated atlas cells of normal/OA, non-enriched libraries with >=2 held-out markers of
their label detected; compartments with >=25 such cells; at most 2,000 cells per compartment, drawn round-robin across
donors. Held-out markers and the seven SNARE genes are not reference genes; nUMI stays the full-library UMI.
Inputs: data/processed/slideseq_<array>.h5ad, data/processed/atlas/atlas_allcells_annotated.h5ad, results/atlas/cell_annotations.tsv.gz.
Outputs (data/processed/atlas/spatial/): <array>_{counts.mtx,genes.tsv,coords.tsv}.gz, reference_strict_{counts.mtx,genes.tsv,meta.tsv}.gz.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, RESOURCES
import gzip
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

from spatial_common import ARRAYS, EXPORTS, HELDOUT_MARKERS, read_array

SEED = 20261001
REFERENCE_DIAGNOSES = ['Normal', 'OA']
MIN_HELDOUT_DETECTED = 2
CELL_MIN_INSTANCE = 25
MAX_CELLS = 2000
# Outcome genes of an earlier analysis round, withheld from the reference genes in the original design.
WITHHELD_MODULE = 'fusion_SNARE'
ANNOTATION_COLUMNS = ['cell_id', 'atlas_cell_type', 'annotation_confident', 'assignment_markers_detected'] + \
                     [f'heldout_detected_{t}' for t in HELDOUT_MARKERS]


def write_mtx(x, path):
    with gzip.open(path, 'wb') as f:
        io.mmwrite(f, sparse.csc_matrix(x).astype(np.int32), field='integer')


def export_array(array):
    a, x, xy, umi = read_array(array)
    pd.DataFrame({'barcode': a.obs_names, 'x_native': xy[:, 0], 'y_native': xy[:, 1], 'total_counts': umi,
                  'author_cell_type_numeric': a.obs.author_cell_type_numeric.astype(str).to_numpy()}
                 ).to_csv(EXPORTS / f'{array}_coords.tsv.gz', sep='\t', index=False)
    pd.DataFrame({'gene': a.var_names}).to_csv(EXPORTS / f'{array}_genes.tsv.gz', sep='\t', index=False)
    write_mtx(x.T, EXPORTS / f'{array}_counts.mtx.gz')


def balanced_sample(cells):
    """Positions of up to MAX_CELLS cells taken round-robin over donors from fixed per-donor shuffles."""
    rng = np.random.default_rng(SEED)
    queues = {d: list(rng.permutation(ix)) for d, ix in cells.groupby('donor_key', observed=True).indices.items()}
    chosen = []
    while len(chosen) < MAX_CELLS and any(queues.values()):
        for d in sorted(queues):
            if queues[d] and len(chosen) < MAX_CELLS:
                chosen.append(queues[d].pop())
    return np.sort(np.asarray(chosen, dtype=int))


def build_reference():
    A = ad.read_h5ad(PROCESSED / 'atlas' / 'atlas_allcells_annotated.h5ad')
    ann = pd.read_csv(RESULTS / 'atlas' / 'cell_annotations.tsv.gz', sep='\t', usecols=ANNOTATION_COLUMNS,
                      index_col='cell_id')
    assert A.obs_names.equals(ann.index)
    # donor_key keeps its atlas category order, which fixes the order of the per-donor shuffles.
    o = A.obs[['donor_key', 'study', 'diagnosis', 'library', 'enrichment', 'total_counts']].join(ann)
    eligible = (o.diagnosis.astype(str).isin(REFERENCE_DIAGNOSES) & o.annotation_confident
                & ~o.enrichment.astype(str).str.contains('enriched', case=False))
    confirmed = np.zeros(len(o), bool)
    for ct, ix in o.loc[eligible].groupby('atlas_cell_type', observed=True).groups.items():
        ids = o.index.get_indexer(ix)
        good = o[f'heldout_detected_{ct}'].to_numpy()[ids] >= MIN_HELDOUT_DETECTED
        if good.sum() >= CELL_MIN_INSTANCE:
            confirmed[ids[good]] = True
    ids = []
    for _, ix in o.loc[confirmed].groupby('atlas_cell_type', observed=True).groups.items():
        cells = o.loc[ix]
        ids.extend(o.index.get_indexer(cells.index[balanced_sample(cells)]))
    ids = np.sort(ids)
    module = json.loads((RESOURCES / 'excluded_module_genes.json').read_text())[WITHHELD_MODULE]
    withheld = set(sum(HELDOUT_MARKERS.values(), []) + module)
    genes = [g for g in A.var_names if g not in withheld]
    meta = o.iloc[ids][['donor_key', 'study', 'diagnosis', 'library', 'atlas_cell_type', 'total_counts', 'assignment_markers_detected']]
    meta = meta.rename(columns={'atlas_cell_type': 'CellType', 'total_counts': 'full_assay_total_counts'})
    meta.insert(0, 'barcode', meta.index)
    meta.to_csv(EXPORTS / 'reference_strict_meta.tsv.gz', sep='\t', index=False)
    write_mtx(A[ids, genes].layers['counts'].T, EXPORTS / 'reference_strict_counts.mtx.gz')
    pd.DataFrame({'gene': genes}).to_csv(EXPORTS / 'reference_strict_genes.tsv.gz', sep='\t', index=False)
    print('reference:', len(meta), 'cells,', len(genes), 'genes,', meta.CellType.value_counts().to_dict())


if __name__ == '__main__':
    for array in ARRAYS:
        export_array(array)
    build_reference()
