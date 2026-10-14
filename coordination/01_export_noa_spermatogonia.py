"""Export full-assay counts of the primary spermatogonia of the external NOA cases (GSE326372).
Feature symbols are upper-cased and harmonized to HGNC symbols; rows that name the same gene are summed.
Inputs: results/atlas/cell_annotations.tsv.gz, data/processed/GSE326372_annotated.h5ad
Outputs: results/coordination/external_inputs/<case>.{counts.mtx.gz,genes.tsv,cells.tsv}, donor_manifest.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import ROOT, PROCESSED, RESULTS, out_dir
import gzip

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import io, sparse

STUDY = 'GSE326372'
# HGNC-verified previous symbols. FAM26F (libraries 1-4) and CALHM6 (library 5) are the same Ensembl
# gene (ENSG00000188820); the two rows never share a cell, so their sum restores one measured gene.
ALIASES = {'C12ORF44': 'ATG101', 'FAM26F': 'CALHM6', 'HIST1H1D': 'H1-3', 'FAM207A': 'SLX9', 'C4ORF3': 'ARLN',
           'RP11-325I22.2': 'LINC01511'}


def symbol_mapper(var_names):
    """Sparse feature-to-gene summation matrix and the gene names (sorted symbols, aliases merged)."""
    upper, inverse = np.unique(np.char.upper(var_names.astype(str)), return_inverse=True)
    genes = list(dict.fromkeys(ALIASES.get(g, g) for g in upper))
    index = {g: j for j, g in enumerate(genes)}
    column = np.array([index[ALIASES.get(g, g)] for g in upper])[inverse]
    mapper = sparse.csr_matrix((np.ones(len(column)), (np.arange(len(column)), column)),
                               shape=(len(column), len(genes)))
    return mapper, genes


def main():
    out = out_dir('coordination', 'external_inputs')
    atlas = pd.read_csv(RESULTS / 'atlas' / 'cell_annotations.tsv.gz', sep='\t', low_memory=False)
    atlas = atlas[atlas.study.eq(STUDY) & atlas.primary_spg].copy()
    assert atlas.qc_analysis_pass.all() and atlas.source_cell_id.is_unique
    with h5py.File(PROCESSED / f'{STUDY}_annotated.h5ad', 'r') as f:
        obs = read_elem(f['obs'])
        mapper, genes = symbol_mapper(read_elem(f['var']).index.to_numpy())
        merged = [np.flatnonzero(mapper[:, j].toarray().ravel()) for j in np.flatnonzero(mapper.sum(0).A1 > 1)]
        counts = sparse_dataset(f['layers/counts'])
        manifest = []
        for donor, cells in atlas.groupby('donor'):
            idx = obs.index.get_indexer(cells.source_cell_id)
            assert (idx >= 0).all()
            x = counts[idx, :].tocsr()
            depth = np.asarray(x.sum(1)).ravel()
            assert np.allclose(depth, cells.full_assay_total_counts, rtol=1e-10, atol=1e-6)
            assert all((x[:, rows] > 0).sum(1).max() <= 1 for rows in merged)
            x = (x @ mapper).tocsr()
            files = {k: out / f'{donor}.{k}' for k in ['counts.mtx.gz', 'genes.tsv', 'cells.tsv']}
            with gzip.open(files['counts.mtx.gz'], 'wb') as fh:
                io.mmwrite(fh, x.T.tocoo(), field='real', precision=17)
            pd.DataFrame({'gene': genes}).to_csv(files['genes.tsv'], sep='\t', index=False)
            cells.assign(state=cells.primary_spg_state, full_assay_counts_recomputed=depth)[
                ['cell_id', 'source_cell_id', 'donor', 'library', 'state', 'full_assay_counts_recomputed']
            ].to_csv(files['cells.tsv'], sep='\t', index=False)
            manifest.append(dict(sample_key=donor, study=STUDY, donor=donor, diagnosis='NOA', n_spg=len(cells),
                                 eligible=True, matrix_file=str(files['counts.mtx.gz'].relative_to(ROOT)),
                                 gene_file=str(files['genes.tsv'].relative_to(ROOT)),
                                 cell_file=str(files['cells.tsv'].relative_to(ROOT))))
    pd.DataFrame(manifest).to_csv(out / 'donor_manifest.tsv', sep='\t', index=False)
    print('NOA cases', len(manifest), 'cells', sum(m['n_spg'] for m in manifest), 'genes', len(genes))


if __name__ == '__main__':
    main()
