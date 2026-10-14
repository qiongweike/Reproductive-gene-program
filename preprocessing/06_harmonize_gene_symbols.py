"""Gene-symbol harmonization for the 66 prior genes, and an ATG101 column for the normal reference.

Missing canonical symbols are matched through the Ensembl IDs of the GSE254315 features; the only alias is
GSE124263 C12orf44 = ATG101. ATG101 counts and ln(1 + CP10k) values are appended to normal_reference.h5ad
from each study's full-assay file, keeping the full-library denominators.
Outputs: results/preprocessing/gene_symbol_harmonization.json, data/processed/normal_reference_geneID_aligned.h5ad.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from utils.paths import PROCESSED, RESOURCES, out_dir
from common import save

GENES = sorted(set(sum(json.loads((RESOURCES / 'excluded_module_genes.json').read_text()).values(), [])))
SOURCES = {'GSE124263': 'normal_spg_GSE124263.h5ad', 'GSE153947': 'normal_spg_GSE153947.h5ad',
           'GSE254315': 'normal_spg_GSE254315.h5ad', 'GSE202647': 'GSE202647_annotated.h5ad',
           'GSE326372': 'GSE326372_annotated.h5ad', 'GSE235321': 'gse235321_spg.h5ad'}


def ensembl_ids(value):
    return set(re.findall(r'ENSG[0-9]+', str(value)))


def find_aliases():
    """{study: {source symbol: canonical symbol}} for prior genes absent under their canonical symbol."""
    variables = {}
    for study, filename in SOURCES.items():
        a = ad.read_h5ad(PROCESSED / filename, backed='r')
        variables[study] = a.var.copy()
        a.file.close()
    expected = {g: ensembl_ids(variables['GSE254315'].loc[g, 'source_gene_ids']) for g in GENES}
    assert len(GENES) == 66 and all(len(ids) == 1 for ids in expected.values())
    aliases = {}
    for study, var in variables.items():
        column = next((c for c in ['source_gene_ids', 'gene_ids'] if c in var), None)
        if column is None:
            assert set(GENES) <= set(var.index), (study, 'symbol-only source lacks a prior gene')
            continue
        lookup = {}
        for symbol, value in var[column].items():
            for gene_id in ensembl_ids(value):
                lookup.setdefault(gene_id, set()).add(symbol)
        for gene in GENES:
            if gene not in var.index:
                matches = sorted(set().union(*(lookup.get(i, set()) for i in expected[gene])))
                assert len(matches) == 1, (study, gene, matches)
                aliases.setdefault(study, {})[matches[0]] = gene
    assert aliases == {'GSE124263': {'C12orf44': 'ATG101'}}, aliases
    return aliases, expected


def append_atg101():
    old = ad.read_h5ad(PROCESSED / 'normal_reference.h5ad')
    assert 'ATG101' not in old.var_names and set(old.layers) == {'counts'}
    counts, x = np.zeros((old.n_obs, 1), dtype=np.float64), np.zeros((old.n_obs, 1), dtype=np.float64)
    for study in old.obs.study.astype(str).unique():
        source = ad.read_h5ad(PROCESSED / SOURCES[study])
        use = np.flatnonzero(old.obs.study.astype(str).eq(study))
        pos = source.obs_names.get_indexer(old.obs_names[use])
        assert (pos >= 0).all()
        j = source.var_names.get_loc('C12orf44' if study == 'GSE124263' else 'ATG101')
        counts[use, 0] = sparse.csr_matrix(source.layers['counts'])[pos, j].toarray().ravel()
        x[use, 0] = sparse.csr_matrix(source.X)[pos, j].toarray().ravel()
    assert np.allclose(x[:, 0], np.log1p(counts[:, 0] * 10000 / old.obs.total_counts.to_numpy()), atol=2e-6, rtol=2e-6)
    var = pd.concat([old.var, pd.DataFrame(index=pd.Index(['ATG101'], name=old.var.index.name))])
    aligned = ad.AnnData(sparse.hstack([old.X, sparse.csr_matrix(x)], format='csr'), obs=old.obs.copy(), var=var)
    aligned.layers['counts'] = sparse.hstack([old.layers['counts'], sparse.csr_matrix(counts)], format='csr')
    save(aligned, 'normal_reference_geneID_aligned.h5ad')


if __name__ == '__main__':
    aliases, expected = find_aliases()
    mapping = dict(scope='Only the 66 prior genes; other assay symbols remain as deposited.',
                   canonical_reference='GSE254315 feature Ensembl IDs, checked against each ID-bearing source',
                   mappings=aliases, stable_ID=next(iter(expected['ATG101'])),
                   limitation='GSE153947 and GSE235321 supply symbols without stable IDs; only exact canonical symbols are checked there.')
    (out_dir('preprocessing') / 'gene_symbol_harmonization.json').write_text(json.dumps(mapping, indent=2) + '\n')
    append_atg101()
