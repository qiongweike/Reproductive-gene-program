"""Gene identities, the shared measured-gene universe and ln(1 + CP10k) expression for clinical scoring.
Verified equivalent features are summed as counts before normalization; the depth is the full-assay total.
"""
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from utils.paths import PROCESSED, RESULTS

# HGNC-verified historical symbols (upper case) and the entity each one denotes.
ALIASES = {'FAM26F': 'CALHM6', 'HIST1H1D': 'H1-3', 'FAM207A': 'SLX9', 'C4ORF3': 'ARLN', 'RP11-325I22.2': 'LINC01511'}
NORMAL_STUDIES = ['GSE124263', 'GSE153947', 'GSE254315']
UNIVERSE_ONLY = {'GSE235321': 'gse235321_spg.h5ad'}  # contributes its gene list to the shared universe only


def canonical(gene):
    """Upper-case symbol with verified aliases resolved; H1-10 is keyed as H1FX."""
    g = str(gene).upper()
    g = ALIASES.get(g, g)
    return 'H1FX' if g == 'H1-10' else g


def symbol_mappings():
    """Per-study symbol harmonization ({study: {assay symbol: canonical symbol}})."""
    return json.loads((RESULTS / 'preprocessing' / 'gene_symbol_harmonization.json').read_text())['mappings']


def measured_genes():
    """Canonical features measured in each study, whether or not any count was observed."""
    files = {s: f'normal_spg_{s}.h5ad' for s in NORMAL_STUDIES} | UNIVERSE_ONLY
    out = {}
    for study, fn in files.items():
        a = ad.read_h5ad(PROCESSED / fn, backed='r')
        out[study] = set(a.var_names)
        a.file.close()
    presence = pd.read_csv(PROCESSED / 'disease_gene_presence_by_study.tsv', sep='\t', index_col=0)
    for study in presence:
        out[study] = set(presence.index[presence[study].eq(1)])
    for study, aliases in symbol_mappings().items():
        for original, target in aliases.items():
            assert original in out[study] and target not in out[study], (study, original, target)
            out[study] = (out[study] - {original}) | {target}
    return {study: {canonical(g) for g in genes} for study, genes in out.items()}


def log_cp10k(counts, names, depth):
    """Sum count columns sharing a canonical name, then ln(1 + counts * 1e4 / depth); returns (genes, cells x genes)."""
    genes = pd.Index(pd.unique(np.asarray(names, dtype=object)))
    x = sparse.csr_matrix(counts, dtype=float)
    if len(genes) < len(names):
        merge = sparse.csr_matrix((np.ones(len(names)), (np.arange(len(names)), genes.get_indexer(names))),
                                  shape=(len(names), len(genes)))
        x = x @ merge
    x = x.multiply(1e4 / depth[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return genes, x
