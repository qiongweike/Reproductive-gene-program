"""Resolve assay features and the fixed gene sets to current HGNC identities.
Several features with one identity: only a unique best-mapped feature is kept; tied features are all flagged ambiguous.
Inputs: data/processed/normal_spg_<study>.h5ad, data/external/hgnc_complete_set.txt, results/programs/MP_genes.tsv,
        resources/candidate_143_genes.txt, resources/excluded_module_genes.json
Outputs: results/knockout/<study>_resolved_identity.tsv, results/knockout/fixed_gene_sets.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, RESOURCES
import json

import anndata as ad
import pandas as pd

from ko_common import KO, N_MP3, STUDY_DONORS, hgnc_resolver

PRIORITY = {'measured_Ensembl': 0, 'approved_symbol': 1, 'unique_alias': 2, 'unique_case_alias': 3,
            'unresolved_preserve_feature': 5, 'ambiguous_alias': 9}

resolve = hgnc_resolver()


def resolve_features(study):
    a = ad.read_h5ad(PROCESSED / f'normal_spg_{study}.h5ad', backed='r')
    genes = a.var_names.astype(str).tolist()
    ensembl = a.var.source_gene_ids.fillna('').astype(str).tolist() if 'source_gene_ids' in a.var else [''] * len(genes)
    a.file.close()
    t = pd.DataFrame([resolve(g, e) for g, e in zip(genes, ensembl)], columns=['canonical', 'hgnc_id', 'mapping'])
    t.insert(0, 'original_gene', genes)
    t.insert(1, 'source_ensembl', ensembl)
    t['canonical_collision'] = t.canonical.duplicated(keep=False)
    t['identity_priority'] = t.mapping.map(PRIORITY)
    by_identity = t.groupby('canonical')
    best = t.identity_priority.eq(by_identity.identity_priority.transform('min'))
    unique_best = best & best.groupby(t.canonical).transform('sum').eq(1)
    t['mapping_ambiguous'] = t.canonical_collision & ~unique_best
    return t


for study in STUDY_DONORS:
    t = resolve_features(study)
    t.to_csv(KO / f'{study}_resolved_identity.tsv', sep='\t', index=False)
    print(study, len(t), 'features;', int(t.hgnc_id.ne('').sum()), 'HGNC;', int(t.mapping_ambiguous.sum()), 'ambiguous')

mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t')
sets = {k: mp.gene[mp.MP.eq(k)].tolist() for k in ['MP_1', 'MP_2', 'MP_3']}
sets['fixed_143'] = (RESOURCES / 'candidate_143_genes.txt').read_text().split()
sets.update({f'module_{k}': v for k, v in json.loads((RESOURCES / 'excluded_module_genes.json').read_text()).items()})
canonical = {k: sorted({resolve(g)[0] for g in v}) for k, v in sets.items()}
assert len(canonical['MP_3']) == N_MP3 and len(canonical['fixed_143']) == 143
pd.DataFrame([(k, g) for k, v in canonical.items() for g in v], columns=['set', 'gene']).to_csv(
    KO / 'fixed_gene_sets.tsv', sep='\t', index=False)
