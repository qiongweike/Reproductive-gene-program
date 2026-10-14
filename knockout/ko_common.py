"""Design constants and helpers shared by the knockout scripts (imported, not run).
HGNC resolution: measured Ensembl ID, then approved symbol, unique previous/alias symbol, unique case-insensitive match.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from utils.paths import EXTERNAL, out_dir

KO = out_dir('knockout')
STUDY_DONORS = {'GSE124263': ['A1', 'A2'], 'GSE153947': ['N2', 'N3']}
SAMPLE_KEYS = [f'{study}__{donor}' for study, donors in STUDY_DONORS.items() for donor in donors]
PAIRED_KEYS = ['GSE124263__A1', 'GSE124263__A2']
CANDIDATE_TARGETS = ['FGFR3', 'IRF2BPL', 'COMMD3', 'STX17']
REFERENCE_TARGETS = ['TCF3']
FOCAL_TARGETS = CANDIDATE_TARGETS + REFERENCE_TARGETS
PAIRED_TARGETS = ['FGFR3', 'COMMD3']
NETWORK_CONFIGS = ['primary', 'seed17', 'down80_seed29']
PAIRED_CONFIGS = ['paired_full_seed1', 'paired_U_seed1']
EARLIER_REFERENCE_READOUT = 'PODXL'  # reference readout of the focal analysis; kept as a forced node and out of control pools
N_CONTROLS = 5
N_MP3 = 49


def hgnc_resolver():
    """Return resolve(symbol, ensembl='') -> (HGNC symbol, HGNC ID, rule); unresolved names are returned unchanged."""
    h = pd.read_csv(EXTERNAL / 'hgnc_complete_set.txt', sep='\t', dtype=str).fillna('')
    h = h[h.status.eq('Approved')]
    approved = dict(zip(h.symbol, h.hgnc_id))
    symbol = dict(zip(h.hgnc_id, h.symbol))
    ens = h[h.ensembl_gene_id.ne('')]
    ens = ens[~ens.ensembl_gene_id.duplicated(keep=False)]
    by_ensembl = dict(zip(ens.ensembl_gene_id, ens.hgnc_id))
    alias, upper = {}, {}
    for z in h.itertuples():
        for a in z.prev_symbol.split('|') + z.alias_symbol.split('|'):
            if a:
                alias.setdefault(a, set()).add(z.hgnc_id)
        for a in [z.symbol] + z.prev_symbol.split('|') + z.alias_symbol.split('|'):
            if a:
                upper.setdefault(a.upper(), set()).add(z.hgnc_id)

    def resolve(gene, ensembl=''):
        if ensembl in by_ensembl:
            return symbol[by_ensembl[ensembl]], by_ensembl[ensembl], 'measured_Ensembl'
        if gene in approved:
            return gene, approved[gene], 'approved_symbol'
        hits = alias.get(gene, set())
        if len(hits) == 1:
            return symbol[next(iter(hits))], next(iter(hits)), 'unique_alias'
        if len(hits) > 1:
            return gene, '', 'ambiguous_alias'
        hits = upper.get(gene.upper(), set())
        if len(hits) == 1:
            return symbol[next(iter(hits))], next(iter(hits)), 'unique_case_alias'
        return gene, '', 'unresolved_preserve_feature'

    return resolve


def fixed_gene_sets():
    """Canonical fixed sets {MP_1, MP_2, MP_3, fixed_143, module_*: sorted genes} written by 01_resolve_identities.py."""
    t = pd.read_csv(KO / 'fixed_gene_sets.tsv', sep='\t', keep_default_na=False)
    return {name: z.gene.tolist() for name, z in t.groupby('set', sort=False)}


def prior_genes(sets):
    """Fixed-set genes, focal targets and the earlier reference readout: forced network nodes, never controls."""
    return set(sum(sets.values(), [])) | set(FOCAL_TARGETS) | {EARLIER_REFERENCE_READOUT}


def all_targets(sets):
    """The 49 MP3 members followed by the focal targets that are not members (primary scan order)."""
    return list(dict.fromkeys(sets['MP_3'] + FOCAL_TARGETS))


def read_tsv(path):
    return pd.read_csv(path, sep='\t')


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA')
