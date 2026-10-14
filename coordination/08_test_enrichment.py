"""GO/Reactome over-representation of MP1-MP3 and the MP3 partitions B1-B4 against the 8,875 genes passing the NMF
input filter in all seven discovery donors: pathways of 10-500 background genes, one-sided hypergeometric test, one BH
family over all tests (zero overlaps included).
Inputs: results/programs/{donor_manifest,MP_genes}.tsv, results/programs/nmf/NMF_target_gene_universe.tsv,
        results/coordination/marker_removal/MP3_branches.tsv, results/coordination/pathways/membership.tsv.gz
Outputs: results/coordination/enrichment/{background_genes.tsv,complete_GO_Reactome.tsv.gz}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

MIN_PATHWAY_GENES = 10
MAX_PATHWAY_GENES = 500
# HGNC-verified previous symbols, mapped before intersecting donors, programs and pathways.
ALIASES = {'C12ORF44': 'ATG101', 'FAM26F': 'CALHM6', 'HIST1H1D': 'H1-3', 'FAM207A': 'SLX9', 'C4ORF3': 'ARLN',
           'RP11-325I22.2': 'LINC01511'}


def canonical(x):
    s = str(x).upper()
    return ALIASES.get(s, s)


def main():
    out = out_dir('coordination', 'enrichment')
    donors = pd.read_csv(RESULTS / 'programs' / 'donor_manifest.tsv', sep='\t').query('eligible').sample_key
    universe = pd.read_csv(RESULTS / 'programs' / 'nmf' / 'NMF_target_gene_universe.tsv', sep='\t')
    universe['gene'] = universe.gene.map(canonical)
    background = set.intersection(*[set(universe.loc[universe.sample_key.eq(k), 'gene']) for k in donors])
    pd.DataFrame({'gene': sorted(background)}).to_csv(out / 'background_genes.tsv', sep='\t', index=False)

    mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t')
    branches = pd.read_csv(RESULTS / 'coordination' / 'marker_removal' / 'MP3_branches.tsv', sep='\t')
    sets = {k: set(v.gene.map(canonical)) for k, v in mp.groupby('MP')}
    sets.update({b: set(v.gene) for b, v in branches.groupby('branch') if b != 'unassigned_coverage'})
    membership = pd.read_csv(RESULTS / 'coordination' / 'pathways' / 'membership.tsv.gz', sep='\t')
    membership['gene'] = membership.gene.map(canonical)

    rows = []
    for pathway, v in membership.groupby('pathway'):
        genes = set(v.gene) & background
        if not MIN_PATHWAY_GENES <= len(genes) <= MAX_PATHWAY_GENES:
            continue
        for key, members in sets.items():
            target = members & background
            overlap = genes & target
            p = hypergeom.sf(len(overlap) - 1, len(background), len(genes), len(target))
            kind = 'whole_fixed_program' if key.startswith('MP_') else 'exploratory_unstable_partition'
            rows.append(dict(gene_set=key, pathway=pathway, database=v['database'].iloc[0], name=v['name'].iloc[0],
                             n_overlap=len(overlap), n_target=len(target), n_pathway=len(genes),
                             n_background=len(background), p=float(p), genes=';'.join(sorted(overlap)),
                             interpretation=kind))
    res = pd.DataFrame(rows)
    res['q_BH_all'] = multipletests(res.p, method='fdr_bh')[1]
    res.to_csv(out / 'complete_GO_Reactome.tsv.gz', sep='\t', index=False)
    hits = res[res.gene_set.str.startswith('MP_') & res.q_BH_all.lt(.05)].groupby('gene_set').size().to_dict()
    print(len(background), 'background genes;', len(res), 'tests; q < 0.05:', hits)


if __name__ == '__main__':
    main()
