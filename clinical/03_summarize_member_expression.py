"""Measured expression and detection of all 49 MP3 members in the GSE202647 undifferentiated donor groups (Fig 4c).
Inputs: data/processed/disease_spg.h5ad, results/programs/MP_genes.tsv, results/atlas/cell_annotations.tsv.gz,
results/clinical/donor_state_scores.tsv (same donor groups as Fig 4b).
Output: results/clinical/expression/MP3_member_expression_U.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, out_dir
import anndata as ad
import numpy as np
import pandas as pd

from assay_genes import canonical, log_cp10k

STUDY = 'GSE202647'
STATE = 'undifferentiated'
MIN_CELLS = 30            # primary donor-state threshold; smaller groups are labelled sensitivity
N_MEMBERS = 49

mp = pd.read_csv(RESULTS / 'programs' / 'MP_genes.tsv', sep='\t').query('MP == "MP_3"')
members = mp.sort_values('rank').gene.tolist()
assert len(set(members)) == N_MEMBERS
a = ad.read_h5ad(PROCESSED / 'disease_spg.h5ad')
ann = pd.read_csv(RESULTS / 'atlas' / 'cell_annotations.tsv.gz', sep='\t',
                  usecols=['study', 'source_cell_id', 'atlas_cell_type'])
atlas_spg = set(ann.loc[ann.study.eq(STUDY) & ann.atlas_cell_type.eq('SPG'), 'source_cell_id'])
use = (a.obs.study.astype(str).eq(STUDY) & a.obs.state.astype(str).eq(STATE) & a.obs_names.isin(atlas_spg)).to_numpy()
obs = a.obs[use]
names = pd.Index([canonical(g) for g in a.var_names])
# Depth stays in its stored precision (float32 in this file), as for the published values.
genes, x = log_cp10k(a.layers['counts'][use], names, pd.to_numeric(obs.total_counts).to_numpy())

rows = []
for donor, z in obs.groupby('donor', observed=True):
    cells = np.flatnonzero(obs.donor.astype(str).eq(str(donor)))
    for order, g in enumerate(members, 1):
        j = genes.get_indexer([g])[0]
        v = x[cells, j].toarray().ravel() if j >= 0 else None
        rows.append(dict(donor=str(donor), diagnosis=str(z.diagnosis.iloc[0]), gene=g, canonical_gene=g,
                         assay_gene=';'.join(a.var_names[names == g]), member_order=order, n_cells=len(cells),
                         qualification='primary' if len(cells) >= MIN_CELLS else 'sensitivity', measured=j >= 0,
                         mean_log1p_CP10k=float(v.mean()) if v is not None else np.nan,
                         detection_fraction=float((v > 0).mean()) if v is not None else np.nan))
table = pd.DataFrame(rows)
assert len(table) == N_MEMBERS * obs.donor.nunique()
scored = pd.read_csv(RESULTS / 'clinical' / 'donor_state_scores.tsv', sep='\t').query(
    'study == @STUDY and scenario == "atlas_SPG_agreement" and state == @STATE')
assert table.groupby('donor').n_cells.first().to_dict() == scored.set_index('donor').n_cells.to_dict()
table.to_csv(out_dir('clinical', 'expression') / 'MP3_member_expression_U.tsv', sep='\t', index=False, na_rep='NA')
print('member expression:', len(table), 'rows,', int((~table.measured).sum()), 'unmeasured')
