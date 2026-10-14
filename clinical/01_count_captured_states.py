"""Captured spermatogonial state composition per GSE202647 donor library (Fig 4a).
Input: results/atlas/cell_annotations.tsv.gz (QC-passing atlas cells classified as SPG).
Output: results/clinical/capture/GSE202647_SPG_states.tsv with every state label, including unknown
and not_primary_SPG (the figure keeps U/T/D and renormalizes).
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir
import pandas as pd

STUDY = 'GSE202647'
GROUP = ['study', 'donor', 'donor_key', 'library', 'enrichment']
IDENTITY_LAYER = 'QC-passing annotation class; not the communication-analysis intersection'

meta = pd.read_csv(RESULTS / 'atlas' / 'cell_annotations.tsv.gz', sep='\t',
                   usecols=['cell_id', *GROUP, 'qc_analysis_pass', 'atlas_cell_type', 'state'])
meta = meta[meta.study.eq(STUDY)]
assert meta.cell_id.is_unique
spg = meta[meta.qc_analysis_pass & meta.atlas_cell_type.eq('SPG')]
states = spg.groupby([*GROUP, 'state'], dropna=False).size().rename('n_cells').reset_index()
states['SPG_denominator'] = states.groupby(['donor_key', 'library'], dropna=False).n_cells.transform('sum')
states['captured_state_fraction'] = states.n_cells / states.SPG_denominator
states['identity_layer'] = IDENTITY_LAYER
states.to_csv(out_dir('clinical', 'capture') / 'GSE202647_SPG_states.tsv', sep='\t', index=False, na_rep='NA')
print('capture states:', len(states), 'rows,', states.donor.nunique(), 'donors')
