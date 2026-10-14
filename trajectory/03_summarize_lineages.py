"""Branch qualification and ten-bin fixed-score means for every Slingshot lineage; largest qualified lineage per donor.
Cells: finite pseudotime, curve weight >=0.5. Qualified: >=100 cells, >=30 U and D cells, marker-consistent root, D median
after U, early/late-marker Spearman rho <0/>0. Ten pseudotime-quantile bins; bins with >=10 cells are eligible.
Inputs: results/trajectory/{lineage_cell_pseudotime.tsv.gz, fixed_cell_scores.tsv.gz}.
Outputs: results/trajectory/{lineage_qualification.tsv, lineage_bin_cells.tsv.gz, trajectory_bin_scores.tsv,
         selected_primary_lineages.tsv (largest qualified lineage per donor; ties by lineage name)}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESOURCES, out_dir
import json

import numpy as np
import pandas as pd
from scipy import stats

MIN_CURVE_WEIGHT = 0.5
MIN_CELLS = 100
MIN_STATE_CELLS = 30
N_BINS = 10
MIN_BIN_CELLS = 10
MODULES = list(json.loads((RESOURCES / 'excluded_module_genes.json').read_text()))
METRICS = ['MP_1', 'MP_2', 'MP_3', 'fixed_143_raw', 'fixed_143_adjusted', *MODULES, 'reference_cell_cycle']
OUT = out_dir('trajectory')


def qualify(x, z):
    """Qualification of one lineage; x = all lineage rows, z = cells passing the pseudotime/weight filter."""
    ns = z.state.value_counts()
    reason = []
    if len(z) < MIN_CELLS: reason.append('fewer_than_100_cells')
    if ns.get('undifferentiated', 0) < MIN_STATE_CELLS: reason.append('fewer_than_30_U')
    if ns.get('differentiating', 0) < MIN_STATE_CELLS: reason.append('fewer_than_30_D')
    if not bool(x.root_marker_ok.iloc[0]): reason.append('root_marker_inconsistent')
    rho_e = stats.spearmanr(z.pseudotime, z.early_marker).statistic if len(z) > 2 else np.nan
    rho_d = stats.spearmanr(z.pseudotime, z.late_marker).statistic if len(z) > 2 else np.nan
    med_u = z.loc[z.state.eq('undifferentiated'), 'pseudotime'].median()
    med_d = z.loc[z.state.eq('differentiating'), 'pseudotime'].median()
    if not med_d > med_u: reason.append('D_not_after_U')
    if not (rho_e < 0 and rho_d > 0): reason.append('marker_gradient_inconsistent')
    return dict(n_cells=len(z), n_U=int(ns.get('undifferentiated', 0)), n_D=int(ns.get('differentiating', 0)),
                rho_early_marker=rho_e, rho_late_marker=rho_d, median_U=med_u, median_D=med_d,
                root_cluster=x.root_cluster.iloc[0], root_marker_ok=x.root_marker_ok.iloc[0], eligible=not reason,
                reason=';'.join(reason) if reason else 'PASS', cluster_path=x.cluster_path.iloc[0])


def main():
    score = pd.read_csv(OUT / 'fixed_cell_scores.tsv.gz', sep='\t')
    assert score.cell_id.is_unique
    pt = pd.read_csv(OUT / 'lineage_cell_pseudotime.tsv.gz', sep='\t')
    pt = pt.merge(score, on=['cell_id', 'sample_key'], validate='many_to_one')
    qual, bins, cells = [], [], []
    for (s, sc, lineage), x in pt.groupby(['sample_key', 'scenario', 'lineage'], sort=False):
        z = x[x.pseudotime.notna() & (x.curve_weight >= MIN_CURVE_WEIGHT)].copy()
        q = qualify(x, z)
        qual.append(dict(sample_key=s, scenario=sc, lineage=lineage, **q))
        if not len(z):
            continue
        span = z.pseudotime.max() - z.pseudotime.min()
        z['time_scaled'] = (z.pseudotime - z.pseudotime.min()) / span if span > 0 else np.nan
        z['trajectory_eligible'] = q['eligible']
        z['bin'] = pd.qcut(z.pseudotime.rank(method='first'), N_BINS, labels=False) + 1
        cells.append(z[['cell_id', 'sample_key', 'scenario', 'lineage', 'pseudotime', 'time_scaled', 'curve_weight', 'bin',
                        'trajectory_eligible']])
        for metric in METRICS:
            for b, g in z.groupby('bin'):
                bins.append(dict(sample_key=s, study=z.study.iloc[0], scenario=sc, lineage=lineage, metric=metric, bin=b,
                                 n_cells=len(g), time_scaled=g.time_scaled.mean(), mean_score=g[metric].mean(),
                                 sd_score=g[metric].std(), eligible=q['eligible'], bin_eligible=len(g) >= MIN_BIN_CELLS))
    q = pd.DataFrame(qual)
    q.to_csv(OUT / 'lineage_qualification.tsv', sep='\t', index=False)
    pd.DataFrame(bins).to_csv(OUT / 'trajectory_bin_scores.tsv', sep='\t', index=False)
    pd.concat(cells).to_csv(OUT / 'lineage_bin_cells.tsv.gz', sep='\t', index=False)
    selected = q[q.eligible].sort_values(['sample_key', 'n_cells', 'lineage'], ascending=[True, False, True])
    selected = selected.drop_duplicates('sample_key')
    selected.to_csv(OUT / 'selected_primary_lineages.tsv', sep='\t', index=False)
    print(f'{q.eligible.sum()} of {len(q)} lineages qualified; selected',
          ', '.join(f'{r.sample_key} {r.lineage} ({r.n_cells} cells)' for r in selected.itertuples()))


if __name__ == '__main__':
    main()
