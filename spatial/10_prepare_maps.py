"""Bead maps with fixed display windows (Fig 6a-d) and per-cluster summaries for the cluster heatmap (Fig 6e).
Window: 40 x 40 median-NN-unit grid squares anchored at floor(min / 40) * 40; squares with >=400 beads and >=5 HC-SPG beads
are eligible; the shown square is the eligible one closest to the median eligible HC-SPG fraction (ties: lower y, then x).
Cluster summaries: non-MP marker profiles (z among clusters within array) and mean reported RCTD weights over HC beads.
Inputs: results/spatial/neighbourhood/<array>_{domain_scores.tsv.gz,BANKSY_annotation.tsv,BANKSY_nonMP_marker_profiles.tsv}.
Outputs (results/spatial/): summaries/<array>_ROI_maps.tsv.gz, summaries/selected_ROIs.tsv,
cluster_bead_counts.tsv, cluster_marker_profiles.tsv, cluster_RCTD_composition.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd

from spatial_common import ARRAYS, HC_CLASSES, NEIGHBOURHOOD, SPATIAL, SUMMARIES, median_nn, save

WINDOW = 40
MIN_BEADS = 400
MIN_HC_SPG = 5
DISPLAY_TYPES = ['SPG', 'Spermatocyte', 'Spermatid', 'Sertoli', 'Peritubular_smooth_muscle', 'Leydig', 'Endothelial', 'Immune']
DENOMINATOR = 'HC singlet or doublet_certain beads in that fixed cluster'


def select_window(z, array):
    nn = median_nn(z[['x_native', 'y_native']].to_numpy())
    z['x_nn'], z['y_nn'] = z.x_native / nn, z.y_native / nn
    origin = {ax: np.floor(z[f'{ax}_nn'].min() / WINDOW) * WINDOW for ax in 'xy'}
    for ax in 'xy':
        z[f'grid_{ax}'] = np.floor((z[f'{ax}_nn'] - origin[ax]) / WINDOW).astype(int)
    g = z.groupby(['grid_y', 'grid_x'], as_index=False).agg(n_beads=('barcode', 'size'), n_HC_SPG=('HC_SPG', 'sum'))
    g['HC_SPG_fraction'] = g.n_HC_SPG / g.n_beads
    g['eligible'] = (g.n_beads >= MIN_BEADS) & (g.n_HC_SPG >= MIN_HC_SPG)
    median = float(g.loc[g.eligible, 'HC_SPG_fraction'].median())
    g['distance_to_eligible_median'] = abs(g.HC_SPG_fraction - median)
    winner = g[g.eligible].sort_values(['distance_to_eligible_median', 'grid_y', 'grid_x']).iloc[0]
    g['selected'] = (g.grid_x == winner.grid_x) & (g.grid_y == winner.grid_y)
    for ax in 'xy':
        g[f'{ax}min_nn'] = origin[ax] + g[f'grid_{ax}'] * WINDOW
        g[f'{ax}max_nn'] = g[f'{ax}min_nn'] + WINDOW
    g['array'], g['median_NN_native'] = array, nn
    z['in_ROI'] = (z.grid_x == winner.grid_x) & (z.grid_y == winner.grid_y)
    return dict(g[g.selected].iloc[0].to_dict(), eligible_fraction_median=median)


def cluster_summaries(z, array):
    clusters, weights = [], []
    for cl, part in z.groupby('domain', sort=True):
        h = part[part.spot_class.isin(HC_CLASSES)]
        common = dict(array=array, cluster=int(cl), cluster_id=f'{array}.C{int(cl):02}')
        clusters.append(dict(common, n_all=len(part), n_HC=len(h), HC_fraction=len(h) / len(part),
                             n_HC_SPG=int(part.HC_SPG.sum())))
        weights += [dict(common, cell_type=t, n_all=len(part), n_HC=len(h), mean_reported_weight=float(h[t].mean()),
                         sum_reported_weight=float(h[t].sum()), denominator=DENOMINATOR) for t in DISPLAY_TYPES]
    return clusters, weights


windows, clusters, weights, markers = [], [], [], []
for array in ARRAYS:
    z = pd.read_csv(NEIGHBOURHOOD / f'{array}_domain_scores.tsv.gz', sep='\t')
    annotation = pd.read_csv(NEIGHBOURHOOD / f'{array}_BANKSY_annotation.tsv', sep='\t')
    z = z.merge(annotation[['domain', 'annotation']], on='domain', validate='many_to_one')
    windows.append(select_window(z, array))
    z.to_csv(SUMMARIES / f'{array}_ROI_maps.tsv.gz', sep='\t', index=False)
    c, w = cluster_summaries(z, array)
    clusters += c
    weights += w
    m = pd.read_csv(NEIGHBOURHOOD / f'{array}_BANKSY_nonMP_marker_profiles.tsv', sep='\t')
    markers.append(m.assign(cluster_id=[f'{array}.C{d:02}' for d in m.domain]))
    print(array, 'window', windows[-1]['n_beads'], 'beads,', len(c), 'clusters')
save(windows, SUMMARIES / 'selected_ROIs.tsv')
pd.concat(markers, ignore_index=True).to_csv(SPATIAL / 'cluster_marker_profiles.tsv', sep='\t', index=False)
save(weights, SPATIAL / 'cluster_RCTD_composition.tsv')
save(clusters, SPATIAL / 'cluster_bead_counts.tsv')
