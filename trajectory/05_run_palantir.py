"""Palantir 1.4.4 on the saved independent-feature PCA of each donor with a label- and marker-consistent root.
Start cell: the U cell of the Slingshot root cluster nearest to the 20-PC centroid of those cells (ties by cell_id).
Diffusion maps: 10 components, 30 neighbours, alpha 0. Palantir: automatic terminal states, 30 neighbours,
min(1,200, donor cells) waypoints, 25 iterations, start cell forced, seed 20261003.
Inputs: results/trajectory/{donor_manifest.tsv, root_clusters.tsv, <donor>.primary.pca.tsv.gz}.
Outputs: results/trajectory/palantir/{donor_start_eligibility.tsv, <donor>.cell_results.tsv.gz}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir

import pandas as pd
import palantir

SEED = 20261003
PCS = [f'PC{i}' for i in range(1, 21)]
MAX_WAYPOINTS = 1200
DIFFUSION = dict(n_components=10, knn=30, alpha=0, seed=SEED, kernel_backend='scanpy')
PALANTIR = dict(terminal_states=None, knn=30, n_jobs=1, scale_components=True, use_early_cell_as_start=True,
                max_iterations=25, seed=SEED)
TRAJ = RESULTS / 'trajectory'
OUT = out_dir('trajectory', 'palantir')


def start_cell(z, root):
    """U cell of the root cluster closest to the centroid of those cells in the 20-PC space."""
    a = z[(z.cluster == root) & z.state.eq('undifferentiated')]
    pc = a[PCS].to_numpy()
    d = ((pc - pc.mean(0)) ** 2).sum(1)
    pick = pd.DataFrame(dict(cell_id=a.cell_id.to_numpy(), squared_distance=d)).sort_values(['squared_distance', 'cell_id']).iloc[0]
    return str(pick.cell_id), float(pick.squared_distance), len(a)


def main():
    manifest = pd.read_csv(TRAJ / 'donor_manifest.tsv', sep='\t')
    roots = pd.read_csv(TRAJ / 'root_clusters.tsv', sep='\t')
    roots = roots[roots.chosen_root].set_index('sample_key')
    starts = []
    for m in manifest.itertuples():
        z = pd.read_csv(TRAJ / f'{m.sample_key}.primary.pca.tsv.gz', sep='\t').set_index('cell_id', drop=False)
        r = roots.loc[m.sample_key]
        eligible = bool(r.root_label_ok and r.root_marker_ok)
        early, dist, n_candidates = start_cell(z, r.cluster) if eligible else ('', None, 0)
        starts.append(dict(sample_key=m.sample_key, study=m.study, donor=m.donor, diagnosis=m.diagnosis, n_cells=len(z),
                           root_cluster=r.cluster, root_label_ok=bool(r.root_label_ok), root_marker_ok=bool(r.root_marker_ok),
                           run_eligible=eligible, root_U_candidates=n_candidates, early_cell=early,
                           root_centroid_squared_distance=dist,
                           reason='PASS' if eligible else ('ROOT_LABEL_INELIGIBLE' if not r.root_label_ok else 'ROOT_MARKER_INELIGIBLE')))
        if not eligible:
            continue
        pcs = z[PCS]
        dm = palantir.utils.run_diffusion_maps(pcs, **DIFFUSION)
        ms = palantir.utils.determine_multiscale_space(dm, n_eigs=None)
        pr = palantir.core.run_palantir(ms, early_cell=early, num_waypoints=min(MAX_WAYPOINTS, len(pcs)), **PALANTIR)
        out = z[['cell_id', 'sample_key', 'state', 'cluster', 'early_marker', 'late_marker']].copy()
        out['pseudotime'] = pr.pseudotime.reindex(z.index)
        out['entropy'] = pr.entropy.reindex(z.index)
        out['is_actual_start'] = out.cell_id.eq(early)
        out.to_csv(OUT / f'{m.sample_key}.cell_results.tsv.gz', sep='\t', index=False)
    starts = pd.DataFrame(starts)
    starts.to_csv(OUT / 'donor_start_eligibility.tsv', sep='\t', index=False)
    print(f'Palantir fitted {starts.run_eligible.sum()} donors; skipped', ', '.join(starts.loc[~starts.run_eligible, 'sample_key']))


if __name__ == '__main__':
    main()
