"""Fixed cell scores and an independent-feature PCA for each of the seven discovery donors.
Scores: MP1-3, 143-gene candidate, cell cycle, nine modules, early/late markers. PCA on log1p(CP10k): 2,000 Seurat HVGs
after removing program, candidate, module, cell-cycle, anchor-marker and MT-/RPL/RPS genes; 20 PCs; four K-means clusters.
Inputs: results/programs/{donor_manifest.tsv, MP_genes.tsv, identity/state_markers.tsv, scores/, input/, nmf/input/,
        cell_scores/cell_scores.tsv.gz}, resources/.
Outputs: results/trajectory/{donor_manifest.tsv, fixed_cell_scores.tsv.gz, <donor>.primary.pca.tsv.gz, <donor>.primary.HVG.tsv}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, RESOURCES, out_dir
import json

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

SEED = 20261003
N_HVG = 2000
N_PCS = 20
N_CLUSTERS = 4
N_INIT = 20
EARLY_MARKERS = ['PIWIL4', 'UTF1', 'FSD1', 'GFRA1']
LATE_MARKERS = ['KIT', 'STRA8', 'DMRTB1', 'SOHLH1', 'SOHLH2']
CYCLE_STATES = ['G1/S', 'G2/M']
NON_FEATURE_PREFIXES = ('MT-', 'RPL', 'RPS')
PROGRAMS = RESULTS / 'programs'
OUT = out_dir('trajectory')


def donor_files(sample):
    return (PROGRAMS / 'nmf/input' / f'{sample}.counts.mtx.gz', PROGRAMS / 'input' / f'{sample}.genes.tsv',
            PROGRAMS / 'nmf/input' / f'{sample}.cells.tsv')


def main():
    manifest = pd.read_csv(PROGRAMS / 'donor_manifest.tsv', sep='\t')
    manifest = manifest[manifest.eligible].reset_index(drop=True)
    assert len(manifest) == 7 and manifest.n_spg.sum() == 8139
    markers = pd.read_csv(PROGRAMS / 'identity/state_markers.tsv', sep='\t')
    cycle = set(markers.loc[markers.state.isin(CYCLE_STATES), 'gene'])
    modules = json.loads((RESOURCES / 'excluded_module_genes.json').read_text())
    mp_genes = set(pd.read_csv(PROGRAMS / 'MP_genes.tsv', sep='\t').gene)
    candidates = (RESOURCES / 'candidate_143_genes.txt').read_text().split()
    excluded = {g.upper() for g in [*cycle, *mp_genes, *candidates, *EARLY_MARKERS, *LATE_MARKERS, *sum(modules.values(), [])]}
    fixed = pd.read_csv(PROGRAMS / 'cell_scores/cell_scores.tsv.gz', sep='\t').set_index('cell_id')

    tables = []
    for di, r in manifest.iterrows():
        counts_file, gene_file, cell_file = donor_files(r.sample_key)
        genes = pd.read_csv(gene_file, sep='\t').gene.astype(str)
        cells = pd.read_csv(cell_file, sep='\t')
        c = io.mmread(counts_file).tocsc().astype(np.float64)
        assert c.shape == (len(genes), len(cells))
        lib = np.asarray(c.sum(0)).ravel()
        assert np.allclose(lib, cells.full_assay_counts_recomputed, rtol=1e-10)
        x = c.T.tocsr().multiply(1e4 / lib[:, None]).tocsr()
        x.data = np.log1p(x.data)
        idx = {g: i for i, g in enumerate(genes)}

        def mean_expression(gs):
            gs = [g for g in gs if g in idx]
            return np.asarray(x[:, [idx[g] for g in gs]].mean(1)).ravel() if gs else np.full(len(cells), np.nan)

        out = cells[['cell_id', 'study', 'donor', 'state', 'sample_key']].copy()
        out['full_library_counts'] = lib
        scores = pd.read_csv(PROGRAMS / 'scores' / f'{r.sample_key}_scores.tsv.gz', sep='\t').set_index('cell_id')
        assert set(scores.index) == set(out.cell_id)
        for col in ['MP_1', 'MP_2', 'MP_3']:
            out[col] = scores.loc[out.cell_id, col].to_numpy()
        for col in ['fixed_143_raw', 'fixed_143_adjusted', 'reference_cell_cycle']:
            out[col] = fixed.loc[out.cell_id, col].to_numpy()
        for name, gs in modules.items():
            out[name] = mean_expression(gs)
        out['early_marker'] = mean_expression(EARLY_MARKERS)
        out['late_marker'] = mean_expression(LATE_MARKERS)
        tables.append(out)

        keep = np.array([g.upper() not in excluded and not g.upper().startswith(NON_FEATURE_PREFIXES) for g in genes])
        a = ad.AnnData(x[:, keep].copy())
        a.var_names = genes[keep].values
        sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor='seurat')
        hv = a[:, a.var.highly_variable].copy()
        sc.pp.scale(hv, max_value=10)
        data = np.asarray(hv.X)
        assert np.isfinite(data).all()
        pcs = PCA(n_components=N_PCS, svd_solver='randomized', random_state=SEED + di).fit_transform(data)
        km = KMeans(n_clusters=N_CLUSTERS, n_init=N_INIT, random_state=SEED + di).fit_predict(pcs)
        tab = out[['cell_id', 'study', 'donor', 'state', 'sample_key', 'early_marker', 'late_marker']].copy()
        tab['cluster'] = ['C' + str(k + 1) for k in km]
        for j in range(N_PCS):
            tab[f'PC{j + 1}'] = pcs[:, j]
        tab.to_csv(OUT / f'{r.sample_key}.primary.pca.tsv.gz', sep='\t', index=False)
        pd.DataFrame({'gene': hv.var_names}).to_csv(OUT / f'{r.sample_key}.primary.HVG.tsv', sep='\t', index=False)

    manifest.to_csv(OUT / 'donor_manifest.tsv', sep='\t', index=False)
    pd.concat(tables).to_csv(OUT / 'fixed_cell_scores.tsv.gz', sep='\t', index=False)
    print(f'Prepared {len(manifest)} donors, {sum(map(len, tables))} cells')


if __name__ == '__main__':
    main()
