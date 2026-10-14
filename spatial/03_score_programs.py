"""Measured meta-program scores per bead and the bead table of the neighbourhood random forests.
Score: mean ln(1 + count / full bead UMI * 1e4) over the measured members, with no smoothing or imputation.
Model beads: all beads retained and not rejected by the original RCTD fit; four quadrants split at the median coordinates.
Inputs: data/processed/slideseq_<array>.h5ad, results/programs/nmf/MP_genes.tsv,
results/spatial/rctd_original/{weights,classes}_<array>_strict_doublet.tsv.gz.
Outputs: results/spatial/program_scores/<array>_measured_programs.tsv.gz, <array>_model_input.tsv.gz.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd

from spatial_common import ARRAYS, RCTD_ORIGINAL, SCORES, TYPES, log_cp10k, map_genes, median_nn, programs, read_array

MIN_MEASURED_FRACTION = 0.8

for array in ARRAYS:
    a, x, xy, umi = read_array(array)
    d = pd.DataFrame(dict(barcode=a.obs_names, array_id=array, x_native=xy[:, 0], y_native=xy[:, 1],
                          total_counts=umi, log_UMI=np.log1p(umi)))
    for name, genes in programs().items():
        ix, measured = map_genes(a, genes)
        assert len(measured) >= MIN_MEASURED_FRACTION * len(genes), (array, name)
        counts = x[:, ix].toarray()
        d[name] = log_cp10k(counts, umi).mean(1)
        d[name + '_detected_n'] = (counts > 0).sum(1)
    d.to_csv(SCORES / f'{array}_measured_programs.tsv.gz', sep='\t', index=False)
    w = pd.read_csv(RCTD_ORIGINAL / f'weights_{array}_strict_doublet.tsv.gz', sep='\t')
    c = pd.read_csv(RCTD_ORIGINAL / f'classes_{array}_strict_doublet.tsv.gz', sep='\t')
    z = d.merge(w[['barcode'] + TYPES], on='barcode', validate='one_to_one').merge(
        c[['barcode', 'spot_class']], on='barcode', validate='one_to_one')
    z = z[z.spot_class.ne('reject')].copy()
    assert np.allclose(z[TYPES].sum(1), 1, atol=1e-7)
    coords = z[['x_native', 'y_native']].to_numpy()
    centre = np.median(coords, axis=0)
    z['quadrant'] = (coords[:, 0] >= centre[0]).astype(int) + 2 * (coords[:, 1] >= centre[1]).astype(int)
    z['NN_median_native'] = median_nn(coords)
    z.to_csv(SCORES / f'{array}_model_input.tsv.gz', sep='\t', index=False)
    print(array, len(d), 'beads,', len(z), 'model beads')
