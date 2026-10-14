"""Fig 6h-j and Supplementary Table 6: held-out R2 gains (percentage points) of the neighbourhood random forests, MP3 and
identity-marker detection in HC-SPG vs other HC beads, per-member MP3 detection (49 entities), and ligand-receptor counts:
ligand+ HC beads with somatic weight > 0, receptor+ HC beads with SPG weight > 0, receivers with a sender among their 30-NN.
Inputs: slideseq h5ad, results/programs/nmf/MP_genes.tsv, neighbourhood/<array>_{domain_scores,buffered_RF_predictions}.tsv.gz.
Outputs (results/spatial/summaries/): <array>_neighbourhood_gain.tsv, member_family_detection.tsv,
MP3_member_bead_detection.tsv, ligand_receptor_bead_counts.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from spatial_common import (ARRAYS, HC_CLASSES, HELDOUT_MARKERS, NEIGHBOURHOOD, SUMMARIES, TYPES, map_genes, programs,
                            read_array, save)

K = 30
# MP3 lists two labels of one entity; only RP11-325I22.2 is measured on the arrays and is reported as LINC01511.
MP3_ALIAS = {'RP11-325I22.2': 'LINC01511'}
N_MP3 = 49
SOMATIC = ['Endothelial', 'Immune', 'Leydig', 'Peritubular_smooth_muscle', 'Sertoli']
# The six ligand-receptor pairs of Fig 6j, fixed in an earlier round: (ligand subunits, receptor subunits).
LR_PAIRS = {'FGF2_FGFR3': (['FGF2'], ['FGFR3']),
            'INHABB_ACVR2B': (['INHA', 'INHBB'], ['ACVR2B']),
            'INHBA_ACVR1B_ACVR2A': (['INHBA'], ['ACVR1B', 'ACVR2A']),
            'INHBA_ACVR1B_ACVR2B': (['INHBA'], ['ACVR1B', 'ACVR2B']),
            'INHBB_ACVR1B_ACVR2A': (['INHBB'], ['ACVR1B', 'ACVR2A']),
            'INHBB_ACVR1B_ACVR2B': (['INHBB'], ['ACVR1B', 'ACVR2B'])}
LR_DENOMINATORS = {'Ligand-positive soma': 'HC beads with reported somatic weight >0',
                   'Receptor-positive SPG': 'HC beads with reported SPG weight >0',
                   'Compatible receiver': 'Receptor-positive SPG-weighted HC beads'}
KEYS = ['array', 'identity_version', 'program', 'fold']


def prediction_gains(array):
    p = pd.read_csv(NEIGHBOURHOOD / f'{array}_buffered_RF_predictions.tsv.gz', sep='\t')
    rows = []
    for key, g in p.groupby(KEYS + ['model']):
        y, pred = g.observed.values, g.prediction.values
        r2 = 1 - float(np.sum((y - pred) ** 2)) / float(np.sum((y - y.mean()) ** 2))
        rows.append(dict(zip(KEYS + ['model'], key), R2=r2))
    gain = pd.DataFrame(rows).pivot(index=KEYS, columns='model', values='R2').reset_index()
    gain['gain_pp'] = 100 * (gain.plus_neighbour_composition - gain.composition_depth_coordinates)
    gain.to_csv(SUMMARIES / f'{array}_neighbourhood_gain.tsv', sep='\t', index=False)
    return gain


def knn_graph(xy):
    nn = cKDTree(xy).query(xy, k=K + 1)[1]
    assert np.array_equal(nn[:, 0], np.arange(len(xy)))
    return sparse.csr_matrix((np.ones(len(xy) * K), (np.repeat(np.arange(len(xy)), K), nn[:, 1:].ravel())), shape=(len(xy), len(xy)))


def coverage(array, a, x, z, mp3):
    hc = z.spot_class.isin(HC_CLASSES).to_numpy()
    spg = z.HC_SPG.to_numpy()
    populations = {'Whole array': np.ones(len(z), bool), 'HC-SPG': spg, 'Other HC': hc & ~spg}
    ix, measured = map_genes(a, mp3)
    assert len({MP3_ALIAS.get(g, g) for g in measured}) == len(measured) == N_MP3
    y = x[:, ix].toarray() > 0
    members = [dict(array=array, population=pop, gene=MP3_ALIAS.get(g, g), n_beads=int(m.sum()), n_detected=int(y[m, j].sum()),
                    detected_fraction=float(y[m, j].mean()), listed_symbol=g, canonical_gene=MP3_ALIAS.get(g, g),
                    raw_symbol_absent=False, raw_assay_symbol=g)
               for pop, m in populations.items() for j, g in enumerate(measured)]
    families = []
    for family, genes in dict({t: g for t, g in HELDOUT_MARKERS.items() if t in TYPES}, MP3=mp3).items():
        ix, measured = map_genes(a, genes)
        detected = x[:, ix].toarray() > 0
        for pop in ['HC-SPG', 'Other HC']:
            m = populations[pop]
            families.append(dict(array=array, population=pop, family=family, n_requested=len({MP3_ALIAS.get(g, g) for g in genes}),
                                 n_measured=len(measured), n_beads=int(m.sum()), n_positive=int(detected[m].any(1).sum()),
                                 fraction_positive=float(detected[m].any(1).mean()), mean_member_detection=float(detected[m].mean()),
                                 genes=';'.join(MP3_ALIAS.get(g, g) for g in measured), n_listed_symbols=len(genes)))
    return members, families


def ligand_receptor(array, a, x, xy, z):
    hc = z.spot_class.isin(HC_CLASSES).to_numpy()
    W = knn_graph(xy)[hc][:, hc]
    xh = x[hc]
    sender = z.loc[hc, SOMATIC].sum(1).to_numpy() > 0
    receiver = z.loc[hc, 'SPG'].to_numpy() > 0
    rows = []
    for name, (ligand, receptor) in LR_PAIRS.items():
        (li, lg), (ri, rg) = map_genes(a, ligand), map_genes(a, receptor)
        assert len(lg) == len(ligand) and len(rg) == len(receptor), name
        lp = np.all(xh[:, li].toarray() > 0, axis=1) & sender
        rp = np.all(xh[:, ri].toarray() > 0, axis=1) & receiver
        compatible = rp & (np.asarray(W.dot(lp.astype(int))).ravel() > 0)
        for (measure, definition), n, den in zip(LR_DENOMINATORS.items(), [lp.sum(), rp.sum(), compatible.sum()],
                                                 [sender.sum(), receiver.sum(), rp.sum()]):
            rows.append(dict(array=array, interaction=name, measure=measure, numerator=int(n), denominator=int(den),
                             fraction=n / den if den else np.nan, denominator_definition=definition))
    return rows


mp3 = programs(sort=True)['MP_3']
members, families, lr = [], [], []
for array in ARRAYS:
    a, x, xy, umi = read_array(array)
    z = pd.read_csv(NEIGHBOURHOOD / f'{array}_domain_scores.tsv.gz', sep='\t').set_index('barcode').loc[a.obs_names]
    assert np.array_equal(z[['x_native', 'y_native']].to_numpy(), xy)
    gain = prediction_gains(array)
    m, f = coverage(array, a, x, z, mp3)
    members += m
    families += f
    lr += ligand_receptor(array, a, x, xy, z)
    hc_spg, other = [next(r for r in f if r['family'] == 'MP3' and r['population'] == p) for p in ['HC-SPG', 'Other HC']]
    mean_gain = gain[gain.program.eq('MP_3')].groupby('identity_version').gain_pp.mean().round(6).to_dict()
    print(f"{array}: MP3 detected in {hc_spg['n_positive']}/{hc_spg['n_beads']} HC-SPG beads; mean member detection "
          f"{hc_spg['mean_member_detection']:.2%} vs {other['mean_member_detection']:.2%} (other HC); MP3 gain (pp) {mean_gain}")
save(families, SUMMARIES / 'member_family_detection.tsv')
save(members, SUMMARIES / 'MP3_member_bead_detection.tsv')
save(lr, SUMMARIES / 'ligand_receptor_bead_counts.tsv')
