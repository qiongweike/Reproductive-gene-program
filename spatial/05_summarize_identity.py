"""Bead identities from the MP-excluded RCTD fit (reported weights; HC = singlet or certain doublet; HC-SPG = HC with SPG
weight >= 0.5) and non-MP marker profiles and annotations of the BANKSY lambda-0.2 clusters (family z >= 0.5, lead >= 0.5).
Inputs: slideseq h5ad; results/spatial/neighbourhood/<array>_{classes,doublet_weights}.tsv.gz,
program_scores/<array>_measured_programs.tsv.gz, domains/<array>_Banksy_domains.tsv.gz.
Outputs (results/spatial/neighbourhood/): <array>_{reported_doublet_weights,bead_annotation,domain_scores}.tsv.gz,
<array>_BANKSY_{nonMP_marker_profiles,annotation}.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd

from spatial_common import ARRAYS, DOMAINS, HC_CLASSES, NEIGHBOURHOOD, SCORES, TYPES, log_cp10k, map_genes, programs, read_array, save

SPG_MIN_WEIGHT = 0.5
CLUSTERING = 'clust_M1_lam0.2_k50_res0.5'
# Interpretation markers per family (also excluded from the BANKSY features); MP members among them are dropped here.
MARKER_FAMILIES = {'SPG': ['UTF1', 'FGFR3', 'PIWIL4', 'GFRA1', 'ZBTB16', 'UCHL1'], 'Sertoli': ['SOX9', 'AMH', 'WT1', 'CLU'],
                   'PTM': ['ACTA2', 'TAGLN', 'MYH11'], 'Leydig': ['STAR', 'CYP17A1', 'INSL3'],
                   'Endothelial': ['VWF', 'PECAM1', 'KDR'], 'Immune': ['PTPRC', 'CD3D', 'LST1'],
                   'Spermatocyte': ['SYCP3', 'DMC1', 'SPO11'], 'Spermatid': ['TNP1', 'TNP2', 'PRM1', 'PRM2']}
MIN_FAMILY_MEMBERS = 2
MIN_Z = 0.5
MIN_MARGIN = 0.5


def reported_weights(array):
    """Reported weights per fitted bead (first type unless rejected; second type only for certain doublets; no renormalization)."""
    c = pd.read_csv(NEIGHBOURHOOD / f'{array}_classes.tsv.gz', sep='\t').set_index('barcode')
    raw = pd.read_csv(NEIGHBOURHOOD / f'{array}_doublet_weights.tsv.gz', sep='\t').set_index('barcode').reindex(c.index)
    assert np.allclose(raw[['first_type', 'second_type']].sum(1), 1, atol=1e-6)
    q = pd.DataFrame(0., index=c.index, columns=TYPES)
    for t in TYPES:
        first = c.spot_class.ne('reject') & c.first_type.eq(t)
        second = c.spot_class.eq('doublet_certain') & c.second_type.eq(t)
        q.loc[first, t] += raw.loc[first, 'first_type']
        q.loc[second, t] += raw.loc[second, 'second_type']
    q['unreported_weight_remainder'] = 1 - q[TYPES].sum(1)
    for col in ['spot_class', 'first_type', 'second_type']:
        q[col] = c[col]
    q.to_csv(NEIGHBOURHOOD / f'{array}_reported_doublet_weights.tsv.gz', sep='\t', index_label='barcode')
    return q


def marker_profiles(array, a, x, umi, domain, mp_members):
    rows = []
    for family, requested in MARKER_FAMILIES.items():
        ix, genes = map_genes(a, [g for g in requested if g not in mp_members])
        vals = log_cp10k(x[:, ix].toarray(), umi)
        for cl, ii in domain.groupby(domain).indices.items():
            rows.append(dict(array=array, domain=int(cl), family=family, n_members=len(genes), genes=';'.join(genes),
                             n_beads=len(ii), mean_expression=float(vals[ii].mean()),
                             detection_fraction=float((vals[ii] > 0).mean())))
    m = pd.DataFrame(rows)
    m['z_among_domains'] = m.groupby('family').mean_expression.transform(lambda s: (s - s.mean()) / s.std(ddof=0))
    return m


def annotate_clusters(array, m, domain, scores, umi, hc_spg):
    rows = []
    for cl, part in m.groupby('domain'):
        pp = part[part.n_members >= MIN_FAMILY_MEMBERS].sort_values('z_among_domains', ascending=False)
        top, margin = pp.iloc[0], pp.iloc[0].z_among_domains - pp.iloc[1].z_among_domains
        resolved = top.z_among_domains >= MIN_Z and margin >= MIN_MARGIN
        ii = np.flatnonzero(domain.to_numpy() == cl)
        rows.append(dict(array=array, domain=int(cl), annotation=top.family if resolved else 'Mixed/unresolved',
                         top_marker_family=top.family, top_marker_z=float(top.z_among_domains), top_margin=float(margin),
                         n_beads=len(ii), MP3_mean=float(scores.MP_3.iloc[ii].mean()),
                         MP3_median=float(scores.MP_3.iloc[ii].median()), logUMI_mean=float(np.log1p(umi[ii]).mean()),
                         HC_SPG_fraction=float(hc_spg[ii].mean()), n_HC_SPG=int(hc_spg[ii].sum())))
    return rows


mp_members = set(sum(programs().values(), []))
for array in ARRAYS:
    a, x, xy, umi = read_array(array)
    q = reported_weights(array).reindex(a.obs_names)
    q['spot_class'] = q.spot_class.fillna('outside_RCTD')
    hc = q.spot_class.isin(HC_CLASSES)
    hc_spg = (hc & (q.SPG >= SPG_MIN_WEIGHT)).to_numpy()
    scores = pd.read_csv(SCORES / f'{array}_measured_programs.tsv.gz', sep='\t').set_index('barcode').reindex(a.obs_names)
    assert np.array_equal(scores.total_counts.to_numpy(), umi)
    meta = pd.DataFrame(dict(barcode=a.obs_names, array=array, x_native=xy[:, 0], y_native=xy[:, 1], total_counts=umi,
                             logUMI=np.log1p(umi), MP_3=scores.MP_3, MP3_detected_n=scores.MP_3_detected_n, HC_SPG=hc_spg,
                             spot_class=q.spot_class, display_type=q.first_type.where(hc, q.spot_class)))
    for t in TYPES:
        meta[t] = q[t].to_numpy()
    save(meta, NEIGHBOURHOOD / f'{array}_bead_annotation.tsv.gz')
    domain = pd.read_csv(DOMAINS / f'{array}_Banksy_domains.tsv.gz', sep='\t').set_index('barcode')[CLUSTERING]
    domain = domain.reindex(a.obs_names)
    assert domain.notna().all()
    m = marker_profiles(array, a, x, umi, domain, mp_members)
    save(m, NEIGHBOURHOOD / f'{array}_BANKSY_nonMP_marker_profiles.tsv')
    save(annotate_clusters(array, m, domain, scores, umi, hc_spg), NEIGHBOURHOOD / f'{array}_BANKSY_annotation.tsv')
    meta['domain'] = domain.to_numpy()
    save(meta, NEIGHBOURHOOD / f'{array}_domain_scores.tsv.gz')
    print(array, int(hc.sum()), 'HC beads,', int(hc_spg.sum()), 'HC-SPG,', domain.nunique(), 'clusters')
