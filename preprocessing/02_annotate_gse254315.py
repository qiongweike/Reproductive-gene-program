"""GSE254315 (25 donors, no public cell labels): QC, Scrublet and transfer of compartment and U/T/D state.

Logistic classifiers trained on GSE124263 author labels assign each cell a compartment and each
spermatogonium a state; spermatogonia need p(SPG) >= 0.60 plus multi-marker support, states p >= 0.60.
Inputs: data/raw/GSE254315/, data/processed/normal_allcells_GSE124263.h5ad.
Outputs: data/processed/normal_age_samples/<GSM>_allcells_qc_annotation.tsv.gz, normal_spg_GSE254315.h5ad.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from utils.paths import PROCESSED
from common import SEED, soft_samples, sample_file, read_mtx, collapse_symbols, quality_and_log, scrublet_calls, save

N_VARIABLE = 2500          # most variable reference genes used by the compartment classifier
C_REGULARIZATION = 0.3
P_MIN = 0.60               # minimum class probability for SPG and state calls
PAN = ['DDX4', 'DAZL', 'MAGEA4']
UNDIFF = ['PIWIL4', 'UTF1', 'FGFR3', 'GFRA1', 'NANOS2', 'NANOS3', 'TSPAN33', 'ZBTB16', 'ID4', 'UCHL1', 'SALL4']
DIFF = ['KIT', 'SOHLH1', 'SOHLH2', 'STRA8', 'DMRT1', 'DMRTB1', 'DMRTC2']
STATE_EXTRA = ['CCND2', 'SPRY1', 'MKI67', 'PCNA']
MARKERS = {'SPG': UNDIFF + DIFF, 'Sertoli': ['SOX9', 'WT1', 'AMH', 'CLDN11', 'FSHR'],
           'Leydig': ['INSL3', 'STAR', 'CYP17A1', 'HSD3B2'], 'PTM': ['ACTA2', 'TAGLN', 'MYH11', 'DCN', 'COL1A1'],
           'Immune': ['PTPRC', 'CD3D', 'CD68', 'LST1', 'C1QA'], 'Endothelial': ['VWF', 'PECAM1', 'KDR', 'EMCN'],
           'Spermatocyte': ['SYCP1', 'DMC1', 'MEIOB', 'SPO11', 'TEX101'],
           'Spermatid': ['PRM1', 'PRM2', 'TNP1', 'TNP2', 'ACRV1']}
AUTHOR_LABELS = {'Spermatogonia (SPG) (37%)': 'SPG', 'Peri-Tubular Myoid (PTM) (23%)': 'PTM', 'EC/Blood (16%)': 'EC_Blood',
                 'Leydig Cells (LC) (13%)': 'Leydig', 'Spermatocytes (SPC) (6%)': 'Spermatocyte',
                 'Macrophage (M) (3%)': 'Macrophage', 'Spermatid (St) (2%)': 'Spermatid'}


def classifier(X, y):
    model = Pipeline([('scale', StandardScaler(with_mean=False)),
                      ('lr', LogisticRegression(C=C_REGULARIZATION, max_iter=500, class_weight='balanced',
                                                solver='lbfgs', random_state=SEED))])
    with threadpool_limits(1):
        model.fit(X, y)
    return model


def train_reference_models():
    """Compartment and state classifiers on basic-QC GSE124263 cells, restricted to GSE254315 genes."""
    ref = ad.read_h5ad(PROCESSED / 'normal_allcells_GSE124263.h5ad')
    assay = pd.read_csv(sample_file(soft_samples('GSE254315')[0], 'features.tsv'), sep='\t', header=None)[1].astype(str)
    ref = ref[:, ref.var_names.isin(assay)].copy()
    ref.obs['broad'] = ref.obs.author_cell_type.map(AUTHOR_LABELS)
    keep = ref.obs.qc_basic_pass & ref.obs.broad.notna() & ((ref.obs.broad != 'SPG') | ref.obs.is_author_spg)
    ref = ref[keep].copy()
    mean = np.asarray(ref.X.mean(axis=0)).ravel()
    var = np.asarray(ref.X.power(2).mean(axis=0)).ravel() - mean ** 2
    allowed = np.flatnonzero(~ref.var_names.str.match(r'^(MT-|RP[SL][0-9])'))
    top = allowed[np.argsort(var[allowed])[-N_VARIABLE:]]
    genes = list(dict.fromkeys(list(ref.var_names[top]) + [g for gs in MARKERS.values() for g in gs if g in ref.var_names]
                               + [g for g in PAN if g in ref.var_names]))
    broad = classifier(ref[:, genes].X.tocsr(), ref.obs.broad.astype(str).values)
    is_spg = (ref.obs.broad == 'SPG').to_numpy()
    state_genes = [g for g in UNDIFF + DIFF + STATE_EXTRA if g in ref.var_names]
    state = classifier(ref[is_spg, state_genes].X, ref.obs.state.astype(str).to_numpy()[is_spg])
    return genes, broad, state_genes, state


def score_markers(a):
    """Standardized marker scores per compartment and detected-marker counts for spermatogonial support."""
    panel = [g for g in dict.fromkeys(PAN + [g for x in MARKERS.values() for g in x]) if g in a.var_names]
    L = a[:, panel].X.toarray()
    D = a.layers['counts'][:, a.var_names.get_indexer(panel)].toarray()
    lookup = {g: i for i, g in enumerate(panel)}
    z = (L - L.mean(axis=0)) / np.maximum(L.std(axis=0), 0.25)
    for name, gs in MARKERS.items():
        a.obs['marker_score_' + name] = z[:, [lookup[g] for g in gs if g in lookup]].mean(axis=1)
    for name, gs in [('pan', PAN), ('germline', ['DDX4', 'DAZL']), ('undiff', UNDIFF), ('diff', DIFF)]:
        a.obs['detected_' + name] = (D[:, [lookup[g] for g in gs if g in lookup]] > 0).sum(axis=1)
    n_spg = a.obs.detected_undiff + a.obs.detected_diff
    a.obs['marker_spg_support'] = (n_spg >= 3) & ((a.obs.detected_germline >= 1) | (n_spg >= 4))
    others = ['marker_score_' + k for k in MARKERS if k != 'SPG']
    a.obs['marker_spg_dominant'] = a.obs.marker_score_SPG >= a.obs[others].max(axis=1)
    return a


def annotate_donor(s, models):
    genes, broad, state_genes, state = models
    donor = s['title']
    features = pd.read_csv(sample_file(s, 'features.tsv'), sep='\t', header=None)
    barcodes = pd.read_csv(sample_file(s, 'barcodes.tsv'), sep='\t', header=None)[0].astype(str)
    X = read_mtx(sample_file(s, 'matrix.mtx'))
    assert X.shape == (len(barcodes), len(features)) and np.allclose(X.data, np.round(X.data))
    X, var = collapse_symbols(X.astype(np.float64), features[1].astype(str), features[0].astype(str))
    obs = pd.DataFrame(index=pd.Index(['GSE254315:' + donor + ':' + b for b in barcodes], name='cell_id'))
    obs['study'] = 'GSE254315'
    obs['donor'] = donor
    obs['diagnosis'] = 'Normal'
    obs['age'] = float(s['characteristics']['age (years)'])
    obs['enrichment'] = 'unfractionated'
    obs['library'] = s['gsm']
    obs['gsm'] = s['gsm']
    obs['barcode'] = barcodes.to_numpy()
    obs['count_type'] = 'integer_UMI'
    obs['clinical_description'] = 'cadaveric organ donor; paper reports offspring or normal sperm production'
    obs['author_cell_type'] = 'not_publicly_provided'
    obs['author_state'] = 'not_publicly_provided'
    obs['is_author_spg'] = False
    a = score_markers(quality_and_log(ad.AnnData(X, obs=obs, var=var), age_qc=True))
    good = (a.obs.qc_paper_pass & a.obs.qc_basic_pass).to_numpy()
    scores, doublet = np.full(a.n_obs, np.nan), np.zeros(a.n_obs, dtype=bool)
    scores[good], doublet[good] = scrublet_calls(a.layers['counts'][good].copy())
    a.obs['doublet_score'] = scores
    a.obs['predicted_doublet'] = doublet
    a.obs['doublet_screen_succeeded'] = True
    with threadpool_limits(1):
        prob = broad.predict_proba(a[:, genes].X)
    a.obs['reference_cell_type'] = broad.classes_[prob.argmax(axis=1)]
    a.obs['reference_probability'] = prob.max(axis=1)
    a.obs['reference_spg_probability'] = prob[:, list(broad.classes_).index('SPG')]
    support = ((a.obs.reference_cell_type == 'SPG') & (a.obs.reference_spg_probability >= P_MIN)
               & a.obs.marker_spg_support & a.obs.marker_spg_dominant)
    a.obs['reference_marker_spg'] = support
    a.obs['normal_training_eligible'] = good & support.to_numpy() & ~doublet
    a.obs['state'] = 'unknown'
    a.obs['state_reference_probability'] = np.nan
    if support.any():
        with threadpool_limits(1):
            sp = state.predict_proba(a[support, state_genes].X)
        pred, conf = state.classes_[sp.argmax(axis=1)], sp.max(axis=1)
        pred[conf < P_MIN] = 'unknown'
        a.obs.loc[support, 'state'] = pred
        a.obs.loc[support, 'state_reference_probability'] = conf
    a.obs.to_csv(PROCESSED / 'normal_age_samples' / f"{s['gsm']}_allcells_qc_annotation.tsv.gz", sep='\t')
    return a[a.obs.normal_training_eligible].copy()


if __name__ == '__main__':
    (PROCESSED / 'normal_age_samples').mkdir(parents=True, exist_ok=True)
    models = train_reference_models()
    spg = [annotate_donor(s, models) for s in soft_samples('GSE254315')]
    save(ad.concat(spg, join='inner', merge='same', uns_merge='same'), 'normal_spg_GSE254315.h5ad')
