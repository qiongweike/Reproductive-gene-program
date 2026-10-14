"""Spermatogonial identity and U/T/D states in GSE202647 and GSE326372, transferred from the GSE124263 reference.

SPG: p(SPG) >= 0.6, >= 3 detected SPG markers (plus DDX4/DAZL or a fourth), SPG or SPG-majority cluster without
somatic conflict; state: p >= 0.6 and margin >= 0.1; GSE202647 NOA3 (no germ cells per the authors) is excluded.
Outputs (data/processed): <study>_annotated.h5ad, disease_allcells.h5ad, disease_spg.h5ad,
disease_gene_presence_by_library.tsv.gz, disease_gene_presence_by_study.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from utils.paths import PROCESSED
from common import SEED, DISEASE_PANELS, DISEASE_STATE_MARKERS, save

STUDIES = ('GSE202647', 'GSE326372')
P_MIN, STATE_MARGIN = 0.6, 0.1
MIN_SPG_MARKERS, CLUSTER_SPG_FRACTION, LATE_FLAG = 3, 0.5, 0.10
AMBIENT_MARKERS = {'PRM1', 'PRM2', 'TNP1', 'TNP2', 'AKAP4', 'ODF1'}       # never used as identity features
EXTRA_LINEAGE_GENES = {'DDX4', 'DAZL', 'GATA1', 'TPSAB1', 'CD79A'}
SOMATIC = {'Sertoli', 'Leydig', 'Peritubular', 'Endothelial', 'Immune'}
REFERENCE_LABELS = {'Spermatogonia (SPG) (37%)': 'SPG', 'Spermatocytes (SPC) (6%)': 'Spermatocyte',
                    'Spermatid (St) (2%)': 'Spermatid', 'Peri-Tubular Myoid (PTM) (23%)': 'Peritubular',
                    'Leydig Cells (LC) (13%)': 'Leydig', 'Macrophage (M) (3%)': 'Immune', 'EC/Blood (16%)': 'EC_Blood'}
GSE202647_AGES = {'OA1': 26, 'OA2': 30, 'NOA1': 32, 'NOA2': 29, 'NOA3': 26}   # eLife 97958 supplementary file 5


def dense_features(a, genes):
    matrix = a[:, genes].X
    return matrix.toarray().astype(np.float64) if sparse.issparse(matrix) else np.asarray(matrix, dtype=np.float64)


def fit_classifier(a, labels, genes):
    model = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=3000, solver='lbfgs',
                                                               class_weight='balanced', random_state=SEED))
    return model.fit(dense_features(a, genes), np.asarray(labels, dtype=str))


def reference_training_sets():
    """GSE124263 basic-QC singlets: compartment labels (author SPG override) and author U/T/D states."""
    ref = ad.read_h5ad(PROCESSED / 'normal_allcells_GSE124263.h5ad')
    broad = ref.obs.author_cell_type.astype(str).map(REFERENCE_LABELS)
    assert broad.notna().all(), 'Unmapped reference label'
    broad.loc[ref.obs.is_author_spg.astype(bool)] = 'SPG'
    uncertain = (broad == 'SPG') & ~ref.obs.is_author_spg.astype(bool)
    singlet = ~ref.obs.predicted_doublet.astype(bool)
    valid = ref.obs.qc_basic_pass.astype(bool) & ~uncertain & singlet
    lineage = ref[valid].copy()
    lineage.obs['training_label'] = broad[valid].to_numpy()
    state = ref[ref.obs.is_author_spg.astype(bool) & ref.obs.qc_basic_pass.astype(bool) & singlet & (ref.obs.state != 'unknown')].copy()
    return lineage, state


def annotate(a, lineage_model, lineage_genes, state_model, state_genes):
    study = str(a.obs.study.iloc[0])
    a.obs['pathology'] = a.obs['pathology'].astype(str)
    probabilities = lineage_model.predict_proba(dense_features(a, lineage_genes))
    a.obs['reference_cell_type'] = lineage_model.classes_[probabilities.argmax(axis=1)]
    a.obs['reference_confidence'] = probabilities.max(axis=1)
    a.obs['reference_SPG_probability'] = probabilities[:, list(lineage_model.classes_).index('SPG')]
    spg_positive = a.obs.reference_cell_type == 'SPG'
    cluster_fraction = spg_positive[a.obs.qc_pass].groupby(a.obs.loc[a.obs.qc_pass, 'leiden'], observed=True).mean()
    a.obs['cluster_reference_SPG_fraction'] = a.obs.leiden.map(cluster_fraction).fillna(0).astype(float)
    reference_counts = pd.crosstab(a.obs.loc[a.obs.qc_pass, 'leiden'], a.obs.loc[a.obs.qc_pass, 'reference_cell_type'])
    reference_majority = reference_counts.idxmax(axis=1)
    reference_fraction = reference_counts.max(axis=1) / reference_counts.sum(axis=1)
    a.obs['cluster_reference_majority'] = a.obs.leiden.map(reference_majority).fillna('QC_excluded').astype(str)
    a.obs['cluster_reference_majority_fraction'] = a.obs.leiden.map(reference_fraction).fillna(0).astype(float)
    germline = [g for g in ('DDX4', 'DAZL') if g in a.var_names]
    a.obs['detected_germline'] = np.asarray((a.layers['counts'][:, a.var_names.get_indexer(germline)] > 0).sum(axis=1)).ravel()
    a.obs['SPG_marker_supported'] = (a.obs.detected_SPG >= MIN_SPG_MARKERS) & ((a.obs.detected_germline >= 1) | (a.obs.detected_SPG >= MIN_SPG_MARKERS + 1))
    a.obs['SPG_cluster_supported'] = (a.obs.cluster_cell_type == 'SPG') | (a.obs.cluster_reference_SPG_fraction >= CLUSTER_SPG_FRACTION)
    a.obs['somatic_cluster_conflict'] = a.obs.cluster_cell_type.isin(SOMATIC)
    a.obs['is_spg_reference_candidate'] = (a.obs.reference_SPG_probability >= P_MIN) & a.obs.SPG_marker_supported
    a.obs['is_spg'] = a.obs.qc_pass & a.obs.is_spg_reference_candidate & a.obs.SPG_cluster_supported & ~a.obs.somatic_cluster_conflict
    a.obs['annotation_disagreement'] = a.obs.is_spg_reference_candidate & (~a.obs.SPG_cluster_supported | a.obs.somatic_cluster_conflict)
    a.obs['late_background_flag'] = a.obs.late_rna_fraction > LATE_FLAG
    # Cell type: marker-defined cluster label, replaced by the reference majority where the cluster agrees
    # (the reference has no Sertoli cells, so Sertoli clusters keep their marker label).
    a.obs['cell_type'] = a.obs.cluster_cell_type.astype(str)
    a.obs['cell_type_basis'] = 'independent_marker_cluster'
    consensus = ((a.obs.cluster_reference_majority_fraction >= 0.5) & (a.obs.reference_confidence >= P_MIN)
                 & (a.obs.cluster_cell_type != 'Sertoli') & a.obs.qc_pass)
    a.obs.loc[consensus, 'cell_type'] = a.obs.loc[consensus, 'cluster_reference_majority'].to_numpy()
    a.obs.loc[consensus, 'cell_type_basis'] = 'normal_reference_majority_plus_disease_cluster'
    a.obs.loc[(a.obs.cell_type == 'SPG') & ~a.obs.is_spg, 'cell_type'] = 'uncertain_SPG'
    a.obs.loc[a.obs.is_spg, 'cell_type'] = 'SPG'
    a.obs.loc[a.obs.is_spg, 'cell_type_basis'] = 'reference_plus_cluster_plus_multiple_markers'
    state_prob = state_model.predict_proba(dense_features(a, state_genes))
    ranked = np.sort(state_prob, axis=1)
    a.obs['reference_state'] = state_model.classes_[state_prob.argmax(axis=1)]
    a.obs['state_confidence'] = state_prob.max(axis=1)
    a.obs['state_margin'] = ranked[:, -1] - ranked[:, -2]
    a.obs['state'] = 'unknown'
    assign = a.obs.is_spg & (a.obs.state_confidence >= P_MIN) & (a.obs.state_margin >= STATE_MARGIN)
    a.obs.loc[assign, 'state'] = a.obs.loc[assign, 'reference_state'].to_numpy()
    a.obs['annotation_basis'] = 'disease_clusters+lineage_markers+independent_normal_author_reference'
    a.obs['age'] = np.nan
    a.obs['qc_author_thresholds'] = False
    a.obs['qc_author_thresholds_available'] = study == 'GSE202647'
    a.obs['donor_identity_caveat'] = 'sample ID used as donor proxy; individual age/pathology mapping unavailable'
    if study == 'GSE202647':
        a.obs['donor_identity_caveat'] = 'five study samples corresponding to two OA and three NOA biopsies'
        a.obs['age'] = a.obs['sample'].map(GSE202647_AGES).astype(float)
        a.obs['clinical_metadata_source'] = 'elife-97958-supp5-v3.docx patient-characteristics table'
        a.obs['AZF_screen'] = 'Normal'
        a.obs['qc_author_thresholds'] = (a.obs.n_genes >= 50) & (a.obs.total_counts <= 20000) & (a.obs.pct_counts_mt <= 30)
        no_germ = a.obs['sample'] == 'NOA3'
        a.obs.loc[no_germ, 'pathology'] = 'author_no_germ_cells'
        a.obs['author_germcell_absence_conflict'] = no_germ & a.obs.is_spg
        a.obs.loc[no_germ, 'is_spg'] = False
        a.obs.loc[no_germ, 'state'] = 'unknown'
        a.obs.loc[no_germ & (a.obs.cell_type == 'SPG'), 'cell_type'] = 'author_germcell_absence_conflict'
    else:
        a.obs['author_germcell_absence_conflict'] = False
        a.obs['clinical_metadata_source'] = 'GEO sample/series only; individual ages unavailable'
        a.obs['AZF_screen'] = 'not_available'
    a.var['gene_ids'] = a.var['gene_ids'].astype(object).fillna('').astype(str)
    return a


def write_gene_presence(allcells, availability):
    """Measured-feature indicators per library and per study (a study requires all of its libraries)."""
    out = availability.astype(int)
    out.index.name = 'gene_symbol'
    out.to_csv(PROCESSED / 'disease_gene_presence_by_library.tsv.gz', sep='\t')
    library_study = allcells.obs[['library', 'study']].drop_duplicates().set_index('library')['study']
    presence = pd.DataFrame(index=allcells.var_names)
    for study in STUDIES:
        presence[study] = availability[library_study.index[library_study == study].tolist()].all(axis=1).astype(int)
    presence.index.name = 'gene_symbol'
    presence.to_csv(PROCESSED / 'disease_gene_presence_by_study.tsv', sep='\t')


def main():
    lineage, state_ref = reference_training_sets()
    studies = [ad.read_h5ad(PROCESSED / (study + '_qc.h5ad')) for study in STUDIES]
    common = set(lineage.var_names)
    for a in studies:
        common &= set(a.var_names)
    lineage_genes = sorted((set(sum(DISEASE_PANELS.values(), [])) | EXTRA_LINEAGE_GENES) & common - AMBIENT_MARKERS)
    state_genes = sorted(set(sum(DISEASE_STATE_MARKERS.values(), [])) & common)
    assert len(lineage_genes) >= 35 and len(state_genes) >= 15
    lineage_model = fit_classifier(lineage, lineage.obs.training_label, lineage_genes)
    state_model = fit_classifier(state_ref, state_ref.obs.state.astype(str), state_genes)
    for a in studies:
        save(annotate(a, lineage_model, lineage_genes, state_model, state_genes), str(a.obs.study.iloc[0]) + '_annotated.h5ad')
    allcells = ad.concat(studies, join='outer', merge='first', fill_value=0)
    allcells.var['gene_ids'] = allcells.var['gene_ids'].astype(object).fillna('').astype(str)
    availability = (pd.concat([a.varm['feature_available_by_library'] for a in studies], axis=1)
                    .reindex(allcells.var_names).fillna(False).astype(bool))
    allcells.varm['feature_available_by_library'] = availability
    save(allcells, 'disease_allcells.h5ad')
    save(allcells[allcells.obs.is_spg].copy(), 'disease_spg.h5ad')
    write_gene_presence(allcells, availability)


if __name__ == '__main__':
    with threadpool_limits(limits=4):
        main()
