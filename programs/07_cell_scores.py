"""Cell-level reference scores in the normal spermatogonial reference: cell-cycle score and fixed 143-gene candidate score.

Scores are mean ln(1 + CP10k) (full-assay depth) over a gene set. The cell-cycle score uses the Bush et al. G1/S and G2/M
markers measured in every assay. The candidate score is also given minus a fixed background: for each candidate gene,
30 detected genes drawn from its expression bin (20 bins of mean expression in the background-training donors),
excluding the candidates and published state, cell-cycle and other cell-type markers.
Inputs:  data/processed/normal_reference_geneID_aligned.h5ad; measured genes from data/processed/normal_spg_<study>.h5ad,
         gse235321_spg.h5ad, disease_gene_presence_by_study.tsv and results/preprocessing/gene_symbol_harmonization.json;
         data/external/prior_gene_lists/{guo2017_clusterA.txt, guo2018_state0_enriched.txt, bush2024_state_markers.tsv};
         resources/candidate_143_genes.txt
Output:  results/programs/cell_scores/cell_scores.tsv.gz
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from utils.paths import EXTERNAL, PROCESSED, RESOURCES, RESULTS, out_dir

SEED = 20261001
N_BINS = 20
CONTROLS_PER_GENE = 30
BACKGROUND_STUDY = "GSE254315"
# GSE254315 donors fixed as the background-training set in an earlier round of the analysis.
BACKGROUND_DONORS = ["Donor10", "Donor11", "Donor12", "Donor13", "Donor15", "Donor18", "Donor19", "Donor20", "Donor30",
                     "Donor34", "Donor4", "Donor5", "Donor9"]
UNDIFFERENTIATED_STATES = ["State 0", "State 0A", "State 0B", "SSC", "SSC-1", "SSC-1A", "SSC-1B", "SSC-1C", "SSC-2"]
CYCLE_STATES = ["G1/S", "G2/M"]
OTHER_CELL_TYPES = ["Endothelial cell", "Leydig cell", "Macrophage", "Myoid cell", "Sertoli cell", "Spermatocyte",
                    "Spermatid", "Spermatozoa"]
PRIORS = EXTERNAL / "prior_gene_lists"


def var_names(path):
    a = ad.read_h5ad(path, backed="r")
    genes = set(a.var_names)
    a.file.close()
    return genes


def measured_genes():
    """Features measured in each assay (whether or not observed), after the recorded symbol harmonization."""
    out = {study: var_names(PROCESSED / f"normal_spg_{study}.h5ad") for study in ["GSE124263", "GSE153947", "GSE254315"]}
    out["GSE235321"] = var_names(PROCESSED / "gse235321_spg.h5ad")
    presence = pd.read_csv(PROCESSED / "disease_gene_presence_by_study.tsv", sep="\t", index_col=0)
    for study in presence:
        out[study] = set(presence.index[presence[study].eq(1)])
    harmonization = json.loads((RESULTS / "preprocessing" / "gene_symbol_harmonization.json").read_text())["mappings"]
    for study, aliases in harmonization.items():
        for original, canonical in aliases.items():
            out[study].remove(original)
            out[study].add(canonical)
    out["GSE235321"].remove("H1-10")  # H1-10 is the current symbol of the candidate H1FX (HGNC:4722)
    out["GSE235321"].add("H1FX")
    return out


def marker_sets():
    """Published state, cell-cycle and other cell-type markers excluded from the background pool."""
    guo = set()
    for name in ["guo2017_clusterA", "guo2018_state0_enriched"]:
        guo |= {z.strip() for z in (PRIORS / f"{name}.txt").read_text().splitlines() if z.strip() and not z.startswith("#")}
    bush = pd.read_csv(PRIORS / "bush2024_state_markers.tsv", sep="\t")
    pick = lambda states: set(bush.gene[bush.state.isin(states)].astype(str))
    return dict(undifferentiated=guo | pick(UNDIFFERENTIATED_STATES), cell_cycle=pick(CYCLE_STATES),
                other_cell_types=pick(OTHER_CELL_TYPES))


def expression(a):
    """ln(1 + CP10k) with the stored full-assay library size."""
    x = sparse.csr_matrix(a.layers["counts"], dtype=np.float64)
    lib = pd.to_numeric(a.obs.total_counts).to_numpy(float)
    assert (lib > 0).all()
    x = x.multiply(1e4 / lib[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x


def background_genes(a, candidates, excluded):
    tr = a[(a.obs.study == BACKGROUND_STUDY) & a.obs.donor.isin(BACKGROUND_DONORS)].copy()
    x = expression(tr)
    mean, det = np.asarray(x.mean(0)).ravel(), np.asarray((x > 0).mean(0)).ravel()
    bins = pd.qcut(pd.Series(mean).rank(method="first"), N_BINS, labels=False).to_numpy()
    rng = np.random.default_rng(SEED)
    background = set()
    for g in candidates:
        j = tr.var_names.get_loc(g)
        pool = np.flatnonzero((bins == bins[j]) & (det > 0) & ~tr.var_names.isin(excluded))
        background.update(tr.var_names[rng.choice(pool, min(CONTROLS_PER_GENE, len(pool)), replace=False)])
    return sorted(background)


def main():
    candidates = (RESOURCES / "candidate_143_genes.txt").read_text().split()
    markers = marker_sets()
    a = ad.read_h5ad(PROCESSED / "normal_reference_geneID_aligned.h5ad")
    common = set.intersection(*measured_genes().values()) & set(a.var_names)
    assert set(candidates) <= common
    background = background_genes(a[:, a.var_names.isin(common)], candidates, set(candidates).union(*markers.values()))

    x = expression(a)
    mean_of = lambda genes: np.asarray(x[:, a.var_names.get_indexer(genes)].mean(1)).ravel()
    d = a.obs[["study", "donor", "diagnosis", "state"]].astype(str)
    d["cell_id"] = a.obs_names.astype(str)
    d["donor_key"] = d.study + ":" + d.donor
    d["score_units"] = "log1p_CP10k_full_assay_denominator"
    d["fixed_143_raw"] = mean_of(candidates)
    d["fixed_background"] = mean_of(background)
    d["fixed_143_adjusted"] = d.fixed_143_raw - d.fixed_background
    d["reference_cell_cycle"] = mean_of(sorted(markers["cell_cycle"] & common))
    cols = ["cell_id", "study", "donor", "donor_key", "diagnosis", "state", "score_units", "fixed_143_raw",
            "fixed_143_adjusted", "reference_cell_cycle", "fixed_background"]
    d[cols].to_csv(out_dir("programs", "cell_scores") / "cell_scores.tsv.gz", sep="\t", index=False)
    print(f"{len(d)} cells; {len(background)} background genes; {len(markers['cell_cycle'] & common)} cell-cycle genes")


if __name__ == "__main__":
    main()
