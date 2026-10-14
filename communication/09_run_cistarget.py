"""pySCENIC modules (0.75/0.90, top 50 targets, top 5/10/50 regulators, >=20 genes, rho > 0.03), cisTarget pruning
(two hg38 v10_clust rankings; rank 1,500, AUC 0.05, NES 3, motif-similarity FDR 0.001) -> 294 regulons, and the
motif support of every regulon edge.
Inputs: scenic/{expression.npy, genes.tsv, cells.tsv.gz, full_GRNBoost_adjacencies.tsv.gz}, data/external/cistarget/
Outputs: results/communication/scenic/ (full_SCENIC_regulons.pkl, motif enrichment, per-motif and per-edge motif support)
"""
import os
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import EXTERNAL, RESULTS
import pickle
import numpy as np
import pandas as pd
np.__dict__.setdefault("object", object)  # pySCENIC 0.12.1 still refers to np.object (removed in NumPy 1.24)
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.prune import df2regulons, prune2df
from pyscenic.utils import modules_from_adjacencies

MODULE_PARAMS = dict(thresholds=(0.75, 0.90), top_n_targets=(50,), top_n_regulators=(5, 10, 50), min_genes=20,
                     rho_threshold=0.03, rho_mask_dropouts=False, keep_only_activating=True)
PRUNE_PARAMS = dict(rank_threshold=1500, auc_threshold=0.05, nes_threshold=3, motif_similarity_fdr=0.001)
MOTIF_TABLE = "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
WORKERS = int(os.environ.get("REGULON_WORKERS", "8"))  # >= 2 (one per ranking database); does not change results
REG = RESULTS / "communication" / "scenic"
CIS = EXTERNAL / "cistarget"


def motif_edge_rows(enrichment):
    """One row per (enriched motif, target gene) behind each regulon."""
    rows = []
    for (tf, motif), row in enrichment.iterrows():
        v = row["Enrichment"]
        context = set(v["Context"])
        regulon = f"{tf}({'-' if 'repressing' in context else '+'})"
        database = ";".join(sorted(c for c in context if c.startswith("hg38_")))
        for target, weight in v["TargetGenes"]:
            rows.append(dict(regulon=regulon, TF=tf, target=target, motif=motif, NES=v["NES"], AUC=v["AUC"],
                             annotation=v["Annotation"], motif_similarity_q=v["MotifSimilarityQvalue"],
                             orthologous_identity=v["OrthologousIdentity"], database=database, importance=weight))
    return pd.DataFrame(rows)


def main():
    genes = pd.read_csv(REG / "genes.tsv", sep="\t").gene.to_numpy()
    cells = pd.read_csv(REG / "cells.tsv.gz", sep="\t")
    expr = pd.DataFrame(np.load(REG / "expression.npy"), index=cells.cell_id, columns=genes)
    adjacency = pd.read_csv(REG / "full_GRNBoost_adjacencies.tsv.gz", sep="\t")
    modules = list(modules_from_adjacencies(adjacency, expr, **MODULE_PARAMS))
    databases = [FeatherRankingDatabase(str(p), p.stem) for p in sorted(CIS.glob("*.rankings.feather"))]
    enrichment = prune2df(databases, modules, str(CIS / MOTIF_TABLE), **PRUNE_PARAMS,
                          client_or_address="custom_multiprocessing", num_workers=WORKERS)
    enrichment.to_pickle(REG / "full_SCENIC_motif_enrichment.pkl")
    regulons = df2regulons(enrichment)
    with open(REG / "full_SCENIC_regulons.pkl", "wb") as f:
        pickle.dump(regulons, f)
    pd.DataFrame([dict(regulon=r.name, TF=r.transcription_factor, target=g, importance=w,
                       context=";".join(sorted(r.context))) for r in regulons for g, w in r.gene2weight.items()]).to_csv(
        REG / "full_SCENIC_motif_supported_edges.tsv", sep="\t", index=False)
    motif_rows = motif_edge_rows(enrichment)
    motif_rows.to_csv(REG / "full_SCENIC_motif_edge_evidence.tsv.gz", sep="\t", index=False)
    summary = motif_rows.groupby(["regulon", "TF", "target"]).agg(
        n_motifs=("motif", "nunique"), max_motif_NES=("NES", "max"),
        databases=("database", lambda a: ";".join(sorted(set(a)))),
        annotation_types=("annotation", lambda a: ";".join(sorted(set(a))))).reset_index()
    summary.to_csv(REG / "full_SCENIC_motif_edge_summary.tsv", sep="\t", index=False)
    print(f"cisTarget: {len(modules)} modules, {len(regulons)} regulons, {len(summary)} motif-supported edges")


if __name__ == "__main__":
    main()
