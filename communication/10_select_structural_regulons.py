"""Motif-supported edges to the fixed program members (canonical names; MP3 = 49 genes) and the structural display set
(Fig 5e): TFs with >=2 non-self MP3 targets, ranked by target count, maximum motif NES and name; the first 12.
Inputs: scenic/full_SCENIC_motif_edge_summary.tsv (step 09),
  results/programs/{MP_genes.tsv, identity/verified_aliases.tsv}
Outputs: results/communication/regulon/{program_motif_edges.tsv, MP3_motif_edges.tsv,
  structural_TF_selection.tsv}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir
import numpy as np
import pandas as pd

MIN_MP3_TARGETS = 2
N_DISPLAY_TFS = 12
OUT = out_dir("communication", "regulon")


def canonical_symbol():
    aliases = pd.read_csv(RESULTS / "programs/identity/verified_aliases.tsv", sep="\t")
    amap = dict(zip(aliases.old_symbol.str.upper(), aliases.approved_symbol))
    return lambda g: amap.get(str(g).upper(), str(g).upper())


def main():
    canon = canonical_symbol()
    mp = pd.read_csv(RESULTS / "programs/MP_genes.tsv", sep="\t")
    motif = pd.read_csv(RESULTS / "communication/scenic/full_SCENIC_motif_edge_summary.tsv", sep="\t")
    motif["assay_target"] = motif.target
    motif["target"] = motif.target.map(canon)
    program_edges = motif.merge(mp[["MP", "gene"]], left_on="target", right_on="gene", how="inner")
    program_edges.to_csv(OUT / "program_motif_edges.tsv", sep="\t", index=False)
    mp3 = program_edges.query('MP == "MP_3" and TF != target').drop_duplicates(["TF", "target"])
    mp3.to_csv(OUT / "MP3_motif_edges.tsv", sep="\t", index=False)
    ranking = (mp3.groupby("TF").agg(n_MP3_targets=("target", "size"), max_NES=("max_motif_NES", "max")).reset_index()
               .query("n_MP3_targets >= @MIN_MP3_TARGETS")
               .sort_values(["n_MP3_targets", "max_NES", "TF"], ascending=[False, False, True]))
    ranking["selected"] = np.arange(len(ranking)) < N_DISPLAY_TFS
    ranking.to_csv(OUT / "structural_TF_selection.tsv", sep="\t", index=False)
    shown = mp3[mp3.TF.isin(ranking.TF[ranking.selected])]
    print(f"{len(mp3)} MP3 motif edges ({mp3.target.nunique()} targets); display: {N_DISPLAY_TFS} TFs, "
          f"{len(shown)} edges, {shown.target.nunique()} targets")


if __name__ == "__main__":
    main()
