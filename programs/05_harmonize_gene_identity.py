"""Harmonize gene identities with five HGNC-verified aliases (symbols are upper-cased) and fix the paper's programs.

MP3's 50 listed symbols contain LINC01511 and RP11-325I22.2, one entity (ENSG00000250584), so MP3 has 49 members.
Count matrices are not modified; the manifests used for scoring point to gene files with the harmonized symbols.
Inputs:  data/external/hgnc_complete_set.txt, data/external/prior_gene_lists/bush2024_state_markers.tsv,
         results/programs/nmf/{MP_genes.tsv, donor_manifest.tsv, NOA_assignment_manifest.tsv, input/<sample>.genes.tsv}
Outputs: results/programs/{MP_genes.tsv, donor_manifest.tsv, NOA_assignment_manifest.tsv, input/<sample>.genes.tsv},
         results/programs/identity/{verified_aliases.tsv, MP_member_symbols.tsv, state_markers.tsv}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd

from utils.paths import EXTERNAL, RESULTS, out_dir

ALIASES = {"FAM26F": "CALHM6", "HIST1H1D": "H1-3", "FAM207A": "SLX9", "C4ORF3": "ARLN", "RP11-325I22.2": "LINC01511"}
MP_SIZES = {"MP_1": 50, "MP_2": 50, "MP_3": 49}

NMF = RESULTS / "programs" / "nmf"
PROGRAMS = out_dir("programs")
IDENTITY = out_dir("programs", "identity")
INPUT = out_dir("programs", "input")


def canon(gene):
    gene = str(gene).upper()
    return ALIASES.get(gene, gene)


def verified_aliases():
    """Each replaced symbol must be a previous or alias symbol of exactly one approved HGNC entry."""
    hgnc = pd.read_csv(EXTERNAL / "hgnc_complete_set.txt", sep="\t", dtype=str).fillna("")
    rows = []
    for old, new in ALIASES.items():
        hit = hgnc[hgnc.symbol.str.upper().eq(new)]
        assert len(hit) == 1, (old, new)
        hit = hit.iloc[0]
        assert old in (hit.prev_symbol + "|" + hit.alias_symbol).upper().split("|"), (old, new)
        rows.append(dict(old_symbol=old, approved_symbol=hit.symbol, canonical_key=new, HGNC=hit.hgnc_id,
                         Ensembl=hit.ensembl_gene_id))
    return pd.DataFrame(rows)


def canonical_manifest(name):
    """Point each exported donor at a canonical gene file (columns: gene, raw_symbol)."""
    m = pd.read_csv(NMF / name, sep="\t")
    for i, r in m.iterrows():
        if pd.isna(r.gene_file):
            continue
        g = pd.read_csv(r.gene_file, sep="\t")
        g["raw_symbol"] = g.gene
        g["gene"] = g.gene.map(canon)
        assert not g.gene.duplicated().any(), r.sample_key
        dest = INPUT / Path(r.gene_file).name
        g.to_csv(dest, sep="\t", index=False)
        m.loc[i, "gene_file"] = str(dest)
    m.to_csv(PROGRAMS / name, sep="\t", index=False)


def main():
    verified_aliases().to_csv(IDENTITY / "verified_aliases.tsv", sep="\t", index=False)
    for name in ["donor_manifest.tsv", "NOA_assignment_manifest.tsv"]:
        canonical_manifest(name)

    mp = pd.read_csv(NMF / "MP_genes.tsv", sep="\t")
    mp["listed_symbol"] = mp.gene
    mp["gene"] = mp.gene.map(canon)
    mp.to_csv(IDENTITY / "MP_member_symbols.tsv", sep="\t", index=False)
    mp = mp.sort_values(["MP", "rank"]).drop_duplicates(["MP", "gene"])
    assert mp.groupby("MP").size().to_dict() == MP_SIZES
    mp.to_csv(PROGRAMS / "MP_genes.tsv", sep="\t", index=False)

    markers = pd.read_csv(EXTERNAL / "prior_gene_lists" / "bush2024_state_markers.tsv", sep="\t")
    markers["gene"] = markers.gene.map(canon)
    markers.drop_duplicates(["state", "gene"]).to_csv(IDENTITY / "state_markers.tsv", sep="\t", index=False)
    print("MP sizes", MP_SIZES)


if __name__ == "__main__":
    main()
