"""Measured expression of the 41 signalling components of Fig 5a: the top ten U-context NicheNet ligands, their
receptors on U spermatogonia, FGFR3, COMMD3, ID1, ID2 and ID4. Non-SPG groups (>=30 cells) are averaged with equal
donor weight; U spermatogonia are shown per donor (including the 28 NOA2 cells).
Inputs: GSE202647_annotated.h5ad, all_cell_eligibility.tsv.gz, nichenet/ (step 05), programs/identity/verified_aliases.tsv
Outputs: results/communication/{signalling_expression.tsv, signalling_gene_identity.tsv}
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import PROCESSED, RESULTS, out_dir
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

N_LIGANDS = 10
FIXED_COMPONENTS = ["FGFR3", "COMMD3", "ID1", "ID2", "ID4"]
N_COMPONENTS = 41
MIN_CELLS = 30
U_CELLS = {"GSE202647_NOA1": 160, "GSE202647_NOA2": 28, "GSE202647_OA1": 137, "GSE202647_OA2": 107}
OUT = out_dir("communication")
NICHENET = RESULTS / "communication/nichenet"


def canonical_symbol():
    aliases = pd.read_csv(RESULTS / "programs/identity/verified_aliases.tsv", sep="\t")
    amap = dict(zip(aliases.old_symbol.str.upper(), aliases.approved_symbol))
    return lambda g: amap.get(str(g).upper(), str(g).upper())


def display_genes(canon):
    best = pd.read_csv(NICHENET / "U_ligand_selection.tsv", sep="\t").head(N_LIGANDS).test_ligand.tolist()
    lr = pd.read_csv(NICHENET / "sender_receptor_evidence.tsv.gz", sep="\t")
    receptors = lr.loc[lr.state.eq("SPG_U") & lr.MP.eq("MP_3") & lr.prior_ligand.isin(best), "to"].drop_duplicates()
    genes = list(dict.fromkeys(map(canon, best + receptors.tolist() + FIXED_COMPONENTS)))
    assert len(genes) == N_COMPONENTS
    return genes


def main():
    canon = canonical_symbol()
    genes = display_genes(canon)
    meta = pd.read_csv(RESULTS / "communication/cellchat_inputs/all_cell_eligibility.tsv.gz",
                       sep="\t", low_memory=False)
    meta = meta[meta.pre_group_eligible]
    rows, identity = [], []
    with h5py.File(PROCESSED / "GSE202647_annotated.h5ad", "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        names = var.index.map(canon)
        assert names.is_unique
        col = {g: i for i, g in enumerate(names)}
        layer = sparse_dataset(f["layers/counts"])
        for g in genes:
            j = col.get(g)
            identity.append(dict(display_gene=g, assay_gene=str(var.index[j]) if j is not None else "",
                                 Ensembl=str(var.gene_ids.iloc[j]) if j is not None else "", measured=j is not None))
        for (donor, group), m in meta.groupby(["donor", "communication_group"]):
            ii = obs.index.get_indexer(m.source_cell_id)
            assert (ii >= 0).all()
            a = layer[ii].tocsr()
            depth = np.asarray(a.sum(1)).ravel()
            assert np.allclose(depth, m.full_assay_total_counts)
            n = len(ii)
            for order, g in enumerate(genes):
                j = col.get(g)
                detection = mean = np.nan
                if j is not None:
                    raw = a[:, j].tocsr()
                    norm = raw.multiply(1e4 / depth[:, None]).tocsr()
                    norm.data = np.log1p(norm.data)
                    detection = float(np.asarray((raw > 0).mean(0)).ravel()[0])
                    mean = float(np.asarray(norm.mean(0)).ravel()[0])
                rows.append(dict(donor=donor, diagnosis=m.diagnosis.iloc[0], group=group, n_cells=n,
                                 eligible_primary=n >= MIN_CELLS, gene=g, gene_order=order, measured=j is not None,
                                 detection=detection, mean_logCP10k=mean,
                                 interpretation="descriptive_measured_expression"))
    z = pd.DataFrame(rows)
    z["mean_logCP10k"] = z.mean_logCP10k.astype("float32")
    send = (z[~z.group.str.startswith("SPG") & z.eligible_primary].groupby(["group", "gene", "gene_order"])
            .agg(detection=("detection", "mean"), mean_logCP10k=("mean_logCP10k", "mean"), n_donors=("donor", "nunique"),
                 n_cells=("n_cells", "sum")).reset_index())
    send["row"] = send.group + " (" + send.n_donors.astype(str) + ")"
    send["unit"] = "equal_weight_qualified_donors"
    send["eligible_primary"] = True
    recv = z[z.group.eq("SPG_U")].copy()
    recv["row"] = recv.donor.str.replace("GSE202647_", "", regex=False) + " / U (" + recv.n_cells.astype(str) + ")"
    recv["unit"] = "one_donor_descriptive_expression"
    recv["n_donors"] = 1
    assert recv.groupby("donor").n_cells.first().to_dict() == U_CELLS
    pd.concat([send, recv], ignore_index=True).to_csv(OUT / "signalling_expression.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(identity).to_csv(OUT / "signalling_gene_identity.tsv", sep="\t", index=False)
    print(f"{len(genes)} components, {z.groupby(['donor', 'group']).ngroups} donor groups")


if __name__ == "__main__":
    main()
