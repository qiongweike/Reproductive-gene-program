"""Final normal spermatogonial reference (9,080 cells from 30 donors) and per-study spermatogonia files.

Eligible: adult Normal or OA donor, author (GSE124263, GSE153947) or transferred (GSE254315) SPG identity,
basic QC, successful whole-library Scrublet and not a predicted doublet. GSE153947 cryptozoospermia
spermatogonia pass the same QC but are kept separately and never enter the reference.
Outputs: data/processed/normal_reference.h5ad (common genes), normal_spg_<study>.h5ad (full assays), crypto_spg.h5ad.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anndata as ad

from utils.paths import PROCESSED
from common import save


def passes_qc(obs):
    return obs.qc_basic_pass & obs.doublet_screen_succeeded & ~obs.predicted_doublet


if __name__ == '__main__':
    a = ad.read_h5ad(PROCESSED / 'normal_allcells_GSE124263.h5ad')
    a = a[a.obs.is_author_spg].copy()
    b = ad.read_h5ad(PROCESSED / 'gse153947_spg_all.h5ad')
    c = ad.read_h5ad(PROCESSED / 'normal_spg_GSE254315.h5ad')
    assert c.obs.normal_training_eligible.all()
    eligible = [x[x.obs.normal_training_eligible].copy() for x in [a, b[b.obs.diagnosis == 'OA'].copy(), c]]
    assert all(passes_qc(x.obs).all() for x in eligible)
    reference = ad.concat(eligible, join='inner', merge='same', uns_merge='same')
    assert not reference.obs.diagnosis.astype(str).str.contains('Crypto|NOA', case=False).any()
    save(reference, 'normal_reference.h5ad')
    crypto = b[(b.obs.diagnosis == 'Crypto') & passes_qc(b.obs)].copy()
    assert not crypto.obs.normal_training_eligible.any()
    save(crypto, 'crypto_spg.h5ad')
    for x in eligible[:2]:
        save(x, f'normal_spg_{x.obs.study.iloc[0]}.h5ad')
