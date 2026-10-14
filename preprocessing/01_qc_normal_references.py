"""QC of the GSE124263 (adult donors A1, A2) and GSE153947 (OA, cryptozoospermia) references.

Cell identities and U/T/D states come from the author annotations (GSE124263 workbook sheets 4 and 5,
GSE153947 cell metadata); GSE153947 Alevin estimated counts stay fractional. Every library is screened
with Scrublet on all of its QC-passing cells, not only spermatogonia. Outputs (data/processed):
normal_allcells_GSE124263.h5ad, gse153947_spg_all.h5ad, GSE153947_<donor>_whole_library_doublet_QC.tsv.gz.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gzip
import xml.etree.ElementTree as ET
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from utils.paths import METADATA, PROCESSED
from common import soft_samples, sample_file, read_mtx, collapse_symbols, quality_and_log, scrublet_calls, save

STATES_GSE124263 = {'SSC1': 'undifferentiated', 'SSC2': 'undifferentiated', 'Transition Cell': 'transition',
                    'Early diff SPG': 'differentiating', 'Diff SPG': 'differentiating'}
STATES_GSE153947 = {'State 0': 'undifferentiated', 'State 0A': 'undifferentiated', 'State 0B': 'undifferentiated',
                    'State 1': 'undifferentiated', 'State 2': 'differentiating', 'State 3': 'differentiating'}


def read_xlsx_sheet(path, number):
    """One worksheet of an .xlsx file as strings, indexed by its first column."""
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        shared = [''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns)]
        rows = []
        for r in ET.fromstring(z.read(f'xl/worksheets/sheet{number}.xml')).findall('m:sheetData/m:row', ns):
            values = {}
            for c in r.findall('m:c', ns):
                v = c.find('m:v', ns)
                text = v.text if v is not None else ''
                column = 0
                for ch in ''.join(x for x in c.get('r', '') if x.isalpha()):
                    column = column * 26 + ord(ch) - 64
                values[column - 1] = shared[int(text)] if c.get('t') == 's' and text else text
            rows.append(values)
    width = max(rows[0]) + 1
    table = [[r.get(i, '') for i in range(width)] for r in rows]
    return pd.DataFrame(table[1:], columns=table[0]).set_index(table[0][0])


def read_alevin(path):
    """Alevin genes x barcodes text table as a barcodes x genes CSR matrix (float64, fractions kept)."""
    columns, values, indptr, genes = [], [], [0], []
    with gzip.open(path, 'rt') as f:
        barcodes = f.readline().rstrip('\r\n').split('\t')
        for line in f:
            if not line.strip():
                continue
            gene, text = line.rstrip('\r\n').split('\t', 1)
            row = np.fromstring(text, sep='\t', dtype=np.float64)
            assert len(row) == len(barcodes) and np.isfinite(row).all() and (row >= 0).all()
            nz = np.flatnonzero(row)
            columns.append(nz)
            values.append(row[nz])
            indptr.append(indptr[-1] + len(nz))
            genes.append(gene)
    X = sp.csr_matrix((np.concatenate(values), np.concatenate(columns), np.array(indptr)),
                      shape=(len(genes), len(barcodes))).T.tocsr()
    assert len(set(barcodes)) == len(barcodes)
    return X, genes, pd.Index(barcodes)


def add_doublets(a, groups):
    """Scrublet on the QC-passing cells of each group (library); unscreened cells keep NaN/False."""
    a.obs['doublet_score'] = np.nan
    a.obs['predicted_doublet'] = False
    a.obs['doublet_screen_succeeded'] = False
    for group in groups:
        sel = np.asarray(group) & a.obs.qc_basic_pass.to_numpy()
        scores, calls = scrublet_calls(a.layers['counts'][sel].copy())
        a.obs.loc[sel, 'doublet_score'] = scores
        a.obs.loc[sel, 'predicted_doublet'] = calls
        a.obs.loc[sel, 'doublet_screen_succeeded'] = True
    return a


def gse124263():
    workbook = METADATA / 'GSE124263_Barcode_Cell_Identity.xlsx'
    cells, spg = read_xlsx_sheet(workbook, 4), read_xlsx_sheet(workbook, 5)
    libraries = []
    for s in soft_samples('GSE124263'):
        if not s['title'].startswith('Adult'):
            continue
        lib = s['characteristics']['sample']
        donor, prefix = lib.split('_')[0], lib.replace('_', '')
        genes = pd.read_csv(sample_file(s, 'genes.tsv'), sep='\t', header=None)
        barcodes = pd.read_csv(sample_file(s, 'barcodes.tsv'), sep='\t', header=None)[0].astype(str)
        X = read_mtx(sample_file(s, 'matrix.mtx'))
        assert X.shape == (len(barcodes), len(genes)) and np.allclose(X.data, np.round(X.data))
        author_ids = pd.Index([prefix + '_' + b.split('-')[0] for b in barcodes])
        keep = author_ids.isin(cells.index)
        assert len(set(author_ids[keep])) == int(keep.sum())
        obs = cells.loc[author_ids[keep]].copy()
        obs['author_cell_type'] = obs['AdultCellAnnotation']
        obs['library'] = lib
        obs['study'] = 'GSE124263'
        obs['donor'] = donor
        obs['diagnosis'] = 'Normal'
        obs['clinical_description'] = 'fertile adult organ donor'
        obs['age'] = float(s['characteristics']['age'].split()[0])
        obs['enrichment'] = 'ITGA6_enriched' if lib.endswith('_I') else 'unfractionated'
        obs['barcode'] = author_ids[keep]
        obs['gsm'] = s['gsm']
        obs['count_type'] = 'integer_UMI'
        obs['author_state'] = spg['TransitionCellIdentity'].reindex(obs.index).fillna('unknown')
        obs['author_SPG_subset'] = spg['SPG_Subsets'].reindex(obs.index).fillna('unknown')
        obs['author_SSC1_subset'] = spg['SSC1_Subsets'].reindex(obs.index).fillna('unknown')
        obs['state'] = obs['author_state'].map(STATES_GSE124263).fillna('unknown')
        obs['is_author_spg'] = obs['author_state'] != 'unknown'
        obs['normal_training_eligible'] = obs['is_author_spg']
        obs.index = 'GSE124263:' + donor + ':' + obs['barcode']
        X, var = collapse_symbols(X[keep].astype(np.float64), genes[1].astype(str), genes[0].astype(str))
        a = quality_and_log(ad.AnnData(X, obs=obs, var=var), author_qc=True)
        a.obs['normal_training_eligible'] &= a.obs.qc_basic_pass
        libraries.append(a)
    a = ad.concat(libraries, join='inner', merge='same', uns_merge='same')
    a = add_doublets(a, [(a.obs.library == lib).to_numpy() for lib in sorted(a.obs.library.unique())])
    a.obs['normal_training_eligible'] = (a.obs.is_author_spg & a.obs.qc_basic_pass
                                         & a.obs.doublet_screen_succeeded & ~a.obs.predicted_doublet)
    save(a, 'normal_allcells_GSE124263.h5ad')


def gse153947():
    meta = pd.read_csv(METADATA / 'GSE153947_Cell_metadata.tsv.gz', sep='\t', index_col=0)
    meta = meta[meta.Spermatogonia_Cluster_identity != 'other'].copy()
    donors = []
    for s in soft_samples('GSE153947'):
        donor = s['title'].replace('Normal ', 'N').replace('Crypto ', 'Cr')
        counts, genes, barcodes = read_alevin(sample_file(s, '_counts.tsv.gz'))
        # Whole-library QC and doublet screen.
        X, var = collapse_symbols(counts, genes)
        library = quality_and_log(ad.AnnData(X, obs=pd.DataFrame(index=barcodes.rename('raw_barcode')), var=var))
        library = add_doublets(library, [np.ones(len(barcodes), dtype=bool)])
        library.obs.to_csv(PROCESSED / f'GSE153947_{donor}_whole_library_doublet_QC.tsv.gz', sep='\t')
        # Author-annotated spermatogonia of this donor.
        obs = meta[meta.Sample_ID == donor].copy()
        rows = barcodes.get_indexer(obs.index.str.rsplit('_', n=1).str[0])
        assert (rows >= 0).all()
        X, var = collapse_symbols(counts[rows], genes)
        obs['study'] = 'GSE153947'
        obs['donor'] = donor
        obs['diagnosis'] = 'OA' if donor.startswith('N') else 'Crypto'
        obs['clinical_description'] = s['characteristics']['diagnosis']
        obs['age'] = float(s['characteristics']['age (years)'])
        obs['enrichment'] = 'unfractionated'
        obs['library'] = s['gsm']
        obs['gsm'] = s['gsm']
        obs['barcode'] = obs.index
        obs['author_cell_type'] = obs.Cluster_identity.astype(str)
        obs['author_state'] = obs.Spermatogonia_Cluster_identity
        obs['state'] = obs.author_state.map(STATES_GSE153947).fillna('unknown')
        obs['is_author_spg'] = True
        obs['count_type'] = 'fractional_Alevin_estimated_UMI'
        obs.index = 'GSE153947:' + donor + ':' + obs['barcode']
        a = quality_and_log(ad.AnnData(X, obs=obs, var=var), author_qc=True)
        a.obs['normal_training_eligible'] = (a.obs.diagnosis == 'OA') & a.obs.qc_basic_pass
        for column in ['doublet_score', 'predicted_doublet', 'doublet_screen_succeeded']:
            a.obs[column] = library.obs[column].to_numpy()[rows]
        a.obs['normal_training_eligible'] = ((a.obs.diagnosis == 'OA') & a.obs.qc_basic_pass
                                             & a.obs.doublet_screen_succeeded & ~a.obs.predicted_doublet)
        donors.append(a)
    save(ad.concat(donors, join='inner', merge='same', uns_merge='same'), 'gse153947_spg_all.h5ad')


if __name__ == '__main__':
    gse124263()
    gse153947()
