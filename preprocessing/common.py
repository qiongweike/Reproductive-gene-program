"""Shared readers, QC helpers and marker panels for the preprocessing scripts."""
import collections
import gzip

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scrublet
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from utils.paths import RAW, METADATA, PROCESSED

SEED = 20260928
MIN_GENES, MAX_MT_PCT = 200, 20                      # basic QC of every reference and disease cell
AGE_QC = dict(min_genes=200, max_genes=5500, min_counts=800, max_counts=17000, max_mt=25)  # GSE254315 paper QC
DOUBLET_RATE, DOUBLET_PCS = 0.05, 30

# Disease compartment panels and U/T/D state markers (GSE202647, GSE326372).
DISEASE_PANELS = {
    'SPG': ['UTF1', 'PIWIL4', 'FGFR3', 'GFRA1', 'ZBTB16', 'UCHL1', 'SALL4', 'LIN28A', 'DMRT1', 'MAGEA4', 'NANOS2', 'ID4', 'KIT', 'SOHLH1', 'SOHLH2'],
    'Spermatocyte': ['SYCP1', 'SYCP2', 'SYCP3', 'DMC1', 'SPO11', 'MEIOB', 'TEX11', 'HORMAD1'],
    'Spermatid': ['ACRV1', 'TNP1', 'TNP2', 'PRM1', 'PRM2', 'SPATA18', 'ODF1', 'AKAP4'],
    'Sertoli': ['SOX9', 'WT1', 'FSHR', 'CLDN11', 'AMH', 'GATA4', 'INHA'],
    'Leydig': ['CYP17A1', 'STAR', 'INSL3', 'HSD3B2', 'CYP11A1', 'LHCGR'],
    'Peritubular': ['ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'DCN', 'COL1A1', 'COL3A1'],
    'Endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN', 'CLDN5', 'RAMP2'],
    'Immune': ['PTPRC', 'CD3D', 'CD3E', 'CD68', 'LYZ', 'C1QA', 'LST1'],
}
DISEASE_STATE_MARKERS = {
    'undifferentiated': ['PIWIL4', 'UTF1', 'GFRA1', 'FGFR3', 'NANOS2', 'ID4', 'TSPAN33', 'PLPPR3', 'ZBTB16', 'ETV5'],
    'transition': ['NANOS3', 'LIN28A', 'SOHLH1', 'SOHLH2', 'SALL4', 'NEUROG3', 'DMRT1', 'UCHL1'],
    'differentiating': ['KIT', 'STRA8', 'SOHLH1', 'SOHLH2', 'DMRT1', 'SALL4', 'DNMT3B', 'MKI67', 'PCNA', 'MCM5'],
}

# Atlas compartment classifier: assignment markers, and withheld markers used only to assess agreement.
ATLAS_ASSIGN = {
    'SPG': ['MAGEA4', 'FGFR3', 'GFRA1', 'UTF1', 'KIT'],
    'Spermatocyte': ['SYCP1', 'SYCP3', 'SPO11', 'MEIOB'],
    'Spermatid': ['PRM1', 'PRM2', 'TNP1', 'ACR'],
    'Sertoli': ['SOX9', 'WT1', 'CLDN11', 'AMH'],
    'Leydig': ['INSL3', 'CYP17A1', 'STAR', 'HSD3B1'],
    'Peritubular_smooth_muscle': ['ACTA2', 'MYH11', 'TAGLN', 'DES'],
    'Endothelial': ['VWF', 'PECAM1', 'CDH5', 'CLDN5'],
    'Immune': ['PTPRC', 'LYZ', 'TYROBP', 'C1QA'],
    'Erythroid': ['HBB', 'HBA1', 'HBA2']}
ATLAS_HELDOUT = {
    'SPG': ['PIWIL4', 'DMRT1', 'UCHL1', 'NANOS3'],
    'Spermatocyte': ['DMC1', 'REC8', 'STAG3', 'PIWIL1', 'HORMAD1', 'MEI1'],
    'Spermatid': ['TNP2', 'IZUMO1', 'SPATA3', 'TSSK6', 'AKAP4', 'ODF1'],
    'Sertoli': ['GATA4', 'INHBB', 'CTSL', 'SHBG', 'FSHR', 'AR'],
    'Leydig': ['CYP11A1', 'DLK1', 'HSD17B3', 'NR5A1', 'LHCGR'],
    'Peritubular_smooth_muscle': ['ACTG2', 'CNN1', 'MYL9', 'WFDC1'],
    'Endothelial': ['EGFL7', 'RAMP2', 'KDR', 'EMCN', 'ESAM'],
    'Immune': ['C1QB', 'AIF1', 'CD68', 'FCER1G', 'LST1'],
    'Erythroid': ['ALAS2', 'AHSP', 'GYPA', 'SLC4A1', 'HBM']}


def soft_samples(study):
    """GSM, title, characteristics and supplementary file names of each sample in a GEO family SOFT file."""
    samples, current = [], None
    with gzip.open(METADATA / f'{study}_family.soft.gz', 'rt') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('^SAMPLE = '):
                current = {'study': study, 'gsm': line.split(' = ', 1)[1], 'files': [], 'characteristics': {}}
                samples.append(current)
            elif current and line.startswith('!Sample_title = '):
                current['title'] = line.split(' = ', 1)[1]
            elif current and line.startswith('!Sample_characteristics_ch1 = '):
                key, value = line.split(' = ', 1)[1].split(': ', 1)
                current['characteristics'][key] = value
            elif current and line.startswith('!Sample_supplementary_file'):
                url = line.split(' = ', 1)[1]
                if url.startswith('ftp://'):
                    current['files'].append(url.rsplit('/', 1)[1])
    return samples


def sample_file(sample, tag):
    """Local path of the sample's supplementary file whose name contains tag."""
    name = next(f for f in sample['files'] if tag in f)
    for path in (RAW / sample['study'] / name, METADATA / name, RAW / name):
        if path.exists():
            return path
    raise FileNotFoundError(name)


def read_mtx(path):
    """10x MatrixMarket file (genes x barcodes) as a barcodes x genes CSR matrix."""
    with gzip.open(path, 'rb') as f, threadpool_limits(1):
        return mmread(f).T.tocsr()


def collapse_symbols(X, symbols, ids=None):
    """Sum features that share a gene symbol; keep the source gene IDs of each symbol."""
    symbols = np.asarray(symbols, dtype=str)
    if ((symbols == '') | (symbols == 'nan')).any():
        raise ValueError('Empty gene symbols')
    unique = pd.Index(pd.unique(symbols))
    if len(unique) != len(symbols):
        mapper = sp.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), unique.get_indexer(symbols))),
                               shape=(len(symbols), len(unique)))
        X = (X @ mapper).tocsr()
    var = pd.DataFrame(index=unique)
    var.index.name = 'gene_symbol'
    if ids is not None:
        groups = collections.defaultdict(list)
        for symbol, gene_id in zip(symbols, ids):
            groups[symbol].append(str(gene_id))
        var['source_gene_ids'] = [';'.join(dict.fromkeys(groups[s])) for s in unique]
    return X, var


def quality_and_log(a, author_qc=False, age_qc=False):
    """Keep raw counts, add full-assay QC metrics and replace X by ln(1 + CP10k)."""
    a.X = a.X.tocsr()
    a.X.eliminate_zeros()
    assert np.isfinite(a.X.data).all() and (a.X.data >= 0).all()
    a.layers['counts'] = a.X.copy()
    totals = np.asarray(a.X.sum(axis=1)).ravel()
    n = np.asarray(a.X.getnnz(axis=1)).ravel()
    mt = np.asarray(a.X[:, a.var_names.str.startswith('MT-')].sum(axis=1)).ravel()
    a.obs['total_counts'] = totals
    a.obs['n_genes'] = n
    a.obs['pct_counts_mt'] = 100 * mt / np.maximum(totals, 1)
    a.obs['qc_basic_pass'] = (n >= MIN_GENES) & (totals > 0) & (a.obs.pct_counts_mt <= MAX_MT_PCT)
    q = AGE_QC
    a.obs['qc_paper_pass'] = ((n >= q['min_genes']) & (n <= q['max_genes']) & (totals >= q['min_counts'])
                              & (totals <= q['max_counts']) & (a.obs.pct_counts_mt <= q['max_mt'])) if age_qc else True
    a.obs['author_qc_retained'] = author_qc
    norm = (sp.diags(10000 / np.maximum(totals, 1)) @ a.X).tocsr()
    norm.data = np.log1p(norm.data)
    a.X = norm
    return a


def scrublet_calls(counts):
    """Scrublet doublet scores and calls for one library (all QC-passing cells)."""
    with threadpool_limits(1):
        screen = scrublet.Scrublet(counts, expected_doublet_rate=DOUBLET_RATE, random_state=SEED)
        scores, calls = screen.scrub_doublets(n_prin_comps=DOUBLET_PCS, verbose=False)
    if calls is None:
        raise RuntimeError('Scrublet found no automatic doublet threshold')
    return scores, calls


def save(a, name):
    """Write an AnnData object to data/processed with gzip compression."""
    a.obs.index.name = 'cell_id'
    assert a.obs_names.is_unique and a.var_names.is_unique
    path = PROCESSED / name
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path, compression='gzip')
    print(f'{name}: {a.n_obs} cells x {a.n_vars} genes', flush=True)
