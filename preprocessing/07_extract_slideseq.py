"""Extract the two normal human Slide-seq arrays (Puck5, Puck6) from the Human.7z release.

Integer bead counts, bead coordinates and the authors' numeric bead cell-type codes are stored unfiltered;
names of the numeric codes are marker-inferred domains. The archive is read with libarchive.
Input: data/raw/Human.7z (github.com/HaiqiChenLab/Slideseq_SSC_Niche).
Outputs: data/processed/slideseq_Puck5.h5ad, data/processed/slideseq_Puck6.h5ad.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ctypes as C
import ctypes.util
import io
from array import array

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from utils.paths import RAW, PROCESSED

ARRAYS = ['Puck5', 'Puck6']
DOMAIN_LABELS = {'1': 'Elongating_elongated_spermatid', '2': 'Round_spermatid', '3': 'SPG_early_like',
                 '4': 'Spermatocyte', '5': 'SPG_late_like', '6': 'Sertoli', '7': 'Leydig',
                 '8': 'Endothelial', '9': 'Macrophage', '10': 'Myoid'}


class Archive:
    """Sequential reader of archive entries through the system libarchive."""

    def __init__(self, path):
        self.lib = C.CDLL(ctypes.util.find_library('archive'))
        self.lib.archive_read_new.restype = C.c_void_p
        for name in ['archive_read_support_filter_all', 'archive_read_support_format_all', 'archive_read_free']:
            getattr(self.lib, name).argtypes = [C.c_void_p]
        self.lib.archive_read_open_filename.argtypes = [C.c_void_p, C.c_char_p, C.c_size_t]
        self.lib.archive_read_next_header.argtypes = [C.c_void_p, C.POINTER(C.c_void_p)]
        self.lib.archive_entry_pathname.argtypes = [C.c_void_p]
        self.lib.archive_entry_pathname.restype = C.c_char_p
        self.lib.archive_read_data.argtypes = [C.c_void_p, C.c_void_p, C.c_size_t]
        self.lib.archive_read_data.restype = C.c_ssize_t
        self.h = self.lib.archive_read_new()
        self.lib.archive_read_support_filter_all(self.h)
        self.lib.archive_read_support_format_all(self.h)
        if self.lib.archive_read_open_filename(self.h, str(path).encode(), 10240) != 0:
            raise RuntimeError('Unable to open archive ' + str(path))

    def entries(self):
        entry = C.c_void_p()
        while True:
            result = self.lib.archive_read_next_header(self.h, C.byref(entry))
            if result == 1:
                break
            if result != 0:
                raise RuntimeError(f'Archive header error {result}')
            yield self.lib.archive_entry_pathname(entry).decode()

    def lines(self):
        buffer, pending = C.create_string_buffer(8 * 1024 * 1024), b''
        while True:
            n = self.lib.archive_read_data(self.h, buffer, len(buffer))
            if n == 0:
                break
            if n < 0:
                raise RuntimeError('Archive decompression error')
            chunks = (pending + buffer.raw[:n]).split(b'\n')
            pending = chunks.pop()
            yield from chunks
        if pending:
            yield pending

    def close(self):
        self.lib.archive_read_free(self.h)


def read_counts(lines):
    """Bead x gene CSV (',barcode,gene1,...') as an int32 CSR matrix."""
    it = iter(lines)
    header = next(it).decode().rstrip('\r').split(',')
    assert header[:2] == ['', 'barcode']
    genes = header[2:]
    assert len(set(genes)) == len(genes)
    values, indices, indptr, barcodes = array('i'), array('i'), array('q', [0]), []
    for line in it:
        if not line:
            continue
        fields = line.rstrip(b'\r').split(b',', 2)
        assert len(fields) == 3
        raw = np.fromstring(fields[2].decode('ascii'), sep=',', dtype=np.float64)
        assert len(raw) == len(genes) and np.isfinite(raw).all() and (raw >= 0).all() and (raw == np.floor(raw)).all()
        assert raw.max(initial=0) <= np.iinfo(np.int32).max
        nz = np.flatnonzero(raw).astype(np.int32)
        values.frombytes(raw[nz].astype(np.int32).tobytes())
        indices.frombytes(nz.tobytes())
        indptr.append(len(values))
        barcodes.append(fields[1].decode())
    assert len(set(barcodes)) == len(barcodes)
    matrix = sp.csr_matrix((np.frombuffer(values, dtype=np.int32).copy(), np.frombuffer(indices, dtype=np.int32).copy(),
                            np.frombuffer(indptr, dtype=np.int64).copy()), shape=(len(barcodes), len(genes)))
    return matrix, barcodes, genes


def build_array(key, matrix, barcodes, genes, coords, labels):
    coords, labels = coords.set_index('barcode'), labels.set_index('barcode')
    assert not coords.index.has_duplicates and not labels.index.has_duplicates
    assert set(barcodes) == set(coords.index) == set(labels.index)
    coords, labels = coords.loc[barcodes], labels.loc[barcodes]
    numeric = labels['max_cell_type'].astype(int).astype(str)
    assert set(numeric) <= set(DOMAIN_LABELS)
    obs = pd.DataFrame(index=pd.Index([key + '__' + b for b in barcodes], name='spatial_id'))
    obs['barcode'] = barcodes
    obs['array_id'] = key
    obs['donor_id'] = 'unknown'
    obs['condition'] = 'normal'
    obs['author_cell_type_numeric'] = numeric.to_numpy()
    obs['domain_label'] = numeric.map(DOMAIN_LABELS).to_numpy()
    obs['author_umi'] = pd.to_numeric(coords['UMI']).to_numpy()
    obs['author_pseudotime'] = pd.to_numeric(coords['Pseudotime'], errors='coerce').to_numpy()
    obs['author_normalized_pseudotime'] = pd.to_numeric(coords['Normalized Pseudotime'], errors='coerce').to_numpy()
    obs['total_counts'] = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
    obs['n_genes_by_counts'] = np.diff(matrix.indptr)
    spatial = coords[['x', 'y']].to_numpy(dtype=np.float64)
    assert np.isfinite(spatial).all()
    a = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame({'symbol': genes}, index=pd.Index(genes, name='gene_symbol')))
    a.layers['counts'] = matrix.copy()
    a.obsm['spatial'] = spatial
    a.write_h5ad(PROCESSED / f'slideseq_{key}.h5ad', compression='gzip')
    print(f'slideseq_{key}.h5ad: {a.n_obs} beads x {a.n_vars} genes', flush=True)


if __name__ == '__main__':
    parts = {key: {} for key in ARRAYS}
    archive = Archive(RAW / 'Human.7z')
    try:
        for name in archive.entries():
            if name.endswith('/'):
                continue
            key = next(k for k in ARRAYS if k in name)
            if 'MappedDGE' in name:
                parts[key]['matrix'] = read_counts(archive.lines())
            else:
                table = pd.read_csv(io.BytesIO(b'\n'.join(archive.lines())), dtype={'barcode': str})
                parts[key]['coords' if 'BeadLocations' in name else 'labels'] = table
    finally:
        archive.close()
    for key, d in parts.items():
        build_array(key, *d['matrix'], d['coords'], d['labels'])
