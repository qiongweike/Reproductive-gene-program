"""Bin means of all 149 program members along every Slingshot lineage, member peak bins and donor-state program means.
Member expression: log1p(CP10k), full-library denominator. Member peak: earliest bin of maximum mean, for members measured
in ten bins of >=10 cells; dynamic range <=1e-12 counts as constant. Program peak: median peak of non-constant members.
Inputs: results/trajectory/{donor_manifest.tsv, lineage_bin_cells.tsv.gz, trajectory_bin_scores.tsv, fixed_cell_scores.tsv.gz},
        results/programs/MP_genes.tsv, donor count matrices.
Outputs: results/trajectory/{member_trajectory_bins.tsv.gz, member_peak_order.tsv, program_member_peak_summary.tsv, MP_donor_state_scores.tsv}.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.paths import RESULTS, out_dir

import numpy as np
import pandas as pd
from scipy import io

N_BINS = 10
CONSTANT_RANGE = 1e-12
MIN_STATE_CELLS = 30
PROGRAMS = RESULTS / 'programs'
OUT = out_dir('trajectory')


def donor_files(sample):
    return (PROGRAMS / 'nmf/input' / f'{sample}.counts.mtx.gz', PROGRAMS / 'input' / f'{sample}.genes.tsv',
            PROGRAMS / 'nmf/input' / f'{sample}.cells.tsv')


def save(rows, name):
    pd.DataFrame(rows).to_csv(OUT / name, sep='\t', index=False, na_rep='NA')


def member_bin_means(manifest, mp, cells, bins):
    rows = []
    for m in manifest.itertuples():
        bb = bins[bins.sample_key.eq(m.sample_key)]
        if not len(bb):
            continue
        counts_file, gene_file, cell_file = donor_files(m.sample_key)
        raw = io.mmread(counts_file).T.tocsr()
        upper = pd.read_csv(gene_file, sep='\t').gene.str.upper()
        dup = upper.duplicated(False)
        gi = {g: i for i, g in enumerate(upper) if not dup.iloc[i]}
        ci = {v: i for i, v in enumerate(pd.read_csv(cell_file, sep='\t').cell_id)}
        lib = np.asarray(raw.sum(1)).ravel()
        selected = {g: raw[:, gi[g]].toarray().ravel() for g in mp.gene if g in gi}
        for br in bb.itertuples():
            cc = cells[cells.sample_key.eq(m.sample_key) & cells.lineage.eq(br.lineage) & cells.bin.eq(br.bin)]
            ix = np.array([ci[x] for x in cc.cell_id])
            assert len(ix) == br.n_cells
            for mm in mp.itertuples():
                v = selected.get(mm.gene)
                v = None if v is None else v[ix]
                rows.append(dict(sample_key=m.sample_key, study=m.study, donor=m.donor, lineage=br.lineage, bin=br.bin,
                                 time_scaled=br.time_scaled, n_cells=len(ix), trajectory_eligible=br.eligible,
                                 bin_eligible=br.bin_eligible, MP=mm.MP, gene=mm.gene, rank=mm.rank,
                                 mean_logCP10k=float(np.log1p(v * 1e4 / lib[ix]).mean()) if v is not None else np.nan,
                                 detection=float(np.mean(v > 0)) if v is not None else np.nan))
    return pd.DataFrame(rows)


def member_peaks(tr):
    peaks = []
    for (sample, lineage, program, gene), x in tr.groupby(['sample_key', 'lineage', 'MP', 'gene'], sort=False):
        valid = len(x) == N_BINS and x.bin_eligible.all() and x.mean_logCP10k.notna().all()
        span = float(x.mean_logCP10k.max() - x.mean_logCP10k.min()) if valid else np.nan
        top = x.mean_logCP10k.max()
        ties = x.loc[x.mean_logCP10k.eq(top), 'bin'].sort_values().tolist() if valid else []
        nonconstant = bool(valid and span > CONSTANT_RANGE)
        peaks.append(dict(sample_key=sample, lineage=lineage, MP=program, gene=gene, peak_bin=min(ties) if ties else np.nan,
                          tied_peak_bins=';'.join(map(str, ties)), n_tied_peak=len(ties), dynamic_range=span,
                          nonconstant=nonconstant, trajectory_eligible=bool(x.trajectory_eligible.iloc[0]),
                          reason='PASS' if nonconstant else ('CONSTANT_MEMBER' if valid else 'MISSING_OR_INSUFFICIENT_BINS')))
    return pd.DataFrame(peaks)


def main():
    manifest = pd.read_csv(OUT / 'donor_manifest.tsv', sep='\t')
    mp = pd.read_csv(PROGRAMS / 'MP_genes.tsv', sep='\t')
    mp.gene = mp.gene.str.upper()
    cells = pd.read_csv(OUT / 'lineage_bin_cells.tsv.gz', sep='\t')
    bins = pd.read_csv(OUT / 'trajectory_bin_scores.tsv', sep='\t')
    bins = bins[bins.metric.eq('MP_2')]

    tr = member_bin_means(manifest, mp, cells, bins)
    save(tr, 'member_trajectory_bins.tsv.gz')
    pk = member_peaks(tr)
    save(pk, 'member_peak_order.tsv')
    summary = pk[pk.nonconstant].groupby(['sample_key', 'lineage', 'MP', 'trajectory_eligible']).agg(
        median_member_peak=('peak_bin', 'median'), n_nonconstant=('gene', 'size')).reset_index()
    save(summary, 'program_member_peak_summary.tsv')

    sc = pd.read_csv(OUT / 'fixed_cell_scores.tsv.gz', sep='\t')
    state = sc.groupby(['sample_key', 'study', 'donor', 'state']).agg(
        n_cells=('cell_id', 'size'), MP_1=('MP_1', 'mean'), MP_2=('MP_2', 'mean'), MP_3=('MP_3', 'mean')).reset_index()
    state['eligible'] = state.n_cells.ge(MIN_STATE_CELLS)
    save(state, 'MP_donor_state_scores.tsv')
    q = summary[summary.trajectory_eligible].pivot_table(index=['sample_key', 'lineage'], columns='MP', values='median_member_peak')
    print('MP1/MP2/MP3 median member peak bins on qualified lineages:',
          '; '.join(f"{s} {l} {'/'.join(f'{v:g}' for v in r)}" for (s, l), r in q.iterrows()))


if __name__ == '__main__':
    main()
