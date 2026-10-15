# A recurrent gene program linking FGFR3 to negative regulators in human early spermatogonia

Analysis code for

A recurrent gene program linking FGFR3 to negative regulators in human early spermatogonia.

Starting from public single-cell RNA-seq and Slide-seq data, the scripts carry out every analysis in the paper and draw Figures 1–7 and Supplementary Figures 1–3.

## Analysis steps

The folders are listed in run order.

| Folder | Analysis |
|---|---|
| `preprocessing` | QC, doublet calls and cell identities for the five single-cell studies; the normal spermatogonial reference (9,080 cells, 30 donors) and its U/T/D states; the adult testis atlas; the two Slide-seq arrays |
| `programs` | Per-donor NMF of spermatogonia in the seven discovery donors (ranks 6–9, ten starts), recurrent meta-programs MP1–MP3, gene-identity harmonization and program scores |
| `coordination` | Within-state coordination of program members against expression-matched gene sets, removal of published state markers, spermatocyte and transcript-background controls, GO/Reactome enrichment |
| `trajectory` | Slingshot and Palantir trajectories on program-independent features, member dynamics along lineages, Hotspot local correlations |
| `clinical` | Captured state composition and fixed-background program scores in OA, NOA and cryptozoospermia donors |
| `communication` | CellChat and NicheNet ligand analysis, GRNBoost2/cisTarget regulons, recovery of motif-supported MP3 edges, regulon activity and dynamics |
| `spatial` | RCTD deconvolution, BANKSY clusters, spatial autocorrelation, neighbourhood analyses and measurement coverage in two normal Slide-seq arrays |
| `knockout` | scTenifoldNet networks and scTenifoldKnk virtual knockouts compared with matched non-member controls |
| `figures` | One script per figure, plus the shared theme |

## Repository layout

```
utils/           paths.py, paths.R: project directories
resources/       candidate_143_genes.txt     fixed 143-gene candidate set
                 excluded_module_genes.json  nine predefined gene modules (66 genes)
preprocessing/   01_qc_normal_references.py … 09_embed_atlas.py
programs/        01_export_donor_matrices.py … 07_cell_scores.py
coordination/    01_export_noa_spermatogonia.py … 08_test_enrichment.py
trajectory/      01_prepare_trajectory_inputs.py … 08_summarize_hotspot.py
clinical/        01_count_captured_states.py … 03_summarize_member_expression.py
communication/   01_export_cellchat_database.R … 14_regulon_dynamics.py
spatial/         01_export_spatial_inputs.py … 10_prepare_maps.py
knockout/        01_resolve_identities.py … 09_compute_response_pca.R
figures/         figure1.R … figure7.R, supplementary_figure1.R … supplementary_figure3.R, theme.R
```

Within each folder the scripts are numbered in run order; files without a number are helpers imported by the others. The candidate set and the nine modules are scored as fixed references and are kept out of trajectory features, background gene pools and knockout controls.

Scripts are run from the repository root. Data, third-party code and outputs live under a project root given by the environment variable `SPG_ROOT` (default: the repository root):

```
$SPG_ROOT/
  data/raw/         GEO supplementary files, one folder per accession, and Human.7z
  data/metadata/    GEO family SOFT files and author annotation tables
  data/external/    reference resources (below)
  vendor/           third-party code at fixed commits (below)
  data/processed/   written by preprocessing
  results/          written by the analysis steps, one folder per step; figures in results/figures/
```

## Data

### Single-cell and spatial data

| Source | Study | Use | Files |
|---|---|---|---|
| GEO GSE124263 | Sohni et al., *Cell Rep* 2019 | Discovery donors A1, A2; atlas; reference for state transfer | 10x files of GSM3526587–GSM3526590 → `data/raw/GSE124263/`; `GSE124263_Barcode_Cell_Identity.xlsx`, `GSE124263_family.soft.gz` → `data/metadata/` |
| GEO GSE153947 | Di Persio et al., *Cell Rep Med* 2021 | Discovery donors N1–N3; cryptozoospermia comparison | `*_counts.tsv.gz` of GSM4658980–GSM4658985 → `data/raw/GSE153947/`; `GSE153947_Cell_metadata.tsv.gz`, `GSE153947_family.soft.gz` → `data/metadata/` |
| GEO GSE254315 | Cui et al., *Nat Aging* 2025 | Discovery donors Donor5, Donor30; normal reference; atlas | 10x files of GSM8037606–GSM8037630 → `data/raw/GSE254315/`; `GSE254315_family.soft.gz` → `data/metadata/` |
| GEO GSE202647 | Wang et al., *eLife* 2025 | OA/NOA comparison; CellChat and NicheNet | 10x files of GSM6127766–GSM6127770 → `data/raw/GSE202647/` |
| GEO GSE326372 | NOA testes, five samples | External NOA cases | 10x files of GSM9630067–GSM9630071 → `data/raw/GSE326372/` |
| GEO GSE235321 | NOA and normal testis | Only its measured gene list, which enters the shared gene universe of the clinical scores | `GSE235321_Info.*.tsv.gz` and `GSE235321_exp.log2TPM.*.tsv.gz` (Normal, NOA1, NOA2) → `data/raw/GSE235321/` |
| [HaiqiChenLab/Slideseq_SSC_Niche](https://github.com/HaiqiChenLab/Slideseq_SSC_Niche) | Chen et al., *Cell Rep* 2021; Rajachandran et al., *Cell Rep* 2023 | Normal human Slide-seq arrays Puck5 and Puck6 | `Human.7z` → `data/raw/` |

Files keep their GEO names. For GSE124263, GSE153947 and GSE254315 the scripts find each sample's files through the GEO family SOFT file.

### Reference resources (`data/external/`)

| File | Source |
|---|---|
| `hgnc_complete_set.txt` | HGNC complete set ([download](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt)), retrieved 5 October 2026 |
| `Reactome/ReactomePathways.gmt` | [Reactome](https://reactome.org/download/current/ReactomePathways.gmt.zip), file dated 21 June 2026 |
| `nichenet/ligand_target_matrix_nsga2r_final.rds`, `nichenet/lr_network_human_21122021.rds` | NicheNet v2 networks, Zenodo [10.5281/zenodo.7074291](https://doi.org/10.5281/zenodo.7074291) |
| `cistarget/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`, `cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | [cisTarget databases](https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/) |
| `cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` | [cisTarget motif annotations](https://resources.aertslab.org/cistarget/motif2tf/) |
| `cistarget/allTFs_hg38.txt` | [cisTarget TF list](https://resources.aertslab.org/cistarget/tf_lists/) |
| `prior_gene_lists/guo2017_clusterA.txt` | Guo et al., *Cell Stem Cell* 2017, Table S3, sheet Gene_Clusters, cluster A (373 genes, one per line) |
| `prior_gene_lists/guo2018_state0_enriched.txt` | Guo et al., *Cell Res* 2018, Table S5, group 1 (490 genes, one per line) |
| `prior_gene_lists/bush2024_state_markers.tsv` | Bush et al., *Cells* 2024, Supplementary Table 4; columns `gene` and `state` (cell type or spermatogonial state) |

GO annotations come from the Bioconductor packages org.Hs.eg.db and GO.db, and the CellChat database from the CellChat package.

### Third-party code (`vendor/`)

Several steps call published analysis code directly. Clone each repository into `vendor/` at the commit used:

| `vendor/` folder | Repository | Commit | Used by |
|---|---|---|---|
| `CCLE_heterogeneity` | [gabrielakinker/CCLE_heterogeneity](https://github.com/gabrielakinker/CCLE_heterogeneity) | `063d17a6d554bedb1df47a95af8f6e0d4fc9ce85` | `programs/02_fit_nmf.R` |
| `3ca` | [tiroshlab/3ca](https://github.com/tiroshlab/3ca) | `21a784103528cdf15008761f9669e2025159d6d2` | `programs/04_meta_programs.R`, `programs/06_score_programs.R` |
| `scalop` | [jlaffy/scalop](https://github.com/jlaffy/scalop) | `d187f6ff89484425f1b8b7be22e2896b1d4869e6` | `programs/06_score_programs.R` |
| `nichenetr` | [saeyslab/nichenetr](https://github.com/saeyslab/nichenetr) | `66f90d5eeafef280b2b2f339b3fd70ffec1781dd` | `communication/05_run_nichenet.R` (gene aliases) |
| `CellChat` | [jinworks/CellChat](https://github.com/jinworks/CellChat) | `75253cd0c9e68410e6e721a6d3a0419a1d7e358f` | `figures/figure5.R` (signalling-role plot) |
| `scTenifoldNet` | [cailab-tamu/scTenifoldNet](https://github.com/cailab-tamu/scTenifoldNet) | `3103a53216c347a3d4c9660d63b03c2ca99b47a1` | `knockout/03_build_networks.R` |
| `scTenifoldKnk` | [cailab-tamu/scTenifoldKnk](https://github.com/cailab-tamu/scTenifoldKnk) | `608fc260283ab24874fe35dc996e1ac15baa9a4b` | `knockout/06_run_knockouts.R`, `knockout/09_compute_response_pca.R` |

```bash
mkdir -p "$SPG_ROOT/vendor" && cd "$SPG_ROOT/vendor"
git clone https://github.com/gabrielakinker/CCLE_heterogeneity && git -C CCLE_heterogeneity checkout 063d17a6d554bedb1df47a95af8f6e0d4fc9ce85
git clone https://github.com/tiroshlab/3ca && git -C 3ca checkout 21a784103528cdf15008761f9669e2025159d6d2
git clone https://github.com/jlaffy/scalop && git -C scalop checkout d187f6ff89484425f1b8b7be22e2896b1d4869e6
git clone https://github.com/saeyslab/nichenetr && git -C nichenetr checkout 66f90d5eeafef280b2b2f339b3fd70ffec1781dd
git clone https://github.com/jinworks/CellChat && git -C CellChat checkout 75253cd0c9e68410e6e721a6d3a0419a1d7e358f
git clone https://github.com/cailab-tamu/scTenifoldNet && git -C scTenifoldNet checkout 3103a53216c347a3d4c9660d63b03c2ca99b47a1
git clone https://github.com/cailab-tamu/scTenifoldKnk && git -C scTenifoldKnk checkout 608fc260283ab24874fe35dc996e1ac15baa9a4b
```

## Software

**Python 3.10**: numpy 1.26.4, pandas 2.3.3, scipy 1.15.3, anndata 0.11.4, scanpy 1.11.5, scikit-learn 1.7.2, statsmodels 0.14.6, h5py 3.15.1, scrublet 0.2.3, threadpoolctl 3.6.0, numba 0.63.1, leidenalg 0.11.0, python-igraph 1.0.0, umap-learn 0.5.9.post2.

Packages used by single steps, installed in separate directories and put on `PYTHONPATH` for those steps:

- hotspotsc 1.1.3: `trajectory/07_run_hotspot.py`, `spatial/06_compute_autocorrelation.py`
- palantir 1.4.4 with jax 0.4.38 and mellon 1.6.1: `trajectory/05_run_palantir.py`
- pyscenic 0.12.1, arboreto 0.1.6, ctxcore 0.2.0, dask and distributed 2023.9.2: `communication/08_run_grnboost2.py`, `09_run_cistarget.py`, `11_refit_mp3_edges.py`, `13_score_regulons.py`

**R 4.4.3**: data.table 1.18.4, Matrix 1.7.6, NMF 0.28, rngtools 1.5.2, Seurat 5.5.1, SeuratObject 5.4.0, slingshot 2.14.0, SingleCellExperiment 1.28.0, CellChat 2.2.0.9001, nichenetr 2.2.1.1, spacexr 2.2.1, Banksy 0.1.6 ([prabhakarlab/Banksy](https://github.com/prabhakarlab/Banksy) at `6130b3d08e19296f05e74dced039b025e7ae718b`, whose original API `spatial/04_cluster_banksy.R` uses), ranger 0.18.0, mistyR 1.18.0, RANN 2.6.2, igraph 2.3.3, AnnotationDbi 1.68.0, org.Hs.eg.db 3.20.0, GO.db 3.20.0, RhpcBLASctl 0.23.42, future 1.75.0, ggplot2 4.0.3, patchwork 1.3.2, ggrepel 0.9.8, ggbeeswarm 0.7.3, ComplexHeatmap 2.22.0, circlize 0.4.18, scales 1.4.0, gtable 0.3.6, pdftools 3.9.0, png 0.1.9.

System libraries: libarchive (used by `preprocessing/07_extract_slideseq.py` to read `Human.7z`), Poppler (for pdftools) and the Liberation Sans font (figures).

## Running the analysis

```bash
export SPG_ROOT=/path/to/project              # holds data/ and vendor/
cd /path/to/this/repository
```

All seeds are set in the scripts. The published results are reproduced with four BLAS threads for `preprocessing` and `communication` and one thread for every other folder, as set in the commands below. Several preprocessing scripts request four BLAS threads, which OpenBLAS cannot supply once it has started with one. `programs/02_fit_nmf.R` and `communication/08`, `11` and `13` set single-threaded BLAS themselves. The Palantir and Hotspot steps were run with

```bash
# trajectory/05_run_palantir.py
PYTHONHASHSEED=20261003 JAX_PLATFORMS=cpu XLA_FLAGS='--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1'
# trajectory/07_run_hotspot.py, spatial/06_compute_autocorrelation.py
PYTHONHASHSEED=20261003 NUMBA_NUM_THREADS=1
```

Each step reads the outputs of earlier ones: `programs` needs `preprocessing`; `coordination`, `clinical` and `spatial` need `preprocessing` and `programs`; `trajectory` needs `programs`; `communication` needs `preprocessing`, `programs` and `trajectory`; `knockout` needs `preprocessing`, `programs` and `coordination`; `figures` needs all of them.

```bash
# Preprocessing
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
for s in preprocessing/0*.py; do python "$s"; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Programs
python  programs/01_export_donor_matrices.py
Rscript programs/02_fit_nmf.R                 # all donors, ranks and starts; one job: <sample_key> [rank [start]]
Rscript programs/03_select_best_fits.R
Rscript programs/04_meta_programs.R
python  programs/05_harmonize_gene_identity.py
Rscript programs/06_score_programs.R
python  programs/07_cell_scores.py

# Coordination
python  coordination/01_export_noa_spermatogonia.py
python  coordination/02_test_member_coordination.py
python  coordination/03_remove_state_markers.py
Rscript coordination/04_export_matched_controls.R
python  coordination/05_score_spermatocytes.py
python  coordination/06_adjust_transcript_background.py
Rscript coordination/07_build_pathway_membership.R
python  coordination/08_test_enrichment.py

# Trajectory
python  trajectory/01_prepare_trajectory_inputs.py
Rscript trajectory/02_fit_slingshot.R
python  trajectory/03_summarize_lineages.py
python  trajectory/04_export_member_dynamics.py
python  trajectory/05_run_palantir.py          # palantir on PYTHONPATH
python  trajectory/06_summarize_palantir.py
python  trajectory/07_run_hotspot.py           # hotspotsc on PYTHONPATH
python  trajectory/08_summarize_hotspot.py

# Clinical comparisons
for s in clinical/0*.py; do python "$s"; done

# Communication and regulons
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript communication/01_export_cellchat_database.R
python  communication/02_prepare_cellchat_inputs.py
Rscript communication/03_run_cellchat.R        # all four GSE202647 donors; one donor: GSE202647_OA1
python  communication/04_prepare_nichenet_inputs.py
Rscript communication/05_run_nichenet.R
python  communication/06_summarize_signalling_expression.py
python  communication/07_prepare_grn_input.py
python  communication/08_run_grnboost2.py      # 08, 09, 11 and 13 with pyscenic/arboreto on PYTHONPATH
python  communication/09_run_cistarget.py
python  communication/10_select_structural_regulons.py
python  communication/11_refit_mp3_edges.py
python  communication/12_summarize_edge_recovery.py
python  communication/13_score_regulons.py
python  communication/14_regulon_dynamics.py
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Spatial
python  spatial/01_export_spatial_inputs.py
for a in Puck5 Puck6; do Rscript spatial/02_run_rctd.R $a; done   # RCTD with all genes and without MP genes
python  spatial/03_score_programs.py
Rscript spatial/04_cluster_banksy.R
python  spatial/05_summarize_identity.py
python  spatial/06_compute_autocorrelation.py  # hotspotsc on PYTHONPATH
python  spatial/07_analyze_neighbourhood.py
Rscript spatial/08_fit_random_forest.R
python  spatial/09_summarize_measurements.py
python  spatial/10_prepare_maps.py

# Virtual knockouts: 16 donor-configuration networks
ko_jobs() {
  for k in GSE124263__A1 GSE124263__A2 GSE153947__N2 GSE153947__N3; do
    for c in primary seed17 down80_seed29; do echo "$k $c"; done
  done
  for k in GSE124263__A1 GSE124263__A2; do
    for c in paired_full_seed1 paired_U_seed1; do echo "$k $c"; done
  done
}
python  knockout/01_resolve_identities.py
python  knockout/02_prepare_inputs.py
ko_jobs | while read k c; do Rscript knockout/03_build_networks.R $k $c; done
python  knockout/04_match_controls.py
python  knockout/05_match_paired_controls.py
ko_jobs | while read k c; do Rscript knockout/06_run_knockouts.R $k $c; done
python  knockout/07_compute_readouts.py
python  knockout/08_prepare_panels.py
Rscript knockout/09_compute_response_pca.R

# Figures
for f in figure1 figure2 figure3 figure4 figure5 figure6 figure7 \
         supplementary_figure1 supplementary_figure2 supplementary_figure3; do Rscript figures/$f.R; done
```

The NMF fits, RCTD fits, network builds and knockouts are independent jobs and can run in parallel.

## Figures

Each script writes `results/figures/<name>.pdf` and a 300-dpi PNG.

| Figure | Script | Results used |
|---|---|---|
| Figure 1 | `figures/figure1.R` | preprocessing (atlas), programs, coordination, trajectory |
| Figure 2 | `figures/figure2.R` | programs, coordination, trajectory (Hotspot) |
| Figure 3 | `figures/figure3.R` | programs, trajectory |
| Figure 4 | `figures/figure4.R` | programs, coordination, clinical, communication (regulons) |
| Figure 5 | `figures/figure5.R` | communication (CellChat, NicheNet, regulons), NicheNet prior |
| Figure 6 | `figures/figure6.R` | spatial |
| Figure 7 | `figures/figure7.R` | knockout |
| Supplementary Figure 1 | `figures/supplementary_figure1.R` | preprocessing (atlas), coordination |
| Supplementary Figure 2 | `figures/supplementary_figure2.R` | trajectory, communication (regulon dynamics) |
| Supplementary Figure 3 | `figures/supplementary_figure3.R` | knockout |

## Contact

Zhiling Li (stlizhiling@126.com), Reproductive Center, The First Affiliated Hospital of Shantou University Medical College, Shantou, China.
