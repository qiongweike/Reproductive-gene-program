# RCTD (spacexr) in doublet mode on one Slide-seq array with the strict adult reference: the original fit on all
# genes and an MP-excluded fit with every meta-program member removed from spatial and reference counts (full UMI kept).
# Inputs: data/processed/atlas/spatial/<array>_* and reference_strict_*, results/programs/nmf/MP_genes.tsv.
# Outputs: results/spatial/rctd_original/{weights,classes,doublet_weights}_<array>_strict_doublet.tsv.gz and
#          results/spatial/neighbourhood/<array>_{weights_alltypes,classes,doublet_weights}.tsv.gz.
# Usage: Rscript spatial/02_run_rctd.R <array> [original|mp_excluded]   (both variants by default)
source("utils/paths.R")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(spacexr)})

SEEDS <- list(original = c(Puck5 = 20261001, Puck6 = 20261001), mp_excluded = c(Puck5 = 20261004, Puck6 = 20261005))
MAX_CORES <- 4
# Reference type order (RCTD's internal and output order); fixed so it does not depend on the locale's collation.
TYPES <- c("Endothelial", "Immune", "Leydig", "Peritubular_smooth_muscle", "Sertoli", "Spermatid", "Spermatocyte", "SPG")
CELL_MIN_INSTANCE <- 25
IN <- file.path(PROCESSED, "atlas", "spatial")

args <- commandArgs(trailingOnly = TRUE)
array <- args[1]
variants <- if (length(args) > 1) args[2] else names(SEEDS)
stopifnot(array %in% c("Puck5", "Puck6"), all(variants %in% names(SEEDS)))

read_counts <- function(stem) {
  m <- as(readMM(gzfile(file.path(IN, paste0(stem, "_counts.mtx.gz")))), "CsparseMatrix")
  rownames(m) <- fread(file.path(IN, paste0(stem, "_genes.tsv.gz")))$gene
  m
}

co <- fread(file.path(IN, paste0(array, "_coords.tsv.gz")))
spatial_counts <- read_counts(array); colnames(spatial_counts) <- co$barcode
stopifnot(all(colSums(spatial_counts) == co$total_counts))
coords <- data.frame(x = co$x_native, y = co$y_native, row.names = co$barcode)
numi <- setNames(co$total_counts, co$barcode)
meta <- fread(file.path(IN, "reference_strict_meta.tsv.gz"))
reference_counts <- read_counts("reference_strict"); colnames(reference_counts) <- meta$barcode
stopifnot(setequal(meta$CellType, TYPES))
cell_types <- setNames(factor(meta$CellType, levels = TYPES), meta$barcode)
reference_numi <- setNames(meta$full_assay_total_counts, meta$barcode)
stopifnot(all(table(cell_types) >= CELL_MIN_INSTANCE))
mp_genes <- toupper(fread(file.path(RESULTS, "programs", "nmf", "MP_genes.tsv"))$gene)

fit_rctd <- function(drop_genes, seed) {
  sc <- spatial_counts[!toupper(rownames(spatial_counts)) %in% drop_genes, , drop = FALSE]
  rc <- reference_counts[!toupper(rownames(reference_counts)) %in% drop_genes, , drop = FALSE]
  set.seed(seed)
  rctd <- create.RCTD(SpatialRNA(coords, sc, numi), Reference(rc, cell_types, reference_numi),
                      max_cores = MAX_CORES, CELL_MIN_INSTANCE = CELL_MIN_INSTANCE, keep_reference = TRUE)
  stopifnot(!any(toupper(c(rctd@internal_vars$gene_list_reg, rctd@internal_vars$gene_list_bulk)) %in% drop_genes))
  run.RCTD(rctd, doublet_mode = "doublet")
}

write_fit <- function(fit, paths) {
  w <- as.matrix(normalize_weights(fit@results$weights))
  stopifnot(all(is.finite(w)), all(w >= 0), max(abs(rowSums(w) - 1)) < 1e-6)
  fwrite(as.data.table(w, keep.rownames = "barcode"), paths[1], sep = "\t")
  fwrite(as.data.table(fit@results$results_df, keep.rownames = "barcode"), paths[2], sep = "\t")
  fwrite(as.data.table(as.matrix(fit@results$weights_doublet), keep.rownames = "barcode"), paths[3], sep = "\t")
}

for (variant in variants) {
  if (variant == "original") {
    fit <- fit_rctd(character(0), SEEDS$original[[array]])
    paths <- file.path(out_dir("spatial", "rctd_original"),
                       paste0(c("weights_", "classes_", "doublet_weights_"), array, "_strict_doublet.tsv.gz"))
  } else {
    fit <- fit_rctd(mp_genes, SEEDS$mp_excluded[[array]])
    paths <- file.path(out_dir("spatial", "neighbourhood"),
                       paste0(array, c("_weights_alltypes", "_classes", "_doublet_weights"), ".tsv.gz"))
  }
  write_fit(fit, paths)
  classes <- table(fit@results$results_df$spot_class)
  message(array, " ", variant, ": ", paste(names(classes), classes, collapse = ", "))
  rm(fit); gc()
}
