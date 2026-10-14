# Per-donor NMF of spermatogonia: Kinker et al. preprocessing, then snmf/r at ranks 6-9 with ten random starts per rank
# (7 donors x 4 ranks x 10 starts = 280 fits). Each start uses its own L'Ecuyer stream from rngtools::RNGseq().
# Input:  results/programs/nmf/donor_manifest.tsv and the exported count matrices (01_export_donor_matrices.py)
# Output: results/programs/nmf/starts/<donor>/K<rank>/start<NN>_fit.rds (one NMFfit per start)
# Usage:  Rscript programs/02_fit_nmf.R [donor [rank [start]]]
#         Without arguments all 280 fits run serially (large donors take hours per fit); run several processes with
#         donor/rank arguments to parallelize. Every fit depends only on its own stream, so results do not depend on the split.
source("utils/paths.R")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(NMF)})
RhpcBLASctl::blas_set_num_threads(1)  # the original fits used single-threaded BLAS; more threads change W and H at ~1e-13

RANKS <- 6:9
N_STARTS <- 10
SEED_BASE <- 20261001
METHOD <- "snmf/r"

args <- commandArgs(trailingOnly = TRUE)
nmf_dir <- file.path(RESULTS, "programs", "nmf")
manifest <- fread(file.path(nmf_dir, "donor_manifest.tsv"))[eligible == TRUE]

# Streams for donor i and rank k come from seed SEED_BASE + 1000 i + k. rngtools::RNGseq() seeds differently under a
# Mersenne-Twister and an L'Ecuyer-CMRG generator; in the original run (all 28 sets drawn in sequence in one fresh
# session) the first set (first donor, rank 6) was drawn under Mersenne-Twister and all later sets under L'Ecuyer-CMRG.
streams <- list()
for (i in seq_len(nrow(manifest))) for (k in RANKS) {
  RNGkind(if (i == 1 && k == RANKS[1]) "Mersenne-Twister" else "L'Ecuyer-CMRG")
  streams[[paste(manifest$sample_key[i], k)]] <- rngtools::RNGseq(N_STARTS, as.integer(SEED_BASE + 1000 * i + k))
}

# Kinker nmf_programs.R (log2(CPM/10 + 1), genes > 3.5 in > 2% of cells, gene centring, negatives set to zero) is run
# unchanged; its internal nmf() call is intercepted to return the prepared matrix.
kinker <- new.env()
sys.source(file.path(VENDOR, "CCLE_heterogeneity", "nmf_programs.R"), envir = kinker)
kinker$nmf <- function(x, ...) stop(structure(list(message = "", call = NULL, x = x), class = c("nmf_input", "condition")))
prepare_target <- function(d) {
  counts <- as(readMM(gzfile(d$matrix_file)), "CsparseMatrix")
  cpm <- as.matrix(counts %*% Diagonal(x = 1e6 / colSums(counts)))
  dimnames(cpm) <- list(fread(d$gene_file)$gene, fread(d$cell_file)$cell_id)
  tryCatch(kinker$nmf_programs(cpm, is.log = FALSE, rank = RANKS[1], method = METHOD, seed = SEED_BASE),
           nmf_input = function(e) e$x)
}

wanted <- function(value, position) length(args) < position || value == args[position]
for (i in seq_len(nrow(manifest))) {
  d <- manifest[i]
  if (!wanted(d$sample_key, 1)) next
  target <- prepare_target(d)
  for (k in RANKS) {
    if (!wanted(k, 2)) next
    start_dir <- out_dir("programs", "nmf", "starts", d$sample_key, paste0("K", k))
    for (j in seq_len(N_STARTS)) {
      if (!wanted(j, 3)) next
      rngtools::setRNG(streams[[paste(d$sample_key, k)]][[j]])
      fit <- nmf(target, rank = k, method = METHOD, seed = "random", nrun = 1, .pbackend = NA,
                 .options = list(keep.all = TRUE, parallel = FALSE))
      saveRDS(fit, file.path(start_dir, sprintf("start%02d_fit.rds", j)))
      message(d$sample_key, " K", k, " start ", j, ": ", niter(fit), " iterations, objective ", deviance(fit))
    }
  }
}
