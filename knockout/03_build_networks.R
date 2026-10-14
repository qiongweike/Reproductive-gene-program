# Wild-type network of one donor and configuration: ten 500-cell principal-component-regression networks (authors'
# scTenifoldNet pcNet/makeNetworks source) on full-assay CPM, CP tensor denoising and strict direction (installed
# scTenifoldNet/scTenifoldKnk); diagonal removed, rows = outgoing edges. Run for all 16 donor-configuration pairs.
# Usage: Rscript knockout/03_build_networks.R <sample_key> <config>
# Inputs: results/knockout/inputs/<key>/ (inputs_paired/<config>/<key>/ for paired configurations)
# Outputs: results/knockout/networks/<config>/<key>/{WT_final.rds, WT_node_features.tsv, analysis_cells.tsv}
suppressPackageStartupMessages({library(Matrix); library(data.table); library(scTenifoldNet); library(scTenifoldKnk)})
source("utils/paths.R")

CONFIGS <- list(primary = c(seed = 1, fraction = 1), seed17 = c(seed = 17, fraction = 1),
                down80_seed29 = c(seed = 29, fraction = 0.8),
                paired_full_seed1 = c(seed = 1, fraction = 1), paired_U_seed1 = c(seed = 1, fraction = 1))
N_NET <- 10
N_CELLS <- 500
N_COMP <- 3
EDGE_QUANTILE <- 0.9
TD_RANK <- 3
TD_MAX_ERROR <- 1e-5
TD_MAX_ITER <- 1000
TD_DECIMALS <- 3
DIRECTION_LAMBDA <- 0

args <- commandArgs(TRUE)
key <- args[1]
cfg <- args[2]
seed <- CONFIGS[[cfg]][["seed"]]
fraction <- CONFIGS[[cfg]][["fraction"]]
input <- if (startsWith(cfg, "paired_")) file.path(RESULTS, "knockout", "inputs_paired", cfg, key) else
  file.path(RESULTS, "knockout", "inputs", key)
out <- out_dir("knockout", "networks", cfg, key)

net <- new.env(parent = asNamespace("scTenifoldNet"))
for (f in c("pcNet.R", "makeNetworks.R")) sys.source(file.path(VENDOR, "scTenifoldNet", "R", f), net)

genes <- fread(file.path(input, "genes.tsv"))
cells <- fread(file.path(input, "cells.tsv"))
counts <- readMM(gzfile(file.path(input, "counts.mtx.gz")))
rownames(counts) <- genes$gene
colnames(counts) <- cells$cell_id
stopifnot(!anyDuplicated(rownames(counts)), all(colSums(counts) <= cells$full_assay_total_counts + 1e-6))
if (fraction < 1) {
  # Stratified subsample: floor(fraction x n) cells of every state, without replacement.
  set.seed(seed)
  pick <- sort(unlist(lapply(split(seq_len(nrow(cells)), cells$state),
                             function(i) sample(i, max(1, floor(length(i) * fraction)), replace = FALSE))))
  cells <- cells[pick]
  counts <- counts[, pick, drop = FALSE]
}
fwrite(cells[, .(cell_id, state, full_assay_total_counts)], file.path(out, "analysis_cells.tsv"), sep = "\t")

X <- counts %*% Diagonal(x = 1e6 / cells$full_assay_total_counts)
dimnames(X) <- dimnames(counts)
node <- copy(genes)
node[, mean_CP10k := as.numeric(rowMeans(X)) / 100]
node[, detection_fraction := as.numeric(rowMeans(counts > 0))]
L <- X / 100
L@x <- log1p(L@x)
mu <- rowMeans(L)
node[, variance_logCP10k := pmax(as.numeric(rowMeans(L^2)) - mu^2, 0)]

set.seed(seed)
nets <- net$makeNetworks(X, q = EDGE_QUANTILE, nNet = N_NET, nCells = N_CELLS, nComp = N_COMP,
                         scaleScores = TRUE, symmetric = FALSE, nCores = 1)
set.seed(seed)
td <- scTenifoldNet::tensorDecomposition(nets, K = TD_RANK, maxError = TD_MAX_ERROR, maxIter = TD_MAX_ITER,
                                         nDecimal = TD_DECIMALS)
WT <- as.matrix(scTenifoldKnk:::strictDirection(td$X, lambda = DIRECTION_LAMBDA))
diag(WT) <- 0
WT <- t(WT)
stopifnot(identical(rownames(WT), genes$gene), all(is.finite(WT)), max(abs(WT)) > 0)
saveRDS(WT, file.path(out, "WT_final.rds"), compress = FALSE)

node$outdegree <- rowSums(WT != 0)
node$outstrength <- rowSums(abs(WT))
fwrite(node, file.path(out, "WT_node_features.tsv"), sep = "\t")
message(cfg, " ", key, ": ", nrow(WT), " genes, ", nrow(cells), " cells")
