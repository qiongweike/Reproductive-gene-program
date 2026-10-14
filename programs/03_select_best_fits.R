# Keep, for each donor and rank, the start with the lowest objective among its ten ordered starts. Fits that reached the
# 20,000-iteration cap are kept and flagged.
# Input:  results/programs/nmf/starts/<donor>/K<rank>/start<NN>_fit.rds (02_fit_nmf.R)
# Output: results/programs/nmf/donors/<donor>/all_ranks_basis_usage.rds (W and H of the selected fits, Gavish column names),
#         results/programs/nmf/all_start_metrics.tsv, results/programs/nmf/NMF_target_gene_universe.tsv (genes entering NMF)
source("utils/paths.R")
suppressPackageStartupMessages({library(data.table); library(NMF)})

RANKS <- 6:9
N_STARTS <- 10
SEED_BASE <- 20261001
MAX_ITER <- 20000

nmf_dir <- file.path(RESULTS, "programs", "nmf")
manifest <- fread(file.path(nmf_dir, "donor_manifest.tsv"))[eligible == TRUE]
metrics <- list(); universe <- list()
for (i in seq_len(nrow(manifest))) {
  d <- manifest[i]; W <- list(); H <- list()
  for (k in RANKS) {
    files <- file.path(nmf_dir, "starts", d$sample_key, paste0("K", k), sprintf("start%02d_fit.rds", seq_len(N_STARTS)))
    fits <- lapply(files, readRDS)
    objective <- vapply(fits, deviance, numeric(1)); best <- which.min(objective)
    niters <- as.integer(vapply(fits, niter, numeric(1)))
    w <- basis(fits[[best]]); h <- t(coef(fits[[best]]))
    colnames(w) <- colnames(h) <- paste0(d$sample_key, "_rank6_9_nruns10.RDS.", k, ".", seq_len(k))
    W[[length(W) + 1]] <- w; H[[length(H) + 1]] <- h
    metrics[[length(metrics) + 1]] <- data.table(sample_key = d$sample_key, K = k, start = seq_len(N_STARTS),
      rank_seed = as.integer(SEED_BASE + 1000 * i + k), niter = niters, hit_iteration_limit = niters >= MAX_ITER,
      objective = objective, selected_best = seq_len(N_STARTS) == best)
  }
  saveRDS(list(w_basis = do.call(cbind, W), h_coef = do.call(cbind, H), sample_key = d$sample_key),
          file.path(out_dir("programs", "nmf", "donors", d$sample_key), "all_ranks_basis_usage.rds"))
  universe[[i]] <- data.table(sample_key = d$sample_key, study = d$study, donor = d$donor, gene = rownames(W[[1]]))
}
metrics <- rbindlist(metrics)
fwrite(metrics, file.path(nmf_dir, "all_start_metrics.tsv"), sep = "\t")
fwrite(rbindlist(universe), file.path(nmf_dir, "NMF_target_gene_universe.tsv"), sep = "\t")
message(nrow(metrics), " starts; ", sum(metrics$hit_iteration_limit), " at the iteration cap (",
        sum(metrics$hit_iteration_limit & metrics$selected_best), " selected)")
