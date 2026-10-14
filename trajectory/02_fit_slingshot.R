# Root selection and per-donor Slingshot fits on the independent-feature PCA.
# Root: the K-means cluster with the highest U fraction (then U count, then label) among clusters with
# >=10 U cells and >=50% U; its early/late marker means are checked against D cells. No terminal cluster is forced.
# Inputs: results/trajectory/{donor_manifest.tsv, <donor>.primary.pca.tsv.gz}.
# Outputs: results/trajectory/{root_clusters.tsv, <donor>.primary.slingshot.rds, lineage_cell_pseudotime.tsv.gz}.
source("utils/paths.R")
suppressPackageStartupMessages({library(data.table); library(slingshot)})

SEED <- 20261003
N_PCS <- 20
ROOT_MIN_U <- 10
ROOT_MIN_FRACTION_U <- 0.5
STRETCH <- 0L
APPROX_POINTS <- 100L
MAXIT <- 15L
OUT <- out_dir("trajectory")

set.seed(SEED)
manifest <- fread(file.path(OUT, "donor_manifest.tsv"))
roots <- list(); pseudotime <- list(); status <- character()
for (sample in manifest$sample_key) {
  z <- fread(file.path(OUT, paste0(sample, ".primary.pca.tsv.gz")))
  clusters <- z[, .(n_cells = .N, n_U = sum(state == "undifferentiated"), fraction_U = mean(state == "undifferentiated"),
                    early = mean(early_marker), late = mean(late_marker)), by = cluster]
  setorder(clusters, -fraction_U, -n_U, cluster)
  candidates <- clusters[n_U >= ROOT_MIN_U & fraction_U >= ROOT_MIN_FRACTION_U]
  label_ok <- nrow(candidates) > 0
  chosen <- if (label_ok) candidates[1] else clusters[1]
  root <- chosen$cluster
  marker_ok <- chosen$early > mean(z[state == "differentiating", early_marker]) &
    chosen$late < mean(z[state == "differentiating", late_marker])
  clusters[, `:=`(sample_key = sample, scenario = "primary", chosen_root = cluster == root,
                  root_label_ok = label_ok, root_marker_ok = marker_ok)]
  roots[[sample]] <- clusters
  if (!label_ok) { status[sample] <- "no eligible root"; next }

  pcs <- as.matrix(z[, paste0("PC", seq_len(N_PCS)), with = FALSE]); rownames(pcs) <- z$cell_id
  fit <- tryCatch(slingshot(pcs, clusterLabels = z$cluster, start.clus = root, stretch = STRETCH,
                            approx_points = APPROX_POINTS, maxit = MAXIT), error = function(e) e)
  if (inherits(fit, "error")) { status[sample] <- paste("fit failed:", conditionMessage(fit)); next }
  saveRDS(fit, file.path(OUT, paste0(sample, ".primary.slingshot.rds")))
  tm <- slingPseudotime(fit, na = TRUE); wt <- slingCurveWeights(fit); lineages <- slingLineages(fit)
  for (k in seq_len(ncol(tm))) {
    pseudotime[[length(pseudotime) + 1]] <- data.table(
      cell_id = z$cell_id, sample_key = sample, scenario = "primary", lineage = colnames(tm)[k], pseudotime = tm[, k],
      curve_weight = wt[, k], root_cluster = root, root_marker_ok = marker_ok, cluster_path = paste(lineages[[k]], collapse = ";"))
  }
  status[sample] <- paste(ncol(tm), "lineage(s)")
}
fwrite(rbindlist(roots), file.path(OUT, "root_clusters.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(pseudotime), file.path(OUT, "lineage_cell_pseudotime.tsv.gz"), sep = "\t", na = "NA")
message(paste(sub(".*__", "", names(status)), status, sep = ": ", collapse = "; "))
