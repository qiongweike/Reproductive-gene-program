# Per-donor CellChat without unassigned-state spermatogonia: triMean, raw.use = TRUE, population.size = FALSE,
# 100 cell-label permutations, Kh 0.5, n 1, groups >= 30 cells, P < 0.05; seeds from the donor manifest.
# Runs all four donors, or only the donor given as the first argument (e.g. GSE202647_OA1).
# Inputs: results/communication/cellchat_inputs/ (steps 01-02)
# Outputs: results/communication/cellchat/<donor>/ (cellchat.rds, all_tested_LR.tsv.gz, all_SPG_LR.tsv.gz,
#          aggregate_network.tsv, all_pathway_probabilities.tsv, signaling_centrality.tsv)
source("utils/paths.R")
suppressPackageStartupMessages({library(CellChat); library(Matrix); library(data.table)})

MIN_CELLS <- 30
N_PERMUTATIONS <- 100
P_THRESHOLD <- 0.05
IN <- file.path(RESULTS, "communication", "cellchat_inputs")

fit_donor <- function(r) {
  m <- fread(file.path(IN, r$donor, "cells.tsv"))
  x <- as(readMM(file.path(IN, r$donor, "normalized.mtx")), "CsparseMatrix")
  rownames(x) <- fread(file.path(IN, "input_genes.tsv"))$gene
  colnames(x) <- m$cell_id
  stopifnot(ncol(x) == nrow(m), all(table(m$communication_group) >= MIN_CELLS), !any(m$communication_group == "SPG_unknown"))
  meta <- as.data.frame(m[, .(cell_id, labels = communication_group, samples = donor, diagnosis)])
  rownames(meta) <- meta$cell_id
  future::plan("sequential")
  set.seed(r$seed)
  cc <- createCellChat(object = x, meta = meta, group.by = "labels")
  cc@DB <- readRDS(file.path(IN, "protein_database.rds"))
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, thresh.pc = 0, thresh.fc = 0, thresh.p = P_THRESHOLD, only.pos = TRUE,
                                   min.cells = MIN_CELLS)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE, population.size = FALSE, distance.use = FALSE,
                          nboot = N_PERMUTATIONS, seed.use = r$seed, Kh = 0.5, n = 1)
  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc, thresh = P_THRESHOLD)
  cc <- aggregateNet(cc, thresh = P_THRESHOLD)
  netAnalysis_computeCentrality(cc, slot.name = "netP")
}

long_table <- function(a, names) {
  z <- as.data.table(as.table(a))
  setnames(z, names)
  z
}

write_donor <- function(cc, r, d) {
  saveRDS(cc, file.path(d, "cellchat.rds"))
  tab <- merge(long_table(cc@net$prob, c("source", "target", "interaction_name", "probability")),
               long_table(cc@net$pval, c("source", "target", "interaction_name", "p_cell_label_shuffle")),
               by = c("source", "target", "interaction_name"), all = TRUE)
  stopifnot(nrow(tab) == length(cc@net$prob), all(is.finite(tab$probability)), !anyNA(tab$p_cell_label_shuffle))
  tab[, donor := r$donor]
  tab[, diagnosis := r$diagnosis]
  tab[, q_BH_within_donor_all_tested := p.adjust(p_cell_label_shuffle, "BH")]
  tab[, cellchat_internal_supported := probability > 0 & p_cell_label_shuffle < P_THRESHOLD]
  tab <- merge(tab, as.data.table(cc@DB$interaction), by = "interaction_name", all.x = TRUE, sort = FALSE)
  tab[, focus_SPG_link := grepl("^SPG_", source) | grepl("^SPG_", target)]
  fwrite(tab, file.path(d, "all_tested_LR.tsv.gz"), sep = "\t")
  fwrite(tab[focus_SPG_link == TRUE], file.path(d, "all_SPG_LR.tsv.gz"), sep = "\t")
  groups <- names(table(cc@idents))
  net <- CJ(source = groups, target = groups)
  for (s in c("count", "weight")) net <- merge(net, long_table(cc@net[[s]], c("source", "target", s)),
                                               by = c("source", "target"), all = TRUE)
  net[, donor := r$donor]
  fwrite(net, file.path(d, "aggregate_network.tsv"), sep = "\t")
  fwrite(long_table(cc@netP$prob, c("source", "target", "pathway", "probability")),
         file.path(d, "all_pathway_probabilities.tsv"), sep = "\t")
  central <- rbindlist(lapply(names(cc@netP$centr), function(k) {
    z <- as.data.table(cc@netP$centr[[k]])
    stopifnot(nrow(z) == length(groups))
    z[, cell_group := groups]
    z[, pathway := k]
    z
  }), fill = TRUE)
  fwrite(central, file.path(d, "signaling_centrality.tsv"), sep = "\t")
  cat(r$donor, ": ", length(cc@idents), " cells, ", length(groups), " groups, ", dim(cc@net$prob)[3],
      " tested LR pairs, ", sum(tab$cellchat_internal_supported), " supported links\n", sep = "")
}

manifest <- fread(file.path(IN, "donor_manifest.tsv"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) manifest <- manifest[donor %in% args]
stopifnot(nrow(manifest) > 0)
for (i in seq_len(nrow(manifest))) {
  r <- manifest[i]
  write_donor(fit_donor(r), r, out_dir("communication", "cellchat", r$donor))
}
