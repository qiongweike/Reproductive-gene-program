# Virtual knockouts on one wild-type network: zero the target's outgoing row, align WT and KO networks (installed
# scTenifoldNet manifoldAlignment, d = 2) and compute differential regulation (authors' scTenifoldKnk dRegulation source).
# Distances <= sqrt(machine eps) x max|aligned coordinate| are flagged as numerical noise. Finished knockouts are skipped.
# Usage: Rscript knockout/06_run_knockouts.R <sample_key> <config>   (all 16 donor-configuration pairs)
# Inputs: results/knockout/networks/<config>/<key>/WT_final.rds, matching/<config>/perturbation_jobs.tsv, fixed_gene_sets.tsv
# Outputs: results/knockout/perturbations/<config>/<key>/<gene>_diffRegulation.tsv
suppressPackageStartupMessages({library(Matrix); library(data.table); library(scTenifoldNet); library(scTenifoldKnk)})
source("utils/paths.R")

SEEDS <- c(primary = 1, seed17 = 17, down80_seed29 = 29, paired_full_seed1 = 1, paired_U_seed1 = 1)
ALIGNMENT_DIM <- 2

args <- commandArgs(TRUE)
key <- args[1]
cfg <- args[2]
seed <- SEEDS[[cfg]]
ko_dir <- file.path(RESULTS, "knockout")
out <- out_dir("knockout", "perturbations", cfg, key)

WT <- readRDS(file.path(ko_dir, "networks", cfg, key, "WT_final.rds"))
jobs <- fread(file.path(ko_dir, "matching", cfg, "perturbation_jobs.tsv"))[sample_key == key]
mp3 <- fread(file.path(ko_dir, "fixed_gene_sets.tsv"))[set == "MP_3", gene]
stopifnot(nrow(jobs) > 0, length(mp3) == 49, all(jobs$perturbed_gene %in% rownames(WT)))
dr <- new.env(parent = asNamespace("scTenifoldKnk"))
sys.source(file.path(VENDOR, "scTenifoldKnk", "R", "dRegulation.R"), dr)

for (g in jobs$perturbed_gene) {
  path <- file.path(out, paste0(g, "_diffRegulation.tsv"))
  if (file.exists(path)) next
  KO <- WT
  KO[g, ] <- 0
  set.seed(seed)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = ALIGNMENT_DIM, nCores = 1)
  stopifnot(all(is.finite(MA)))
  DR <- as.data.table(dr$dRegulation(MA, gKO = g, empiricalNull = FALSE))
  threshold <- sqrt(.Machine$double.eps) * max(abs(MA))
  DR[, `:=`(config = cfg, sample_key = key, perturbed_gene = g, noise_threshold = threshold)]
  DR[, above_noise := distance > threshold]
  DR[, is_MP3 := gene %in% mp3]
  DR[, rank_score := frank(ifelse(above_noise, distance, 0), ties.method = "average") / .N]
  fwrite(DR, paste0(path, ".part"), sep = "\t")
  file.rename(paste0(path, ".part"), path)
}
message(cfg, " ", key, ": ", nrow(jobs), " knockouts")
