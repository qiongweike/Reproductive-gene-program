# Export the matched control genes drawn by program scoring, with the number of times each was drawn
# (repeated draws are kept, as in scalop::baseScores), so other cells can be scored with the same controls.
# Inputs: results/programs/{donor_manifest,MP_genes}.tsv, results/programs/scores/<donor>_matched_control_genes.rds
# Output: results/coordination/specificity/matched_control_genes.tsv
source("utils/paths.R")
suppressPackageStartupMessages(library(data.table))

PROGRAM_SIZES <- c(MP_1 = 50L, MP_2 = 50L, MP_3 = 49L)

donors <- fread(file.path(RESULTS, "programs", "donor_manifest.tsv"))[eligible == TRUE, sample_key]
members <- fread(file.path(RESULTS, "programs", "MP_genes.tsv"))
programs <- lengths(split(members$gene, members$MP))
stopifnot(identical(programs, PROGRAM_SIZES))
score_dir <- file.path(RESULTS, "programs", "scores")
rows <- list()
for (key in donors) {
  controls <- readRDS(file.path(score_dir, paste0(key, "_matched_control_genes.rds")))
  header <- names(fread(file.path(score_dir, paste0(key, "_scores.tsv.gz")), nrows = 1))
  stopifnot(length(controls) == length(programs), identical(header[-1], names(programs)))
  for (j in seq_along(controls)) {
    tab <- as.data.table(table(controls[[j]]))
    setnames(tab, c("gene", "multiplicity"))
    tab[, `:=`(sample_key = key, MP = names(programs)[j], n_control_draws = length(controls[[j]]))]
    rows[[length(rows) + 1L]] <- tab
  }
}
out <- rbindlist(rows)
fwrite(out, file.path(out_dir("coordination", "specificity"), "matched_control_genes.tsv"), sep = "\t")
cat(length(donors), "donors;", nrow(out), "control gene rows\n")
