# Score MP1-MP3 in every spermatogonium with the Gavish et al. (2023) Score_cells block and the scalop implementation:
# log2(CPM/10 + 1) with full-assay CPM, gene centring within donor, 30 expression bins and 100 matched controls per
# member. Assignment needs >= 25 measured members, a unique maximum score >= 1 and an MP frequency >= 5% among the
# preliminary assignments; every input cell stays in the denominator. Discovery donors also export measured expression.
# Input:  results/programs/{MP_genes.tsv, donor_manifest.tsv, NOA_assignment_manifest.tsv} and the exported count matrices
# Output: results/programs/scores/<sample>_{scores.tsv.gz, matched_control_genes.rds},
#         results/programs/measured_expression/<sample>_assigned_within_donor_centered.tsv.gz, results/programs/MP_*.tsv(.gz)
source("utils/paths.R")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(magrittr)})

SEED_BASE <- 20261001
MIN_GENES <- 25
MIN_SCORE <- 1
MIN_CELL_FRACTION <- 0.05

prog_dir <- file.path(RESULTS, "programs")
mp_table <- fread(file.path(prog_dir, "MP_genes.tsv"))
MPs <- split(mp_table$gene, mp_table$MP)
stopifnot(lengths(MPs)[["MP_3"]] == 49L)
manifest <- rbindlist(list(
  fread(file.path(prog_dir, "donor_manifest.tsv"))[eligible == TRUE][, cohort := "normal_discovery"],
  fread(file.path(prog_dir, "NOA_assignment_manifest.tsv"))[, cohort := "NOA_case_description"]), fill = TRUE)
score_dir <- out_dir("programs", "scores")
expr_dir <- out_dir("programs", "measured_expression")

# scalop sigScores; binmatch() is wrapped to keep each sample's matched control genes.
scalop <- new.env(parent = globalenv())
for (f in c("utils-bin.R", "bin.R", "matrix-ops.R", "sigScores.R")) sys.source(file.path(VENDOR, "scalop", "R", f), envir = scalop)
controls <- new.env()
scalop_binmatch <- scalop$binmatch
scalop$binmatch <- function(...) { z <- scalop_binmatch(...); controls$sets <- c(controls$sets, list(z)); z }

# Gavish MP_distribution.R lines 36-75 (Score_cells). Adaptations: use the scalop copy above, keep a one-row matrix,
# return the unfiltered scores and the preliminary assignments, and leave equal maxima unassigned ("AMBIGUOUS_TIE").
patch <- function(code, from, to) { stopifnot(grepl(from, code, fixed = TRUE)); gsub(from, to, code, fixed = TRUE) }
score_cells_code <- paste(readLines(file.path(VENDOR, "3ca", "ITH_hallmarks", "MPs_distribution", "MP_distribution.R"))[36:75], collapse = "\n") |>
  patch("scalop::sigScores", "scalop$sigScores") |>
  patch("x[-cells_rm , ]", "x[-cells_rm , ,drop=FALSE]") |>
  patch("MP_scores_per_sample <- lapply(MP_scores_per_sample, function(x) remove_cells(x,MinScore))",
        "MP_scores_unfiltered <- MP_scores_per_sample\n  MP_scores_per_sample <- lapply(MP_scores_per_sample, function(x) remove_cells(x,MinScore))") |>
  patch("function(x) apply(x, 1, function(y) colnames(x)[which(y==max(y))] )",
        "function(x) if(nrow(x)==0) setNames(character(),character()) else apply(x, 1, function(y) if(sum(y==max(y))==1) colnames(x)[which.max(y)] else \"AMBIGUOUS_TIE\")") |>
  patch("return(Assign_MP_per_cell_filtered)",
        "return(list(scores=MP_scores_unfiltered,above_threshold=Assign_MP_per_cell,retained=Assign_MP_per_cell_filtered))")
eval(parse(text = score_cells_code))

assignments <- list(); counts_all <- list(); coverage <- list(); expr_manifest <- list(); expr_order <- list()
for (i in seq_len(nrow(manifest))) {
  d <- manifest[i]
  cells <- fread(d$cell_file)
  genes <- fread(d$gene_file)$gene
  out <- data.table(cell_id = cells$cell_id, sample_key = d$sample_key, study = d$study, donor = d$donor, cohort = d$cohort,
                    assignment = "Unassigned", reason = "maximum_score_below_1", best_score = NA_real_, above_threshold_MP = NA_character_)
  cv <- data.table(MP = names(MPs), n_genes = lengths(MPs), n_available = vapply(MPs, function(g) sum(g %in% genes), integer(1)),
                   sample_key = d$sample_key, study = d$study)
  cv[, score_eligible := n_available >= MIN_GENES]
  coverage[[i]] <- cv
  usable <- lapply(MPs[cv[score_eligible == TRUE]$MP], intersect, y = genes)
  stopifnot(length(usable) > 0)

  counts <- as(readMM(gzfile(d$matrix_file)), "CsparseMatrix")
  dimnames(counts) <- list(genes, cells$cell_id)
  depth <- as.numeric(cells$total_counts)  # full-assay library size, also for the NOA scoring subset
  stopifnot(all(colSums(counts) <= depth + 1e-6))
  cpm <- as.matrix(counts %*% Diagonal(x = 1e6 / depth))
  dimnames(cpm) <- dimnames(counts)
  set.seed(as.integer(SEED_BASE + 400000 + i))
  controls$sets <- list()
  z <- Score_cells(setNames(list(cpm), d$sample_key), cell_type = "SPG", MinGenes = MIN_GENES, MinScore = MIN_SCORE,
                   MinCells = MIN_CELL_FRACTION, MP_list = list(SPG = usable))
  fwrite(as.data.table(z$scores[[1]], keep.rownames = "cell_id"), file.path(score_dir, paste0(d$sample_key, "_scores.tsv.gz")), sep = "\t")
  saveRDS(controls$sets, file.path(score_dir, paste0(d$sample_key, "_matched_control_genes.rds")))

  out[, best_score := apply(as.matrix(z$scores[[1]]), 1, max)]
  prelim <- z$above_threshold[[1]]; retained <- z$retained[[1]]
  if (length(prelim))
    out[match(names(prelim), cell_id), `:=`(above_threshold_MP = as.character(prelim), reason = "MP_below_5pct_of_preliminary_assigned_cells")]
  retained <- retained[retained != "AMBIGUOUS_TIE"]
  if (length(retained))
    out[match(names(retained), cell_id), `:=`(assignment = as.character(retained), reason = "assigned")]
  out[above_threshold_MP == "AMBIGUOUS_TIE", `:=`(assignment = "Unassigned", reason = "equal_maximum_ambiguous")]

  if (d$cohort == "normal_discovery") {
    # MP_distribution.R lines 143-154: measured log2(CPM/10 + 1) of all MP members, centred among assigned cells,
    # cells ordered by program frequency.
    display_genes <- intersect(unique(unlist(MPs)), genes)
    expr <- log2(cpm[display_genes, , drop = FALSE] / 10 + 1)
    assigned <- out[assignment != "Unassigned"]
    freq <- assigned[, .N, by = assignment][order(-N, assignment)]
    assigned[, program_order := match(assignment, freq$assignment)]
    setorder(assigned, program_order, cell_id)
    centered <- expr[, assigned$cell_id, drop = FALSE]
    centered <- centered - rowMeans(centered)
    f <- file.path(expr_dir, paste0(d$sample_key, "_assigned_within_donor_centered.tsv.gz"))
    fwrite(as.data.table(centered, keep.rownames = "gene"), f, sep = "\t")
    assigned[, column_order := seq_len(.N)]
    expr_order[[i]] <- assigned
    expr_manifest[[i]] <- data.table(sample_key = d$sample_key, study = d$study, donor = d$donor, n_genes = nrow(expr),
                                     n_cells = ncol(expr), n_assigned_cells = nrow(assigned), assigned_centered_expression_file = f)
  }
  assignments[[i]] <- out
  counts_all[[i]] <- out[, .N, by = .(assignment, reason)][, `:=`(sample_key = d$sample_key, study = d$study, donor = d$donor,
    cohort = d$cohort, total_input_cells = nrow(out), fraction_all_input_cells = N / nrow(out))]
}
fwrite(rbindlist(assignments), file.path(prog_dir, "MP_cell_assignments.tsv.gz"), sep = "\t")
fwrite(rbindlist(counts_all), file.path(prog_dir, "MP_assignment_counts.tsv"), sep = "\t")
fwrite(rbindlist(coverage), file.path(prog_dir, "MP_scoring_feature_coverage.tsv"), sep = "\t")
fwrite(rbindlist(expr_manifest), file.path(prog_dir, "MP_measured_expression_manifest.tsv"), sep = "\t")
fwrite(rbindlist(expr_order), file.path(prog_dir, "MP_expression_cell_order.tsv.gz"), sep = "\t")
message(nrow(manifest), " samples scored; ", sum(rbindlist(assignments)$assignment != "Unassigned"), " cells assigned")
