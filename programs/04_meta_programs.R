# Robust NMF programs and meta-programs with the published code of Gavish et al. (2023), whose robustness filter is that
# of Kinker et al. (2020): top 50 genes per program; intra_min 35, intra_max 10, inter_min 10; MP seeding and extension at
# intersection 10 with seeds needing more than five neighbours. An MP is qualified with >= 3 donors and >= 2 studies.
# Run twice: with the NMF gene names (fixed MP1-MP3) and with five HGNC aliases harmonized across studies (program
# similarity used for the Fig 2a distances).
# Input:  results/programs/nmf/{donor_manifest.tsv, all_start_metrics.tsv, donors/<donor>/all_ranks_basis_usage.rds}
# Output: results/programs/nmf/{candidate_programs, candidate_top50_genes, robust_program_similarity, MP_program_membership,
#         MP_genes, MP_qualification, robust_program_plot_order}.tsv, results/programs/meta_alias/canonical/robust_program_similarity.tsv
source("utils/paths.R")
suppressPackageStartupMessages(library(data.table))

TOP_GENES <- 50
MIN_DONORS <- 3
MIN_STUDIES <- 2
ALIASES <- c("FAM26F" = "CALHM6", "HIST1H1D" = "H1-3", "FAM207A" = "SLX9", "C4ORF3" = "ARLN", "RP11-325I22.2" = "LINC01511")

nmf_dir <- file.path(RESULTS, "programs", "nmf")
manifest <- fread(file.path(nmf_dir, "donor_manifest.tsv"))[eligible == TRUE]
W <- lapply(manifest$sample_key, function(s) readRDS(file.path(nmf_dir, "donors", s, "all_ranks_basis_usage.rds"))$w_basis)
names(W) <- paste0(manifest$sample_key, "_rank6_9_nruns10.RDS")

# Vendor code. Guard patches keep matrix shape when a subset has one column and stop the loops on exhausted inputs.
gavish <- file.path(VENDOR, "3ca", "ITH_hallmarks", "Generating_MPs")
patch <- function(code, from, to) { stopifnot(grepl(from, code, fixed = TRUE)); gsub(from, to, code, fixed = TRUE) }
robust_code <- paste(readLines(file.path(gavish, "robust_nmf_programs.R")), collapse = "\n") |>
  patch("nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]",
        "nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min,drop=FALSE]") |>
  patch("inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]",
        "inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect)),drop=FALSE]") |>
  patch("nmf_sel_unlist <- do.call(cbind, nmf_sel)",
        "nmf_sel_unlist <- do.call(cbind, nmf_sel)\n  if(is.null(dim(nmf_sel_unlist)) || ncol(nmf_sel_unlist)<2) return(character())") |>
  patch("b <- sort(apply(a, 2, max), decreasing = T)",
        "if(nrow(a)==0 || ncol(a)==0) next\n    b <- sort(apply(a, 2, max), decreasing = T)")
mp_code <- readLines(file.path(gavish, "Generate_Meta_Programs.R"))
cluster_code <- paste(mp_code[64:134], collapse = "\n") |>
  patch("while (Sorted_intersection[1]>Min_group_size)",
        "while (length(Sorted_intersection)>0 && is.finite(Sorted_intersection[1]) && Sorted_intersection[1]>Min_group_size)") |>
  patch("while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster)",
        "while (length(Intersection_with_Genes_MP)>0 && is.finite(Intersection_with_Genes_MP[1]) && Intersection_with_Genes_MP[1] >= Min_intersect_cluster)") |>
  patch("nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]",
        "nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs)),drop=FALSE]") |>
  patch("nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]",
        "nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs)),drop=FALSE]") |>
  patch("nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]",
        "nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ),drop=FALSE ]")

run_gavish <- function(W) {
  env <- new.env(parent = globalenv())
  env$Genes_nmf_w_basis <- W
  eval(parse(text = robust_code), envir = env)
  eval(parse(text = mp_code[34:47]), envir = env)  # thresholds, top-50 genes, robust and non-redundant programs
  selected <- colnames(env$nmf_programs)
  eval(parse(text = mp_code[49:55]), envir = env)  # top-50 overlaps ordered by average-linkage clustering
  capture.output(eval(parse(text = cluster_code), envir = env))  # iterative MP clustering
  list(selected = selected, similarity = env$nmf_intersect_original, MPs = env$MP_list, clusters = env$Cluster_list)
}
similarity_table <- function(sim) {
  long <- as.data.table(as.table(sim))
  setnames(long, c("program1", "program2", "overlap_count"))
  long[, jaccard_percent := 100 * overlap_count / (100 - overlap_count)]
}

res <- run_gavish(W)

fits <- fread(file.path(nmf_dir, "all_start_metrics.tsv"))[selected_best == TRUE,
  .(sample_key, K, selected_fit_id = sprintf("%s::K%d::start%02d", sample_key, K, start), selected_start = start,
    selected_fit_niter = niter, selected_fit_hit_iteration_limit = hit_iteration_limit)]
programs <- rbindlist(lapply(seq_along(W), function(i) data.table(program = colnames(W[[i]]),
  sample_key = manifest$sample_key[i], study = manifest$study[i], donor = manifest$donor[i],
  K = as.integer(sub(".*RDS\\.([0-9]+)\\..*", "\\1", colnames(W[[i]]))))))
programs <- cbind(programs, fits[match(paste(programs$sample_key, programs$K), paste(sample_key, K)), !c("sample_key", "K")])
programs[, author_robust_selected := program %in% res$selected]
top50 <- rbindlist(lapply(W, function(w) rbindlist(lapply(colnames(w), function(p) {
  v <- w[, p]; g <- names(sort(v, decreasing = TRUE))[1:TOP_GENES]
  data.table(program = p, gene = g, rank = seq_along(g), nmf_basis_score = as.numeric(v[g]))
}))))

membership <- rbindlist(lapply(seq_along(res$MPs), function(i)
  cbind(data.table(MP = names(res$MPs)[i]), programs[match(res$clusters[[i]], program)])))
mp_genes <- rbindlist(lapply(names(res$MPs), function(mp) data.table(MP = mp, gene = res$MPs[[mp]], rank = seq_along(res$MPs[[mp]]))))
qualification <- membership[, .(n_programs = .N, n_donors = uniqueN(sample_key), n_studies = uniqueN(study),
  donors = paste(sort(unique(sample_key)), collapse = ";"), studies = paste(sort(unique(study)), collapse = ";"),
  n_member_programs_from_cap_fits = sum(selected_fit_hit_iteration_limit),
  cap_fit_ids = paste(sort(unique(selected_fit_id[selected_fit_hit_iteration_limit])), collapse = ";")), by = MP]
qualification[, eligible_main := n_donors >= MIN_DONORS & n_studies >= MIN_STUDIES]
qualification[, qualification_reason := fifelse(eligible_main, "at_least_3_donors_and_2_studies",
                                                 "insufficient_distinct_donor_or_study_replication")]
setcolorder(qualification, c("MP", "n_programs", "n_donors", "n_studies", "donors", "studies", "eligible_main", "qualification_reason"))

clustered <- unlist(res$clusters, use.names = FALSE)
plot_order <- programs[match(c(clustered, setdiff(rownames(res$similarity), clustered)), program)]
plot_order[, `:=`(plot_order = seq_len(.N), MP = membership$MP[match(program, membership$program)])]
plot_order[, MP_eligible_main := qualification$eligible_main[match(MP, qualification$MP)]]

out <- function(x, name) fwrite(x, file.path(nmf_dir, name), sep = "\t")
out(programs, "candidate_programs.tsv"); out(top50, "candidate_top50_genes.tsv")
out(similarity_table(res$similarity), "robust_program_similarity.tsv")
out(membership, "MP_program_membership.tsv"); out(mp_genes, "MP_genes.tsv")
out(qualification, "MP_qualification.tsv"); out(plot_order, "robust_program_plot_order.tsv")

# Same procedure after harmonizing the five aliases in every donor (FAM26F/CALHM6 and RP11-325I22.2/LINC01511 differ
# between studies). The same robust programs are selected; only their pairwise overlaps can change.
harmonize <- function(w) {
  hit <- toupper(rownames(w)) %in% names(ALIASES)
  rownames(w)[hit] <- ALIASES[toupper(rownames(w))[hit]]
  stopifnot(!anyDuplicated(toupper(rownames(w))))
  w
}
canonical <- run_gavish(lapply(W, harmonize))
stopifnot(setequal(canonical$selected, res$selected))
fwrite(similarity_table(canonical$similarity), file.path(out_dir("programs", "meta_alias", "canonical"), "robust_program_similarity.tsv"), sep = "\t")
message(length(res$selected), " robust programs; ", paste(names(res$MPs), lengths(res$clusters), sep = ": ", collapse = ", "),
        " contributing programs; qualified: ", paste(qualification[eligible_main == TRUE]$MP, collapse = ", "))
