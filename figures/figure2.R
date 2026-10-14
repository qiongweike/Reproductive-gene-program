# Figure 2: recurrent programs and coordination beyond established state markers.
# Inputs: alias-harmonized program similarity (results/programs/meta_alias), program plot order (results/programs/nmf),
# program genes, cell order and measured expression (results/programs), coordination,
# marker-removal and enrichment tables (results/coordination) and Hotspot summaries
# (results/trajectory/hotspot). Output: results/figures/figure2.pdf and figure2.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({library(data.table); library(scales)})

PROGRAMS <- file.path(RESULTS, "programs")
META <- file.path(PROGRAMS, "meta_alias/canonical")
HEATMAP_SEED <- 20261001
HEATMAP_GENES <- 10          # top-ranked genes shown per program
HEATMAP_CELLS <- 150         # cells shown per donor at most
GO_Q <- .05                  # MP3 shows every GO term below this global BH q
REDUNDANT_JACCARD <- .6      # MP1/MP2 terms overlapping an earlier term more than this are skipped
TERMS_PER_LIBRARY <- 3       # MP1/MP2 terms shown per library

MP_COLOURS <- c(MP_1 = "#368A88", MP_2 = "#7763A6", MP_3 = "#C47745")
MP_LABELS <- c(MP_1 = "MP1", MP_2 = "MP2", MP_3 = "MP3")
STUDY_COLOURS <- c(GSE124263 = "#4477AA", GSE153947 = "#228833", GSE254315 = "#AA3377")
DONOR_COLOURS <- c(A1 = "#4477AA", A2 = "#66AADD", N1 = "#228833", N2 = "#44AA99", N3 = "#999933",
                   Donor30 = "#AA3377", Donor5 = "#CC6677")
STATES <- c(undifferentiated = "U", transition = "T", differentiating = "D")
no_axis_lines <- theme(axis.line = element_blank(), axis.ticks = element_blank())

## a  average-linkage tree of the 47 robust programs (distance = 50 - top-50 overlap)
programs <- fread(file.path(PROGRAMS, "nmf/robust_program_plot_order.tsv"))
programs[, display := paste0(fifelse(donor == "Donor30", "D30", fifelse(donor == "Donor5", "D5", donor)), ".K", K, ".",
                             sub(".*\\.", "", program), fifelse(selected_fit_hit_iteration_limit, "*", ""))]
sim <- fread(file.path(META, "robust_program_similarity.tsv"))
overlap <- matrix(NA_real_, nrow(programs), nrow(programs), dimnames = list(programs$program, programs$program))
overlap[cbind(match(sim$program1, programs$program), match(sim$program2, programs$program))] <- sim$overlap_count
stopifnot(all(is.finite(overlap)), all(diag(overlap) == 50))
tree <- hclust(as.dist(50 - overlap), method = "average")

leaf_x <- setNames(seq_along(tree$order), tree$order)
node_x <- node_y <- numeric(nrow(tree$merge)); edges <- list()
for (i in seq_len(nrow(tree$merge))) {
  ch <- tree$merge[i, ]
  cx <- vapply(ch, function(j) if (j < 0) unname(leaf_x[as.character(-j)]) else node_x[j], numeric(1))
  cy <- vapply(ch, function(j) if (j < 0) 0 else node_y[j], numeric(1))
  node_x[i] <- mean(cx); node_y[i] <- tree$height[i]
  edges[[i]] <- data.table(x = c(cx, cx[1]), xend = c(cx, cx[2]), y = c(cy, node_y[i]), yend = rep(node_y[i], 3))
}
edges <- rbindlist(edges)
leaves <- programs[tree$order, .(x = seq_len(.N), display, program = fifelse(is.na(MP) | !nzchar(MP), "Unassigned", MP))]
pa <- ggplot() +
  geom_segment(data = edges, aes(x = x, xend = xend, y = y, yend = yend), linewidth = .28, colour = INK) +
  geom_point(data = leaves, aes(x = x, y = 0, colour = program), size = 1.3) +
  scale_x_continuous(breaks = leaves$x, labels = leaves$display, expand = expansion(add = .6)) +
  scale_y_continuous(breaks = c(0, 25, 50), limits = c(0, 50), expand = expansion(mult = c(.02, .04))) +
  scale_colour_manual(values = c(MP_COLOURS, Unassigned = "#8A9298"), labels = c(MP_LABELS, Unassigned = "Unassigned"),
                      name = "Program") +
  labs(x = NULL, y = "50 − overlap") + theme_spg() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6.5),
        axis.ticks.x = element_blank(), axis.line.x = element_blank())

## b  measured expression of the top genes of each program in at most 150 assigned cells per donor
genes <- fread(file.path(PROGRAMS, "MP_genes.tsv"))[MP %in% names(MP_LABELS) & rank <= HEATMAP_GENES]
gene_order <- unique(genes$gene)
expression_files <- fread(file.path(PROGRAMS, "MP_measured_expression_manifest.tsv"))
cells <- fread(file.path(PROGRAMS, "MP_expression_cell_order.tsv.gz"))[
  cohort == "normal_discovery" & assignment %in% names(MP_LABELS)]
cells <- cells[sample_key %in% expression_files$sample_key]
set.seed(HEATMAP_SEED)
shown <- cells[, .SD[sort(sample(seq_len(.N), min(.N, HEATMAP_CELLS)))], by = sample_key]
setorder(shown, sample_key, program_order, column_order)
stopifnot(!anyDuplicated(shown$cell_id))
shown[, display_column := seq_len(.N)]
expr <- matrix(NA_real_, length(gene_order), nrow(shown), dimnames = list(gene_order, shown$cell_id))
for (sk in unique(shown$sample_key)) {
  f <- basename(expression_files[sample_key == sk, assigned_centered_expression_file])
  x <- fread(file.path(PROGRAMS, "measured_expression", f))
  take <- shown[sample_key == sk, cell_id]
  v <- as.matrix(x[, ..take]); rownames(v) <- x[[1]]
  g <- intersect(gene_order, rownames(v)); expr[g, take] <- v[g, take, drop = FALSE]
}
heat <- setnames(as.data.table(as.table(expr)), c("gene", "cell_id", "value"))
heat <- merge(heat, shown[, .(cell_id, display_column, sample_key)], by = "cell_id")
heat <- merge(heat, genes[, .(gene, MP)], by = "gene", allow.cartesian = TRUE)
heat[, Gene := factor(gene, levels = rev(gene_order))]
donor_blocks <- shown[, .(mid = mean(display_column), end = max(display_column)), by = sample_key]
gene_labels <- setNames(gene_order, gene_order)
first_genes <- genes[order(rank), .(gene = gene[1]), by = MP]
gene_labels[first_genes$gene] <- paste0(MP_LABELS[first_genes$MP], ": ", first_genes$gene)
row_dividers <- length(gene_order) - head(cumsum(genes[, .N, by = MP]$N), -1) + .5
pb <- ggplot(heat, aes(display_column, Gene, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = "#426F9C", mid = "white", high = "#BC594C", limits = c(-4, 4), oob = squish,
                       na.value = "#B5B5B5", breaks = c(-4, 0, 4), labels = c("<= -4", "0", ">= 4"),
                       name = "Within-donor centered\nlog2(CP100K + 1)") +
  geom_vline(xintercept = head(donor_blocks$end, -1) + .5, colour = "white", linewidth = .4) +
  geom_hline(yintercept = row_dividers, colour = "white", linewidth = .6) +
  scale_y_discrete(labels = gene_labels, expand = c(0, 0)) +
  scale_x_continuous(breaks = donor_blocks$mid, labels = sub(".*__", "", donor_blocks$sample_key), expand = c(0, 0)) +
  labs(x = "Normal / OA donor", y = NULL) + theme_spg() + no_axis_lines +
  theme(axis.text.y = element_text(size = 6.5)) + guides(fill = colourbar_guide())

## c  within-state coordination against matched gene sets
coh <- fread(file.path(RESULTS, "coordination/coherence/coherence_results.tsv"))[
  model == "technical_cycle" & level == "primary" & module %in% names(MP_LABELS) & is.finite(mean_pairwise_rho)]
coh <- rbind(coh[, .(module, rho = mean_pairwise_rho, distribution = "Observed")],
             coh[, .(module, rho = null_mean, distribution = "Matched null")])
pc <- ggplot(coh, aes(rho, colour = module, linetype = distribution)) + stat_ecdf(linewidth = .5) +
  scale_colour_manual(values = MP_COLOURS, labels = MP_LABELS, name = "Program") +
  labs(x = "Mean residual ρ", y = "Cumulative fraction", linetype = "Reference") + theme_spg()

## d  donor-mean Hotspot local pair C, whole donor versus within state
hs <- fread(file.path(RESULTS, "trajectory/hotspot/donor_residual_comparison.tsv"))
hs[, `:=`(donor = sub(".*__", "", sample_key),
          scope = factor(scope, levels = c("whole_donor", "within_state"), labels = c("All", "State")),
          model = factor(model, levels = c("bernoulli", "danb"), labels = c("Detection", "Counts")))]
pd <- ggplot(hs, aes(scope, mean_Hotspot_C_all_available, group = sample_key, colour = donor)) +
  geom_hline(yintercept = 0, colour = "#BCC5CA", linewidth = .25) +
  geom_line(linewidth = .4, na.rm = TRUE) + geom_point(size = 1.4, na.rm = TRUE) +
  facet_wrap(~model, nrow = 1) + scale_colour_manual(values = DONOR_COLOURS, name = "Donor") +
  labs(x = NULL, y = "Mean local pair C") + theme_spg()

## e  fraction of measurable donors with positive within-state local C for all MP3 member pairs
members <- fread(file.path(PROGRAMS, "MP_genes.tsv"))[MP == "MP_3"][order(rank), gene]
pairs <- fread(file.path(RESULTS, "trajectory/hotspot/MP3_pair_consensus.tsv"))[
  scope == "within_state" & model == "bernoulli"]
pairs[, fraction_positive := fifelse(covered_ge4_donors_ge2_studies, n_positive_donors / n_donors, NA_real_)]
pairs <- rbind(pairs[, .(gene1, gene2, fraction_positive)], pairs[, .(gene1 = gene2, gene2 = gene1, fraction_positive)])
pairs <- merge(CJ(gene1 = members, gene2 = members), pairs, by = c("gene1", "gene2"), all.x = TRUE)
stopifnot(nrow(pairs) == length(members)^2)
pairs[, `:=`(gene1 = factor(gene1, levels = members), gene2 = factor(gene2, levels = rev(members)))]
ticks <- seq(1, length(members), 8)
pe <- ggplot(pairs, aes(gene1, gene2, fill = fraction_positive)) + geom_tile() +
  scale_fill_gradient(low = "#F4F2E8", high = "#245B78", limits = c(0, 1), breaks = c(0, .5, 1), na.value = "#D9D9D9",
                      name = "Donor fraction\nC > 0") +
  scale_x_discrete(breaks = members[ticks], labels = ticks) + scale_y_discrete(breaks = members[ticks], labels = ticks) +
  labs(x = "Member index", y = "Member index") + theme_spg() + no_axis_lines + theme(aspect.ratio = 1) +
  guides(fill = colourbar_guide())

## f  MP3 residual coordination before and after removing published state markers
mr <- dcast(fread(file.path(RESULTS, "coordination/marker_removal/known_marker_sensitivity.tsv")),
            sample_key + study + state ~ sensitivity, value.var = "mean_rho")
stopifnot(nrow(mr) == 14)
lim <- range(c(0, mr$fixed_MP3, mr$remove_published_state_markers)); lim[2] <- ceiling(lim[2] * 100) / 100
pf <- ggplot(mr, aes(fixed_MP3, remove_published_state_markers, colour = study, shape = state)) +
  geom_abline(slope = 1, intercept = 0, colour = "#A8B1B9", linewidth = .35) + geom_point(size = 1.9, stroke = .35) +
  scale_colour_manual(values = STUDY_COLOURS, name = "Study") +
  scale_shape_manual(values = c(undifferentiated = 16, transition = 17, differentiating = 15), labels = STATES,
                     name = "State") +
  coord_fixed(xlim = lim, ylim = lim) + labs(x = "Complete MP3: residual ρ", y = "State genes removed: ρ") + theme_spg()

## g  GO/Reactome over-representation: all significant MP3 GO terms; for MP1 and MP2 the
##    first three non-redundant terms per library
ora <- fread(file.path(RESULTS, "coordination/enrichment/complete_GO_Reactome.tsv.gz"))[gene_set %in% names(MP_LABELS)]
ora[, library := fifelse(grepl("^GO", database), "GO", "Reactome")]
setorder(ora, gene_set, library, q_BH_all, -n_overlap, pathway)
non_redundant <- function(genes) {
  used <- list(); take <- integer()
  for (i in seq_along(genes)) {
    g <- strsplit(genes[i], ";", fixed = TRUE)[[1]]
    jaccard <- vapply(used, function(u) length(intersect(u, g)) / max(1, length(union(u, g))), numeric(1))
    if (any(jaccard > REDUNDANT_JACCARD)) next
    take <- c(take, i); used[[length(used) + 1]] <- g
    if (length(take) == TERMS_PER_LIBRARY) break
  }
  take
}
mp3_terms <- ora[gene_set == "MP_3" & library == "GO" & q_BH_all < GO_Q]
stopifnot(nrow(mp3_terms) == 7)
terms <- rbind(ora[gene_set != "MP_3", .SD[non_redundant(genes)], by = .(gene_set, library)], mp3_terms,
               use.names = TRUE)
setorder(terms, gene_set, library, q_BH_all, -n_overlap, pathway)
terms[, `:=`(score = -log10(pmax(q_BH_all, 1e-300)),
             key = factor(paste(gene_set, pathway), levels = rev(paste(gene_set, pathway))),
             label = vapply(name, function(x) paste(strwrap(x, width = 55), collapse = "\n"), character(1)))]
pg <- ggplot(terms, aes(score, key)) +
  geom_segment(aes(x = 0, xend = score, yend = key), colour = "#BECAD0", linewidth = .35) +
  geom_point(aes(size = n_overlap, colour = library)) +
  facet_grid(gene_set ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(MP_LABELS)) +
  scale_y_discrete(labels = setNames(terms$label, terms$key)) +
  scale_colour_manual(values = c(GO = "#4477AA", Reactome = "#AA6677"), name = "Library") +
  scale_size_area(max_size = 2.7, name = "Overlap") +
  labs(x = "−log10(global BH q)", y = NULL) + theme_spg() +
  theme(axis.text.y = element_text(size = 6.5), strip.text.y = element_text(angle = 0))

compose_figure(list(pa, pb, pc, pd, pe, pf, pg), rows = list(1, 2, 3:4, 5:6, 7), heights = c(37, 69, 57, 65, 76),
               widths = list(1, 1, c(1, 1.05), c(1, 1.05), 1), file = file.path(FIGURES, "figure2"))
cat("Figure 2 written:", nrow(programs), "programs,", nrow(shown), "heatmap cells,", nrow(terms), "enrichment terms\n")
