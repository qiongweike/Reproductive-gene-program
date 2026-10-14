# Figure 5: candidate ligand inputs and regulatory associations of MP3 (panels a-i).
# Inputs (results/communication/): measured signalling expression, NicheNet ligand
# activities with the NicheNet prior matrix, per-donor CellChat models, fixed regulon
# motif edges, non-MP activity correlations, edge recovery, regulon peak offsets and
# AUCell summaries. Output: results/figures/figure5.pdf and figure5.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({
  library(data.table); library(scales); library(ggrepel); library(CellChat)
})

COMM     <- file.path(RESULTS, "communication")
NICHENET <- file.path(COMM, "nichenet")
REGULON  <- file.path(COMM, "regulon")
CELLCHAT <- file.path(COMM, "cellchat")
PRIOR    <- file.path(EXTERNAL, "nichenet", "ligand_target_matrix_nsga2r_final.rds")

N_LIGANDS   <- 10                                   # top U-context ligands shown in b and c
N_MP3       <- 49
CASE_DONORS <- c("OA1", "OA2", "NOA1", "NOA2")      # GSE202647; d shows OA1 on limits shared by all four
SENSITIVITY <- c(fixed_MP3 = "Fixed", remove_ID1_ID2_ID4 = "− IDs", remove_published_state_markers = "− State")
STATE_CODE  <- c(undifferentiated = "U", transition = "T", differentiating = "D")
REPEL_SEED  <- 20261004
# g: column order of the published panel; donor/seed models, then the three leave-one-study-out models.
RECOVERY_COLUMNS <- c("N1 / 03", "N1 / 17", "N2 / 17", "N2 / 03", "N3 / 17", "N3 / 03", "Donor30 / 03", "Donor30 / 17",
                      "Donor5 / 03", "Donor5 / 17", "A2 / 03", "A2 / 17", "A1 / 17", "A1 / 03",
                      "−124263 / 03", "−153947 / 03", "−254315 / 03")

DONOR_COLOURS  <- c(A1 = "#4477AA", A2 = "#66AADD", N1 = "#228833", N2 = "#44AA99", N3 = "#999933",
                    Donor30 = "#AA3377", Donor5 = "#CC6677")
LIGAND_COLOURS <- c("#4477AA", "#CC6677", "#228833", "#AA3377", "#66AABB", "#999933", "#88CCEE",
                    "#DDCC77", "#882255", "#44AA99")   # assigned to the ligands in alphabetical order
GROUP_COLOURS  <- c(Spermatocyte = "#69A9AF", Spermatid = "#89B67E", Sertoli = "#CE9A45", Leydig = "#B45E65",
                    Peritubular_smooth_muscle = "#8D7AAB", Endothelial = "#6B84B7", Immune = "#8C785D",
                    Erythroid = "#CE8982", SPG_U = "#48799B", SPG_T = "#D3A249", SPG_D = "#459386")
GROUP_LABELS   <- c(Spermatocyte = "SPC", Spermatid = "SPT", Sertoli = "Sertoli", Leydig = "Leydig",
                    Peritubular_smooth_muscle = "PTM", Endothelial = "EC", Immune = "Imm", Erythroid = "Ery",
                    SPG_U = "SPG-U", SPG_T = "SPG-T", SPG_D = "SPG-D")

tile_theme <- theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6.5),
                    axis.text.y = element_text(size = 6.5))

# netAnalysis_signalingRole_scatter from the pinned CellChat source, evaluated in the package namespace.
cellchat_role_scatter <- function() {
  code <- as.list(parse(file.path(VENDOR, "CellChat", "R", "analysis.R")))
  def <- Filter(function(x) is.call(x) && as.character(x[[1]]) %in% c("<-", "=") &&
                  identical(x[[2]], as.name("netAnalysis_signalingRole_scatter")), code)
  stopifnot(length(def) == 1)
  env <- new.env(parent = asNamespace("CellChat"))
  eval(def[[1]], env)
  env$netAnalysis_signalingRole_scatter
}

tfs <- fread(file.path(REGULON, "structural_TF_selection.tsv"))[selected == TRUE, TF]
stopifnot(length(tfs) == 12L)
P <- list()

# a: measured expression of the signalling components; U donor rows last, in donor order.
expr <- fread(file.path(COMM, "signalling_expression.tsv"))
expr[, row := gsub("Peritubular_smooth_muscle", "PTM", row, fixed = TRUE)]
u_rows <- unique(expr[group == "SPG_U", row])
u_rows <- u_rows[order(match(sub(" / .*", "", u_rows), CASE_DONORS))]
expr[, gene := factor(gene, levels = fread(file.path(COMM, "signalling_gene_identity.tsv"))$display_gene)]
expr[, row := factor(row, levels = rev(c(setdiff(unique(row), u_rows), u_rows)))]
P$a <- ggplot(expr, aes(gene, row)) +
  geom_point(aes(size = detection, colour = mean_logCP10k), stroke = 0) +
  scale_size_area(max_size = 2.05, limits = c(0, 1), breaks = c(.1, .5, 1), labels = percent_format(),
                  name = "Detected\ncells") +
  scale_colour_gradient(low = "#DEE8EC", high = "#27566E", name = "Mean log1p\nCP10k") +
  labs(x = NULL, y = NULL) + theme_spg() + tile_theme +
  guides(colour = colourbar_guide(), size = guide_legend(order = 2))

# b: NicheNet prior regulatory potential, displayed ligands x all MP3 members.
ligands <- head(fread(file.path(NICHENET, "U_ligand_selection.tsv"))$test_ligand, N_LIGANDS)
members <- fread(file.path(NICHENET, "fixed_program_entity_mapping.tsv"))[MP == "MP_3"][order(rank)]
members <- members[!duplicated(ensembl)]
stopifnot(nrow(members) == N_MP3, !anyDuplicated(members$prior_target))
prior <- t(readRDS(PRIOR)[members$prior_target, rev(ligands), drop = FALSE])
colnames(prior) <- members$approved_symbol
potential <- data.table(ligand = factor(rownames(prior)[row(prior)], levels = rownames(prior)),
                        target = factor(colnames(prior)[col(prior)], levels = colnames(prior)),
                        potential = as.vector(prior))
P$b <- ggplot(potential, aes(target, ligand, fill = potential)) +
  geom_tile(colour = "white", linewidth = .5) +
  scale_fill_gradient(low = "#F4F2F2", high = "#865170", name = "Prior\npotential") +
  labs(x = NULL, y = NULL) + theme_spg() + tile_theme + guides(fill = colourbar_guide())

# c: median ligand rank across eligible U contexts, fixed targets and two removal sensitivities.
act <- fread(file.path(NICHENET, "ligand_activities.tsv"))
act <- act[state == "SPG_U" & MP == "MP_3" & test_ligand %in% ligands]
act[, sensitivity := "fixed_MP3"]
act <- rbind(act, fread(file.path(NICHENET, "ligand_sensitivity.tsv"))[test_ligand %in% ligands], fill = TRUE)
act[, sensitivity := factor(sensitivity, levels = names(SENSITIVITY))]
ranks <- act[, .(median_rank = median(rank), n_contexts = .N), by = .(test_ligand, sensitivity)]
fixed <- ranks[sensitivity == "fixed_MP3"]
P$c <- ggplot(ranks, aes(sensitivity, median_rank, group = test_ligand, colour = test_ligand)) +
  geom_line(linewidth = .35) + geom_point(size = 1.25) +
  scale_y_reverse() + scale_x_discrete(labels = SENSITIVITY) +
  scale_colour_manual(values = setNames(LIGAND_COLOURS, sort(ligands)),
                      labels = setNames(paste0(fixed$test_ligand, " (", fixed$n_contexts, ")"), fixed$test_ligand),
                      name = "Ligand") +
  labs(x = NULL, y = "Median rank") + theme_spg() +
  guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = .35, size = 1.5)))

# d: OA1 CellChat outgoing and incoming strengths; axis and size limits shared by the four donors.
role_scatter <- cellchat_role_scatter()
roles <- lapply(setNames(CASE_DONORS, CASE_DONORS), function(d) {
  cc <- readRDS(file.path(CELLCHAT, paste0("GSE202647_", d), "cellchat.rds"))
  suppressMessages(role_scatter(cc, color.use = GROUP_COLOURS[levels(cc@idents)], do.label = FALSE))$data
})
all_roles <- rbindlist(roles)
P$d <- ggplot(roles$OA1, aes(x, y)) +
  geom_point(aes(size = Count, colour = labels)) +
  scale_size_continuous(range = c(1.3, 3.4), limits = range(all_roles$Count), breaks = c(250, 750, 1250),
                        name = "Incident\nLR pairs") +
  scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = "Cell group") +
  coord_cartesian(xlim = c(0, max(all_roles$x) * 1.07), ylim = c(0, max(all_roles$y) * 1.07)) +
  labs(x = "Outgoing strength", y = "Incoming strength") + theme_spg() +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.1)), size = guide_legend(order = 2))

# e: motif-supported edges from the twelve fixed factors to MP3 targets (self-edges excluded).
edges <- fread(file.path(REGULON, "MP3_motif_edges.tsv"))[TF %in% tfs & TF != target]
edges <- edges[, .(n_motifs = max(n_motifs), max_motif_NES = max(max_motif_NES)), by = .(TF, target)]
stopifnot(nrow(edges) == 65L, uniqueN(edges$target) == 29L)
targets <- sort(unique(edges$target))
adjacency <- merge(CJ(TF = tfs, target = targets), edges, by = c("TF", "target"), all.x = TRUE)
adjacency[, `:=`(TF = factor(TF, levels = rev(tfs)), target = factor(target, levels = targets))]
P$e <- ggplot(adjacency, aes(target, TF)) +
  geom_tile(fill = "#EEF1F2", colour = "white", linewidth = .3) +
  geom_point(data = adjacency[is.finite(n_motifs)], aes(size = n_motifs, colour = max_motif_NES), shape = 15) +
  scale_colour_gradient(low = "#92B5B2", high = "#215C69", name = "Motif\nNES") +
  scale_size_area(max_size = 2.25, breaks = c(1, 3, 5), name = "Motifs") +
  labs(x = NULL, y = NULL) + theme_spg() + tile_theme +
  guides(colour = colourbar_guide(), size = guide_legend(order = 2))

# f: donor-mean partial correlation with MP3, MP targets removed vs MP + 143-gene set removed.
rho <- fread(file.path(REGULON, "nonMP_vs_extended_associations.tsv"))[program == "MP_3"]
rho <- rho[is.finite(mean_rho_nonMP) & is.finite(mean_rho_extended)]
rho[, label := ifelse(structurally_selected_TF, TF, "")]
P$f <- ggplot(rho, aes(mean_rho_nonMP, mean_rho_extended)) +
  geom_abline(slope = 1, intercept = 0, colour = "#BAC3C6", linewidth = .35) +
  geom_point(colour = "#608C99", size = .9, alpha = .7) +
  geom_text_repel(aes(label = label), size = 6.5 / .pt, seed = REPEL_SEED, max.overlaps = Inf,
                  min.segment.length = 0, box.padding = .12) +
  labs(x = "ρ: MP removed", y = "ρ: MP + 143 removed") + theme_spg()

# g: recovery of eligible fixed edges in donor-seed and leave-study-out regulon models.
recovery <- fread(file.path(COMM, "regulon_edges", "MP3_edges_by_model.tsv"))[TF %in% tfs]
recovery <- recovery[, .(fraction = if (any(eligible)) sum(edge_selected & eligible) / sum(eligible) else NA_real_),
                     by = .(TF, block, scope, key, seed)]
recovery[, context := paste0(ifelse(scope == "donor", sub(".*__", "", key), paste0("−", sub("GSE", "", key))),
                             " / ", sub("202610", "", as.character(seed)))]
stopifnot(setequal(recovery$context, RECOVERY_COLUMNS))
recovery[, `:=`(TF = factor(TF, levels = rev(tfs)), context = factor(context, levels = RECOVERY_COLUMNS))]
setorder(recovery, context)   # tiles drawn column by column, in the published order
P$g <- ggplot(recovery, aes(context, TF, fill = fraction)) +
  geom_tile(colour = "white", linewidth = .25) +
  scale_fill_gradient(low = "#EDF2EE", high = "#3D776B", limits = c(0, 1), na.value = "#B9BCBE",
                      labels = percent_format(), name = "Recovered\nedges") +
  labs(x = NULL, y = NULL) + theme_spg() + tile_theme + guides(fill = colourbar_guide())

# h: regulon peak minus MP3 score peak in the selected trajectories; bars span tied maximum bins.
peaks <- fread(file.path(COMM, "regulon_peak_offsets.tsv"))
peaks[, TF := factor(TF, levels = rev(tfs))]
P$h <- ggplot(peaks, aes(peak_offset, TF, colour = donor)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#B8C1C7", linewidth = .3) +
  geom_errorbar(aes(xmin = offset_min, xmax = offset_max), orientation = "y", width = .1, linewidth = .35,
                position = position_dodge(width = .5)) +
  geom_point(size = 1.2, position = position_dodge(width = .5)) +
  scale_colour_manual(values = DONOR_COLOURS, name = "Donor") +
  labs(x = "Regulon − MP3 peak (bins)", y = NULL) + theme_spg() + theme(axis.text.y = element_text(size = 6.5))

# i: non-MP AUCell and target detection by donor and state; crosses mark groups without cells.
auc <- fread(file.path(REGULON, "nonMP_AUCell_by_state.tsv"))[TF %in% tfs]
row_levels <- rev(unlist(lapply(names(DONOR_COLOURS), function(d) paste(d, STATE_CODE, sep = " / "))))
auc[, `:=`(TF = factor(TF, levels = tfs),
           row = factor(paste(sub(".*__", "", sample_key), STATE_CODE[state], sep = " / "), levels = row_levels))]
P$i <- ggplot(auc, aes(TF, row)) +
  geom_point(data = auc[!is.finite(mean_AUCell)], aes(shape = "No cells"), size = 1, colour = "#A9AFB4") +
  geom_point(aes(size = mean_detected_target_fraction, colour = mean_AUCell), na.rm = TRUE) +
  scale_shape_manual(values = c("No cells" = 4), name = NULL) +
  scale_size_area(max_size = 2, limits = c(0, 1), breaks = c(.1, .3, .6), labels = percent_format(),
                  name = "Detected\ntargets") +
  scale_colour_gradient(low = "#DCE6E8", high = "#2E657A", name = "Non-MP\nAUCell") +
  scale_y_discrete(limits = row_levels, drop = FALSE) +
  labs(x = NULL, y = NULL) + theme_spg() + tile_theme +
  guides(colour = colourbar_guide(), size = guide_legend(order = 2), shape = guide_legend(order = 3))

compose_figure(P, rows = list(1, 2, c(3, 4), c(5, 6), c(7, 8, 9)), heights = c(54, 49, 55, 60, 80),
               widths = list(1, 1, c(1, 1.15), c(1.55, 1), c(1.2, .9, 1.1)), file = file.path(FIGURES, "figure5"))
cat("figure5:", length(P), "panels ->", file.path(FIGURES, "figure5.{pdf,png}"), "\n")
