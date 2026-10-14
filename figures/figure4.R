# Figure 4: program measurements in residual spermatogonia from infertile men.
# Inputs: results/clinical/ (capture states, donor-state scores, cohort effects, measured MP3 members),
#   results/programs/MP_assignment_counts.tsv,
#   results/coordination/coherence/external_NOA/coherence_results.tsv and
#   results/communication/regulon/ (projection_U_AUCell.tsv, structural_TF_selection.tsv).
# Output: results/figures/figure4.pdf and figure4.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages(library(data.table))

CLINICAL <- file.path(RESULTS, "clinical")
DONORS <- c("OA1", "OA2", "NOA1", "NOA2")             # GSE202647 donors, display order
MIN_CELLS <- 30                                       # primary donor-state cell minimum
MP_LABELS <- c(MP_1 = "MP1", MP_2 = "MP2", MP_3 = "MP3")
MP_COLOURS <- c(MP_1 = "#368A88", MP_2 = "#7763A6", MP_3 = "#C47745")
STATE_COLOURS <- c(undifferentiated = "#48799B", transition = "#D3A249", differentiating = "#459386")
STATE_LABELS <- c(undifferentiated = "U", transition = "T", differentiating = "D")
DIAGNOSIS_COLOURS <- c(OA = "#48799B", NOA = "#B75F50", Crypto = "#A386B7", Normal = "#459386")
TILE <- theme(axis.line = element_blank(), axis.ticks = element_blank())
VERTICAL_X <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6.5))

short_donor <- function(x) sub("GSE202647_", "", x)
# Row labels such as "OA1 (137)" from the U-cell counts, OA1 on top.
donor_rows <- function(donor, n_cells) {
  label <- paste0(short_donor(donor), " (", n_cells, ")")
  factor(label, levels = rev(unique(label[order(match(short_donor(donor), DONORS))])))
}

# a: captured U/T/D composition per donor among state-labelled primary spermatogonia.
z <- fread(file.path(CLINICAL, "capture", "GSE202647_SPG_states.tsv"))[state %in% names(STATE_COLOURS)]
z[, fraction := n_cells / sum(n_cells), by = donor]
z[, state := factor(state, levels = names(STATE_COLOURS))]
p_a <- ggplot(z, aes(short_donor(donor), fraction, fill = state)) + geom_col(width = .65) +
  scale_fill_manual(values = STATE_COLOURS, labels = STATE_LABELS, name = "State") +
  scale_x_discrete(limits = DONORS) + scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Captured SPG fraction") + theme_spg() +   # a keeps 6.5-pt ticks and 2.5-mm legend keys
  theme(axis.text = element_text(size = 6.5), legend.key.height = unit(2.5, "mm"), legend.key.width = unit(2.5, "mm")) +
  guides(fill = guide_legend(ncol = 1))

# b: background-adjusted program scores in U cells (filled: primary 30-cell criterion met).
z <- fread(file.path(CLINICAL, "donor_state_scores.tsv"))[
  study == "GSE202647" & state == "undifferentiated" & scenario == "atlas_SPG_agreement"]
z[, eligible := n_cells >= MIN_CELLS]
z <- melt(z, id.vars = c("donor", "diagnosis", "eligible"), measure.vars = paste0(names(MP_LABELS), "_adjusted"),
          variable.name = "MP", value.name = "score")
z[, MP := sub("_adjusted", "", MP)]
p_b <- ggplot(z, aes(short_donor(donor), score, colour = diagnosis, shape = eligible)) + geom_point(size = 2) +
  facet_wrap(~MP, nrow = 1, labeller = as_labeller(MP_LABELS)) +
  scale_colour_manual(values = DIAGNOSIS_COLOURS, name = "Diagnosis") +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), labels = c("10–29", "≥30"), name = "U cells") +
  scale_x_discrete(limits = DONORS) + labs(x = NULL, y = "Score − background") + theme_spg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

# c: all 49 MP3 members measured in the same U-cell groups.
z <- fread(file.path(CLINICAL, "expression", "MP3_member_expression_U.tsv"))
setorder(z, member_order)
z[, gene := factor(canonical_gene, levels = unique(canonical_gene))]
z[, row := donor_rows(donor, n_cells)]
p_c <- ggplot(z, aes(gene, row)) + geom_point(aes(size = detection_fraction, colour = mean_log1p_CP10k), stroke = 0) +
  scale_size_area(max_size = 2.5, limits = c(0, 1), breaks = c(.1, .5, .9), labels = scales::percent_format(), name = "Detected") +
  scale_colour_gradient(low = "#DAE7E9", high = "#286C83", name = "Mean ln\n(1+CP10k)") +
  labs(x = NULL, y = NULL) + theme_spg() + TILE + VERTICAL_X + guides(colour = colourbar_guide())

# d: cryptozoospermia versus OA donor-mean differences with 95% Welch intervals (GSE153947).
z <- fread(file.path(CLINICAL, "all_effects.tsv"))[study == "GSE153947" & scenario == "primary" & min_cells == MIN_CELLS &
                                                   grepl("_adjusted$", metric) & state %in% c("undifferentiated", "differentiating")]
z[, MP := sub("_adjusted", "", metric)]
z <- z[order(match(state, names(STATE_LABELS)), MP)]
z[, label := paste0(MP_LABELS[MP], " / ", STATE_LABELS[state])]
z[, label := factor(label, levels = rev(label))]
p_d <- ggplot(z, aes(effect, label, colour = state)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, colour = "#999999") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = .12, linewidth = .4) + geom_point(size = 1.8) +
  scale_colour_manual(values = STATE_COLOURS, labels = STATE_LABELS, name = "State") +
  labs(x = "Crypto − OA (95% Welch CI)", y = NULL) + theme_spg()

# e: fixed program assignments in four NOA cases; unassigned cells stay in the denominator.
z <- fread(file.path(RESULTS, "programs", "MP_assignment_counts.tsv"))[cohort == "NOA_case_description"]
z <- z[, .(N = sum(N), total_input_cells = unique(total_input_cells)), by = .(sample_key, donor, assignment)]
stopifnot(z[, sum(N) == unique(total_input_cells), by = sample_key]$V1)
z[, fraction := N / total_input_cells]
z[, assignment := factor(assignment, levels = c(names(MP_COLOURS), "Unassigned"))]
z[, case := paste0("Case ", sub(".*testis_", "", donor), "\n(", format(total_input_cells, big.mark = ",", trim = TRUE), ")")]
p_e <- ggplot(z, aes(case, fraction, fill = assignment)) + geom_col(width = .65) +
  scale_fill_manual(values = c(MP_COLOURS, Unassigned = "#D6DADE"), labels = c(MP_LABELS, Unassigned = "Unassigned"),
                    drop = FALSE, name = "Program") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Program fraction") + theme_spg()

# f: residual member coordination above matched-module means in the NOA cases.
z <- fread(file.path(RESULTS, "coordination", "coherence", "external_NOA", "coherence_results.tsv"))[
  model == "technical_cycle" & level == "primary" & module %in% names(MP_LABELS)]
z[, case := sub("GSE326372_testis_", "Case ", donor)]
p_f <- ggplot(z, aes(state, excess_rho, colour = module, group = module)) +
  geom_hline(yintercept = 0, colour = "#AAAAAA", linewidth = .3) +
  geom_line(linewidth = .4, na.rm = TRUE) + geom_point(size = 1.5, na.rm = TRUE) + facet_wrap(~case, nrow = 1) +
  scale_x_discrete(limits = names(STATE_LABELS), labels = STATE_LABELS) +
  scale_colour_manual(values = MP_COLOURS, labels = MP_LABELS, name = "Program") +
  labs(x = "State", y = "Mean ρ − matched null") + theme_spg()

# g: projected non-MP AUCell of the twelve structural regulons in the four U populations.
tfs <- fread(file.path(RESULTS, "communication", "regulon", "structural_TF_selection.tsv"))[selected == TRUE, TF]
z <- fread(file.path(RESULTS, "communication", "regulon", "projection_U_AUCell.tsv"))
z[, TF := factor(TF, levels = tfs)]
z[, row := donor_rows(sample_key, n_cells)]
p_g <- ggplot(z, aes(TF, row, fill = mean_AUCell)) + geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradient(low = "#EDF2F0", high = "#416F7A", na.value = "#B5B8BB", name = "Non-MP\nAUCell") +
  labs(x = NULL, y = NULL) + theme_spg() + TILE + VERTICAL_X + guides(fill = colourbar_guide())

panels <- list(p_a, p_b, p_c, p_d, p_e, p_f, p_g)
compose_figure(panels, rows = list(1:2, 3, 4:5, 6, 7), heights = c(52, 44, 57, 45, 38),
               widths = list(c(.85, 1.4), 1, c(1.05, 1), 1, 1), file = file.path(FIGURES, "figure4"))
cat("figure4: 7 panels;", length(tfs), "regulons;", uniqueN(z$sample_key), "GSE202647 donors\n")
