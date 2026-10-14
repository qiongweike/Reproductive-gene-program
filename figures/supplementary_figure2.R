# Supplementary Figure 2: trajectories of the two further selected donors and regulon dynamics.
# Inputs: results/trajectory/<donor>.primary.pca.tsv.gz and .primary.slingshot.rds (A2, N2),
#   results/communication/regulon_trajectory_curves.tsv and
#   results/communication/regulon/structural_TF_selection.tsv.
# Output: results/figures/supplementary_figure2.pdf and supplementary_figure2.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({library(data.table); library(slingshot)})

TRAJ <- file.path(RESULTS, "trajectory")
TRAJECTORY_DONORS <- c("GSE124263__A2", "GSE153947__N2")   # selected donors besides A1 (Fig. 3a)
STATE_COLOURS <- c(undifferentiated = "#48799B", transition = "#D3A249", differentiating = "#459386")
STATE_LABELS <- c(undifferentiated = "U", transition = "T", differentiating = "D")
DONOR_COLOURS <- c(A1 = "#4477AA", A2 = "#66AADD", N1 = "#228833", N2 = "#44AA99", N3 = "#999933",
                   Donor30 = "#AA3377", Donor5 = "#CC6677")

# a, b: independent-feature PCA with states and the Slingshot cluster-centre MST.
trajectory_panel <- function(sample_key) {
  cells <- fread(file.path(TRAJ, paste0(sample_key, ".primary.pca.tsv.gz")))[state != "unknown"]
  fit <- SlingshotDataSet(readRDS(file.path(TRAJ, paste0(sample_key, ".primary.slingshot.rds"))))
  X <- reducedDim(fit)[, 1:2]; weights <- slingClusterLabels(fit); mst <- slingMST(fit)
  centres <- t(vapply(rownames(mst), function(k) apply(X, 2, weighted.mean, w = weights[, k]), numeric(2)))
  ij <- which(upper.tri(mst) & mst == 1, arr.ind = TRUE)
  edges <- data.table(x = centres[ij[, 1], 1], y = centres[ij[, 1], 2], xend = centres[ij[, 2], 1], yend = centres[ij[, 2], 2])
  ggplot(cells, aes(PC1, PC2)) + geom_point(aes(colour = state), size = .25, alpha = .65, stroke = 0) +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "#26313C", linewidth = .4) +
    scale_colour_manual(values = STATE_COLOURS, breaks = names(STATE_COLOURS), labels = STATE_LABELS, name = "State") +
    labs(x = paste0(sub(".*__", "", sample_key), " · PC1"), y = "PC2") + theme_spg() +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.1, alpha = 1)))
}

# c: mean non-MP AUCell of the twelve structural regulons per pseudotime bin, all eligible lineages.
tfs <- fread(file.path(RESULTS, "communication", "regulon", "structural_TF_selection.tsv"))[selected == TRUE, TF]
curves <- fread(file.path(RESULTS, "communication", "regulon_trajectory_curves.tsv"))
stopifnot(setequal(curves$TF, tfs))
curves[, TF := factor(TF, levels = tfs)]
curves[, series := paste(donor, lineage)]
curves[, branch := sub("Lineage", "", lineage)]
p_c <- ggplot(curves, aes(bin, mean_AUCell, colour = donor, group = series, linetype = branch)) +
  geom_line(linewidth = .4) + geom_point(size = .5) + facet_wrap(~TF, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = DONOR_COLOURS, name = "Donor") +
  scale_linetype_manual(values = c("1" = "solid", "3" = "dashed"), name = "Lineage") +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  labs(x = "Within-lineage quantile bin", y = "Non-MP AUCell") + theme_spg() + theme(panel.spacing = unit(3, "mm")) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2))

panels <- c(lapply(TRAJECTORY_DONORS, trajectory_panel), list(p_c))
compose_figure(panels, rows = list(1:2, 3), heights = c(62, 108), file = file.path(FIGURES, "supplementary_figure2"))
cat("supplementary_figure2: 3 panels;", uniqueN(curves$series), "lineages;", length(tfs), "regulons\n")
