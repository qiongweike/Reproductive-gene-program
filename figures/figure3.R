# Figure 3: donor-resolved program dynamics along independent-feature Slingshot lineages.
# Inputs: results/trajectory/ (donor manifest, selected primary lineages, PCA and Slingshot fit of
#   the largest donor, bin scores, member trajectories and peak bins, Palantir peak order,
#   donor-state means) and results/programs/MP_genes.tsv.
# Output: results/figures/figure3.pdf and figure3.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({library(data.table); library(slingshot); library(ggbeeswarm)})

TRAJ <- file.path(RESULTS, "trajectory")
N_BINS <- 10
Z_LIMITS <- c(-2.5, 2.5)
MP_LABELS <- c(MP_1 = "MP1", MP_2 = "MP2", MP_3 = "MP3")
MP_COLOURS <- c(MP_1 = "#368A88", MP_2 = "#7763A6", MP_3 = "#C47745")
DONOR_COLOURS <- c(A1 = "#4477AA", A2 = "#66AADD", N1 = "#228833", N2 = "#44AA99", N3 = "#999933",
                   Donor30 = "#AA3377", Donor5 = "#CC6677")
STATE_COLOURS <- c(undifferentiated = "#48799B", transition = "#D3A249", differentiating = "#459386")
STATE_LABELS <- c(undifferentiated = "U", transition = "T", differentiating = "D")
METHOD_COLOURS <- c(Slingshot = "#345F80", Palantir = "#A8536F")
MODULE_LABELS <- c(initiation = "Initiation", formation = "Formation", fusion_SNARE = "SNARE fusion",
                   fusion_tether = "Fusion tether", acidification = "Acidification", degradation = "Degradation",
                   iron_handling = "Iron handling", antioxidant = "Antioxidant", lipid_susceptibility = "Lipid response")

SMALL_TICKS <- theme(axis.text = element_text(size = 6.5))
POINT_KEYS <- guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.2, alpha = 1)))
TILE <- theme(axis.ticks = element_blank(), axis.line = element_blank())
z_fill <- function(name, ...) scale_fill_gradient2(low = "#355F8D", mid = "#F7F6F2", high = "#B85B46", midpoint = 0,
                                                   limits = Z_LIMITS, oob = scales::squish, name = name, ...)

read_traj <- function(file) fread(file.path(TRAJ, file))
donor_of <- function(sample_key) sub(".*__", "", sample_key)
zscore <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  if (sd(x, na.rm = TRUE) < 1e-12) return(0)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

manifest <- read_traj("donor_manifest.tsv")
selected <- read_traj("selected_primary_lineages.tsv")   # largest eligible primary lineage per donor
keep_selected <- function(x) merge(x, selected[, .(sample_key, lineage)], by = c("sample_key", "lineage"), sort = FALSE)
bins <- keep_selected(read_traj("trajectory_bin_scores.tsv")[scenario == "primary"])[eligible == TRUE & bin_eligible == TRUE]
members <- keep_selected(read_traj("member_trajectory_bins.tsv.gz")[trajectory_eligible == TRUE & bin_eligible == TRUE])
peaks <- keep_selected(read_traj("member_peak_order.tsv"))[trajectory_eligible == TRUE & nonconstant == TRUE]
mp_genes <- fread(file.path(RESULTS, "programs", "MP_genes.tsv"))
samples <- sort(unique(members$sample_key))
column_of <- function(sample_key, bin) (match(sample_key, samples) - 1) * N_BINS + bin

# a: PCA of the largest donor with states, cluster-centre MST and the U-rich root cluster.
largest <- manifest[order(-n_spg, sample_key), sample_key][1]
cells <- read_traj(paste0(largest, ".primary.pca.tsv.gz"))[state != "unknown"]
fit <- SlingshotDataSet(readRDS(file.path(TRAJ, paste0(largest, ".primary.slingshot.rds"))))
X <- reducedDim(fit)[, 1:2]; weights <- slingClusterLabels(fit); mst <- slingMST(fit)
centres <- t(vapply(rownames(mst), function(k) apply(X, 2, weighted.mean, w = weights[, k]), numeric(2)))
ij <- which(upper.tri(mst) & mst == 1, arr.ind = TRUE)
edges <- data.table(x = centres[ij[, 1], 1], y = centres[ij[, 1], 2], xend = centres[ij[, 2], 1], yend = centres[ij[, 2], 2])
root <- data.table(PC1 = centres[slingParams(fit)$start.clus, 1], PC2 = centres[slingParams(fit)$start.clus, 2])
p_a <- ggplot(cells, aes(PC1, PC2)) +
  geom_point(aes(colour = state), size = .28, alpha = .6, stroke = 0) +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "#26313C", linewidth = .4) +
  geom_point(data = root, aes(PC1, PC2, shape = "U-rich root"), fill = "white", size = 2.2, stroke = .5, inherit.aes = FALSE) +
  scale_shape_manual(values = c("U-rich root" = 23), name = NULL) +
  scale_colour_manual(values = STATE_COLOURS, breaks = names(STATE_COLOURS), labels = STATE_LABELS, name = "State") +
  labs(x = "PC 1", y = "PC 2") + theme_spg() + SMALL_TICKS +
  guides(shape = guide_legend(order = 1), colour = guide_legend(ncol = 1, override.aes = list(size = 1.2, alpha = 1), order = 2))

# b: observed program-score means per pseudotime bin.
z <- bins[metric %in% names(MP_LABELS), .(mean_score = mean(mean_score)), by = .(sample_key, metric, bin)]
z[, donor := donor_of(sample_key)]
z[, metric := factor(metric, levels = names(MP_LABELS), labels = MP_LABELS)]
p_b <- ggplot(z, aes(bin, mean_score, colour = donor, group = donor)) +
  geom_line(linewidth = .55) + geom_point(size = .65) + facet_wrap(~metric, nrow = 1) +
  scale_colour_manual(values = DONOR_COLOURS, name = "Donor") + scale_x_continuous(breaks = c(1, N_BINS)) +
  labs(x = "Within-lineage quantile bin", y = "Fixed program score") + theme_spg() + SMALL_TICKS + POINT_KEYS

# c, d: all MP2 and MP3 members, z-scored within donor and gene, ordered by median peak bin.
member_heatmap <- function(mp) {
  z <- members[MP == mp, .(mean_logCP10k = mean(mean_logCP10k)), by = .(sample_key, gene, bin)]
  z[, scaled := zscore(mean_logCP10k), by = .(sample_key, gene)]
  ord <- peaks[MP == mp, .(peak = median(peak_bin)), by = gene][order(peak, gene), gene]
  ord <- c(ord, setdiff(mp_genes[MP == mp, gene], ord))
  stopifnot(length(ord) == mp_genes[MP == mp, .N])
  z[, gene := factor(gene, levels = ord)]
  z[, column := column_of(sample_key, bin)]
  ggplot(z, aes(gene, column, fill = scaled)) + geom_tile() +
    geom_hline(yintercept = seq(N_BINS + .5, length(samples) * N_BINS - 1, N_BINS), colour = "white", linewidth = .7) +
    z_fill("Within-donor\ngene z score", na.value = "#D0D0D0") +
    scale_y_reverse(breaks = seq(N_BINS / 2 + .5, length(samples) * N_BINS, N_BINS), labels = donor_of(samples), expand = c(0, 0)) +
    labs(x = NULL, y = "Early to late bins") + theme_spg() + TILE +
    theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = .5)) + SMALL_TICKS +
    guides(fill = colourbar_guide())
}
p_c <- member_heatmap("MP_2")
p_d <- member_heatmap("MP_3")

# e: peak bins of non-constant member profiles; jitter only along the categorical program axis.
z <- copy(peaks)[, donor := donor_of(sample_key)][, MP := factor(MP, levels = names(MP_COLOURS))]
p_e <- ggplot(z, aes(MP, peak_bin, colour = MP)) +
  geom_quasirandom(method = "maxout", size = .6, width = .28) +
  scale_colour_manual(values = MP_COLOURS, labels = MP_LABELS, name = "Program") +
  scale_x_discrete(labels = MP_LABELS) + scale_y_continuous(breaks = c(1, 4, 7, 10), limits = c(.5, 10.5)) +
  facet_wrap(~donor, nrow = 1) + coord_flip() +
  labs(x = NULL, y = "Bin of maximum member expression") + theme_spg() + SMALL_TICKS + POINT_KEYS

# f: MP3 minus MP2 median member-peak bin, Slingshot (selected lineage) versus Palantir (whole donor).
po <- read_traj("palantir/MP2_MP3_peak_order.tsv")[trajectory_eligible == TRUE]
z <- rbind(po[method == "Palantir" & scope == "global_all_cells"],
           keep_selected(po[method == "Slingshot" & scope == "paired_primary_shared_cells"]), use.names = TRUE)
z <- z[, .(sample_key, method, peak_difference = MP3_minus_MP2_peak)]
stopifnot(!anyDuplicated(z[, .(sample_key, method)]))
z <- merge(CJ(sample_key = manifest$sample_key, method = names(METHOD_COLOURS)), z, by = c("sample_key", "method"), all.x = TRUE)
z[, donor := factor(donor_of(sample_key), levels = rev(manifest$donor))]
z[, donor_position := as.numeric(donor) + ifelse(method == "Slingshot", .13, -.13)]
empty <- z[, .(n_finite = sum(is.finite(peak_difference))), by = donor][n_finite == 0]
p_f <- ggplot(z, aes(peak_difference, donor_position, shape = method, colour = method)) +
  geom_text(data = empty, aes(x = .8, y = as.numeric(donor), label = "NA"), inherit.aes = FALSE,
            size = 7 / .pt, family = FONT, colour = "#777777") +
  geom_vline(xintercept = 0, linetype = 2, colour = "#BCBCBC", linewidth = .3) +
  geom_line(aes(group = sample_key), colour = "#B4B4B4", linewidth = .4, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_shape_manual(values = c(Slingshot = 16, Palantir = 1), name = "Method") +
  scale_colour_manual(values = METHOD_COLOURS, name = "Method") +
  scale_x_continuous(breaks = seq(floor(min(z$peak_difference, na.rm = TRUE) / 3) * 3,
                                  ceiling(max(z$peak_difference, na.rm = TRUE) / 3) * 3, 3),
                     limits = range(c(0, z$peak_difference), na.rm = TRUE) + c(-.6, .6)) +
  scale_y_continuous(breaks = seq_along(levels(z$donor)), labels = levels(z$donor), limits = c(.55, nlevels(z$donor) + .45)) +
  labs(x = "MP3 minus MP2 median member-peak bin", y = NULL) + theme_spg() + SMALL_TICKS + POINT_KEYS

# g: program means in donor-state groups with at least 30 cells.
z <- read_traj("MP_donor_state_scores.tsv")[eligible == TRUE & state %in% names(STATE_LABELS)]
z <- melt(z, id.vars = c("sample_key", "donor", "state"), measure.vars = names(MP_LABELS), variable.name = "MP", value.name = "mean_score")
z[, state := factor(state, levels = names(STATE_LABELS), labels = STATE_LABELS)]
p_g <- ggplot(z, aes(state, mean_score, group = sample_key, colour = donor)) +
  geom_line(linewidth = .4, alpha = .8) + geom_point(size = 1.15) +
  facet_wrap(~MP, nrow = 1, labeller = as_labeller(MP_LABELS)) +
  scale_colour_manual(values = DONOR_COLOURS, name = "Donor") +
  labs(x = "Observed state (at least 30 cells)", y = "Mean fixed program score") + theme_spg() + SMALL_TICKS + POINT_KEYS

# h: nine predefined autophagy, iron, antioxidant and lipid modules, z-scored within donor.
z <- bins[metric %in% names(MODULE_LABELS), .(value = mean(mean_score)), by = .(sample_key, metric, bin)]
stopifnot(uniqueN(z$metric) == length(MODULE_LABELS))
z[, scaled := zscore(value), by = .(sample_key, metric)]
z[, column := column_of(sample_key, bin)]
z[, metric := factor(metric, levels = rev(names(MODULE_LABELS)))]
p_h <- ggplot(z, aes(column, metric, fill = scaled)) + geom_tile() +
  geom_vline(xintercept = seq(N_BINS + .5, length(samples) * N_BINS - 1, N_BINS), colour = "white", linewidth = .7) +
  z_fill("Within-donor\nmodule z score") + scale_y_discrete(labels = MODULE_LABELS) +
  scale_x_continuous(breaks = seq(N_BINS / 2 + .5, length(samples) * N_BINS, N_BINS), labels = donor_of(samples), expand = c(0, 0)) +
  labs(x = "Early to late bins within donor", y = NULL) + theme_spg() + TILE + SMALL_TICKS + guides(fill = colourbar_guide())

panels <- list(p_a, p_b, p_c, p_d, p_e, p_f, p_g, p_h)
compose_figure(panels, rows = list(1:2, 3, 4, 5:6, 7:8), heights = c(51, 53, 53, 53, 51),
               widths = list(c(1, 1.25), 1, 1, c(1.25, 1), c(1, 1.1)), file = file.path(FIGURES, "figure3"))
cat("figure3: 8 panels;", nrow(peaks), "member peak profiles;", length(samples), "selected lineages\n")
