# Figure 6: normal human spatial context and measurement coverage of MP3 (panels a-j).
# Inputs (results/spatial/): bead maps with BANKSY clusters and the fixed ROI windows,
# cluster marker profiles and RCTD weight summaries, member Hotspot consistency,
# Sertoli neighbour effects, random-forest gains, marker-family detection and
# ligand-receptor bead counts. Output: results/figures/figure6.pdf and figure6.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({
  library(data.table); library(scales); library(ggrepel); library(ComplexHeatmap); library(circlize)
})

SPATIAL       <- file.path(RESULTS, "spatial")
SUMMARIES     <- file.path(SPATIAL, "summaries")
NEIGHBOURHOOD <- file.path(SPATIAL, "neighbourhood")

PAGE_WIDTH  <- 183                                  # mm
PAGE_PAD    <- 1.5
TAG_BAND    <- 3                                    # mm above each panel kept for its letter
ROWS        <- list(c(1, 2), c(3, 4), 5, c(6, 7), c(8, 9), 10)
ROW_HEIGHTS <- c(49, 42, 70, 44, 44, 51)
ROW_WIDTHS  <- list(c(1, 1), c(1, 1), 1, c(1.12, 1), c(1, 1), 1)

ARRAYS        <- c(Puck5 = "Array 1", Puck6 = "Array 2")
ARRAY_COLOURS <- c(Puck5 = "#357E89", Puck6 = "#C3794E")
MP3_PALETTE   <- c("#E3E6E7", "#89B9B1", "#367D87", "#123F5F")
MAP_BREAKS    <- c(100, 200, 300, 400)
ROI_STEP      <- 20
LABELLED_MEMBERS <- c("FGFR3", "COMMD3", "IRF2BPL")
REPEL_SEED    <- 20261004
CELL_TYPES    <- c(SPG = "SPG", Spermatocyte = "Spermatocyte", Spermatid = "Spermatid", Sertoli = "Sertoli",
                   Peritubular_smooth_muscle = "PTM", Leydig = "Leydig", Endothelial = "Endothelial",
                   Immune = "Immune")
COVERAGE_FAMILIES <- c(MP3 = "MP3", SPG = "SPG", Spermatocyte = "Spermatocyte", Spermatid = "Spermatid",
                       Sertoli = "Sertoli", Leydig = "Leydig", Peritubular_smooth_muscle = "PTM",
                       Endothelial = "Endothelial", Immune = "Immune")
LR_PAIRS      <- c(FGF2_FGFR3 = "FGF2–FGFR3", INHABB_ACVR2B = "(INHA+INHBB)–ACVR2B",
                   INHBA_ACVR1B_ACVR2A = "INHBA–(ACVR1B+ACVR2A)", INHBA_ACVR1B_ACVR2B = "INHBA–(ACVR1B+ACVR2B)",
                   INHBB_ACVR1B_ACVR2A = "INHBB–(ACVR1B+ACVR2A)", INHBB_ACVR1B_ACVR2B = "INHBB–(ACVR1B+ACVR2B)")
LR_MEASURES   <- c("Ligand-positive soma" = "Ligand+", "Receptor-positive SPG" = "Receptor+",
                   "Compatible receiver" = "Neighbor+")
LR_OFFSET     <- .18                                # array 1 above, array 2 below each pair

# Fig. 6 panels use 6.5-pt tick and legend text and compact legend keys.
spatial_theme <- function() theme_spg() + theme(
  axis.text = element_text(size = 6.5, colour = INK), axis.title = element_text(size = 7),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5),
  legend.key.height = unit(2.4, "mm"), legend.key.width = unit(2.4, "mm"), legend.box.spacing = unit(1, "mm"),
  legend.spacing.y = unit(.6, "mm"), plot.margin = margin(2, 2, 2, 2, "pt"), plot.tag = element_blank())
array_colour <- function() scale_colour_manual(values = ARRAY_COLOURS, labels = ARRAYS, name = "Array")
read_summary <- function(suffix) {
  rbindlist(lapply(names(ARRAYS), function(a) fread(file.path(SUMMARIES, paste0(a, suffix)))))
}
# ComplexHeatmap measures text on pdf(NULL), which cannot resolve Liberation Sans; measure on Cairo instead.
assignInNamespace("dev.null", function(...) cairo_pdf(tempfile(fileext = ".pdf"), family = FONT, ...),
                  ns = "ComplexHeatmap")

maps <- read_summary("_ROI_maps.tsv.gz")
roi <- fread(file.path(SUMMARIES, "selected_ROIs.tsv"))
mp3_max <- max(maps$MP_3)
stopifnot(min(maps$MP_3) == 0)

# a, b: BANKSY RNA clusters on every bead; the box marks the fixed display window.
cluster_map <- function(a) {
  z <- maps[array == a]
  ids <- sort(unique(z$domain))
  z[, cluster := factor(domain, levels = ids, labels = sprintf("C%02d", ids))]
  ggplot(z, aes(x_nn, y_nn, colour = cluster)) +
    geom_point(size = .13) +
    geom_rect(data = roi[array == a], aes(xmin = xmin_nn, xmax = xmax_nn, ymin = ymin_nn, ymax = ymax_nn),
              inherit.aes = FALSE, fill = NA, colour = "#263740", linewidth = .4) +
    coord_fixed(expand = FALSE) +
    scale_x_continuous(breaks = MAP_BREAKS) + scale_y_continuous(breaks = MAP_BREAKS) +
    scale_colour_manual(values = setNames(hcl.colors(length(ids), palette = "Dark 3"), levels(z$cluster)),
                        name = paste0(ARRAYS[a], "\nRNA cluster"), drop = FALSE) +
    labs(x = "x / median NN", y = "y / median NN") +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.4))) + spatial_theme()
}

# c, d: measured MP3 score in the windows on one shared scale; rings mark high-confidence SPG beads.
mp3_window <- function(a) {
  z <- maps[array == a & in_ROI == TRUE]
  r <- roi[array == a]
  ggplot(z, aes(x_nn, y_nn)) +
    geom_point(aes(colour = MP_3), size = .7, shape = 16, stroke = 0) +
    geom_point(data = z[HC_SPG == TRUE], aes(shape = "HC-SPG"), size = 1.1, stroke = .22, colour = "#16272E") +
    scale_colour_gradientn(colours = MP3_PALETTE, limits = c(0, mp3_max),
                           breaks = c(0, round(mp3_max / 2, 2), round(mp3_max, 2)), name = "MP3\nMean ln(1+CP10k)") +
    scale_shape_manual(values = c("HC-SPG" = 1), name = "ROI overlay") +
    scale_x_continuous(limits = c(r$xmin_nn, r$xmax_nn), breaks = seq(r$xmin_nn, r$xmax_nn, ROI_STEP),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(r$ymin_nn, r$ymax_nn), breaks = seq(r$ymin_nn, r$ymax_nn, ROI_STEP),
                       expand = c(0, 0)) +
    coord_fixed() +
    labs(x = paste0(ARRAYS[a], " · x / median NN"), y = "y / median NN") +
    guides(shape = guide_legend(order = 1, override.aes = list(size = 1.2, colour = "#16272E")),
           colour = colourbar_guide(order = 2)) +
    spatial_theme()
}

P <- list(a = cluster_map("Puck5"), b = cluster_map("Puck6"), c = mp3_window("Puck5"), d = mp3_window("Puck6"))

# e: bead counts, within-array marker-family Z and mean reported RCTD weights among HC beads, per cluster.
counts  <- setorder(fread(file.path(SPATIAL, "cluster_bead_counts.tsv")), array, cluster)
markers <- fread(file.path(SPATIAL, "cluster_marker_profiles.tsv"))
weights <- fread(file.path(SPATIAL, "cluster_RCTD_composition.tsv"))
stopifnot(nrow(counts) == 27L, nrow(markers) == 216L, nrow(weights) == 216L)
cluster_matrix <- function(d, key, value, keys) {
  m <- as.matrix(dcast(d, as.formula(paste(key, "~ cluster_id")), value.var = value), rownames = key)
  m <- m[keys, counts$cluster_id]
  stopifnot(!anyNA(m))
  colnames(m) <- sub("Puck[56]\\.", "", colnames(m))
  m
}
marker_z <- cluster_matrix(markers, "family", "z_among_domains", CELL_TYPES)
weight <- cluster_matrix(weights, "cell_type", "mean_reported_weight", names(CELL_TYPES))
member_n <- markers[, .(n = unique(n_members)), by = family]
rownames(marker_z) <- paste0(CELL_TYPES, " (", member_n$n[match(CELL_TYPES, member_n$family)], ")")
rownames(weight) <- CELL_TYPES
z_max <- max(abs(marker_z))
array_split <- factor(ARRAYS[counts$array], levels = ARRAYS)
small <- gpar(fontsize = 6.5)
named <- gpar(fontsize = 6.5, fontfamily = FONT)
count_bars <- function(n, fill, limit, tick) {
  anno_barplot(n, gp = gpar(fill = fill, col = NA), height = unit(6, "mm"), ylim = c(0, limit),
               axis_param = list(at = c(0, tick), labels = c("0", paste0(tick / 1000, "k")), gp = small))
}
bars <- HeatmapAnnotation(`All n` = count_bars(counts$n_all, "#7F909A", 11000, 10000),
                          `HC n` = count_bars(counts$n_HC, "#3C737F", 7000, 6000),
                          gap = unit(2, "mm"), annotation_name_gp = small, annotation_name_side = "left")
# ComplexHeatmap raises 13-mm colour bars to its own minimum height and warns; the published panel has the same bars.
legend_param <- list(title_gp = small, labels_gp = small, legend_height = unit(13, "mm"))
cluster_heatmap <- function(m, name, colours, ...) {
  Heatmap(m, name = name, col = colours, cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE,
          column_split = array_split, row_names_side = "left", row_names_gp = named, height = unit(18, "mm"),
          width = unit(126, "mm"), column_gap = unit(3, "mm"), ...)
}
marker_map <- cluster_heatmap(marker_z, "Marker Z", colorRamp2(c(-z_max, 0, z_max), c("#4678A0", "white", "#B96447")),
                              column_title_gp = gpar(fontsize = 7), show_column_names = FALSE, top_annotation = bars,
                              heatmap_legend_param = legend_param)
weight_map <- cluster_heatmap(weight, "HC mean\nweight", colorRamp2(c(0, .5, 1), c("#F0F5F4", "#8BB9B0", "#1D5868")),
                              column_title = NULL, column_names_gp = named, column_names_rot = 90,
                              heatmap_legend_param = c(legend_param, list(at = c(0, .5, 1))))
# grid.grabExpr() returns to the device that was current, so one is open while the heatmap is captured.
cairo_pdf(tempfile(fileext = ".pdf"), family = FONT)
P$e <- grid.grabExpr(
  draw(marker_map %v% weight_map, heatmap_legend_side = "right", merge_legends = TRUE,
       align_heatmap_legend = "heatmap_center", ht_gap = unit(3, "mm"), padding = unit(c(1, 1, 1, 1), "mm")),
  width = (PAGE_WIDTH - 2 * PAGE_PAD) / 25.4, height = (ROW_HEIGHTS[3] - TAG_BAND) / 25.4,
  device = function(width, height) cairo_pdf(tempfile(fileext = ".pdf"), width, height, family = FONT))
invisible(dev.off())

# f: cross-array Hotspot autocorrelation of the 49 members, whole arrays and HC-SPG beads.
consistency <- fread(file.path(NEIGHBOURHOOD, "MP3_spatial_member_consistency.tsv"))
consistency <- consistency[!is.na(Puck5_C) & !is.na(Puck6_C)]
stopifnot(all(consistency[, .N, by = population]$N == 49L))
consistency[, population := factor(population, levels = c("whole_array", "HC_SPG"),
                                   labels = c("Whole array", "HC-SPG"))]
P$f <- ggplot(consistency, aes(Puck5_C, Puck6_C)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#B5BDC1", linewidth = .25) +
  geom_point(size = .7, alpha = .85, colour = "#AA7858") +
  geom_text_repel(data = consistency[gene %in% LABELLED_MEMBERS], aes(label = gene), seed = REPEL_SEED,
                  size = 6.5 / .pt, max.overlaps = Inf, min.segment.length = 0, box.padding = unit(.2, "lines"),
                  point.padding = .08, segment.size = .2) +
  facet_wrap(~population, nrow = 1, scales = "free") +
  scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 3) +
  labs(x = "Array 1 C", y = "Array 2 C") + spatial_theme()

# g: out-of-fold Sertoli-neighbour weight per SD of MP3 residual in four annuli, 95% block-bootstrap CI.
neighbour <- fread(file.path(NEIGHBOURHOOD, "neighbour_effects.tsv"))
neighbour <- neighbour[type == "Sertoli" & method == "OOF_continuous_secondary"]
neighbour[, ring := factor(mid, levels = c(1, 3, 5, 7), labels = c("0–2", "2–4", "4–6", "6–8"))]
dodge <- position_dodge(.15)
P$g <- ggplot(neighbour, aes(ring, effect, colour = array, group = array)) +
  geom_hline(yintercept = 0, colour = "#ABB4BA", linewidth = .25) +
  geom_errorbar(aes(ymin = low, ymax = high), width = .09, linewidth = .35, position = dodge) +
  geom_line(linewidth = .4, position = dodge) + geom_point(size = 1.25, position = dodge) +
  array_colour() + scale_y_continuous(n.breaks = 3) +
  labs(x = "Distance / median NN", y = "Sertoli weight / SD MP3\n(OOF; 95% block CI)") + spatial_theme()

# h: neighbourhood random-forest gain per held-out fold, MP-excluded RCTD inputs.
gain <- read_summary("_neighbourhood_gain.tsv")
gain <- gain[program %in% c("MP_1", "MP_2", "MP_3") & identity_version == "MP_excluded_RCTD"]
gain[, `:=`(program = factor(program, levels = c("MP_1", "MP_2", "MP_3"), labels = c("MP1", "MP2", "MP3")),
            fold = factor(fold, levels = 0:3, labels = 1:4))]
P$h <- ggplot(gain, aes(program, gain_pp, colour = array, shape = fold)) +
  geom_hline(yintercept = 0, colour = "#ADB6BB", linewidth = .25) +
  geom_point(size = 1.5, position = position_dodge(width = .65)) +
  array_colour() + scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 18), name = "Fold") +
  scale_y_continuous(n.breaks = 3) + labs(x = NULL, y = "ΔR² (percentage points)") +
  spatial_theme() + theme(legend.spacing.y = unit(.3, "mm"))

# i: equally weighted member detection for MP3 and the identity-marker families; (n) measured members.
coverage <- fread(file.path(SUMMARIES, "member_family_detection.tsv"))
n_measured <- unique(coverage[, .(family, n_measured)])
family_labels <- setNames(paste0(COVERAGE_FAMILIES[n_measured$family], " (", n_measured$n_measured, ")"),
                          n_measured$family)
family_order <- rev(names(COVERAGE_FAMILIES))
coverage[, family := factor(family, levels = family_order, labels = family_labels[family_order])]
P$i <- ggplot(coverage, aes(mean_member_detection, family, colour = array, shape = population)) +
  geom_point(size = 1.3, position = position_dodge(.5)) +
  array_colour() + scale_shape_manual(values = c("HC-SPG" = 16, "Other HC" = 1), name = "Population") +
  scale_x_continuous(labels = percent_format(accuracy = 1), n.breaks = 3) +
  labs(x = "Mean member detection", y = NULL) + spatial_theme()

# j: six ligand-receptor pairs, each with its own n/N denominator for senders, receivers and neighbours.
lr <- fread(file.path(SUMMARIES, "ligand_receptor_bead_counts.tsv"))
stopifnot(nrow(lr) == 36L, all(lr$denominator > 0))
lr[, `:=`(y = match(interaction, rev(names(LR_PAIRS))) + ifelse(array == "Puck5", LR_OFFSET, -LR_OFFSET),
          measure = factor(LR_MEASURES[measure], levels = LR_MEASURES), label = paste0(numerator, "/", denominator))]
P$j <- ggplot(lr, aes(numerator, y, fill = array)) +
  geom_col(orientation = "y", width = .3, position = "identity") +
  geom_text(aes(label = label), hjust = -.13, size = 6.5 / .pt, show.legend = FALSE) +
  facet_wrap(~measure, nrow = 1, scales = "free_x") +
  scale_y_continuous(breaks = seq_along(LR_PAIRS), labels = rev(unname(LR_PAIRS)), expand = expansion(add = .55)) +
  scale_x_continuous(n.breaks = 3, expand = expansion(mult = c(.01, .72))) +
  scale_fill_manual(values = ARRAY_COLOURS, breaks = names(ARRAYS), labels = ARRAYS, name = "Array") +
  labs(x = "Measured beads", y = NULL) + spatial_theme()

## Page: rows of fixed height as in compose_figure(), but with this figure's 1.5-mm padding and
## 8-pt letters; panels are drawn as built (no clean_panel), and panel e is already a grob.
compose_rows <- function(panels, rows, heights, widths, file, width = PAGE_WIDTH, pad = PAGE_PAD) {
  height <- sum(heights) + 2 * pad
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  metrics <- paste0(file, "_metrics.pdf")
  cairo_pdf(metrics, width / 25.4, height / 25.4, family = FONT)
  parts <- list(); top <- pad
  for (i in seq_along(rows)) {
    w <- (width - 2 * pad) * widths[[i]] / sum(widths[[i]]); left <- pad
    for (j in seq_along(rows[[i]])) {
      k <- rows[[i]][j]
      g <- if (inherits(panels[[k]], "grob")) panels[[k]] else ggplotGrob(panels[[k]])
      parts <- c(parts, list(
        grobTree(g, vp = viewport(x = unit(left, "mm"), y = unit(height - top - TAG_BAND, "mm"),
                                  width = unit(w[j], "mm"), height = unit(heights[i] - TAG_BAND, "mm"),
                                  just = c("left", "top"), clip = "off")),
        textGrob(letters[k], x = unit(left + .25, "mm"), y = unit(height - top - .4, "mm"),
                 just = c("left", "top"), gp = gpar(fontfamily = FONT, fontsize = 8, fontface = "bold", col = INK))))
      left <- left + w[j]
    }
    top <- top + heights[i]
  }
  dev.off(); unlink(metrics)
  page <- do.call(grobTree, parts)
  cairo_pdf(paste0(file, ".pdf"), width / 25.4, height / 25.4, family = FONT)
  grid.newpage(); grid.draw(page); dev.off()
  pdf_to_png(file)
  invisible(page)
}

compose_rows(P, ROWS, ROW_HEIGHTS, ROW_WIDTHS, file.path(FIGURES, "figure6"))
cat("figure6:", length(P), "panels ->", file.path(FIGURES, "figure6.{pdf,png}"), "\n")
