# Figure 7: matched virtual-knockout responses in the four reference donor networks.
#   a MP3 response excess of the 49 members; b-g FGFR3 and COMMD3 against their five
#   matched controls; h response-profile PCA of the focal knockouts in donor A1.
# Inputs: results/knockout/panels/ (member, program, function and gene response tables; PCA coordinates and variance).
# Output: results/figures/figure7.pdf and figure7.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({library(data.table); library(ggrepel); library(ComplexHeatmap); library(circlize)})

PANELS <- file.path(RESULTS, "knockout", "panels")
DONORS <- c("A1", "A2", "N2", "N3")
FOCUS <- c("FGFR3", "COMMD3")
PROGRAMS <- c(MP_1 = "MP1", MP_2 = "MP2", MP_3 = "MP3", fixed_143 = "143-gene")
FUNCTIONS <- c(   # fixed display functions of panel d, top to bottom
  "GO:1902850" = "Mitotic microtubules", "GO:0034502" = "Chromosome targeting",
  "R-HSA-141424" = "Kinetochore signalling", "R-HSA-2467813" = "Sister separation",
  "R-HSA-69620" = "Cell-cycle checkpoints", "GO:0000776" = "Kinetochore",
  "GO:0051304" = "Chromosome separation", "GO:0033047" = "Sister segregation regulation",
  "GO:0007088" = "Nuclear division regulation", "GO:0050673" = "Epithelial proliferation",
  "GO:0007091" = "Metaphase–anaphase", "GO:0051303" = "Chromosome localization",
  "GO:0000070" = "Mitotic sister segregation", "R-HSA-69206" = "G1/S transition")
EXTERNAL_GENES <- c("GAN", "IGFBP6", "MFHAS1", "ZWINT", "PAGE1", "B2M", "KLF6", "MTCH1",
                    "CDC45", "GINS2", "HOXC9", "TMED9")   # fixed display genes of panel e
IRON_GENES <- c("FTH1", "FTL", "GPX4", "NCOA4", "TFRC")
PAIRED_FUNCTIONS <- c(FGFR3 = "GO:0008283", COMMD3 = "GO:0000075")   # proliferation, checkpoint
PCA_DONOR <- "GSE124263__A1"
PCA_MODE <- "Author max-abs"
KO_COLOURS <- c(FGFR3 = "#B65738", COMMD3 = "#367B88", IRF2BPL = "#7863A5", STX17 = "#6C7953",
                TCF3 = "#565656", Control = "#BAC3C7")
NEGATIVE <- "#367C9E"; POSITIVE <- "#BE6038"; UNAVAILABLE <- "#E9E9E9"
TEXT_PT <- 6.5
JITTER_SEED <- 1                  # point jitter in b
LABEL_SEEDS <- c(f = 1, h = 3)    # label repulsion in f and h
# Panel boxes (mm from the top-left corner of the page): f spans the height of g and h.
BOXES <- data.table(x = c(1.2, 1.2, 100.8, 1.2, 105.2, 1.2, 82.2, 82.2), y = c(1.2, 52.2, 52.2, 102.2, 102.2, 177.2, 177.2, 221.2),
                    w = c(180.6, 99.6, 81, 104, 76.6, 81, 99.6, 99.6), h = c(51, 50, 50, 75, 75, 91, 44, 47))
PAGE_HEIGHT <- 269.4

read_panel <- function(name) fread(file.path(PANELS, name), na.strings = "NA")
round_up <- function(x) ceiling(max(abs(x), na.rm = TRUE) * 10) / 10

## Devices: the published page measured its panels on three devices, kept here so the
## page is reproduced exactly. Panels are built while R's pdf() device is current (its
## fallback metrics for this font set the ggplot text extents; the font warnings are
## muffled), heatmaps are drawn on Cairo capture devices, and ComplexHeatmap's own
## scratch device is Cairo instead of pdf(NULL).
assignInNamespace("dev.null", function(...) cairo_pdf(tempfile(fileext = ".pdf"), family = FONT, ...), "ComplexHeatmap")
globalCallingHandlers(warning = function(w) if (grepl(FONT, conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
capture_device <- function(width, height, ...) cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height, family = FONT)
pdf(NULL)

theme_ko <- function() theme_spg() + theme(
  axis.text = element_text(size = TEXT_PT), axis.title = element_text(size = TEXT_PT),
  legend.text = element_text(size = TEXT_PT), legend.title = element_text(size = TEXT_PT),
  strip.text = element_text(size = TEXT_PT), plot.margin = margin(2, 2, 2, 2, "pt"))
ggplot_panel <- function(p) ggplotGrob(p + theme_ko())

## Heatmaps (a, c, d, e) ----------------------------------------------------
text_gp <- gpar(fontsize = TEXT_PT, fontfamily = FONT)
legend_style <- list(title_gp = gpar(fontfamily = FONT, fontsize = TEXT_PT, fontface = "plain"), labels_gp = text_gp,
                     grid_width = unit(2, "mm"), legend_height = unit(16, "mm"))
dot_legend <- function() Legend(labels = "Above all 5", type = "points", background = "white", pch = 16,
                                size = unit(.75, "mm"), labels_gp = text_gp, grid_height = unit(3, "mm"))
na_legend <- function() Legend(labels = "Unavailable", legend_gp = gpar(fill = UNAVAILABLE, col = NA), labels_gp = text_gp,
                               grid_height = unit(3, "mm"), grid_width = unit(2, "mm"))
scale_legend <- function(title, limit) c(list(title = title, at = c(-limit, 0, limit)), legend_style)
diverging <- function(limit) colorRamp2(c(-limit, 0, limit), c(NEGATIVE, "white", POSITIVE))

# Draws a heatmap with its legends on the right into a grob of w x h mm; `legends` is
# evaluated inside the capture so its text is measured there.
heatmap_panel <- function(ht, w, h, legends = list())
  grid.grabExpr(draw(ht, newpage = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right",
                     align_heatmap_legend = "heatmap_center", align_annotation_legend = "heatmap_center",
                     heatmap_legend_list = legends, padding = unit(c(1, 1, 1, 1), "mm")),
                width = w / 25.4, height = h / 25.4, device = capture_device)

# Target / donor columns of the focal heatmaps (c, d, e).
FOCAL_COLUMNS <- paste(rep(FOCUS, each = length(DONORS)), DONORS, sep = " / ")
focal_matrix <- function(z, key, keys, labels, values, empty = NA_real_) {
  m <- matrix(empty, length(keys), length(FOCAL_COLUMNS), dimnames = list(unname(labels), FOCAL_COLUMNS))
  m[cbind(match(z[[key]], keys), match(paste(z$target, z$donor, sep = " / "), FOCAL_COLUMNS))] <- values
  m
}
# hit: logical matrix of responses above all five controls, drawn as dots.
focal_heatmap <- function(m, name, colours, legend, hit = NULL)
  Heatmap(m, col = colours, name = name, na_col = UNAVAILABLE, cluster_rows = FALSE, cluster_columns = FALSE,
          rect_gp = gpar(col = "white", lwd = .6), column_split = factor(rep(FOCUS, each = length(DONORS)), levels = FOCUS),
          cluster_column_slices = FALSE, column_gap = unit(2, "mm"), column_labels = rep(DONORS, length(FOCUS)),
          column_names_rot = 0, column_names_gp = text_gp, column_title_gp = text_gp, row_names_side = "left",
          row_names_gp = text_gp, heatmap_legend_param = legend,
          cell_fun = if (!is.null(hit)) function(j, i, x, y, w, h, fill) if (hit[i, j]) grid.points(x, y, pch = 16, size = unit(.65, "mm")))
above_all <- function(z) z$comparison_eligible %in% TRUE & !is.na(z$delta_vs_control_max) & z$delta_vs_control_max > 0

members <- read_panel("MP3_member_responses.tsv")
member_order <- read_panel("MP3_member_order.tsv")
programs <- read_panel("focal_program_responses.tsv")[target %in% FOCUS]
functions <- read_panel("function_responses.tsv")[target %in% FOCUS]
excess_limit <- round_up(c(members$delta_vs_control_median, programs$delta_vs_control_median,
                           functions$delta_vs_control_median))
excess_legend <- scale_legend("Excess\nrank-biserial", excess_limit)

a_matrix <- matrix(NA_real_, length(DONORS), nrow(member_order), dimnames = list(DONORS, member_order$target))
a_matrix[cbind(match(members$donor, DONORS), match(members$target, member_order$target))] <- members$delta_vs_control_median
donor_bars <- HeatmapAnnotation(
  Donors = anno_barplot(member_order$n_estimable, ylim = c(0, 4), border = FALSE, gp = gpar(fill = "#718996", col = NA),
                        axis_param = list(at = c(0, 4), labels = c("0", "4"), gp = text_gp)),
  height = unit(6, "mm"), annotation_name_gp = text_gp, annotation_name_side = "left")
pa <- Heatmap(a_matrix, name = "Excess", col = diverging(excess_limit), na_col = UNAVAILABLE, cluster_rows = FALSE,
              cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = .3), row_names_side = "left",
              row_names_gp = text_gp, column_names_gp = text_gp, column_names_rot = 90,
              column_split = factor(ifelse(member_order$n_estimable == 4, "Complete", "Partial"), levels = c("Complete", "Partial")),
              cluster_column_slices = FALSE, column_title = NULL, column_gap = unit(2, "mm"), top_annotation = donor_bars,
              heatmap_legend_param = excess_legend)
grob_a <- heatmap_panel(pa, 180.6, 48, list(na_legend()))

c_matrix <- focal_matrix(programs, "pathway", names(PROGRAMS), PROGRAMS, programs$delta_vs_control_median)
c_hits <- focal_matrix(programs, "pathway", names(PROGRAMS), PROGRAMS, above_all(programs), empty = FALSE)
grob_c <- heatmap_panel(focal_heatmap(c_matrix, "Excess", diverging(excess_limit), excess_legend, c_hits), 81, 47,
                        c(list(dot_legend()), if (anyNA(c_matrix)) list(na_legend())))

stopifnot(setequal(functions$pathway, names(FUNCTIONS)))
d_matrix <- focal_matrix(functions, "pathway", names(FUNCTIONS), FUNCTIONS, functions$delta_vs_control_median)
d_hits <- focal_matrix(functions, "pathway", names(FUNCTIONS), FUNCTIONS, above_all(functions), empty = FALSE)
grob_d <- heatmap_panel(focal_heatmap(d_matrix, "Excess", diverging(excess_limit), excess_legend, d_hits), 104, 72,
                        c(list(dot_legend()), if (anyNA(d_matrix)) list(na_legend())))

genes <- read_panel("member_donor_responses.tsv")[target %in% FOCUS]
stopifnot(setequal(genes$gene, EXTERNAL_GENES))
e_matrix <- focal_matrix(genes, "gene", EXTERNAL_GENES, EXTERNAL_GENES, genes$delta_rank_percentile)
percentile_limit <- round_up(e_matrix)
grob_e <- heatmap_panel(focal_heatmap(e_matrix, "Percentile", diverging(percentile_limit),
                                      scale_legend("Percentile\nexcess", percentile_limit)),
                        76.6, 72, if (anyNA(e_matrix)) list(na_legend()))

## b: focal and matched-control MP3 responses -------------------------------
points <- read_panel("focal_and_control_responses.tsv")[target %in% FOCUS]
points[, `:=`(donor = factor(donor, DONORS), target = factor(target, FOCUS),
              kind = ifelse(probe_role == "matched_control", "Matched", "Candidate"))]
pb <- ggplot(points, aes(donor, rank_biserial, shape = kind, colour = kind)) +
  geom_point(position = position_jitter(width = .13, height = 0, seed = JITTER_SEED), size = 1.6, stroke = .35, fill = "white") +
  facet_grid(. ~ target) +
  scale_shape_manual("KO", values = c(Candidate = 18, Matched = 21)) +
  scale_colour_manual("KO", values = c(Candidate = "#B65738", Matched = "#7E8E95")) +
  labs(x = NULL, y = "MP3 rank-biserial")

## f: program-external gene responses in A1 ---------------------------------
responses <- read_panel("donor_gene_responses.tsv")[donor == "A1" & !is.na(FGFR3) & !is.na(COMMD3)]
iron <- responses[gene %in% IRON_GENES]
pf <- ggplot(responses, aes(FGFR3, COMMD3)) +
  geom_hline(yintercept = 0, colour = "#CDD3D6", linewidth = .2) +
  geom_vline(xintercept = 0, colour = "#CDD3D6", linewidth = .2) +
  geom_point(colour = "#728A94", alpha = .32, size = .36) +
  geom_point(data = iron, shape = 21, fill = "white", colour = "#B65738", size = 1.25) +
  geom_text_repel(data = iron, aes(label = gene), size = TEXT_PT / .pt, family = FONT, seed = LABEL_SEEDS[["f"]],
                  min.segment.length = 0, max.overlaps = Inf, box.padding = .2) +
  labs(x = "FGFR3 percentile excess", y = "COMMD3 percentile excess")

## g: same-node full-SPG versus U-only networks -----------------------------
paired <- read_panel("paired_function_responses.tsv")[pathway == PAIRED_FUNCTIONS[target]]
rows <- c("FGFR3 / A1", "FGFR3 / A2", "COMMD3 / A1", "COMMD3 / A2")
paired[, `:=`(condition = ifelse(config == "paired_full_seed1", "Full SPG", "U only"),
              row = factor(paste(target, donor, sep = " / "), rev(rows)))]
paired[, y := as.numeric(row) + ifelse(condition == "Full SPG", .13, -.13)]
pg <- ggplot(paired, aes(delta_vs_control_median, y, colour = condition, shape = condition)) +
  geom_vline(xintercept = 0, colour = "#BEBEBE", linewidth = .2) +
  geom_segment(aes(x = delta_vs_control_max, xend = delta_vs_control_min, yend = y), linewidth = .3) +
  geom_point(size = 1.4) +
  scale_y_continuous(breaks = 1:4, labels = rev(rows), expand = expansion(add = .4)) +
  scale_colour_manual("Network\nmin–max", values = c("Full SPG" = "#B65738", "U only" = "#367B88")) +
  scale_shape_manual("Network\nmin–max", values = c("Full SPG" = 16, "U only" = 17)) +
  labs(x = "Excess rank-biserial", y = NULL)

## h: response-profile PCA ---------------------------------------------------
pca <- read_panel("PCA_coordinates.tsv")[sample_key == PCA_DONOR & mode == PCA_MODE]
variance <- read_panel("PCA_explained_variance.tsv")[sample_key == PCA_DONOR & mode == PCA_MODE]
pca[, display := ifelse(role == "Control", "Control", target)]
ph <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(colour = display, shape = role), size = 1.8, stroke = .4) +
  geom_text_repel(data = pca[role != "Control"], aes(label = target, colour = display), size = TEXT_PT / .pt,
                  family = FONT, seed = LABEL_SEEDS[["h"]], min.segment.length = 0, max.overlaps = Inf, box.padding = .25,
                  show.legend = FALSE) +
  scale_colour_manual(values = KO_COLOURS, guide = "none") +
  scale_shape_manual("KO", values = c(Candidate = 18, Reference = 17, Control = 1)) +
  labs(x = sprintf("PC1 (%.1f%%)", variance[PC == 1, percent]), y = sprintf("PC2 (%.1f%%)", variance[PC == 2, percent]))

grobs <- list(grob_a, ggplot_panel(pb), grob_c, grob_d, grob_e, ggplot_panel(pf), ggplot_panel(pg), ggplot_panel(ph))
invisible(dev.off())

## Page: panels in fixed boxes, each below a 3-mm band holding its 8-pt letter -----
compose_boxes <- function(grobs, boxes, file, width = 183, height = PAGE_HEIGHT) {
  parts <- unlist(lapply(seq_along(grobs), function(k) { b <- boxes[k]; list(
    grobTree(grobs[[k]], vp = viewport(x = unit(b$x, "mm"), y = unit(height - b$y - 3, "mm"), width = unit(b$w, "mm"),
                                       height = unit(b$h - 3, "mm"), just = c("left", "top"), clip = "off")),
    textGrob(letters[k], x = unit(b$x + .2, "mm"), y = unit(height - b$y - .2, "mm"), just = c("left", "top"),
             gp = gpar(fontfamily = FONT, fontsize = 8, fontface = "bold"))) }), recursive = FALSE)
  page <- do.call(grobTree, parts)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(paste0(file, ".pdf"), width / 25.4, height / 25.4, family = FONT); grid.newpage(); grid.draw(page); invisible(dev.off())
  pdf_to_png(file)
}
compose_boxes(grobs, BOXES, file.path(FIGURES, "figure7"))
cat("Figure 7 written to", FIGURES, "\n")
