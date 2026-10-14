# Figure 1: donor-resolved reference of human spermatogonial identity and state.
# Inputs: atlas embeddings, identity-marker summaries and withheld-marker agreement
# (results/atlas), the discovery-cell index (results/programs/nmf), fixed cell scores
# (results/trajectory) and paired SPG/SPC donor scores (results/coordination/specificity).
# Output: results/figures/figure1.pdf and figure1.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages({library(data.table); library(scales)})

ATLAS <- file.path(RESULTS, "atlas")
MIN_AGREEMENT_CELLS <- 30            # minimum cells behind a donor-cell-type point in e
JITTER_SEED <- 20261004
SEURAT_STROKE <- 0.600075815011372   # point stroke of Seurat::DimPlot, which drew the embeddings

CELL_TYPES <- c(SPG = "SPG", Spermatocyte = "Spermatocyte", Spermatid = "Spermatid", Sertoli = "Sertoli",
                Leydig = "Leydig", Peritubular_smooth_muscle = "Peritubular / muscle", Endothelial = "Endothelial",
                Immune = "Immune", Erythroid = "Erythroid candidate")
CELL_COLOURS <- c(SPG = "#376F91", Spermatocyte = "#69A9AF", Spermatid = "#89B67E", Sertoli = "#CE9A45",
                  Leydig = "#B45E65", Peritubular_smooth_muscle = "#8D7AAB", Endothelial = "#6B84B7",
                  Immune = "#8C785D", Erythroid = "#CE8982")
SHORT_TYPES <- c(Endothelial = "Endo.", Erythroid = "Eryth.", Immune = "Immune", Leydig = "Leydig",
                 Peritubular_smooth_muscle = "PTM", Sertoli = "Sertoli", SPG = "SPG", Spermatocyte = "SPC",
                 Spermatid = "Spermatid")
STATES <- c(undifferentiated = "U", transition = "T", differentiating = "D")
STATE_COLOURS <- c(undifferentiated = "#48799B", transition = "#D3A249", differentiating = "#459386")
STUDY_COLOURS <- c(GSE124263 = "#4477AA", GSE153947 = "#228833", GSE254315 = "#AA3377",
                   GSE202647 = "#B08A20", GSE326372 = "#CC6677")
MP_LABELS <- c(MP_1 = "MP1", MP_2 = "MP2", MP_3 = "MP3")

# Panels a-e keep the 6.5-pt tick labels and compact legends of the atlas panels.
theme_atlas <- function() theme_spg() + theme(
  axis.text = element_text(size = 6.5), panel.spacing = unit(1.4, "mm"),
  legend.key.height = unit(2.5, "mm"), legend.key.width = unit(2.5, "mm"),
  legend.spacing = unit(1.5, "mm"), legend.spacing.y = unit(1.5, "mm"))
no_axis_text <- theme(axis.text = element_blank(), axis.ticks = element_blank())
no_axis_lines <- theme(axis.line = element_blank(), axis.ticks = element_blank())

# UMAP coordinates in the drawing order of Seurat::DimPlot(shuffle = TRUE, seed = 1).
read_umap <- function(stem, columns) {
  u <- fread(file.path(ATLAS, paste0(stem, "_umap.tsv.gz")), select = c(columns, "UMAP_1", "UMAP_2"))
  set.seed(1)
  u[sample(nrow(u))]
}

umap_panel <- function(d, group, values, labels, name, size) {
  ggplot(d, aes(UMAP_1, UMAP_2, colour = .data[[group]])) +
    geom_point(size = size, stroke = SEURAT_STROKE) +
    scale_colour_manual(values = values, labels = labels, name = name) +
    coord_fixed() + theme_atlas() + no_axis_text +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.5)))
}

## a, b  reference embeddings (compartment-assigned and state-labelled cells)
atlas <- read_umap("adult_reference", "cell_type")[cell_type != "Uncertain"]
atlas[, cell_type := factor(cell_type, levels = names(CELL_TYPES))]
pa <- umap_panel(atlas, "cell_type", CELL_COLOURS, CELL_TYPES, "Cell type", .12)

spg <- read_umap("normal_spg", "state")[state != "unknown"]
spg[, state := factor(state, levels = names(STATES))]
pb <- umap_panel(spg, "state", STATE_COLOURS, STATES, "State", .35)

## c  state composition of the seven discovery donors
discovery <- fread(file.path(RESULTS, "programs/donor_manifest.tsv"))[eligible == TRUE, sample_key]
cells <- rbindlist(lapply(discovery, function(k)
  fread(file.path(RESULTS, "programs/nmf/input", paste0(k, ".cells.tsv")), select = c("sample_key", "study", "state"))))[
  state != "unknown"]
comp <- cells[, .(n_cells = .N), by = .(sample_key, study, state)]
comp[, `:=`(fraction = n_cells / sum(n_cells), donor = sub(".*__", "", sample_key)), by = sample_key]
comp[, state := factor(state, levels = names(STATES))]
pc <- ggplot(comp, aes(donor, fraction, fill = state)) + geom_col(width = .72) +
  scale_fill_manual(values = STATE_COLOURS, labels = STATES, name = "State") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "SPG fraction") + theme_atlas() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) + guides(fill = guide_legend(ncol = 1))

## d  identity markers: the first two measured assignment and withheld markers per cell type
manifest <- fread(file.path(ATLAS, "identity_marker_manifest.tsv"))
shown <- manifest[measured == TRUE & cell_type %in% names(CELL_TYPES), head(.SD, 2), by = .(cell_type, role)]
shown <- shown[order(match(cell_type, names(CELL_TYPES)), role)]
withheld <- manifest[role == "heldout", gene]
markers <- fread(file.path(ATLAS, "identity_marker_donor_equal_summary.tsv"))[
  gene %in% shown$gene & cell_type %in% names(CELL_TYPES)]
markers[, `:=`(gene = factor(gene, levels = unique(shown$gene)), cell_type = factor(cell_type, levels = names(CELL_TYPES)))]
pd <- ggplot(markers, aes(gene, cell_type)) +
  geom_point(aes(size = detection_fraction, colour = log1p(mean_CP10k)), stroke = 0) +
  scale_size_area(max_size = 3.5, limits = c(0, 1), breaks = c(.1, .5, .9), labels = percent_format(), name = "Detected") +
  scale_colour_gradient(low = "#E9EFF2", high = "#38708D", name = "log1p(mean\nCP10k)") +
  scale_x_discrete(labels = function(g) ifelse(g %in% withheld, paste0(g, "*"), g)) +
  scale_y_discrete(labels = CELL_TYPES) + labs(x = NULL, y = NULL) + theme_atlas() + no_axis_lines +
  theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = .5), axis.text.y = element_text(size = 7)) +
  guides(colour = colourbar_guide(), size = guide_legend(ncol = 1))

## e  withheld-marker agreement per donor and cell type
agree <- fread(file.path(ATLAS, "annotation_heldout_diagnostics.tsv"))[
  n_cells >= MIN_AGREEMENT_CELLS & atlas_cell_type %in% names(SHORT_TYPES)]
agree[, atlas_cell_type := factor(atlas_cell_type, levels = names(SHORT_TYPES))]
pe <- ggplot(agree, aes(atlas_cell_type, heldout_type_agreement, colour = study)) +
  geom_point(position = position_jitter(width = .16, height = 0, seed = JITTER_SEED), size = 1.4, alpha = .8) +
  scale_colour_manual(values = STUDY_COLOURS, name = "Study") + scale_x_discrete(labels = SHORT_TYPES) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(x = NULL, y = "Withheld-marker agreement") + theme_atlas() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) + guides(colour = guide_legend(ncol = 1))

## f  cell-cycle score by state, with donor-state medians
cycle <- fread(file.path(RESULTS, "trajectory/fixed_cell_scores.tsv.gz"))[state != "unknown"]
cycle[, state := factor(state, levels = names(STATES))]
medians <- cycle[, .(median_cycle = median(reference_cell_cycle)), by = .(sample_key, state, study)]
pf <- ggplot(cycle, aes(state, reference_cell_cycle, fill = state)) +
  geom_violin(scale = "width", linewidth = .3, alpha = .6) +
  geom_point(data = medians, aes(y = median_cycle, colour = study),
             position = position_jitter(width = .1, height = 0, seed = JITTER_SEED), size = 1.5) +
  scale_fill_manual(values = STATE_COLOURS, guide = "none") + scale_colour_manual(values = STUDY_COLOURS, name = "Study") +
  scale_x_discrete(labels = STATES) + labs(x = "State", y = "Cell-cycle score") + theme_spg()

## g  paired SPG and spermatocyte program scores; donors offset along the categorical axis
paired <- fread(file.path(RESULTS, "coordination/specificity/MP_paired_donor_values.tsv"))[eligible == TRUE]
paired[, stratum := ifelse(identity_stratum == "author_label_auxiliary", "Author labels", "Strict labels")]
paired <- melt(paired, id.vars = c("sample_key", "study", "MP", "stratum"), measure.vars = c("SPG_mean", "SPC_mean"),
               variable.name = "cell", value.name = "score")
donors <- sort(unique(paired$sample_key))
offset <- setNames(seq(-.15, .15, length.out = length(donors)), donors)
paired[, position := as.integer(cell) + offset[sample_key]]
pg <- ggplot(paired, aes(position, score, group = sample_key, colour = study, linetype = stratum)) +
  geom_line(linewidth = .4) + geom_point(size = 1.7) +
  facet_wrap(~MP, nrow = 1, labeller = as_labeller(MP_LABELS)) +
  scale_colour_manual(values = STUDY_COLOURS, name = "Study") +
  scale_x_continuous(breaks = 1:2, labels = c("SPG", "SPC"), limits = c(.73, 2.27)) +
  labs(x = NULL, y = "SPG-centred score", linetype = "Identity") + theme_spg()

## Page: Fig. 1 is not a row layout (b and c stack beside a), so panels are placed in
## boxes (mm from the top-left corner) with the conventions of compose_figure();
## the embedding boxes shrink to the fixed aspect ratio of their coordinates.
BOXES <- data.table(x = c(0, 105, 105, 0, 0, 109, 0), y = c(0, 0, 26, 64, 112, 112, 164),
                    w = c(104, 78, 78, 183, 108, 74, 183), h = c(63, 26, 37, 47, 51, 51, 44),
                    fixed_aspect = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))

compose_boxes <- function(panels, boxes, file, width = 183, height = 210.4) {
  pad <- 1.2; body <- width - 2 * pad
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  boxes <- copy(boxes)[, `:=`(x = pad + x * body / width, w = w * body / width, y = y + pad)]
  metrics <- paste0(file, "_metrics.pdf")
  cairo_pdf(metrics, width / 25.4, height / 25.4, family = FONT)
  grobs <- lapply(panels, function(p) ggplotGrob(clean_panel(p)))
  for (k in which(boxes$fixed_aspect)) {
    pp <- ggplot_build(panels[[k]])$layout$panel_params[[1]]
    fixed_w <- convertWidth(sum(grobs[[k]]$widths), "mm", valueOnly = TRUE)
    fixed_h <- convertHeight(sum(grobs[[k]]$heights), "mm", valueOnly = TRUE)
    boxes[k, h := min(h, 3 + fixed_h + (w - fixed_w) * diff(pp$y.range) / diff(pp$x.range))]
  }
  dev.off(); unlink(metrics)
  parts <- unlist(lapply(seq_along(panels), function(k) { b <- boxes[k]; list(
    grobTree(grobs[[k]], vp = viewport(x = unit(b$x, "mm"), y = unit(height - b$y - 3, "mm"), width = unit(b$w, "mm"),
                                       height = unit(b$h - 3, "mm"), just = c("left", "top"), clip = "off")),
    textGrob(letters[k], x = unit(b$x + .3, "mm"), y = unit(height - b$y - .3, "mm"), just = c("left", "top"),
             gp = gpar(fontfamily = FONT, fontsize = 9, fontface = "bold", col = INK))) }), recursive = FALSE)
  page <- do.call(grobTree, parts)
  cairo_pdf(paste0(file, ".pdf"), width / 25.4, height / 25.4, family = FONT); grid.newpage(); grid.draw(page); dev.off()
  pdf_to_png(file)
  invisible(page)
}

compose_boxes(list(pa, pb, pc, pd, pe, pf, pg), BOXES, file.path(FIGURES, "figure1"))
cat("Figure 1 written:", nrow(atlas), "atlas cells,", nrow(spg), "SPG,", nrow(cycle), "discovery cells\n")
