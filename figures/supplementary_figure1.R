# Supplementary Figure 1: reference composition and robustness of MP3 coordination.
# Inputs: atlas embeddings (results/atlas) and donor-level changes in MP3 residual
# coordination after adding transcript covariates (results/coordination/specificity).
# Output: results/figures/supplementary_figure1.pdf and supplementary_figure1.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages(library(data.table))

STUDY_COLOURS <- c(GSE124263 = "#4477AA", GSE153947 = "#228833",   # also the legend order in a, b
                   GSE254315 = "#AA3377", GSE202647 = "#B08A20")
COVARIATES <- c(plus_late = "Late RNA", plus_somatic = "Somatic RNA", plus_late_somatic = "Both")

# The study-coloured embeddings keep their 6.5-pt base theme without axes.
theme_embedding <- function() theme_classic(base_size = 6.5, base_family = FONT) + theme(
  text = element_text(colour = "#25303B"), axis.title = element_text(size = 7),
  axis.title.x = element_text(margin = margin(t = 2)), axis.title.y = element_text(margin = margin(r = 2)),
  axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
  axis.ticks.length = unit(.8, "mm"), legend.key.size = unit(2.5, "mm"), legend.spacing = unit(.5, "mm"),
  legend.margin = margin(0, 0, 0, 0))

# UMAP coordinates in the drawing order of Seurat::DimPlot(shuffle = TRUE, seed = 1).
read_umap <- function(stem, columns) {
  u <- fread(file.path(RESULTS, "atlas", paste0(stem, "_umap.tsv.gz")), select = c(columns, "UMAP_1", "UMAP_2"))
  set.seed(1)
  u[sample(nrow(u))]
}

study_panel <- function(d) {
  ggplot(d, aes(UMAP_1, UMAP_2, colour = study)) + geom_point(size = .13, alpha = .65, stroke = 0) +
    scale_colour_manual(values = STUDY_COLOURS, limits = intersect(names(STUDY_COLOURS), d$study), name = NULL) +
    coord_fixed() + labs(x = "UMAP 1", y = "UMAP 2") + theme_embedding() +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.6, alpha = 1)))
}

## a, b  the cells of Fig. 1a,b coloured by source study
pa <- study_panel(read_umap("adult_reference", c("study", "cell_type"))[cell_type != "Uncertain"])
pb <- study_panel(read_umap("normal_spg", c("study", "state"))[state != "unknown"])

## c  change in donor-mean MP3 residual coordination with added transcript covariates
delta <- fread(file.path(RESULTS, "coordination/specificity/background_donor_differences.tsv"))[
  MP == "MP_3" & model %in% names(COVARIATES)]
delta[, `:=`(model = factor(model, levels = names(COVARIATES)), study = sub("__.*", "", sample_key))]
pc <- ggplot(delta, aes(model, difference, group = sample_key, colour = study)) +
  geom_hline(yintercept = 0, colour = "#A3AAAE", linewidth = .3) +
  geom_line(linewidth = .35, alpha = .5) + geom_point(size = 1.6) +
  scale_colour_manual(values = STUDY_COLOURS, name = "Study") + scale_x_discrete(labels = COVARIATES) +
  labs(x = "Additional transcript covariates", y = "Δ MP3 residual ρ") + theme_spg() +
  theme(axis.text = element_text(size = 6.5)) + guides(colour = guide_legend(ncol = 1))

compose_figure(list(pa, pb, pc), rows = list(1:2, 3), heights = c(66, 52),
               file = file.path(FIGURES, "supplementary_figure1"))
cat("Supplementary Figure 1 written:", nrow(pa$data), "atlas cells,", nrow(pb$data), "SPG,",
    nrow(delta), "donor contrasts\n")
