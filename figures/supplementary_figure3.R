# Supplementary Figure 3: MP3 response excess of the five focal knockouts over their matched-control
# median in the four donor networks, built under the primary, seed-17 and 80%-cell configurations.
# Input: results/knockout/panels/configuration_sensitivity.tsv.
# Output: results/figures/supplementary_figure3.pdf and supplementary_figure3.png.
source("utils/paths.R")
source("figures/theme.R")
suppressPackageStartupMessages(library(data.table))

CONFIGS <- c(primary = "Primary", seed17 = "Seed 17", down80_seed29 = "80% cells")
DONORS <- c("A1", "A2", "N2", "N3")
KO_COLOURS <- c(FGFR3 = "#B65738", COMMD3 = "#367B88", IRF2BPL = "#7863A5", STX17 = "#6C7953", TCF3 = "#565656")
WIDTH <- 183; ROW_HEIGHT <- 70; PAD <- 1.2   # mm; one row of a composed page
PAGE_HEIGHT <- ROW_HEIGHT + 2 * PAD

sens <- fread(file.path(RESULTS, "knockout", "panels", "configuration_sensitivity.tsv"), na.strings = "NA")
sens[, `:=`(config = factor(config, names(CONFIGS), CONFIGS), donor = factor(sub(".*__", "", sample_key), DONORS))]
# The four combinations without an eligible comparison (no excess) are left out of points and lines.
p <- ggplot(sens, aes(config, delta_vs_control_median, colour = target, group = target)) +
  geom_hline(yintercept = 0, colour = "#CCCCCC", linewidth = .2) +
  geom_line(linewidth = .3, na.rm = TRUE) + geom_point(size = 1.3, na.rm = TRUE) +
  facet_wrap(~donor, nrow = 1) +
  scale_colour_manual("Knockout", values = KO_COLOURS) +
  labs(x = "Network sampling", y = "MP3 excess over control median") +
  theme_spg() + theme(axis.text.x = element_text(angle = 35, hjust = 1))

# A single unlettered panel, placed as a row of a composed page (measured on a Cairo device,
# below the 3-mm band that would hold a panel letter).
cairo_pdf(tempfile(fileext = ".pdf"), WIDTH / 25.4, PAGE_HEIGHT / 25.4, family = FONT)
page <- grobTree(ggplotGrob(clean_panel(p)),
                 vp = viewport(x = unit(PAD, "mm"), y = unit(PAGE_HEIGHT - PAD - 3, "mm"), width = unit(WIDTH - 2 * PAD, "mm"),
                               height = unit(ROW_HEIGHT - 3, "mm"), just = c("left", "top"), clip = "off"))
invisible(dev.off())
out <- file.path(FIGURES, "supplementary_figure3")
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(paste0(out, ".pdf"), WIDTH / 25.4, PAGE_HEIGHT / 25.4, family = FONT); grid.newpage(); grid.draw(page); invisible(dev.off())
pdf_to_png(out)
cat("Supplementary Figure 3 written to", FIGURES, "\n")
