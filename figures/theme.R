# Shared figure theme and page composition: 7-pt Liberation Sans, 183-mm pages,
# panels placed in rows of fixed height with bold panel letters.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(gtable)
})

FONT <- "Liberation Sans"
INK <- "#293942"

theme_spg <- function() {
  theme_classic(base_size = 7, base_family = FONT) + theme(
    text = element_text(colour = INK), axis.text = element_text(size = 7, colour = INK),
    axis.title = element_text(size = 7), axis.line = element_line(linewidth = .25),
    axis.ticks = element_line(linewidth = .25), strip.background = element_blank(),
    strip.text = element_text(size = 7), panel.spacing = unit(1.5, "mm"),
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt"), legend.position = "right", legend.justification = "center",
    legend.box.just = "center", legend.box = "vertical", legend.direction = "vertical",
    legend.title = element_text(size = 7), legend.text = element_text(size = 7),
    legend.key.height = unit(2.6, "mm"), legend.key.width = unit(2.6, "mm"),
    legend.box.spacing = unit(1.3, "mm"), legend.spacing.y = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0))
}

# Vertical colour bar used for continuous scales.
colourbar_guide <- function(order = 1) {
  guide_colourbar(direction = "vertical", barheight = unit(13, "mm"), barwidth = unit(2, "mm"),
                  title.position = "top", order = order)
}

# Final per-panel settings applied when a page is composed.
clean_panel <- function(p) {
  p <- p + labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL) + theme(
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    plot.tag = element_blank(), text = element_text(family = FONT), legend.position = "right",
    legend.justification = "center", legend.box.just = "center", legend.box = "vertical",
    legend.direction = "vertical", legend.title = element_text(size = 7),
    legend.text = element_text(size = 7), legend.box.spacing = unit(1.3, "mm"),
    plot.margin = margin(2, 2, 2, 2, "pt"))
  for (a in c("colour", "fill", "size", "shape", "linetype", "alpha")) {
    s <- p$scales$get_scales(a)
    if (!is.null(s) && is.character(s$name)) s$name <- sub(" \\([a-z,]+\\)$", "", s$name)
    if (is.character(p$labels[[a]])) p$labels[[a]] <- sub(" \\([a-z,]+\\)$", "", p$labels[[a]])
  }
  p
}

# Write <file>.png by rasterizing <file>.pdf at 300 dpi, as for the published figures.
pdf_to_png <- function(file) {
  bitmap <- pdftools::pdf_render_page(paste0(file, ".pdf"), page = 1, dpi = 300, numeric = TRUE)
  png::writePNG(bitmap, paste0(file, ".png"), dpi = 300)
  invisible(paste0(file, ".png"))
}

# Compose panels into one page. `rows` lists panel indices per row, `heights`
# gives row heights in mm, `widths` relative widths within each row. Panel k
# is lettered letters[k] unless tags = FALSE. Writes <file>.pdf and <file>.png.
compose_figure <- function(panels, rows, heights, file, widths = NULL, width = 183, tags = TRUE) {
  stopifnot(sum(lengths(rows)) == length(panels), length(rows) == length(heights))
  pad <- 1.2; body <- width - 2 * pad; height <- sum(heights) + 2 * pad
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  metrics <- paste0(file, "_metrics.pdf")
  cairo_pdf(metrics, width / 25.4, height / 25.4, family = FONT)   # text metrics for grob layout
  parts <- list(); top <- pad
  for (i in seq_along(rows)) {
    w <- if (is.null(widths)) rep(1, length(rows[[i]])) else widths[[i]]
    w <- body * w / sum(w); left <- pad
    for (j in seq_along(rows[[i]])) {
      k <- rows[[i]][j]; p <- clean_panel(panels[[k]])
      g <- if (inherits(p, "patchwork")) patchworkGrob(p + plot_layout(guides = "keep")) else ggplotGrob(p)
      parts <- c(parts, list(
        grobTree(g, vp = viewport(x = unit(left, "mm"), y = unit(height - top - 3, "mm"),
                                  width = unit(w[j], "mm"), height = unit(heights[i] - 3, "mm"),
                                  just = c("left", "top"), clip = "off"))))
      if (tags) parts <- c(parts, list(
        textGrob(letters[k], x = unit(left + .3, "mm"), y = unit(height - top - .2, "mm"),
                 just = c("left", "top"), gp = gpar(fontfamily = FONT, fontsize = 9, fontface = "bold", col = INK))))
      left <- left + w[j]
    }
    top <- top + heights[i]
  }
  dev.off(); unlink(metrics)
  page <- do.call(grobTree, parts)
  cairo_pdf(paste0(file, ".pdf"), width / 25.4, height / 25.4, family = FONT); grid.newpage(); grid.draw(page); dev.off()
  pdf_to_png(file)
  invisible(page)
}
