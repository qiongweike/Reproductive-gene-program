# Response-profile PCA per donor (Fig. 7h): the authors' negative-control expressions from scTenifoldKnk
# inst/manuscript/NEGCONTROL/negControls.R are evaluated on the response-Z matrix: max-absolute scaling of every
# knockout profile, then prcomp (centred, unscaled) of the knockouts.
# Inputs: results/knockout/panels/<key>_response_Z.tsv.gz, PCA_profile_eligibility.tsv
# Outputs: results/knockout/panels/PCA_coordinates.tsv, PCA_explained_variance.tsv
suppressPackageStartupMessages(library(data.table))
source("utils/paths.R")

MODE <- "Author max-abs"
panels <- out_dir("knockout", "panels")

author <- as.list(parse(file.path(VENDOR, "scTenifoldKnk", "inst", "manuscript", "NEGCONTROL", "negControls.R")))
author_expression <- function(pattern) {
  hit <- Filter(function(z) grepl(pattern, paste(deparse(z), collapse = " "), fixed = TRUE), author)
  stopifnot(length(hit) == 1)
  hit[[1]]
}
normalize_expr <- author_expression("apply(abs(tProfile), 2, max)")
pca_expr <- author_expression("PCA <- as.data.frame(prcomp(t(tProfile))$x)")

meta <- fread(file.path(panels, "PCA_profile_eligibility.tsv"))
coords <- list()
variance <- list()
for (key in unique(meta$sample_key)) {
  z <- as.data.frame(fread(file.path(panels, paste0(key, "_response_Z.tsv.gz"))))
  rownames(z) <- z$gene
  z$gene <- NULL
  en <- new.env()
  en$tProfile <- as.matrix(z)
  stopifnot(all(is.finite(en$tProfile)), all(apply(abs(en$tProfile), 2, max) > 0))
  eval(normalize_expr, en)
  eval(pca_expr, en)
  sdev <- prcomp(t(en$tProfile), center = TRUE, scale. = FALSE)$sdev
  pct <- 100 * sdev^2 / sum(sdev^2)
  mm <- as.data.frame(meta[meta$sample_key == key & meta$eligible, ])
  co <- data.frame(target = rownames(en$PCA), en$PCA[, 1:2], sample_key = key, mode = MODE)
  coords[[key]] <- merge(co, mm[, c("target", "donor", "role", "associated_candidates")], by = "target", sort = FALSE)
  variance[[key]] <- data.frame(sample_key = key, mode = MODE, PC = seq_along(pct), percent = pct)
  message(key, ": ", ncol(en$tProfile), " profiles, PC1/PC2 ", round(pct[1], 1), "/", round(pct[2], 1), "%")
}
fwrite(rbindlist(coords, fill = TRUE), file.path(panels, "PCA_coordinates.tsv"), sep = "\t")
fwrite(rbindlist(variance), file.path(panels, "PCA_explained_variance.tsv"), sep = "\t")
