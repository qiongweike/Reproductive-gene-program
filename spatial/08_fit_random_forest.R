# Buffered neighbourhood random forests (ranger) for the measured MP scores, with original or MP-excluded RCTD weights.
# Beads: not rejected by either RCTD fit. Each geometric quadrant is held out in turn; training beads within 5 median-NN
# units of it are dropped. Baseline: log UMI, 8 all-type weights, two kernel masses, x, y (13 features); augmented adds
# Gaussian-weighted neighbour compositions over the 50 nearest beads at 2 and 5 median-NN units (MISTy weight function).
# Inputs: results/spatial/program_scores/<array>_model_input.tsv.gz, neighbourhood/<array>_{classes,weights_alltypes}.tsv.gz.
# Output: results/spatial/neighbourhood/<array>_buffered_RF_predictions.tsv.gz (held-out bead predictions).
source("utils/paths.R")
suppressPackageStartupMessages({library(data.table); library(ranger); library(RANN); library(mistyR)})

TYPES <- c("Endothelial", "Immune", "Leydig", "Peritubular_smooth_muscle", "Sertoli", "Spermatid", "Spermatocyte", "SPG")
TARGETS <- c("MP_1", "MP_2", "MP_3")
SCALES <- c(2, 5)          # Gaussian kernel scales in median-NN units
K_VIEW <- 50
BUFFER <- 5
NUM_TREES <- 100
SEED <- 20261003           # ranger seed = SEED + held-out fold
COVARIATES <- c("log_UMI", TYPES, "kernel_mass_2", "kernel_mass_5", "x_native", "y_native")
VIEWS <- c(paste0("near_", TYPES), paste0("broad_", TYPES))
gaussian_weight <- getFromNamespace("get_weight", "mistyR")
OUT <- out_dir("spatial", "neighbourhood")

for (array in c("Puck5", "Puck6")) {
  model_input <- fread(file.path(RESULTS, "spatial", "program_scores", paste0(array, "_model_input.tsv.gz")))
  classes <- fread(file.path(OUT, paste0(array, "_classes.tsv.gz")))
  excluded_weights <- fread(file.path(OUT, paste0(array, "_weights_alltypes.tsv.gz")))
  ids <- model_input$barcode[model_input$barcode %in% classes[spot_class != "reject", barcode]]
  base <- as.data.frame(model_input[match(ids, barcode)])
  excluded_weights <- as.data.frame(excluded_weights[match(ids, barcode)])
  stopifnot(identical(base$barcode, excluded_weights$barcode))
  position <- base[, c("x_native", "y_native")]
  mnn <- unique(model_input$NN_median_native)
  stopifnot(length(mnn) == 1)
  kn <- RANN::nn2(position, k = K_VIEW + 1, eps = 0)
  stopifnot(all(kn$nn.idx[, 1] == seq_len(nrow(base))))
  neighbours <- kn$nn.idx[, -1]; distances <- kn$nn.dists[, -1]
  predictions <- list()
  for (version in c("original_RCTD", "MP_excluded_RCTD")) {
    d <- base
    if (version == "MP_excluded_RCTD") d[, TYPES] <- excluded_weights[, TYPES]
    for (scale in SCALES) {
      w <- gaussian_weight("gaussian", distances, scale * mnn, 0)
      mass <- rowSums(w)
      d[[paste0("kernel_mass_", scale)]] <- mass
      for (t in TYPES) d[[paste0(if (scale == 2) "near_" else "broad_", t)]] <-
        rowSums(matrix(d[[t]][neighbours], nrow = nrow(d)) * w) / mass
    }
    for (target in TARGETS) for (fold in 0:3) {
      test <- which(d$quadrant == fold); train <- which(d$quadrant != fold)
      train <- train[RANN::nn2(data = position[test, ], query = position[train, ], k = 1)$nn.dists[, 1] > BUFFER * mnn]
      stopifnot(length(train) > 1000, length(test) > 1000)
      for (model in c("composition_depth_coordinates", "plus_neighbour_composition")) {
        cols <- if (model == "composition_depth_coordinates") COVARIATES else c(COVARIATES, VIEWS)
        fit <- ranger(x = d[train, cols], y = d[train, target], num.trees = NUM_TREES, importance = "none",
                      num.threads = 1, seed = SEED + fold)
        predictions[[length(predictions) + 1]] <- data.table(
          barcode = d$barcode[test], array = array, identity_version = version, program = target, fold = fold, model = model,
          observed = d[test, target], prediction = predict(fit, data = d[test, cols])$predictions)
      }
    }
  }
  fwrite(rbindlist(predictions), file.path(OUT, paste0(array, "_buffered_RF_predictions.tsv.gz")), sep = "\t")
  message(array, ": ", length(predictions), " models on ", nrow(base), " beads")
}
