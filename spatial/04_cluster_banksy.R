# BANKSY RNA clusters (legacy BANKSY 0.1.6 API) at lambda 0 and 0.2; Fig 6a,b use lambda 0.2.
# Features: 2,000 Seurat VST genes detected in >=20 beads, excluding all meta-program members and the held-out
# interpretation markers; expression is counts / full bead UMI * median UMI.
# Arrays run in one process in fixed order: the PCA (irlba) draws from the R random stream seeded once.
# Inputs: data/processed/atlas/spatial/<array>_{counts.mtx.gz,genes.tsv.gz,coords.tsv.gz}, results/programs/nmf/MP_genes.tsv.
# Outputs: results/spatial/domains/<array>_Banksy_domains.tsv.gz.
source("utils/paths.R")
suppressPackageStartupMessages({library(Banksy); library(Seurat); library(Matrix); library(data.table)})

SEED <- 20261003
N_FEATURES <- 2000
MIN_BEADS <- 20
K_GEOM <- c(15, 30)
LAMBDA <- c(0, 0.2)
NPCS <- 20
K_EXPRESSION <- 50
RESOLUTION <- 0.5
MARKERS <- list(SPG = c("UTF1", "FGFR3", "PIWIL4", "GFRA1", "ZBTB16", "UCHL1"), Sertoli = c("SOX9", "AMH", "WT1", "CLU"),
                PTM = c("ACTA2", "TAGLN", "MYH11"), Leydig = c("STAR", "CYP17A1", "INSL3"),
                Endothelial = c("VWF", "PECAM1", "KDR"), Immune = c("PTPRC", "CD3D", "LST1"),
                Spermatocyte = c("SYCP3", "DMC1", "SPO11"), Spermatid = c("TNP1", "TNP2", "PRM1", "PRM2"))
IN <- file.path(PROCESSED, "atlas", "spatial")
OUT <- out_dir("spatial", "domains")

set.seed(SEED)
exclude <- toupper(unique(c(fread(file.path(RESULTS, "programs", "nmf", "MP_genes.tsv"))$gene, unlist(MARKERS))))
for (array in c("Puck5", "Puck6")) {
  x <- as(readMM(gzfile(file.path(IN, paste0(array, "_counts.mtx.gz")))), "dgCMatrix")
  genes <- fread(file.path(IN, paste0(array, "_genes.tsv.gz")))$gene
  meta <- as.data.frame(fread(file.path(IN, paste0(array, "_coords.tsv.gz"))))
  dimnames(x) <- list(genes, meta$barcode)
  total <- colSums(x)
  stopifnot(all(total == meta$total_counts))
  eligible <- !(toupper(genes) %in% exclude) & Matrix::rowSums(x > 0) >= MIN_BEADS
  ids <- sprintf("gene%06d", seq_len(sum(eligible)))   # syntactic names for Seurat; mapped back below
  xx <- x[eligible, ]; rownames(xx) <- ids
  seu <- FindVariableFeatures(CreateSeuratObject(counts = xx, min.cells = 0, min.features = 0),
                              selection.method = "vst", nfeatures = N_FEATURES, verbose = FALSE)
  selected <- genes[eligible][match(VariableFeatures(seu), ids)]
  stopifnot(length(selected) == N_FEATURES, !anyNA(selected))
  vals <- sweep(as.matrix(x[selected, ]), 2, total, "/") * median(total)
  coords <- meta[, c("x_native", "y_native")]; rownames(coords) <- meta$barcode
  bank <- BanksyObject(own.expr = vals, cell.locs = coords, meta.data = data.frame(row.names = meta$barcode))
  bank <- ComputeBanksy(bank, k_geom = K_GEOM, M = 1, verbose = TRUE)
  bank <- ScaleBanksy(bank)
  bank <- RunBanksyPCA(bank, lambda = LAMBDA, npcs = NPCS, M = 1)
  bank <- ClusterBanksy(bank, lambda = LAMBDA, M = 1, npcs = NPCS, k.neighbors = K_EXPRESSION,
                        resolution = RESOLUTION, seed = SEED, num.cores = 1)
  cl <- as.data.frame(bank@meta.data)
  fwrite(data.table(barcode = colnames(bank@own.expr), cl), file.path(OUT, paste0(array, "_Banksy_domains.tsv.gz")), sep = "\t")
  message(array, ": ", paste(sapply(cl[grep("^clust", names(cl))], function(v) length(unique(v))), collapse = "/"),
          " clusters at lambda ", paste(LAMBDA, collapse = "/"))
  rm(x, xx, seu, bank, vals); gc()
}
