# Export the CellChatDB.human protein-signalling subset used by every per-donor CellChat model:
# subsetDB() defaults keep Secreted Signaling, ECM-Receptor and Cell-Cell Contact and drop non-protein signalling.
# Outputs: results/communication/cellchat_inputs/{protein_database.rds, database_interactions.tsv,
#          database_genes.tsv}
source("utils/paths.R")
suppressPackageStartupMessages(library(CellChat))

out <- out_dir("communication", "cellchat_inputs")
db <- subsetDB(CellChatDB.human)
saveRDS(db, file.path(out, "protein_database.rds"))
write.table(data.frame(gene = extractGene(db)), file.path(out, "database_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(db$interaction, file.path(out, "database_interactions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("CellChat protein database:", nrow(db$interaction), "interactions,", length(extractGene(db)), "genes\n")
