# GO and Reactome gene membership with upper-case symbols, for the enrichment tests (08) and knockout readouts.
# GO: org.Hs.eg.db GOALL (BP, CC and MF with ancestor propagation, no annotation-code filter), names from GO.db.
# Reactome: human (R-HSA) pathways of the downloaded ReactomePathways.gmt.
# Inputs: org.Hs.eg.db 3.20.0, GO.db 3.20.0, data/external/Reactome/ReactomePathways.gmt
# Output: results/coordination/pathways/membership.tsv.gz (pathway, gene, database, name)
source("utils/paths.R")
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(org.Hs.eg.db); library(GO.db)})

go <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                          keytype = "ENTREZID", columns = c("SYMBOL", "GOALL", "ONTOLOGYALL")))
go <- unique(go[!is.na(GOALL) & !is.na(SYMBOL),
                .(pathway = GOALL, gene = toupper(SYMBOL), database = paste0("GO_", ONTOLOGYALL))])
term <- vapply(as.list(GO.db::GOTERM), AnnotationDbi::Term, character(1))
go[, name := unname(term[pathway])]

gmt <- strsplit(readLines(file.path(EXTERNAL, "Reactome", "ReactomePathways.gmt")), "\t")
reactome <- rbindlist(lapply(gmt, function(x) if (startsWith(x[2], "R-HSA-"))
  data.table(pathway = x[2], gene = toupper(x[-c(1, 2)]), database = "Reactome", name = x[1])))

membership <- unique(rbind(go, reactome))
fwrite(membership, file.path(out_dir("coordination", "pathways"), "membership.tsv.gz"), sep = "\t")
cat(sprintf("%d pathways; %d membership rows (org.Hs.eg.db %s, GO.db %s)\n", uniqueN(membership$pathway),
            nrow(membership), packageVersion("org.Hs.eg.db"), packageVersion("GO.db")))
