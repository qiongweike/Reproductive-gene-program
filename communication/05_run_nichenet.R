# NicheNet ligand activities (v2 human prior; predict_ligand_activities(single = TRUE), ranked by Pearson) for the fixed
# MP3 targets (MP2 comparator) in U/T/D spermatogonia per donor. Prior names map to measured genes via Ensembl/HGNC
# identity (ambiguous ones dropped); >=10% detection for receiver genes, sender ligands and receptors; >=15 targets.
# Sensitivities: the three U contexts without ID1/ID2/ID4 or without published state markers.
# Inputs: step 04 outputs, results/programs/nmf/MP_genes.tsv, data/external/ (NicheNet v2 prior, HGNC, state markers),
#         VENDOR/nichenetr/data/geneinfo_alias_human.rda
# Outputs: results/communication/nichenet/ (ligand activities, sensitivities, U ligand selection, mapping tables)
source("utils/paths.R")
suppressPackageStartupMessages({library(nichenetr); library(data.table)})

MIN_DETECTION <- 0.1
MIN_TARGETS <- 15
WORKERS <- 4L
STATES <- c("SPG_U", "SPG_T", "SPG_D")
PROGRAMS <- c("MP_3", "MP_2")
ID_GENES <- c("ID1", "ID2", "ID4")
IN <- file.path(RESULTS, "communication", "nichenet_inputs")
O <- out_dir("communication", "nichenet")
NN <- file.path(EXTERNAL, "nichenet")
save_tsv <- function(z, fname, ...) fwrite(z, file.path(O, fname), sep = "\t", ...)

lt <- readRDS(file.path(NN, "ligand_target_matrix_nsga2r_final.rds"))
lr <- unique(as.data.table(readRDS(file.path(NN, "lr_network_human_21122021.rds")))[, .(from, to)])

# ---- Assay identities: Ensembl -> unique approved HGNC; unmatched old Ensembl IDs keep an unambiguous HGNC symbol ----
load(file.path(VENDOR, "nichenetr", "data", "geneinfo_alias_human.rda"))
alias <- unique(as.data.table(geneinfo_alias_human)[!is.na(entrez), .(symbol, entrez = as.character(entrez), alias)])
ae <- alias[, .(n_id = uniqueN(entrez), entrez = first(entrez)), by = alias][n_id == 1]
v <- fread(file.path(O, "assay_gene_identity.tsv"))
v[, assay_order := .I]
v <- merge(v, ae[, .(alias, entrez)], by.x = "assay_symbol", by.y = "alias", all.x = TRUE, sort = FALSE)
setorder(v, assay_order)
stopifnot(!anyDuplicated(v$ensembl), !anyDuplicated(v$assay_symbol))
hgnc <- fread(file.path(EXTERNAL, "hgnc_complete_set.txt"), na.strings = "", colClasses = "character")[status == "Approved"]
hu <- hgnc[!is.na(ensembl_gene_id)][!duplicated(ensembl_gene_id) & !duplicated(ensembl_gene_id, fromLast = TRUE)]
v[, `:=`(hgnc_id = hu$hgnc_id[match(ensembl, hu$ensembl_gene_id)],
         approved_symbol = hu$symbol[match(ensembl, hu$ensembl_gene_id)])]
ha <- rbindlist(lapply(c("symbol", "prev_symbol", "alias_symbol"), function(field) {
  zz <- hgnc[!is.na(get(field)), .(hgnc_id, names = get(field))]
  zz <- zz[, .(symbol = unlist(strsplit(names, "|", fixed = TRUE))), by = hgnc_id]
  zz[, priority := match(field, c("symbol", "prev_symbol", "alias_symbol"))]
  zz
}))
ha <- ha[ha[, .I[priority == min(priority)], by = symbol]$V1]  # an approved symbol outranks an alias of another gene
ha <- ha[, .(n_id = uniqueN(hgnc_id), hgnc_id = first(hgnc_id)), by = symbol][n_id == 1]
vi <- which(is.na(v$hgnc_id))
hi <- ha$hgnc_id[match(v$assay_symbol[vi], ha$symbol)]
ok <- !is.na(hi) & !hi %in% v$hgnc_id & !duplicated(hi) & !duplicated(hi, fromLast = TRUE)
v[vi[ok], `:=`(hgnc_id = hi[ok], approved_symbol = hgnc$symbol[match(hi[ok], hgnc$hgnc_id)])]
stopifnot(!anyDuplicated(v$hgnc_id[!is.na(v$hgnc_id)]))
vh <- v[!is.na(hgnc_id)][!duplicated(hgnc_id) & !duplicated(hgnc_id, fromLast = TRUE)]
save_tsv(v, "assay_HGNC_identity.tsv", na = "NA")
unique_entrez <- v[!is.na(entrez), .(n_assay = .N, ensembl = first(ensembl)), by = entrez][n_assay == 1]

# Map names to measured assay genes. Exact symbol first, then case restore, HGNC approved/previous/alias symbols and
# NicheNet aliases resolving to one measured Ensembl ID; a measured HGNC identity overrides legacy names, and
# symbol-versus-Ensembl HGNC conflicts stay unmapped.
resolve <- function(names) {
  z <- data.table(symbol = unique(names))
  z[, assay_i := match(symbol, v$assay_symbol)]
  z[, mapping := "exact_assay_symbol"]
  ii <- which(is.na(z$assay_i)); ci <- match(toupper(z$symbol[ii]), v$symbol_upper); good <- !is.na(ci)
  z[ii[good], `:=`(assay_i = ci[good], mapping = "unique_case_restore")]
  ii <- which(is.na(z$assay_i)); hid <- ha$hgnc_id[match(z$symbol[ii], ha$symbol)]
  ci <- match(vh$ensembl[match(hid, vh$hgnc_id)], v$ensembl); good <- !is.na(ci)
  z[ii[good], `:=`(assay_i = ci[good], mapping = "HGNC_approved_previous_alias_to_assay_Ensembl")]
  ii <- which(is.na(z$assay_i)); ent <- ae$entrez[match(z$symbol[ii], ae$alias)]
  ci <- match(unique_entrez$ensembl[match(ent, unique_entrez$entrez)], v$ensembl); good <- !is.na(ci)
  z[ii[good], `:=`(assay_i = ci[good], mapping = "author_alias_unique_assay_Ensembl")]
  z[symbol == "RP11-325I22.2", `:=`(assay_i = match("ENSG00000250584", v$ensembl), mapping = "HGNC51200_ENSG00000250584")]
  query_hgnc <- ha$hgnc_id[match(z$symbol, ha$symbol)]
  preferred_i <- match(vh$ensembl[match(query_hgnc, vh$hgnc_id)], v$ensembl)
  remap <- which(!is.na(preferred_i) & (is.na(z$assay_i) | preferred_i != z$assay_i))
  z[remap, `:=`(assay_i = preferred_i[remap], mapping = "HGNC_measured_Ensembl_precedes_legacy_assay_name")]
  assay_hgnc <- v$hgnc_id[z$assay_i]
  conflict <- which(!is.na(query_hgnc) & !is.na(assay_hgnc) & query_hgnc != assay_hgnc)
  z[conflict, `:=`(assay_i = NA_integer_, mapping = "HGNC_symbol_vs_assay_Ensembl_conflict")]
  z[is.na(assay_i) & mapping != "HGNC_symbol_vs_assay_Ensembl_conflict", mapping := "unresolved"]
  z[, `:=`(ensembl = v$ensembl[assay_i], assay_symbol = v$assay_symbol[assay_i], hgnc_id = v$hgnc_id[assay_i],
           approved_symbol = v$approved_symbol[assay_i], query_hgnc_id = query_hgnc)]
  z
}

# One prior row per measured gene: exact symbol > case variant > current approved symbol > other alias; ties dropped.
select_prior <- function(names, kind) {
  z <- resolve(names)
  z[, kind := kind]
  z[, priority := fifelse(symbol == assay_symbol, 0L, fifelse(toupper(symbol) == toupper(assay_symbol), 1L, 3L))]
  z[priority == 3L & symbol == approved_symbol, priority := 2L]
  z[, selected := FALSE]
  z[, selection_reason := fifelse(is.na(ensembl), "unresolved_identity", "alternate_same_entity")]
  best <- z[!is.na(ensembl), {ix <- .I[priority == min(priority)]; .(index = ix, n_best = length(ix))}, by = ensembl]
  z[best[n_best == 1]$index, `:=`(selected = TRUE, selection_reason = "unique_preferred_prior_row")]
  z[best[n_best > 1]$index, selection_reason := "ambiguous_prior_aliases_not_selected"]
  z
}

prior_targets <- select_prior(rownames(lt), "target")
prior_ligands <- select_prior(colnames(lt), "ligand")
save_tsv(rbindlist(list(prior_targets, prior_ligands)), "prior_identity_mapping.tsv", na = "NA")
mapt <- prior_targets[selected == TRUE]
mapl <- prior_ligands[selected == TRUE]
stopifnot(!anyDuplicated(mapt$hgnc_id[!is.na(mapt$hgnc_id)]), !anyDuplicated(mapl$hgnc_id[!is.na(mapl$hgnc_id)]))

expr <- fread(file.path(IN, "donor_group_expression.tsv.gz"))
ix <- match(expr$gene, v$symbol_upper)
stopifnot(!anyNA(ix))
expr[, `:=`(legacy_gene = gene, gene = v$assay_symbol[ix], ensembl = v$ensembl[ix])]
mp <- fread(file.path(RESULTS, "programs", "nmf", "MP_genes.tsv"))
mp <- merge(mp, resolve(mp$gene), by.x = "gene", by.y = "symbol", all.x = TRUE, sort = FALSE)
mp[, prior_target := mapt$symbol[match(ensembl, mapt$ensembl)]]
save_tsv(mp, "fixed_program_entity_mapping.tsv", na = "NA")
stopifnot(uniqueN(mp[MP == "MP_3"]$ensembl) == 49, !anyNA(mp[MP == "MP_3"]$ensembl))
lr_id <- resolve(unique(c(lr$from, lr$to)))
lr[, `:=`(ligand_entity = lr_id$ensembl[match(from, lr_id$symbol)],
          receptor_entity = lr_id$ensembl[match(to, lr_id$symbol)])]
lr[, prior_ligand := mapl$symbol[match(ligand_entity, mapl$ensembl)]]
save_tsv(lr, "LR_entity_mapping.tsv", na = "NA")

# single = TRUE scores each ligand on its own, so the ligand list can be split across workers.
predict_activities <- function(geneset, background, ligands) {
  chunks <- split(ligands, ceiling(seq_along(ligands) / ceiling(length(ligands) / WORKERS)))
  ans <- parallel::mclapply(chunks, function(pp) predict_ligand_activities(
    geneset = geneset, background_expressed_genes = background, ligand_target_matrix = lt, potential_ligands = pp,
    single = TRUE), mc.cores = WORKERS, mc.preschedule = TRUE)
  stopifnot(!any(vapply(ans, inherits, logical(1), "try-error")))
  ans <- rbindlist(ans)[match(ligands, test_ligand)]
  setorder(ans, -pearson)
  ans
}

# ---- Main contexts: donor x receiver state x program ----
cov <- list(); acts <- list(); edges <- list(); targets <- list(); bgrows <- list(); targetrows <- list()
for (d in unique(expr$donor)) for (st in STATES) {
  rec <- expr[donor == d & group == st]
  if (!nrow(rec)) next
  eg <- rec[detection >= MIN_DETECTION]$ensembl
  bt <- mapt[ensembl %in% eg]
  bg <- bt$symbol
  send <- expr[donor == d & !grepl("^SPG_", group) & detection >= MIN_DETECTION]
  elr <- merge(send, lr, by.x = "ensembl", by.y = "ligand_entity", allow.cartesian = TRUE)
  elr <- elr[receptor_entity %in% eg & !is.na(prior_ligand)]
  potential <- colnames(lt)[colnames(lt) %in% unique(elr$prior_ligand)]
  bgrows[[length(bgrows) + 1]] <- bt[, .(donor = d, state = st, assay_symbol, ensembl, prior_target = symbol)]
  for (m in PROGRAMS) {
    gs <- unique(mp[MP == m & ensembl %in% eg]$prior_target)
    gs <- gs[!is.na(gs) & gs %in% bg]
    cov[[length(cov) + 1]] <- data.table(donor = d, state = st, MP = m, n_cells = unique(rec$n_cells), n_target = length(gs),
      background = length(bg), n_ligands = length(potential),
      status = ifelse(length(gs) >= MIN_TARGETS && length(potential) > 0, "estimated", "ineligible"))
    if (length(gs) < MIN_TARGETS || !length(potential)) next
    a <- predict_activities(gs, bg, potential)
    a[, `:=`(donor = d, state = st, MP = m, rank = seq_len(.N))]
    acts[[length(acts) + 1]] <- a
    ee <- copy(elr)
    ee[, `:=`(state = st, MP = m)]
    receptors <- rec[, .(ensembl, receptor_detection = detection, receptor_expression = mean_logCP10k)]
    edges[[length(edges) + 1]] <- merge(ee, receptors, by.x = "receptor_entity", by.y = "ensembl")
    tl <- as.data.table(as.table(lt[gs, potential, drop = FALSE]))
    setnames(tl, c("target", "ligand", "regulatory_potential"))
    tl[, `:=`(donor = d, state = st, MP = m)]
    targets[[length(targets) + 1]] <- tl
    targetrows[[length(targetrows) + 1]] <- data.table(donor = d, state = st, MP = m, target = gs)
  }
}
act <- rbindlist(acts)
save_tsv(act, "ligand_activities.tsv")
save_tsv(rbindlist(cov), "nichenet_eligibility.tsv")
save_tsv(rbindlist(edges, fill = TRUE), "sender_receptor_evidence.tsv.gz")
save_tsv(rbindlist(targets), "all_ligand_target_potential.tsv.gz")
save_tsv(rbindlist(targetrows), "all_targets.tsv")
save_tsv(rbindlist(bgrows), "all_background_genes.tsv.gz")

# ---- Sensitivities: U contexts, MP3 without ID genes or without published state markers ----
pr <- fread(file.path(EXTERNAL, "prior_gene_lists", "bush2024_state_markers.tsv"))
marker <- resolve(unique(pr[!state %in% c("G1/S", "G2/M")]$gene))
save_tsv(marker, "sensitivity_marker_entity_mapping.tsv", na = "NA")
sens <- list(); cover <- list()
for (d in unique(act$donor)) {
  eg <- expr[donor == d & group == "SPG_U" & detection >= MIN_DETECTION]$ensembl
  bg <- mapt[ensembl %in% eg]$symbol
  pot <- act[donor == d & state == "SPG_U" & MP == "MP_3"]$test_ligand
  for (s in c("remove_published_state_markers", "remove_ID1_ID2_ID4")) {
    rem <- if (s == "remove_ID1_ID2_ID4") resolve(ID_GENES)$ensembl else marker$ensembl
    gs <- unique(mp[MP == "MP_3" & ensembl %in% eg & !ensembl %in% rem]$prior_target)
    gs <- gs[!is.na(gs) & gs %in% bg]
    cover[[length(cover) + 1]] <- data.table(donor = d, sensitivity = s, n_genes = length(gs),
                                             genes = paste(gs, collapse = ";"))
    a <- predict_activities(gs, bg, pot)
    a[, `:=`(donor = d, sensitivity = s, rank = seq_len(.N))]
    sens[[length(sens) + 1]] <- a
  }
}
save_tsv(rbindlist(sens), "ligand_sensitivity.tsv")
save_tsv(rbindlist(cover), "ligand_sensitivity_members.tsv")

sel <- act[state == "SPG_U" & MP == "MP_3", .(median_pearson = median(pearson), n_contexts = .N, median_rank = median(rank)),
           by = test_ligand][order(-median_pearson, test_ligand)]
save_tsv(sel, "U_ligand_selection.tsv")
cat("NicheNet:", uniqueN(act[, .(donor, state, MP)]), "main contexts,", length(sens), "sensitivity contexts; top U ligands:",
    paste(head(sel$test_ligand, 10), collapse = ", "), "\n")
