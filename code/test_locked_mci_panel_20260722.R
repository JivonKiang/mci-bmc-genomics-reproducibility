suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(edgeR)
})

root_dir <- "E:/20241004_MCI"
panel_source_stamp <- "20260722_144351"
panel_path <- file.path(root_dir, "20260709", paste0("mci_centered_panel_compression_", panel_source_stamp), "provisional_panels_by_size.csv")
panel_all <- fread(panel_path, data.table = FALSE, check.names = FALSE)
panel <- panel_all %>% filter(panel_size == 12) %>% arrange(trajectory_class, gene)
panel$gene <- toupper(trimws(panel$gene))
panel <- panel %>% mutate(score_weight = ifelse(trajectory_class == "MCI_peak", 0.5 / sum(trajectory_class == "MCI_peak"), -0.5 / sum(trajectory_class == "MCI_trough")))
if (nrow(panel) != 12 || anyDuplicated(panel$gene)) stop("Locked 12-feature panel is not valid")

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(root_dir, "20260709", paste0("locked_mci_panel_tests_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(panel, file.path(out_dir, "locked_12_feature_panel.csv"), row.names = FALSE)

zscore_genes <- function(mat) {
  z <- t(apply(mat, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
  rownames(z) <- rownames(mat); colnames(z) <- colnames(mat); z
}

score_matrix <- function(mat, meta, dataset) {
  mapped <- intersect(panel$gene, rownames(mat))
  use <- panel %>% filter(gene %in% mapped)
  z <- zscore_genes(mat[use$gene, , drop = FALSE])
  weights <- use$score_weight; names(weights) <- use$gene
  score <- as.numeric(crossprod(weights, z[use$gene, , drop = FALSE])); names(score) <- colnames(mat)
  meta <- meta[match(names(score), meta$sample), , drop = FALSE]
  meta$mci_centered_score <- score; meta$dataset <- dataset
  list(meta = meta, mapped = mapped, missing = setdiff(panel$gene, mapped))
}

fit_contrasts <- function(meta, stage_var = "stage") {
  meta[[stage_var]] <- factor(meta[[stage_var]])
  fit <- lm(reformulate(c("age_num", "sex_factor", stage_var), response = "mci_centered_score"), data = meta)
  b <- coef(fit); v <- vcov(fit); cn <- names(b)
  make_one <- function(label, w) {
    ww <- setNames(rep(0, length(cn)), cn)
    for (nm in names(w)) if (nm %in% cn) ww[nm] <- w[[nm]]
    est <- sum(ww * b); se <- sqrt(as.numeric(t(ww) %*% v %*% ww)); tval <- est / se
    data.frame(dataset = unique(meta$dataset), contrast = label, estimate = est, se = se, df = df.residual(fit), p_value = 2 * pt(abs(tval), df = df.residual(fit), lower.tail = FALSE), n = nrow(meta), stringsAsFactors = FALSE)
  }
  lev <- levels(meta[[stage_var]])
  if (all(c("HC", "MCI", "AD") %in% lev)) {
    bind_rows(make_one("MCI_vs_HC", c(stageMCI = 1)), make_one("AD_vs_MCI", c(stageAD = 1, stageMCI = -1)), make_one("AD_vs_HC", c(stageAD = 1)))
  } else if (all(c("S-MCI", "P-MCI") %in% lev)) {
    bind_rows(make_one("P-MCI_vs_S-MCI", c(`stageP-MCI` = 1)))
  } else stop("Unexpected stage levels")
}

read_array <- function(path, design_path, dataset) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  genes <- toupper(trimws(as.character(x[[1]]))); x <- x[, -1, drop = FALSE]
  design <- fread(design_path, data.table = FALSE, check.names = FALSE); colnames(design)[1] <- "sample"
  sample_ids <- intersect(colnames(x), design$sample); meta <- design[match(sample_ids, design$sample), , drop = FALSE]
  keep <- meta$included.in.case..control.study == "yes" & meta$status %in% c("CTL", "MCI", "AD")
  sample_ids <- sample_ids[keep]; meta <- meta[keep, , drop = FALSE]
  mat <- as.matrix(x[, sample_ids, drop = FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- genes
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]; mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  meta$stage <- factor(ifelse(meta$status == "CTL", "HC", meta$status), levels = c("HC", "MCI", "AD")); meta$age_num <- suppressWarnings(as.numeric(meta$age)); meta$sex_factor <- factor(meta$gender)
  keep_meta <- complete.cases(meta[, c("stage", "age_num", "sex_factor")]); list(mat = mat[, meta$sample[keep_meta], drop = FALSE], meta = meta[keep_meta, , drop = FALSE], dataset = dataset)
}

read_gse249477 <- function() {
  count_path <- file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE249477_count_clean_full.gz")
  x <- fread(count_path, data.table = FALSE, check.names = FALSE); total_cols <- grep(" - Total counts$", colnames(x), value = TRUE)
  sample_key <- sub(" \\(GE\\) - Total counts$", "", total_cols); mat <- as.matrix(x[, total_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
  genes <- toupper(trimws(as.character(x$Name))); keep <- !is.na(genes) & genes != ""; mat <- mat[keep, , drop = FALSE]; genes <- genes[keep]
  complete_rows <- rowSums(is.na(mat)) == 0; mat <- mat[complete_rows, , drop = FALSE]; genes <- genes[complete_rows]; rownames(mat) <- genes; mat <- rowsum(mat, group = rownames(mat), reorder = FALSE); colnames(mat) <- sample_key
  audit <- read.csv(file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE249477_sample_label_audit.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  audit$sample_key <- vapply(regmatches(audit$title, gregexpr("DK[0-9]+_[0-9]+", audit$title)), function(z) if (length(z)) z[[1]] else NA_character_, character(1))
  audit$stage <- factor(audit$group, levels = c("HC", "MCI", "AD")); audit$age_num <- suppressWarnings(as.numeric(audit$age)); audit$sex_factor <- factor(audit$sex)
  meta <- audit %>% transmute(sample = sample_key, stage, age_num, sex_factor) %>% filter(sample %in% colnames(mat)); mat <- mat[, meta$sample, drop = FALSE]
  list(mat = edgeR::cpm(edgeR::DGEList(counts = round(mat)), log = TRUE, prior.count = 1), meta = meta, dataset = "GSE249477")
}

read_gse282742 <- function() {
  count_path <- file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE282742_Expected_count.txt.gz")
  x <- fread(count_path, data.table = FALSE, check.names = FALSE); sample_cols <- setdiff(colnames(x), "gene_id"); ids <- sub("\\..*$", "", as.character(x$gene_id))
  ann <- fread(file.path(root_dir, "20250816 revise", "GEO", "Human.GRCh38.p13.annot.tsv.gz"), data.table = FALSE, select = c("EnsemblGeneID", "Symbol")); ann <- ann[ann$EnsemblGeneID != "" & ann$Symbol != "", , drop = FALSE]
  symbols <- toupper(ann$Symbol[match(ids, ann$EnsemblGeneID)]); keep <- !is.na(symbols) & symbols != ""; mat <- as.matrix(x[keep, sample_cols, drop = FALSE]); storage.mode(mat) <- "numeric"; complete_rows <- rowSums(is.na(mat)) == 0
  mat <- mat[complete_rows, , drop = FALSE]; symbols <- symbols[keep][complete_rows]; rownames(mat) <- symbols; mat <- rowsum(mat, group = rownames(mat), reorder = FALSE); colnames(mat) <- sample_cols
  meta0 <- read.csv(file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE282742_sample_subject_audit.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  meta <- meta0 %>% transmute(sample = vgh_id, group, subject_id, age_num = as.numeric(age), sex_factor = factor(sex)) %>% filter(sample %in% colnames(mat), group %in% c("P-MCI", "S-MCI", "AD")); meta$stage <- factor(meta$group, levels = c("S-MCI", "P-MCI", "AD")); mat <- mat[, meta$sample, drop = FALSE]
  list(mat = edgeR::cpm(edgeR::DGEList(counts = round(mat)), log = TRUE, prior.count = 1), meta = meta, dataset = "GSE282742")
}

disc <- read_array(file.path(root_dir, "20260709", "data", "GSE63061_data_matrix.csv"), file.path(root_dir, "20260709", "data", "GSE63061_design_matrix.csv"), "GSE63061")
ext <- read_gse249477(); prog <- read_gse282742()
scored_disc <- score_matrix(disc$mat, disc$meta, "GSE63061")
scored_ext <- score_matrix(ext$mat, ext$meta, "GSE249477")
scored_prog <- score_matrix(prog$mat, prog$meta, "GSE282742")

write.csv(scored_disc$meta, file.path(out_dir, "GSE63061_locked_panel_scores.csv"), row.names = FALSE)
write.csv(scored_ext$meta, file.path(out_dir, "GSE249477_locked_panel_scores.csv"), row.names = FALSE)
write.csv(scored_prog$meta, file.path(out_dir, "GSE282742_locked_panel_scores.csv"), row.names = FALSE)
write.csv(bind_rows(fit_contrasts(scored_disc$meta), fit_contrasts(scored_ext$meta)), file.path(out_dir, "locked_panel_cross_sectional_contrasts.csv"), row.names = FALSE)
write.csv(bind_rows(
  data.frame(dataset = "GSE63061", n_panel = nrow(panel), n_mapped = length(scored_disc$mapped), n_missing = length(scored_disc$missing), missing_genes = paste(scored_disc$missing, collapse = ";")),
  data.frame(dataset = "GSE249477", n_panel = nrow(panel), n_mapped = length(scored_ext$mapped), n_missing = length(scored_ext$missing), missing_genes = paste(scored_ext$missing, collapse = ";")),
  data.frame(dataset = "GSE282742", n_panel = nrow(panel), n_mapped = length(scored_prog$mapped), n_missing = length(scored_prog$missing), missing_genes = paste(scored_prog$missing, collapse = ";"))
), file.path(out_dir, "locked_panel_mapping_audit.csv"), row.names = FALSE)

subject_prog <- scored_prog$meta %>% filter(group %in% c("P-MCI", "S-MCI")) %>% group_by(subject_id, group) %>% summarise(sample = first(sample), age_num = mean(age_num, na.rm = TRUE), sex_factor = first(sex_factor), mci_centered_score = mean(mci_centered_score, na.rm = TRUE), dataset = first(dataset), .groups = "drop") %>% mutate(stage = factor(group, levels = c("S-MCI", "P-MCI")))
progression <- fit_contrasts(subject_prog)
progression$mapping_fraction <- length(scored_prog$mapped) / nrow(panel)
progression$evaluable <- progression$mapping_fraction >= 0.80
write.csv(subject_prog, file.path(out_dir, "GSE282742_locked_panel_subject_scores.csv"), row.names = FALSE)
write.csv(progression, file.path(out_dir, "GSE282742_locked_panel_progression_contrast.csv"), row.names = FALSE)

writeLines(c(
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "The 12-feature panel was locked from GSE63060 discovery-only bootstrap output before testing.",
  paste0("GSE63061 mapping: ", length(scored_disc$mapped), "/", nrow(panel), "; GSE249477 mapping: ", length(scored_ext$mapped), "/", nrow(panel), "; GSE282742 mapping: ", length(scored_prog$mapped), "/", nrow(panel), "."),
  "GSE282742 progression is considered evaluable only when the locked panel mapping fraction is at least 0.80 and subject-level P-MCI/S-MCI labels are available.",
  "No genes were reselected after opening any test cohort. Results remain research-stage and do not establish clinical utility."
), file.path(out_dir, "locked_panel_QA_notes.txt"))
cat("Wrote locked panel tests to", out_dir, "\n")
