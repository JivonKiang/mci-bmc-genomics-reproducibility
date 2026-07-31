suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(edgeR)
})

root_dir <- "E:/20241004_MCI"
stamp <- "20260722_140059"
out_dir <- file.path(root_dir, "20260709", paste0("mci_centered_score_external_tests_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_path <- file.path(root_dir, "20260709", paste0("mci_centered_score_denovo_", stamp), "mci_centered_candidate_registry_weighted.csv")
candidates <- read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
candidates$gene <- toupper(trimws(candidates$gene))

zscore_genes <- function(mat) {
  z <- t(apply(mat, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

score_matrix <- function(mat, meta, candidates, dataset) {
  mapped <- intersect(candidates$gene, rownames(mat))
  use <- candidates %>% filter(gene %in% mapped)
  z <- zscore_genes(mat[use$gene, , drop = FALSE])
  weights <- use$score_weight
  names(weights) <- use$gene
  score <- as.numeric(crossprod(weights, z[use$gene, , drop = FALSE]))
  names(score) <- colnames(mat)
  meta <- meta[match(names(score), meta$sample), , drop = FALSE]
  meta$mci_centered_score <- score
  meta$dataset <- dataset
  list(meta = meta, mapped = mapped, missing = setdiff(candidates$gene, mapped))
}

fit_contrasts <- function(meta, stage_var = "stage", contrasts = c("MCI_vs_HC", "AD_vs_MCI", "AD_vs_HC")) {
  meta[[stage_var]] <- factor(meta[[stage_var]])
  fit <- lm(reformulate(c("age_num", "sex_factor", stage_var), response = "mci_centered_score"), data = meta)
  b <- coef(fit); v <- vcov(fit); cn <- names(b)
  stage_levels <- levels(meta[[stage_var]])
  make_one <- function(label, w) {
    ww <- setNames(rep(0, length(cn)), cn)
    for (nm in names(w)) if (nm %in% cn) ww[nm] <- w[[nm]]
    est <- sum(ww * b); se <- sqrt(as.numeric(t(ww) %*% v %*% ww)); tval <- est / se
    data.frame(dataset = unique(meta$dataset), contrast = label, estimate = est, se = se, df = df.residual(fit), p_value = 2 * pt(abs(tval), df = df.residual(fit), lower.tail = FALSE), n = nrow(meta), stringsAsFactors = FALSE)
  }
  out <- list()
  if (all(c("HC", "MCI", "AD") %in% stage_levels)) {
    out <- list(make_one("MCI_vs_HC", c(stageMCI = 1)), make_one("AD_vs_MCI", c(stageAD = 1, stageMCI = -1)), make_one("AD_vs_HC", c(stageAD = 1)))
  } else if (all(c("P-MCI", "S-MCI") %in% stage_levels)) {
    out <- list(make_one("P-MCI_vs_S-MCI", c(`stageP-MCI` = 1)))
  }
  bind_rows(out)
}

read_gse249477 <- function() {
  count_path <- file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE249477_count_clean_full.gz")
  x <- fread(count_path, data.table = FALSE, check.names = FALSE)
  total_cols <- grep(" - Total counts$", colnames(x), value = TRUE)
  sample_key <- sub(" \\(GE\\) - Total counts$", "", total_cols)
  mat <- as.matrix(x[, total_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
  genes <- toupper(trimws(as.character(x$Name)))
  keep <- !is.na(genes) & genes != ""
  mat <- mat[keep, , drop = FALSE]; genes <- genes[keep]
  complete_rows <- rowSums(is.na(mat)) == 0
  mat <- mat[complete_rows, , drop = FALSE]; genes <- genes[complete_rows]
  rownames(mat) <- genes; mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  colnames(mat) <- sample_key
  audit <- read.csv(file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE249477_sample_label_audit.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  audit$sample_key <- vapply(regmatches(audit$title, gregexpr("DK[0-9]+_[0-9]+", audit$title)), function(z) if (length(z)) z[[1]] else NA_character_, character(1))
  audit$stage <- factor(audit$group, levels = c("HC", "MCI", "AD"))
  audit$age_num <- suppressWarnings(as.numeric(audit$age)); audit$sex_factor <- factor(audit$sex)
  meta <- audit %>% transmute(sample = sample_key, stage, age_num, sex_factor) %>% filter(sample %in% colnames(mat))
  mat <- mat[, meta$sample, drop = FALSE]
  list(mat = edgeR::cpm(edgeR::DGEList(counts = round(mat)), log = TRUE, prior.count = 1), meta = meta)
}

read_gse282742 <- function() {
  count_path <- file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE282742_Expected_count.txt.gz")
  x <- fread(count_path, data.table = FALSE, check.names = FALSE)
  sample_cols <- setdiff(colnames(x), "gene_id")
  ids <- sub("\\..*$", "", as.character(x$gene_id))
  ann <- fread(file.path(root_dir, "20250816 revise", "GEO", "Human.GRCh38.p13.annot.tsv.gz"), data.table = FALSE, select = c("EnsemblGeneID", "Symbol"))
  ann <- ann[ann$EnsemblGeneID != "" & ann$Symbol != "", , drop = FALSE]
  symbols <- toupper(ann$Symbol[match(ids, ann$EnsemblGeneID)])
  keep <- !is.na(symbols) & symbols != ""
  mat <- as.matrix(x[keep, sample_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
  complete_rows <- rowSums(is.na(mat)) == 0
  mat <- mat[complete_rows, , drop = FALSE]; symbols <- symbols[keep][complete_rows]
  rownames(mat) <- symbols; mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  colnames(mat) <- sample_cols
  meta0 <- read.csv(file.path(root_dir, "20260709", "external_public_data_20260721_151749", "GSE282742_sample_subject_audit.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  meta <- meta0 %>% transmute(sample = vgh_id, group, subject_id, age_num = as.numeric(age), sex_factor = factor(sex)) %>% filter(sample %in% colnames(mat), group %in% c("P-MCI", "S-MCI", "AD"))
  meta$stage <- factor(meta$group, levels = c("S-MCI", "P-MCI", "AD"))
  mat <- mat[, meta$sample, drop = FALSE]
  list(mat = edgeR::cpm(edgeR::DGEList(counts = round(mat)), log = TRUE, prior.count = 1), meta = meta)
}

ext <- read_gse249477()
score_ext <- score_matrix(ext$mat, ext$meta, candidates, "GSE249477")
write.csv(score_ext$meta, file.path(out_dir, "GSE249477_sample_mci_centered_scores.csv"), row.names = FALSE)
write.csv(data.frame(dataset = "GSE249477", n_candidates = nrow(candidates), n_mapped = length(score_ext$mapped), n_missing = length(score_ext$missing), missing_genes = paste(score_ext$missing, collapse = ";")), file.path(out_dir, "GSE249477_mapping_audit.csv"), row.names = FALSE)

prog <- read_gse282742()
score_prog <- score_matrix(prog$mat, prog$meta, candidates, "GSE282742")
write.csv(score_prog$meta, file.path(out_dir, "GSE282742_sample_mci_centered_scores.csv"), row.names = FALSE)
write.csv(data.frame(dataset = "GSE282742", n_candidates = nrow(candidates), n_mapped = length(score_prog$mapped), n_missing = length(score_prog$missing), missing_genes = paste(score_prog$missing, collapse = ";")), file.path(out_dir, "GSE282742_mapping_audit.csv"), row.names = FALSE)

external_contrasts <- fit_contrasts(score_ext$meta)
subject_scores <- score_prog$meta %>% filter(group %in% c("P-MCI", "S-MCI")) %>% group_by(subject_id, group) %>% summarise(sample = first(sample), age_num = mean(age_num, na.rm = TRUE), sex_factor = first(sex_factor), mci_centered_score = mean(mci_centered_score, na.rm = TRUE), dataset = first(dataset), .groups = "drop") %>% mutate(stage = factor(group, levels = c("S-MCI", "P-MCI")))
progression_contrasts <- fit_contrasts(subject_scores, stage_var = "stage")
progression_mapping_fraction <- length(score_prog$mapped) / nrow(candidates)
progression_contrasts$mapping_fraction <- progression_mapping_fraction
progression_contrasts$evaluable <- progression_mapping_fraction >= 0.80
write.csv(external_contrasts, file.path(out_dir, "GSE249477_mci_centered_score_contrasts.csv"), row.names = FALSE)
write.csv(subject_scores, file.path(out_dir, "GSE282742_subject_level_mci_centered_scores.csv"), row.names = FALSE)
write.csv(progression_contrasts, file.path(out_dir, "GSE282742_subject_level_mci_centered_contrasts.csv"), row.names = FALSE)

writeLines(c(
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "GSE249477: fixed 95-gene discovery score tested on independent cross-sectional count data.",
  paste0("GSE282742: fixed 95-gene discovery score tested on processed expected counts with subject-level aggregation; mapping fraction = ", sprintf("%.3f", progression_mapping_fraction), "."),
  "GSE282742 is not evaluable as a definitive progression validation when the mapping fraction is below 0.80.",
  "These tests are transportability and progression-concordance analyses, not clinical prediction validation.",
  "The score remains a candidate molecular staging signature until a fully independent cohort and clinical increment analysis are available."
), file.path(out_dir, "external_score_QA_notes.txt"))

cat("Wrote external score tests to", out_dir, "\n")
