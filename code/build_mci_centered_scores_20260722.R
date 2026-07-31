suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

root_dir <- "E:/20241004_MCI"
stamp <- "20260722_140059"
out_dir <- file.path(root_dir, "20260709", paste0("mci_centered_score_denovo_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_path <- file.path(root_dir, "20260709", "mci_stage_marker_reclassification_20260722_140059", "strict_mci_centered_candidates.csv")
registry <- read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
registry <- registry %>% rename(
  estimate_MCI_vs_HC = estimate__MCI_vs_HC,
  estimate_AD_vs_MCI = estimate__AD_vs_MCI,
  fdr_MCI_vs_HC = fdr__MCI_vs_HC,
  fdr_AD_vs_MCI = fdr__AD_vs_MCI
)
registry$gene <- toupper(trimws(registry$gene))
registry <- registry %>% filter(!is.na(gene), gene != "")

mci_candidates <- registry %>% mutate(mci_centered = TRUE) %>% distinct(gene, trajectory_class, estimate_MCI_vs_HC, estimate_AD_vs_MCI, fdr_MCI_vs_HC, fdr_AD_vs_MCI, .keep_all = TRUE)
if (!nrow(mci_candidates)) stop("No strict MCI-centered candidates in the discovery registry")

# Equalise the contribution of the MCI-peak and MCI-trough groups.
n_peak <- sum(mci_candidates$trajectory_class == "MCI_peak")
n_trough <- sum(mci_candidates$trajectory_class == "MCI_trough")
mci_candidates$score_weight <- ifelse(mci_candidates$trajectory_class == "MCI_peak", 0.5 / n_peak, -0.5 / n_trough)
mci_candidates$max_transition_fdr <- pmax(mci_candidates$fdr_MCI_vs_HC, mci_candidates$fdr_AD_vs_MCI)
mci_candidates$effect_amplitude <- abs(mci_candidates$estimate_MCI_vs_HC) + abs(mci_candidates$estimate_AD_vs_MCI)
mci_candidates$stability_rank <- rank(-log10(pmax(mci_candidates$max_transition_fdr, .Machine$double.xmin)) * mci_candidates$effect_amplitude, ties.method = "first")
mci_candidates <- mci_candidates %>% arrange(stability_rank)
write.csv(mci_candidates, file.path(out_dir, "mci_centered_candidate_registry_weighted.csv"), row.names = FALSE)

read_array <- function(path, design_path, dataset) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  gene_col <- colnames(x)[1]
  genes <- toupper(trimws(as.character(x[[gene_col]])))
  x <- x[, -1, drop = FALSE]
  design <- fread(design_path, data.table = FALSE, check.names = FALSE)
  colnames(design)[1] <- "sample"
  sample_ids <- intersect(colnames(x), design$sample)
  meta <- design[match(sample_ids, design$sample), , drop = FALSE]
  keep <- meta$included.in.case..control.study == "yes" & meta$status %in% c("CTL", "MCI", "AD")
  sample_ids <- sample_ids[keep]
  meta <- meta[keep, , drop = FALSE]
  mat <- as.matrix(x[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  meta$stage <- factor(ifelse(meta$status == "CTL", "HC", meta$status), levels = c("HC", "MCI", "AD"))
  meta$age_num <- suppressWarnings(as.numeric(meta$age))
  meta$sex_factor <- factor(meta$gender)
  keep_meta <- complete.cases(meta[, c("stage", "age_num", "sex_factor")])
  list(dataset = dataset, expr = mat[, meta$sample[keep_meta], drop = FALSE], meta = meta[keep_meta, , drop = FALSE])
}

zscore_genes <- function(mat) {
  z <- t(apply(mat, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

score_dataset <- function(obj, candidates) {
  mapped <- intersect(candidates$gene, rownames(obj$expr))
  missing <- setdiff(candidates$gene, mapped)
  use <- candidates %>% filter(gene %in% mapped)
  z <- zscore_genes(obj$expr[mapped, , drop = FALSE])
  weights <- use$score_weight[match(rownames(z), use$gene)]
  score <- as.numeric(crossprod(weights, z[use$gene, , drop = FALSE]))
  names(score) <- colnames(obj$expr)
  meta <- obj$meta[match(names(score), obj$meta$sample), , drop = FALSE]
  meta$mci_centered_score <- score
  meta$stage_num <- as.numeric(meta$stage) - 1
  meta$dataset <- obj$dataset
  list(meta = meta, mapped = mapped, missing = missing)
}

contrast_fit <- function(meta, score_name = "mci_centered_score") {
  fit <- lm(reformulate(c("age_num", "sex_factor", "stage"), response = score_name), data = meta)
  b <- coef(fit)
  v <- vcov(fit)
  cn <- names(b)
  contrast_one <- function(label, weights) {
    w <- setNames(rep(0, length(cn)), cn)
    for (nm in names(weights)) if (nm %in% cn) w[nm] <- weights[[nm]]
    est <- sum(w * b)
    se <- sqrt(as.numeric(t(w) %*% v %*% w))
    tval <- est / se
    data.frame(dataset = unique(meta$dataset), score = score_name, contrast = label, estimate = est, se = se, df = df.residual(fit), p_value = 2 * pt(abs(tval), df = df.residual(fit), lower.tail = FALSE), n = nrow(meta), stringsAsFactors = FALSE)
  }
  rbind(
    contrast_one("MCI_vs_HC", c(stageMCI = 1)),
    contrast_one("AD_vs_MCI", c(stageAD = 1, stageMCI = -1)),
    contrast_one("AD_vs_HC", c(stageAD = 1))
  )
}

score_objects <- list(
  GSE63060 = read_array(file.path(root_dir, "20250816 revise", "GEO", "data", "GSE63060_data_matrix.csv"), file.path(root_dir, "20260709", "data", "GSE63060_design_matrix.csv"), "GSE63060"),
  GSE63061 = read_array(file.path(root_dir, "20260709", "data", "GSE63061_data_matrix.csv"), file.path(root_dir, "20260709", "data", "GSE63061_design_matrix.csv"), "GSE63061")
)

scored <- lapply(score_objects, score_dataset, candidates = mci_candidates)
sample_scores <- bind_rows(lapply(scored, function(x) x$meta))
write.csv(sample_scores, file.path(out_dir, "sample_level_mci_centered_scores.csv"), row.names = FALSE)

mapping_audit <- bind_rows(lapply(names(scored), function(nm) {
  x <- scored[[nm]]
  data.frame(dataset = nm, n_candidates = nrow(mci_candidates), n_mapped = length(x$mapped), n_missing = length(x$missing), missing_genes = paste(x$missing, collapse = ";"), stringsAsFactors = FALSE)
}))
write.csv(mapping_audit, file.path(out_dir, "score_mapping_audit.csv"), row.names = FALSE)

contrasts <- bind_rows(lapply(scored, function(x) contrast_fit(x$meta)))
contrasts$significant_nominal <- contrasts$p_value < 0.05
write.csv(contrasts, file.path(out_dir, "mci_centered_score_contrasts.csv"), row.names = FALSE)

stage_summary <- sample_scores %>% group_by(dataset, stage) %>% summarise(n = n(), mean_score = mean(mci_centered_score), sd_score = sd(mci_centered_score), median_score = median(mci_centered_score), .groups = "drop")
write.csv(stage_summary, file.path(out_dir, "mci_centered_score_stage_summary.csv"), row.names = FALSE)

# A provisional minimal panel is a reproducibility input, not a validated biomarker panel.
provisional_panel <- mci_candidates %>% filter(stability_rank <= 20) %>% arrange(stability_rank)
write.csv(provisional_panel, file.path(out_dir, "provisional_top20_mci_centered_panel.csv"), row.names = FALSE)

writeLines(c(
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Score definition: within-cohort gene z-scores; equal total weight for MCI-peak and MCI-trough groups; discovery directions frozen before scoring GSE63061.",
  paste0("De novo strict MCI-centered candidates: ", nrow(mci_candidates), " (MCI-peak ", n_peak, "; MCI-trough ", n_trough, ")."),
  "Interpretation status: candidate molecular staging score, not a validated clinical biomarker.",
  "The provisional top-20 panel is for the next locked modelling step only; it is not a final gene list."
), file.path(out_dir, "score_analysis_QA_notes.txt"))

cat("Wrote MCI-centered score analysis to", out_dir, "\n")
