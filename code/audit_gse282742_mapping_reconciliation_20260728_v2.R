options(stringsAsFactors = FALSE)

ROOT <- "E:/20241004_MCI/20260709"
PUBLIC <- file.path(ROOT, "external_public_data_20260721_151749")
OUT <- file.path(ROOT, "gse282742_mapping_reconciliation_20260728_v2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

annotation_path <- "E:/20241004_MCI/20250816 revise/GEO/Human.GRCh38.p13.annot.tsv.gz"
count_path <- file.path(PUBLIC, "GSE282742_Expected_count.full.gz")
sample_path <- file.path(PUBLIC, "GSE282742_sample_subject_audit.csv")
programme_path <- file.path(ROOT, "figures_mci_manuscript_rebuild_20260728_105919", "source_figure2_95_gene_programme.csv")
panel_path <- file.path(ROOT, "mci_development_optimized_panel_20260722_161919", "development_locked_12_gene_panel.csv")

annotation <- read.delim(gzfile(annotation_path), sep = "\t", quote = "", check.names = FALSE,
                         colClasses = "character", na.strings = c("", "NA"))
annotation$Symbol <- ifelse(is.na(annotation$Symbol), "", trimws(annotation$Symbol))
annotation$Synonyms <- ifelse(is.na(annotation$Synonyms), "", annotation$Synonyms)
annotation$EnsemblGeneID <- ifelse(is.na(annotation$EnsemblGeneID), "", trimws(annotation$EnsemblGeneID))

exact <- list()
aliases <- list()
for (i in which(annotation$EnsemblGeneID != "")) {
  symbol <- toupper(annotation$Symbol[i])
  ens <- annotation$EnsemblGeneID[i]
  if (nzchar(symbol)) exact[[symbol]] <- unique(c(exact[[symbol]], ens))
  tokens <- unlist(strsplit(annotation$Synonyms[i], "\\|"))
  tokens <- toupper(trimws(tokens[nzchar(trimws(tokens))]))
  for (token in tokens) aliases[[token]] <- unique(c(aliases[[token]], ens))
}

resolve_symbol <- function(symbol) {
  key <- toupper(trimws(symbol))
  hits <- exact[[key]]
  if (is.null(hits) || length(hits) == 0) hits <- aliases[[key]]
  if (is.null(hits) || length(hits) != 1) return(NA_character_)
  hits[[1]]
}

counts <- read.delim(gzfile(count_path), sep = "\t", quote = "", check.names = FALSE,
                     row.names = 1, na.strings = c("", "NA"))
rownames(counts) <- sub("\\.[0-9]+$", "", rownames(counts))
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
counts[] <- lapply(counts, as.numeric)
sample_audit <- read.csv(sample_path, check.names = FALSE)
sample_audit <- sample_audit[match(colnames(counts), sample_audit$vgh_id), , drop = FALSE]
if (any(is.na(sample_audit$vgh_id))) stop("Count/sample audit mismatch")

log_cpm <- log2(sweep(as.matrix(counts), 2, colSums(counts), "/") * 1e6 + 1)

programme <- read.csv(programme_path, check.names = FALSE)
programme$ensembl_id_alias_aware <- vapply(programme$gene, resolve_symbol, character(1))
programme$mapped_in_matrix <- !is.na(programme$ensembl_id_alias_aware) & programme$ensembl_id_alias_aware %in% rownames(counts)
programme$mapping_status <- ifelse(programme$mapped_in_matrix, "mapped_alias_aware", "unresolved_or_not_in_matrix")
write.csv(programme[, c("gene", "trajectory_class", "ensembl_id_alias_aware", "mapping_status")],
          file.path(OUT, "GSE282742_95_gene_mapping_alias_aware.csv"), row.names = FALSE, na = "")

panel <- read.csv(panel_path, check.names = FALSE)
panel$ensembl_id_alias_aware <- vapply(panel$gene, resolve_symbol, character(1))
panel$mapped_in_matrix <- !is.na(panel$ensembl_id_alias_aware) & panel$ensembl_id_alias_aware %in% rownames(counts)
panel$mapping_status <- ifelse(panel$mapped_in_matrix, "mapped_alias_aware", "unresolved_or_not_in_matrix")
write.csv(panel[, c("gene", "trajectory_class", "ensembl_id_alias_aware", "mapping_status")],
          file.path(OUT, "GSE282742_12_panel_mapping_alias_aware.csv"), row.names = FALSE, na = "")

score_for <- function(registry, renormalize = TRUE) {
  keep <- registry$mapped_in_matrix
  dat <- registry[keep, , drop = FALSE]
  Z <- t(apply(log_cpm[dat$ensembl_id_alias_aware, , drop = FALSE], 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
  direction <- ifelse(dat$trajectory_class == "MCI_peak", 1, -1)
  weights <- direction / length(unique(registry$gene))
  if (renormalize) weights <- weights / sum(abs(weights))
  as.numeric(weights %*% Z)
}

score_table <- data.frame(sample_audit[, c("vgh_id", "sample", "group", "subject_id", "age", "sex")],
                          score = score_for(panel, TRUE))
write.csv(score_table, file.path(OUT, "GSE282742_12_panel_alias_aware_sample_scores.csv"), row.names = FALSE, na = "")

ord <- order(score_table$subject_id, ifelse(is.na(score_table$age), Inf, score_table$age), score_table$vgh_id)
base <- score_table[ord, , drop = FALSE]
base <- base[!duplicated(base$subject_id), , drop = FALSE]
ps <- base[base$group %in% c("P-MCI", "S-MCI"), , drop = FALSE]
ps$P_MCI <- as.integer(ps$group == "P-MCI")
ps$age_filled <- ifelse(is.na(ps$age), median(ps$age, na.rm = TRUE), ps$age)
ps$age_c <- ps$age_filled - mean(ps$age_filled)
ps$male <- as.integer(ps$sex == "M")
model <- lm(score ~ age_c + male + P_MCI, data = ps)
coef_row <- summary(model)$coefficients["P_MCI", , drop = FALSE]
ci <- confint(model, "P_MCI", level = 0.95)
contrast <- data.frame(dataset = "GSE282742", contrast = "P_MCI_vs_S_MCI",
                       mapping = paste0(sum(panel$mapped_in_matrix), "/", nrow(panel)),
                       mapping_fraction = mean(panel$mapped_in_matrix),
                       estimate = coef_row[1, "Estimate"], se = coef_row[1, "Std. Error"],
                       ci_low = ci[1], ci_high = ci[2], p_value = coef_row[1, "Pr(>|t|)"],
                       n_subjects = nrow(ps), n_P_MCI = sum(ps$P_MCI == 1), n_S_MCI = sum(ps$P_MCI == 0))
write.csv(contrast, file.path(OUT, "GSE282742_12_panel_alias_aware_progression_contrast.csv"), row.names = FALSE, na = "")

split_by_subject <- split(score_table, score_table$subject_id)
paired <- do.call(rbind, lapply(names(split_by_subject), function(subject) {
  dat <- split_by_subject[[subject]]
  p <- dat[dat$group == "P-MCI", , drop = FALSE]
  a <- dat[dat$group == "AD", , drop = FALSE]
  if (nrow(p) == 0 || nrow(a) == 0) return(NULL)
  p <- p[order(ifelse(is.na(p$age), Inf, p$age), p$vgh_id), , drop = FALSE][1, ]
  a <- a[order(ifelse(is.na(a$age), Inf, a$age), a$vgh_id), , drop = FALSE]
  a <- a[nrow(a), ]
  data.frame(subject_id = subject, P_MCI_sample = p$vgh_id, AD_sample = a$vgh_id,
             P_MCI_score = p$score, AD_score = a$score,
             delta_AD_minus_P_MCI = a$score - p$score, stringsAsFactors = FALSE)
}))
paired_test <- t.test(paired$AD_score, paired$P_MCI_score, paired = TRUE)
paired_summary <- data.frame(dataset = "GSE282742", n_transition_subjects = nrow(paired),
                             mean_delta_AD_minus_P_MCI = mean(paired$delta_AD_minus_P_MCI),
                             p_value = paired_test$p.value)
write.csv(paired, file.path(OUT, "GSE282742_12_panel_alias_aware_paired_transitions.csv"), row.names = FALSE, na = "")
write.csv(paired_summary, file.path(OUT, "GSE282742_12_panel_alias_aware_paired_summary.csv"), row.names = FALSE, na = "")

group_counts <- table(sample_audit$group)
audit_lines <- c(
  "GSE282742 mapping reconciliation audit v2",
  "Generated: 2026-07-28",
  paste0("n_samples=", nrow(sample_audit)),
  paste0("sample_group_counts=P-MCI:", unname(group_counts["P-MCI"]), ";S-MCI:", unname(group_counts["S-MCI"]), ";AD:", unname(group_counts["AD"])),
  paste0("n_unique_subjects=", length(unique(sample_audit$subject_id))),
  paste0("n_subjects_with_repeated_sample_rows=", sum(table(sample_audit$subject_id) > 1)),
  paste0("n_repeated_sample_rows=", sum(table(sample_audit$subject_id)[table(sample_audit$subject_id) > 1])),
  paste0("n_baseline_P_or_S_subjects=", length(unique(sample_audit$subject_id[sample_audit$group %in% c("P-MCI", "S-MCI")]))),
  paste0("n_P_MCI_to_AD_paired_subjects=", nrow(paired)),
  paste0("n_95_programme_mapped_alias_aware=", sum(programme$mapped_in_matrix), "/", nrow(programme)),
  paste0("n_12_panel_mapped_alias_aware=", sum(panel$mapped_in_matrix), "/", nrow(panel)),
  "old_1_of_12_audit_status=SUPERSEDED; it used the older strict lookup table.",
  "old_15_of_95_audit_status=SUPERSEDED; it was produced by an incompatible preliminary mapping table and is not used for the locked score.",
  paste0("mapping_only_sanity_check_not_primary_estimate=", contrast$estimate, ";p=", contrast$p_value),
  paste0("mapping_only_sanity_check_not_primary_paired_delta=", paired_summary$mean_delta_AD_minus_P_MCI, ";p=", paired_summary$p_value),
  "primary_locked_score_result=estimate -0.255; 95% CI -0.637 to 0.127; HC3 P=0.191; paired mean delta -0.108; paired P=0.558 (see mci_key_analyses_20260723_091500)."
)
writeLines(audit_lines, file.path(OUT, "GSE282742_mapping_reconciliation_audit_v2.txt"), useBytes = TRUE)
cat(paste(audit_lines, collapse = "\n"), "\n")
