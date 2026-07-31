suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readr)
  library(dplyr)
})

lib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", "4.3")
.libPaths(c(lib, .libPaths()))

base_dir <- "E:/20241004_MCI/20260709"
input_dir <- file.path(base_dir, "mci_spi1_standard_mr_final_20260727_115216")
input_file <- file.path(input_dir, "all_harmonised_rows_combined.csv")
out_dir <- file.path(base_dir, paste0("mci_spi1_twosamplemr_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read_csv(input_file, show_col_types = FALSE)
exposure_raw <- raw %>%
  transmute(
    SNP = rsid,
    beta = beta_exp,
    se = se_exp,
    effect_allele = effect_allele,
    other_allele = other_allele,
    pval = pval_exp,
    id.exposure = "SPI1::eQTLGen",
    exposure = "SPI1 eQTLGen cis-eQTL"
  ) %>% distinct()

outcome_raw <- raw %>%
  transmute(
    SNP = rsid,
    beta = beta_out,
    se = se_out,
    effect_allele = outcome_alt,
    other_allele = outcome_ref,
    pval = pval_out,
    eaf = as.numeric(outcome_af_alt),
    id.outcome = outcome,
    outcome = outcome_class
  ) %>% distinct()

exposure_dat <- format_data(
  exposure_raw, type = "exposure", snp_col = "SNP", beta_col = "beta",
  se_col = "se", effect_allele_col = "effect_allele",
  other_allele_col = "other_allele", pval_col = "pval",
  id_col = "id.exposure", phenotype_col = "exposure"
)
outcome_dat <- format_data(
  outcome_raw, type = "outcome", snp_col = "SNP", beta_col = "beta",
  se_col = "se", effect_allele_col = "effect_allele",
  other_allele_col = "other_allele", pval_col = "pval", eaf_col = "eaf",
  id_col = "id.outcome", phenotype_col = "outcome"
)

harmonised <- harmonise_data(exposure_dat, outcome_dat, action = 2)
write_csv(harmonised, file.path(out_dir, "twosamplemr_harmonised.csv"))

methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_wald_ratio")
mr_results <- mr(harmonised, method_list = methods) %>%
  mutate(
    OR = exp(b),
    OR_lci95 = exp(b - 1.96 * se),
    OR_uci95 = exp(b + 1.96 * se),
    p_bonferroni_5_outcomes = NA_real_
  )
ivw_rows <- mr_results$method == "Inverse variance weighted"
mr_results$p_bonferroni_5_outcomes[ivw_rows] <- p.adjust(mr_results$pval[ivw_rows], method = "bonferroni")
write_csv(mr_results, file.path(out_dir, "twosamplemr_results.csv"))

heterogeneity <- mr_heterogeneity(harmonised, method_list = c("mr_ivw", "mr_egger_regression"))
write_csv(heterogeneity, file.path(out_dir, "twosamplemr_heterogeneity.csv"))

pleiotropy <- mr_pleiotropy_test(harmonised)
write_csv(pleiotropy, file.path(out_dir, "twosamplemr_egger_intercept.csv"))

loo <- mr_leaveoneout(harmonised, method = mr_ivw)
write_csv(loo, file.path(out_dir, "twosamplemr_leaveoneout.csv"))

writeLines(c(
  paste0("# TwoSampleMR 0.7.9 verification (", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ")"),
  "",
  paste0("Input: ", normalizePath(input_file, winslash = "/", mustWork = FALSE)),
  "Exposure: SPI1 eQTLGen; three independent cis-eQTL instruments.",
  "Outcome rows were already allele-harmonised by the frozen local pipeline and were re-harmonised with TwoSampleMR action=2 for verification.",
  "IVW is the primary multi-SNP estimator. Weighted median, MR-Egger, heterogeneity, intercept and leave-one-out are sensitivity analyses.",
  "Single-SNP Wald ratios are retained as supplementary estimates and are not interpreted as robust causal biomarker evidence.",
  "",
  paste0("TwoSampleMR version: ", as.character(packageVersion("TwoSampleMR"))),
  paste0("Output: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
), file.path(out_dir, "README.md"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Wrote", out_dir, "\n")
