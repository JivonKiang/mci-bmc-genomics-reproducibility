#!/usr/bin/env Rscript

# Rebuild the single-cell virtual-knockout screen from the current locked
# 12-gene MCI panel. This script deliberately does not read any legacy
# RPS27L/EEF1G/RPL17 result object.

local_lib <- "E:/20241004_MCI/20260709/Rlib_20260729"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(scTenifoldKnk)
})

project_root <- "E:/20241004_MCI"
data_root <- file.path(project_root, "20250816 revise/scRNA/GSE134578_RAW")
panel_root <- file.path(project_root, "20260709/mci_development_optimized_panel_20260722_161919")
default_out <- file.path(project_root, "20260709/mci_locked_panel_scTenifoldKnk_20260729_093000")
out_root <- Sys.getenv("MCI_KNK_OUT_ROOT", unset = default_out)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

seurat_path <- file.path(data_root, "annotated.rds")
panel_path <- file.path(panel_root, "development_locked_12_gene_panel.csv")
stopifnot(file.exists(seurat_path), file.exists(panel_path))

panel <- read.csv(panel_path, stringsAsFactors = FALSE, check.names = FALSE)
required_panel_cols <- c("gene", "trajectory_class")
stopifnot(all(required_panel_cols %in% colnames(panel)))
panel$gene <- as.character(panel$gene)

sc <- readRDS(seurat_path)
meta <- sc[[]]
stopifnot(all(c("disease_group", "sample_id", "integrated_annotation") %in% colnames(meta)))
counts <- SeuratObject::LayerData(sc, assay = "RNA", layer = "counts")
gene_names <- rownames(counts)
screen_mode <- tolower(Sys.getenv("MCI_KNK_MODE", unset = "screen"))
default_max_features <- if (screen_mode == "high") 3000L else 1500L
max_features <- as.integer(Sys.getenv("MCI_KNK_MAX_FEATURES", unset = as.character(default_max_features)))
default_nnet <- if (screen_mode == "high") 10L else 3L
default_ncells <- if (screen_mode == "high") 300L else 100L
default_lsize <- if (screen_mode == "high") 1000L else 500L
default_max_iter <- if (screen_mode == "high") 1000L else 300L

map_symbol <- function(g) {
  hit <- gene_names[toupper(gene_names) == toupper(g)]
  if (length(hit) == 1L) hit else NA_character_
}
panel$scRNA_symbol <- vapply(panel$gene, map_symbol, character(1))
panel$mapping_status <- ifelse(is.na(panel$scRNA_symbol), "not_found_in_scRNA", "exact_symbol")
panel$legacy_excluded <- TRUE
write.csv(panel, file.path(out_root, "panel_gene_mapping.csv"), row.names = FALSE, na = "")

writeLines(c(
  "MCI locked-panel scTenifoldKnk screen",
  "Generated: 2026-07-29",
  "",
  "This run was rebuilt from the current development_locked_12_gene_panel.csv.",
  "Legacy RPS27L, EEF1G and RPL17 perturbation files were not read and are excluded from all outputs.",
  "",
  "Primary context: GSE134578 annotated Seurat object, disease-stratified CD8 T-cell states.",
  "The network perturbation is computational and exploratory; pooled cells are not donor-level replicates.",
  "MCI has two donors in this object, so no clinical or causal claim is made from this screen.",
  "PHF15 and LAT1-3TM were not present as exact symbols in the matrix and were not replaced by aliases.",
  "",
  paste0("Method: scTenifoldKnk 1.0.3, mode=", screen_mode,
         "; nNet=", Sys.getenv("MCI_KNK_NNET", unset = as.character(default_nnet)),
         "; nCells=", Sys.getenv("MCI_KNK_NCELLS", unset = as.character(default_ncells)),
         "; max features=", max_features, "; mapped panel genes are forced into the input."),
  "Target-specific outputs include the target gene, context, cell count, donor count, mapped symbol, effect table and adjusted P values."
), file.path(out_root, "README.md"))

panel_mapped <- panel$scRNA_symbol[!is.na(panel$scRNA_symbol)]
variable_genes <- Seurat::VariableFeatures(sc)
if (length(variable_genes) == 0L) stop("No stored variable features in annotated Seurat object")
variable_genes <- intersect(variable_genes, gene_names)
variable_genes <- head(variable_genes, max_features)
analysis_genes <- unique(c(variable_genes, panel_mapped))

cd8_states <- c(
  "Central memory CD8 T cells",
  "Effector memory CD8 T cells",
  "MAIT cells",
  "Naive CD8 T cells",
  "Terminal effector CD8 T cells"
)

make_stratified_index <- function(index, max_cells = 3000L, seed = 20260729L) {
  if (length(index) <= max_cells) return(index)
  set.seed(seed)
  by_sample <- split(index, meta$sample_id[index])
  prop <- vapply(by_sample, length, numeric(1)) / length(index)
  take <- pmax(1L, floor(prop * max_cells))
  while (sum(take) > max_cells) {
    j <- which.max(take)
    if (take[j] > 1L) take[j] <- take[j] - 1L else break
  }
  while (sum(take) < max_cells) {
    room <- which(take < vapply(by_sample, length, integer(1)))
    if (!length(room)) break
    take[room[1]] <- take[room[1]] + 1L
  }
  unlist(Map(function(x, n) sample(x, n), by_sample, take), use.names = FALSE)
}

contexts <- list(
  MCI_CD8_states = which(meta$disease_group == "MCI" & meta$integrated_annotation %in% cd8_states),
  HC_CD8_states = which(meta$disease_group == "HC" & meta$integrated_annotation %in% cd8_states),
  AD_CD8_states = which(meta$disease_group == "AD" & meta$integrated_annotation %in% cd8_states)
)

target_filter <- Sys.getenv("MCI_KNK_TARGETS", unset = "")
targets <- panel_mapped
if (nzchar(target_filter)) targets <- intersect(targets, trimws(strsplit(target_filter, ",", fixed = TRUE)[[1]]))
context_filter <- Sys.getenv("MCI_KNK_CONTEXTS", unset = "")
if (nzchar(context_filter)) contexts <- contexts[intersect(names(contexts), trimws(strsplit(context_filter, ",", fixed = TRUE)[[1]]))]

knk_params <- list(
  qc = TRUE,
  qc_mtThreshold = 0.10,
  qc_minLSize = as.integer(Sys.getenv("MCI_KNK_MIN_LSIZE", unset = as.character(default_lsize))),
  qc_minCells = 25,
  nc_lambda = 0,
  nc_nNet = as.integer(Sys.getenv("MCI_KNK_NNET", unset = as.character(default_nnet))),
  nc_nCells = as.integer(Sys.getenv("MCI_KNK_NCELLS", unset = as.character(default_ncells))),
  nc_nComp = 3,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q = 0.90,
  td_K = 3,
  td_maxIter = as.integer(Sys.getenv("MCI_KNK_MAX_ITER", unset = as.character(default_max_iter))),
  td_maxError = 1e-5,
  td_nDecimal = 3,
  ma_nDim = 2,
  nCores = max(1L, min(4L, parallel::detectCores(logical = FALSE)))
)

manifest <- list()
all_effects <- list()
all_summary <- list()

standardize_result <- function(result, target, context_name, n_cells, n_donors, sampled_cells) {
  d <- as.data.frame(result$diffRegulation, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(d)) d$gene <- rownames(d)
  d$target_gene <- target
  d$context <- context_name
  d$n_input_cells <- n_cells
  d$n_sampled_cells <- sampled_cells
  d$n_donors <- n_donors
  d$significant <- !is.na(d$p.adj) & d$p.adj < 0.05 & abs(d$Z) >= 2
  d$log2FC <- log2(pmax(as.numeric(d$FC), .Machine$double.xmin))
  d[order(d$p.adj, -abs(d$Z)), , drop = FALSE]
}

for (context_name in names(contexts)) {
  raw_index <- contexts[[context_name]]
  if (length(raw_index) == 0L) next
  selected_index <- make_stratified_index(raw_index, max_cells = 3000L,
                                           seed = 20260729L + match(context_name, names(contexts)))
  context_dir <- file.path(out_root, context_name)
  dir.create(context_dir, recursive = TRUE, showWarnings = FALSE)
  context_counts <- counts[analysis_genes, selected_index, drop = FALSE]
  sample_ids <- unique(meta$sample_id[selected_index])
  context_info <- data.frame(
    context = context_name,
    disease_group = unique(meta$disease_group[selected_index]),
    input_cells = length(raw_index),
    sampled_cells = length(selected_index),
    donors = length(sample_ids),
    sample_ids = paste(sample_ids, collapse = ";"),
    stringsAsFactors = FALSE
  )
  write.csv(context_info, file.path(context_dir, "context_metadata.csv"), row.names = FALSE)

  for (target in targets) {
    target_dir <- file.path(context_dir, target)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    seed <- 20260729L + match(target, targets) + 100L * match(context_name, names(contexts))
    set.seed(seed)
    target_ok <- target %in% rownames(context_counts)
    manifest_row <- data.frame(
      context = context_name,
      target_gene = target,
      mapped_symbol = target,
      input_cells = length(raw_index),
      sampled_cells = length(selected_index),
      donors = length(sample_ids),
      seed = seed,
      status = ifelse(target_ok, "pending", "target_not_in_matrix"),
      error = "",
      stringsAsFactors = FALSE
    )
    if (!target_ok) {
      manifest[[length(manifest) + 1L]] <- manifest_row
      next
    }
    result <- tryCatch(
      do.call(scTenifoldKnk::scTenifoldKnk, c(list(countMatrix = context_counts, gKO = target), knk_params)),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      manifest_row$status <- "error"
      manifest_row$error <- conditionMessage(result)
      manifest[[length(manifest) + 1L]] <- manifest_row
      writeLines(conditionMessage(result), file.path(target_dir, "error.txt"))
      next
    }
    saveRDS(result, file.path(target_dir, paste0(target, "_scTenifoldKnk_result.rds")))
    effects <- standardize_result(result, target, context_name, length(raw_index), length(sample_ids), length(selected_index))
    write.csv(effects, file.path(target_dir, paste0(target, "_diffRegulation.csv")), row.names = FALSE)
    sig <- effects[effects$significant, , drop = FALSE]
    write.csv(sig, file.path(target_dir, paste0(target, "_significant.csv")), row.names = FALSE)
    all_effects[[length(all_effects) + 1L]] <- effects
    all_summary[[length(all_summary) + 1L]] <- data.frame(
      context = context_name,
      target_gene = target,
      input_cells = length(raw_index),
      sampled_cells = length(selected_index),
      donors = length(sample_ids),
      significant_genes = nrow(sig),
      top_perturbed_gene = if (nrow(effects)) effects$gene[1] else NA_character_,
      top_adjusted_p = if (nrow(effects)) effects$p.adj[1] else NA_real_,
      top_Z = if (nrow(effects)) effects$Z[1] else NA_real_,
      status = "completed",
      stringsAsFactors = FALSE
    )
    manifest_row$status <- "completed"
    manifest[[length(manifest) + 1L]] <- manifest_row
    rm(result, effects, sig)
    gc(verbose = FALSE)
  }
}

manifest_df <- if (length(manifest)) do.call(rbind, manifest) else data.frame()
write.csv(manifest_df, file.path(out_root, "run_manifest.csv"), row.names = FALSE)
summary_df <- if (length(all_summary)) do.call(rbind, all_summary) else data.frame()
write.csv(summary_df, file.path(out_root, "panel_virtual_knockout_summary.csv"), row.names = FALSE)
effects_df <- if (length(all_effects)) do.call(rbind, all_effects) else data.frame()
write.csv(effects_df, file.path(out_root, "panel_virtual_knockout_all_genes.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
writeLines(c(
  "Completed current-panel scTenifoldKnk screen.",
  paste("Contexts:", paste(names(contexts), collapse = ", ")),
  paste("Targets:", paste(targets, collapse = ", ")),
  paste("Output:", out_root)
), file.path(out_root, "completion.txt"))

message("Completed locked-panel scTenifoldKnk screen: ", out_root)
