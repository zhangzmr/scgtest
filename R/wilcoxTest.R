#' Pairwise Wilcoxon Rank-Sum Tests for Single-Cell Data
#'
#' Performs pairwise Wilcoxon rank-sum tests between groups and calculates
#' log fold changes using t-tests for single-cell expression data.
#'
#' @param sce A SingleCellExperiment object containing expression data
#' @param group_vector A factor or character vector specifying group assignments for each cell
#' @param pairs Character vector specifying which pairwise comparisons to extract (e.g., "Group1--Group2")
#' @param ... Additional arguments passed to `scran::pairwiseWilcox` and `scran::pairwiseTTests`
#'
#' @return A named list where each element corresponds to a pairwise comparison and contains:
#'   genes, logFC, pvalue, and FDR
#'
#' @examples
#' \dontrun{
#' comparisons <- c("CD4--CD8", "CD4--NK")
#' results <- wilcoxTest(sce_object, group_vector = cell_types, pairs = comparisons)
#' }
#'
#' @export
wilcoxTest <- function(sce, group_vector, pairs, ...) {
  g1 <- sapply(strsplit(pairs, "--"), `[`, 1)
  g2 <- sapply(strsplit(pairs, "--"), `[`, 2)
  unique_groups <- unique(c(g1, g2))
  message("Perform Wilcoxon test...")
  marker_wil <- scran::pairwiseWilcox(
    sce, group_vector, assay.type = "logcounts",
    restrict = unique_groups, ...
  )
  all_pairs <- paste0(marker_wil$pairs$first, "--", marker_wil$pairs$second)
  marker_wil <- marker_wil$statistics
  names(marker_wil) <- all_pairs
  mucl_wil <- marker_wil[pairs]

  message("Perform t test for LogFC calculation...")
  marker_t <- scran::pairwiseTTests(
    sce, group_vector, assay.type = "logcounts",
    restrict = unique_groups, ...
  )
  all_pairs <- paste0(marker_t$pairs$first, "--", marker_t$pairs$second)
  marker_t <- marker_t$statistics
  names(marker_t) <- all_pairs
  mucl_t <- marker_t[pairs]

  message("Making result tables...")
  combined <- list()
  for (p in pairs) {
    m1 <- mucl_wil[[p]]
    m2 <- mucl_t[[p]]
    top <- data.frame(
      genes = rownames(m1), logFC = m2$logFC,
      pvalue = m1$p.value, FDR = m1$FDR
    )
    top <- top[order(top$FDR), ]
    combined[[p]] <- top
  }
  combined
}
