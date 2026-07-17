#' @title Convert Single Cell to Bulk-like Data with Filtering
#' @description Convert single-cell data to bulk-like data by summing or averaging
#'   counts across groups, with optional gene filtering.
#' @param sce A SingleCellExperiment object containing single-cell data.
#' @param group A vector or factor indicating the grouping of cells.
#' @param min_counts Minimum count threshold for filtering genes. Default is 5.
#' @param min_sample Minimum number of samples that must exceed `min_counts`. Default is 2.
#' @param filter Logical indicating whether to apply filtering. Default is TRUE.
#' @param average Logical indicating whether to average counts instead of summing. Default is FALSE.
#' @param return_se Logical; if TRUE return the aggregated SummarizedExperiment/SCE, otherwise a list. Default is FALSE.
#' @param ... Additional arguments passed to `scuttle::sumCountsAcrossCells`.
#' @return If `return_se` is TRUE, the aggregated object; otherwise a list with `counts` and `metadata`.
#' @export
convertBulkFilter <- function(sce, group, min_counts = 5, min_sample = 2,
                              filter = TRUE, average = FALSE,
                              return_se = FALSE, ...) {
  summed <- scuttle::sumCountsAcrossCells(sce, ids = group, average = average, ...)
  if (isFALSE(average)) {
    message("Calculate by Sum...")
    counts <- assay(summed, "sum")
  } else {
    message("Calculate by Average...")
    counts <- assay(summed, "average")
  }

  if (isTRUE(filter)) {
    message("Running filtering...")
    keep <- apply(counts, 1, function(x) length(x[x > min_counts]) >= min_sample)
    print(table(keep))
    counts <- counts[keep, ]
  } else {
    message("No filtering performed..")
  }

  if (isTRUE(return_se)) {
    return(summed)
  }
  list(counts = counts, metadata = data.frame(colData(summed)))
}

#' @title Filter significant DE results
#' @description Keep rows with FDR below a threshold.
#' @param data A data frame with an `FDR` column.
#' @param pvalue FDR cutoff (default 0.05).
#' @return A filtered data frame.
#' @export
getTop <- function(data, pvalue = 0.05) {
  data[data$FDR < pvalue, ]
}
