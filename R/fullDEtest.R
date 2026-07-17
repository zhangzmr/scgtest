#' Full Differential Expression Testing with Optional RUV Correction
#'
#' Performs differential expression testing using edgeR with optional RUVSeq
#' normalization for unwanted variation removal.
#'
#' @param counts A numeric matrix of raw count data with genes as rows and samples as columns
#' @param metadata A data frame containing sample metadata and group information
#' @param design A design matrix for the linear model
#' @param group Column name in metadata specifying the grouping variable (default: "group")
#' @param contrasts Character string specifying the contrast to test (default: "T1D-HC")
#' @param method Character string specifying the edgeR method: "lrt" or "qlf" (default: "lrt")
#' @param filter Logical indicating whether to filter lowly expressed genes (default: FALSE)
#' @param normalize Character string specifying normalization method (default: "TMM")
#' @param RUV Logical indicating whether to apply RUV correction (default: TRUE)
#' @param k Integer specifying the number of factors for RUV correction (default: 1)
#' @param nGenes Integer specifying the number of top genes to exclude when finding empirical genes for RUV (default: 5000)
#' @param model_out Logical indicating whether to return the fitted model instead of results (default: FALSE)
#'
#' @return If `model_out = FALSE`: a list with `top` and `topSign`.
#'   If `model_out = TRUE`: the fitted edgeR model object.
#'
#' @examples
#' \dontrun{
#' result <- fullDEtest(
#'   counts = count_matrix, metadata = sample_info,
#'   design = design_matrix, group = "condition"
#' )
#' }
#'
#' @export
fullDEtest <- function(counts, metadata, design, group = "group",
                       contrasts = "T1D-HC", method = "lrt", filter = FALSE,
                       normalize = "TMM", RUV = TRUE, k = 1, nGenes = 5000,
                       model_out = FALSE) {
  if (isTRUE(RUV)) {
    message("Running RUV correction...")
    set2 <- RUVcorrection(
      counts = counts, metadata = metadata, design = design,
      group = group, contrasts = contrasts, k = k, nGenes = nGenes
    )

    message("Add RUV factors to design matrix...")
    design <- cbind(design, set2$W)
  } else {
    message("NO RUV correction...")
  }

  message("Running differential expression test with edgeR...")
  edgeRtest(
    counts = counts, metadata = metadata, group = group,
    design = design, contrasts = contrasts, method = method,
    normalize = normalize, filter = filter, model_out = model_out
  )
}
