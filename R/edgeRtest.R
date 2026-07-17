#' EdgeR Differential Expression Testing
#'
#' Performs differential expression analysis using edgeR with support for
#' likelihood ratio tests (LRT) and quasi-likelihood F-tests (QLF).
#'
#' @param counts A numeric matrix of raw count data with genes as rows and samples as columns
#' @param metadata A data frame containing sample metadata and group information
#' @param group Column name in metadata specifying the grouping variable (default: "group")
#' @param design A design matrix for the linear model
#' @param contrasts Character string specifying the contrast to test (default: "T1D-HC")
#' @param method Character string specifying the edgeR method: "lrt" or "qlf" (default: "lrt")
#' @param normalize Character string specifying normalization method (default: "TMM")
#' @param filter Logical indicating whether to filter lowly expressed genes using filterByExpr (default: TRUE)
#' @param model_out Logical indicating whether to return the fitted model instead of results (default: FALSE)
#'
#' @return If `model_out = FALSE`: a list with `top` and `topSign`.
#'   If `model_out = TRUE`: the fitted edgeR model object.
#'
#' @examples
#' \dontrun{
#' result <- edgeRtest(
#'   counts = count_matrix, metadata = sample_info,
#'   design = design_matrix, contrasts = "Treatment-Control"
#' )
#' }
#'
#' @export
edgeRtest <- function(counts, metadata, group = "group", design,
                      contrasts = "T1D-HC", method = "lrt",
                      normalize = "TMM", filter = TRUE, model_out = FALSE) {
  message("Create a DGE object...")
  dge <- edgeR::DGEList(counts = counts, samples = metadata, group = metadata[[group]])

  if (isTRUE(filter)) {
    message("Further filtering...")
    keep <- edgeR::filterByExpr(dge, group = metadata[[group]])
    print(table(keep))
    dge <- dge[keep, , keep.lib.sizes = FALSE]
  } else {
    message("No gene filtering...")
  }

  dge <- edgeR::calcNormFactors(dge, method = normalize)
  dge <- edgeR::estimateDisp(dge, design)

  contrast <- limma::makeContrasts(contrasts = contrasts, levels = design)
  print(contrast)
  if (method == "lrt") {
    message(paste0("edgeR method is ", method))
    fit <- edgeR::glmFit(dge, design)
    res_fit <- edgeR::glmLRT(fit, contrast = contrast)
  } else if (method == "qlf") {
    message(paste0("edgeR method is ", method))
    fit <- edgeR::glmQLFit(dge, design)
    res_fit <- edgeR::glmQLFTest(fit, contrast = contrast)
  } else {
    stop("method must be \"lrt\" or \"qlf\".")
  }

  if (isTRUE(model_out)) {
    message("Output fitted model..")
    return(res_fit)
  }

  message("Getting top genes..")
  top <- data.frame(edgeR::topTags(res_fit, adjust.method = "BH", n = NULL))
  topSign <- top[top$FDR < 0.05, ]
  list(top = top, topSign = topSign)
}
