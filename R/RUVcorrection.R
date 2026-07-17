#' RUV Correction for RNA-seq Data
#'
#' Performs Remove Unwanted Variation (RUV) correction using RUVSeq.
#' Identifies empirical genes and applies RUVg normalization to remove
#' unwanted technical variation from count data.
#'
#' @param counts A numeric matrix of raw count data with genes as rows and samples as columns
#' @param metadata A data frame containing sample metadata and group information
#' @param design A design matrix for the linear model used to identify empirical genes
#' @param group Column name in metadata specifying the grouping variable (default: "group")
#' @param contrasts Character string specifying the contrast to test for empirical gene selection (default: "T1D-HC")
#' @param k Integer specifying the number of factors of unwanted variation to remove (default: 1)
#' @param nGenes Integer specifying the number of top differential genes to exclude when selecting empirical genes (default: 5000)
#'
#' @return A RUVSeq object containing the normalized count matrix and estimated factors of unwanted variation
#'
#' @examples
#' \dontrun{
#' ruv_result <- RUVcorrection(
#'   counts = count_matrix, metadata = sample_info,
#'   design = design_matrix, group = "treatment"
#' )
#' }
#'
#' @export
RUVcorrection <- function(counts, metadata, design, group = "group",
                          contrasts = "T1D-HC", k = 1, nGenes = 5000) {
  message("Finding empirical genes using edgeR...")
  y <- edgeR::DGEList(counts = counts, samples = metadata, group = metadata[[group]])
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateGLMCommonDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)
  fit <- edgeR::glmFit(y, design)

  message(paste0("contrast group is ", contrasts))
  contrast <- limma::makeContrasts(contrasts = contrasts, levels = design)
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  top <- data.frame(edgeR::topTags(lrt, n = nrow(y)))
  message("Perform RUV correction...")
  empirical <- rownames(counts)[which(!(rownames(counts) %in% rownames(top)[1:nGenes]))]
  set2 <- RUVSeq::RUVg(counts, empirical, k = k)
  return(set2)
}
