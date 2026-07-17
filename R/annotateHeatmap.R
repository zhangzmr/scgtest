#' Annotate Heatmap with Gene Expression Data
#'
#' Creates an annotated heatmap using dittoSeq for visualizing gene expression
#' patterns across samples with optional additional column annotations.
#'
#' @param bulk A SingleCellExperiment or similar object containing expression data
#' @param genes Character vector of gene names to include in the heatmap
#' @param colunm Optional data frame or matrix with additional column annotations to add to the object
#' @param ... Additional arguments passed to dittoHeatmap function
#'
#' @return A heatmap object from dittoHeatmap
#'
#' @examples
#' \dontrun{
#' heatmap_plot <- annotateHeatmap(bulk_data, genes = c("CD3D", "CD4", "CD8A"))
#' }
#'
#' @export
annotateHeatmap <- function(bulk, genes, colunm = NULL, ...) {
  intergene <- intersect(genes, rownames(bulk))
  if (!is.null(colunm)) {
    colData(bulk) <- cbind(colData(bulk), colunm)
  }
  dittoSeq::dittoHeatmap(bulk, intergene, ...)
}
