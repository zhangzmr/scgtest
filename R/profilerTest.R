#' KEGG Pathway Analysis Using clusterProfiler
#'
#' Performs KEGG pathway enrichment analysis on differential expression results
#' using either over-representation analysis or gene set enrichment analysis (GSEA).
#'
#' @param top A data frame containing differential expression results with gene symbols as row names
#' @param logfc Character string specifying the column name containing log fold changes (default: "logFC")
#' @param method Character string specifying analysis method: "enrich" or "gse" (default: "gse")
#' @param lfc Numeric threshold for log fold change when using "enrich" method (default: 0.5)
#' @param cutoff Numeric p-value cutoff for significance (default: 0.05)
#' @param ... Additional arguments passed to gseKEGG or enrichKEGG
#'
#' @return A data frame containing pathway enrichment results
#'
#' @examples
#' \dontrun{
#' pathway_results <- profilerTest(de_results)
#' pathway_results <- profilerTest(de_results, method = "enrich", lfc = 1.0)
#' }
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit
#' @export
profilerTest <- function(top, logfc = "logFC", method = "gse",
                         lfc = 0.5, cutoff = 0.05, ...) {
  marker <- top
  marker <- marker[order(-marker[[logfc]]), ]
  message("Covert gene symbol to ENTREID...")
  marker$ensid <- mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(marker), column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"
  )
  marker <- na.omit(marker)

  genelist <- marker[[logfc]]
  names(genelist) <- marker$ensid

  if (method == "enrich") {
    message("clusterProfiler enrich test...")
    gene <- names(genelist)[abs(genelist) > lfc]
    kk <- clusterProfiler::enrichKEGG(
      gene = gene, organism = "hsa", pvalueCutoff = cutoff, ...
    )
    return(kk@result)
  }

  if (method == "gse") {
    message("clusterProfiler gse test...")
    kk2 <- clusterProfiler::gseKEGG(
      geneList = genelist, organism = "hsa",
      pvalueCutoff = cutoff, verbose = TRUE, ...
    )
    return(kk2@result)
  }

  stop("method must be \"enrich\" or \"gse\".")
}
