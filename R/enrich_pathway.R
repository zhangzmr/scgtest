#' Pathway Enrichment Analysis
#'
#' Performs gene set enrichment analysis (GSEA) using Gene Ontology (GO) or
#' KEGG pathway databases on differential expression results.
#'
#' @param marker A data frame containing differential expression results with columns:
#'   genes (gene symbols) and logFC (log fold changes)
#' @param database Character string specifying the database: "GO" or "KEGG"
#'
#' @return A data frame containing enrichment results from clusterProfiler
#'
#' @examples
#' \dontrun{
#' go_results <- enrich_pathway(de_results, database = "GO")
#' kegg_results <- enrich_pathway(de_results, database = "KEGG")
#' }
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit
#' @export
enrich_pathway <- function(marker, database) {
  marker <- marker[order(-marker$logFC), ]
  marker$ensid <- mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = marker$genes, column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"
  )
  marker <- na.omit(marker)
  genelist <- marker$logFC
  names(genelist) <- marker$ensid
  if (database == "GO") {
    kk2 <- clusterProfiler::gseGO(
      geneList = genelist, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      pvalueCutoff = 0.05, verbose = TRUE
    )
  } else if (database == "KEGG") {
    kk2 <- clusterProfiler::gseKEGG(
      geneList = genelist, organism = "hsa",
      pvalueCutoff = 0.05, verbose = TRUE
    )
  } else {
    stop("database must be \"GO\" or \"KEGG\".")
  }
  kk2@result
}
