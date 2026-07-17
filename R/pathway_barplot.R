#' Create Pathway Enrichment Bar Plot
#'
#' Generates a horizontal bar plot for pathway enrichment results showing
#' enrichment scores colored by adjusted p-values.
#'
#' @param res A data frame containing pathway enrichment results with columns:
#'   Description, enrichmentScore, and p.adjust
#' @param pathway Character string for the y-axis label describing the pathway database
#' @param title Character string for the plot title
#'
#' @return A ggplot object representing the pathway enrichment bar plot
#'
#' @examples
#' \dontrun{
#' kegg_plot <- pathway_barplot(
#'   kegg_results, pathway = "KEGG Pathways",
#'   title = "KEGG Enrichment Analysis"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_classic scale_fill_viridis_c labs theme element_text ggtitle
#' @importFrom stringr str_to_title str_replace_all
#' @export
pathway_barplot <- function(res, pathway, title) {
  temp <- res[res$enrichmentScore > 0, ]
  temp <- temp[order(temp$p.adjust), ]

  temp2 <- res[res$enrichmentScore < 0, ]
  temp2 <- temp2[order(temp2$p.adjust), ]
  res <- rbind(temp, temp2)

  res$Description <- str_to_title(res$Description)
  res$Description <- str_replace_all(res$Description, "Of", "of")
  res$Description <- factor(res$Description, levels = rev(res$Description))

  ggplot(res, aes(x = enrichmentScore, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_viridis_c(direction = -1) +
    labs(x = "Enrichment Score", y = pathway, fill = "Adjusted P-value") +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    ) +
    ggtitle(title)
}
