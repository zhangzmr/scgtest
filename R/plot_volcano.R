#' Create Volcano Plot for Differential Expression Results
#'
#' Generates an enhanced volcano plot to visualize differential expression results
#' with labeled top up- and down-regulated genes.
#'
#' @param top A data frame containing differential expression results with columns:
#'   genes, logFC, and FDR
#' @param title Character string for the plot title
#' @param n Integer specifying the number of top up- and down-regulated genes to label (default: 20)
#'
#' @return A ggplot object representing the volcano plot
#'
#' @examples
#' \dontrun{
#' volcano <- plot_volcano(de_results, title = "Treatment vs Control")
#' }
#'
#' @importFrom ggplot2 theme_classic ggtitle theme element_text element_blank xlim ylim
#' @export
plot_volcano <- function(top, title, n = 20) {
  top$FDR[top$FDR %in% 0] <- 1e-300
  top$FDR[top$FDR < 1e-300] <- 1e-300
  topSign <- top[top$FDR < 0.05, ]
  up <- topSign[topSign$logFC > 0, ]
  up <- up[order(-up$logFC), ]
  down <- topSign[topSign$logFC < 0, ]
  down <- down[order(down$logFC), ]
  top20 <- c(up$genes[seq_len(min(n, nrow(up)))],
             down$genes[seq_len(min(n, nrow(down)))])

  EnhancedVolcano::EnhancedVolcano(
    top, lab = top$genes, x = "logFC", y = "FDR",
    selectLab = top20, xlab = bquote(~Log[2] ~ "fold change"), pCutoff = 5e-2,
    FCcutoff = 0.2, pointSize = 1.0, labSize = 4, labCol = "black",
    labFace = "bold", boxedLabels = TRUE, colAlpha = 4 / 5,
    legendPosition = "right", legendLabSize = 14, legendIconSize = 4.0,
    drawConnectors = TRUE, widthConnectors = 1.0, subtitle = NULL,
    colConnectors = "black", max.overlaps = 60
  ) +
    theme_classic() +
    ggtitle(title) +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      panel.grid = element_blank(),
      title = element_text(face = "bold"),
      legend.position = "none"
    ) +
    xlim(min(top$logFC) - 0.3, max(top$logFC) + 0.3) +
    ylim(0, max(-log10(top$FDR)) + 5)
}
