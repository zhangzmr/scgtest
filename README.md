# scgtest

R package for differential expression testing on single-cell and pseudobulk RNA-seq data: edgeR DE tests, optional RUV correction, Wilcoxon pairwise tests, pathway enrichment, and visualization helpers.

## Installation

```r
# install.packages("remotes")
remotes::install_github("zhangzmr/scgtest")
```

For local development from a clone:

```r
# install.packages("devtools")
devtools::install("path/to/scgtest")
```

## Example workflow

```r
library(scgtest)

# bulk <- convertBulkFilter(sce, group = sce$celltype)
# result <- fullDEtest(
#   counts = bulk$counts,
#   metadata = bulk$metadata,
#   design = model.matrix(~ 0 + group, data = bulk$metadata),
#   group = "group",
#   contrasts = "Treatment-Control",
#   RUV = TRUE
# )
# plot_volcano(result$top, title = "Treatment vs Control")
# kegg <- enrich_pathway(result$topSign, database = "KEGG")
# pathway_barplot(kegg, pathway = "KEGG", title = "Pathway enrichment")
```

## License

MIT © Mengrui Zhang
