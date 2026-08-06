# scgtest

R package for differential expression testing on single-cell and pseudobulk RNA-seq data: edgeR DE tests, optional RUV correction, Wilcoxon pairwise tests, pathway enrichment, and visualization helpers.

## Citation

If you use scgtest, please cite:

Xiaohao Wu, Mengrui Zhang, Jae-Seung Moon, Eun Kyung Song, Joshua C. Abrams, Laura S. van Dam, Orr Sharpe, Marie Feser, Laurie Moss, Peggy P. Ho, Melanie H. Smith, Laura T. Donlin, Jane H. Buckner, Eddie A. James, Gary S. Firestein, Yuko Okamoto, Tobiaz V. Lanz, Eric Meffre, V. Michael Holers, Kevin D. Deane, William H. Robinson. Senescent Activated Naive B Cells Promote Anti-Citrullinated Antigen T Cell Responses and the Transition to Clinical Rheumatoid Arthritis. bioRxiv. 2025. doi:[10.1101/2025.10.14.682430](https://doi.org/10.1101/2025.10.14.682430).

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
