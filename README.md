# corinet

<!-- badges: start -->
[![R-CMD-check](https://github.com/lmuir16/corinet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lmuir16/corinet/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/lmuir16/corinet/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/lmuir16/corinet/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**corinet** is a lightweight R package for building **COR**relation **NET**works 
to analyze gene co-expression in bulk RNA-seq data. Starting from a
`SummarizedExperiment` object, corinet supports the full workflow from
data preparation through network construction, module detection,
visualization, and export.

## Repository Structure

```
corinet/
├── r-package/    # R package source
└── nextflow/     # Nextflow pipeline for large-scale analysis
└── results/      # Final output files from example Nextflow run
```
Both the `r-package`  and `nextflow` subdirectories contain their own dedicated README files. 

- [r-package/README.md](r-package/README.md)
- [nextflow/README.md](nextflow/README.md)

## Features

- Select top variable genes for network construction
- Compute pairwise Pearson or Spearman gene-gene correlations
- Threshold correlations into a binary adjacency matrix
- Detect gene modules using Louvain or walktrap community detection
- Identify hub genes by degree and betweenness centrality
- Visualize networks as correlation heatmaps or force-directed graphs
- Export results to TSV files
- Command line interface via Rapp

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "recount3"))

pak::pak("lmuir16/corinet")
```

## Documentation

Full documentation and a getting-started vignette are available at the
[corinet pkgdown site](https://lmuir16.github.io/corinet/).

## License

MIT © Liane Muir
