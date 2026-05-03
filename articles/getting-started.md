# Getting Started with corinet

## Introduction

*corinet* provides a lightweight pipeline for building and analyzing
gene co-expression networks from bulk RNA-seq data. Starting from a
`SummarizedExperiment` object, the package supports the full workflow
from data preparation through network construction, module detection,
visualization, and export.

Let’s see what the complete pipeline looks like…

## Installation

``` r

# Install from GitHub
devtools::install_github("lmuir16/corinet")
```

## Data Preparation

The package includes `example_se`, a `SummarizedExperiment` containing
GTEx skeletal muscle bulk RNA-seq data for 2000 genes across 200
samples.

``` r

library(corinet)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

data(example_se)
example_se
#> class: RangedSummarizedExperiment 
#> dim: 2000 200 
#> metadata(8): time_created recount3_version ... annotation recount3_url
#> assays(1): counts
#> rownames(2000): ENSG00000143632.14 ENSG00000198804.2 ...
#>   ENSG00000128578.9 ENSG00000174903.15
#> rowData names(10): source type ... havana_gene tag
#> colnames(200): GTEX-S3XE-2026-SM-3K2B5.1 GTEX-1MJIX-2126-SM-E9J3A.1 ...
#>   GTEX-1E2YA-0526-SM-7P8RY.1 GTEX-S3LF-0526-SM-EZ6LK.1
#> colData names(2): gtex.age gtex.sex
```

Before building a network, we select the most variable genes and prepare
a log-normalized expression matrix. Working with the top variable genes
keeps pairwise correlation computation tractable and focuses the network
on biologically informative signal.

``` r

# Select top 500 most variable genes
se_top <- top_variable_features(example_se, n = 500)

# Optionally subset to a random sample of donors for faster computation
set.seed(42)
se_sub <- select_random_samples(se_top, n = 200)

# Log2-normalize counts and remap Ensembl IDs to HGNC gene symbols
mat <- prepare_expression_matrix(se_sub)

dim(mat)
#> [1] 500 200
mat[1:5, 1:3]
#>         GTEX-1LGOU-2126-SM-DIPFG.1 GTEX-1RQEC-0126-SM-E9J2H.1
#> MT-CO1                    21.18882                   20.48683
#> ACTA1                     19.84421                   20.83015
#> MT-RNR2                   18.70732                   18.04020
#> TTN                       18.82295                   19.24694
#> MYH7                      19.03442                   18.98329
#>         GTEX-ZYFC-0526-SM-5GIDF.1
#> MT-CO1                   21.15531
#> ACTA1                    20.70645
#> MT-RNR2                  18.40462
#> TTN                      18.97047
#> MYH7                     19.80240
```

## Computing Gene-Gene Correlations

With the expression matrix prepared, we compute pairwise Pearson
correlations across all gene pairs. This produces a symmetric genes x
genes matrix where each entry reflects the co-expression strength
between two genes.

``` r

cor_mat <- pairwise_corr(mat, cor_method = "pearson")

dim(cor_mat)
#> [1] 500 500
cor_mat[1:5, 1:5]
#>              MT-CO1       ACTA1     MT-RNR2         TTN        MYH7
#> MT-CO1   1.00000000  0.38569940  0.02169463 -0.11184901  0.35278624
#> ACTA1    0.38569940  1.00000000 -0.21346768  0.08683613  0.76805953
#> MT-RNR2  0.02169463 -0.21346768  1.00000000 -0.18679021 -0.09860481
#> TTN     -0.11184901  0.08683613 -0.18679021  1.00000000  0.12158525
#> MYH7     0.35278624  0.76805953 -0.09860481  0.12158525  1.00000000
```

## Building the Adjacency Matrix

We threshold the correlation matrix into a binary adjacency matrix. Gene
pairs with absolute correlation at or above `cor_threshold` are
connected by an edge. Choosing an appropriate threshold is important —
too high produces a sparse disconnected network, too low produces a
dense uninformative one. A threshold of 0.7 is a common starting point
for bulk RNA-seq co-expression networks.

``` r

adj_mat <- build_adjacency_mat(cor_mat, cor_threshold = 0.7)

# Proportion of possible edges retained
sum(adj_mat) / (nrow(adj_mat) * (nrow(adj_mat) - 1))
#> [1] 0.01356313
```

## Constructing the Network

The adjacency matrix is passed to
[`build_network()`](../reference/build_network.md), which constructs an
igraph object and computes per-gene degree and betweenness centrality.

- **Degree** — the number of other genes a gene is correlated with above
  the threshold. High-degree genes are candidates for hub genes.
- **Betweenness centrality** — how often a gene lies on the shortest
  path between other gene pairs. High betweenness genes act as bridges
  between co-expression modules.

``` r

net <- build_network(adj_mat)

# Preview degree distribution
summary(net$degree)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.000   0.000   3.000   6.768  10.000  52.000
```

## Detecting Gene Modules

Gene modules are groups of highly interconnected genes that likely share
biological function. We use the Louvain algorithm for community
detection, which optimizes modularity — a measure of how well-separated
the modules are.

``` r

gene_mods <- detect_network_modules(net, community_method = "louvain")
#> Modularity: 0.617
#> Modules found: 161

# Modularity score — higher values indicate clearer module structure
gene_mods$modularity
#> [1] 0.6170515

# Number of modules detected
length(unique(gene_mods$modules))
#> [1] 161

# Preview module assignments
head(gene_mods$module_df)
#>            gene module degree
#> MT-CO1   MT-CO1     M1      3
#> ACTA1     ACTA1     M2     15
#> MT-RNR2 MT-RNR2     M3      1
#> TTN         TTN     M4      5
#> MYH7       MYH7     M2     14
#> MT-ND4   MT-ND4     M1     10
```

## Summarizing the Network

[`summarize_network()`](../reference/summarize_network.md) returns a
data frame of global network statistics, including total edge count,
mean degree, hub gene count, and modularity.

``` r

net_summary <- summarize_network(net, gene_mods$modules, gene_mods$modularity)
#>          metric     value
#> 1       n_genes  500.0000
#> 2       n_edges 1692.0000
#> 3   mean_degree    6.7700
#> 4 median_degree    3.0000
#> 5  n_hubs_deg20   45.0000
#> 6     n_modules  161.0000
#> 7    modularity    0.6171
```

## Identifying Hub Genes

Hub genes are the most highly connected genes in the network. They often
correspond to master regulators or genes central to the biological
processes captured by the co-expression structure.

``` r

hub_genes <- get_hub_genes(net, modules = gene_mods$modules, n = 20)
head(hub_genes)
#>              gene degree betweenness module
#> PHC2         PHC2     52    3164.134     M6
#> IL6R         IL6R     48    1873.185     M6
#> TNKS1BP1 TNKS1BP1     45    2967.964     M6
#> NAP1L1     NAP1L1     45    7714.047     M5
#> LDB3         LDB3     44    4912.569    M10
#> KIAA0368 KIAA0368     41    1830.652     M6
```

## Visualization

### Correlation Heatmap

The correlation heatmap shows pairwise gene-gene correlations with
hierarchical clustering, revealing block structure that corresponds to
co-expression modules.

``` r

plot_corr_map(cor_mat, cor_method = "pearson")
```

![](getting-started_files/figure-html/heatmap-1.png)

### Force-Directed Network Graph

The network plot shows the top hub genes as nodes sized by degree and
colored by module, with the top 20 hub genes labeled by gene symbol.

``` r

plot_network(net,
             modules        = gene_mods$modules,
             n_top          = 80,
             cor_threshold  = 0.7,
             community_method = "louvain")
```

![](getting-started_files/figure-html/network-plot-1.png)

Note that the number of labeled hub genes can be controlled with
`n_label` to reduce overlap in dense networks:

``` r

# Reduce label density for cleaner visualization
plot_network(net,
             modules          = gene_mods$modules,
             n_top            = 80,
             n_label          = 10,
             cor_threshold    = 0.7,
             community_method = "louvain")
```

## Exporting Results

All results can be exported to TSV files for downstream analysis or
reporting. Six files are written to a temporary directory by default.

``` r

exported_files <- net_export(
  cor_mat,
  net,
  net_summary,
  gene_mods$modules,
  gene_mods$module_df,
  cor_threshold = 0.7
)

exported_files
#> [1] "gene_correlations.tsv"  "hub_genes.tsv"          "module_assignments.tsv"
#> [4] "module_summary.tsv"     "network_summary.tsv"    "node_stats.tsv"
```

The exported files are:

- `gene_correlations.tsv` — all gene pairs above the correlation
  threshold
- `network_summary.tsv` — global network statistics
- `hub_genes.tsv` — top 20 hub genes ranked by degree
- `module_assignments.tsv` — module membership for all genes
- `node_stats.tsv` — degree, betweenness, and module for every gene
- `module_summary.tsv` — per-module size, mean degree, and top hub gene

## Command Line Interface

For non-interactive use, corinet provides a CLI via Rapp. After
installing the package, install the CLI with:

``` r

Rapp::install_pkg_cli_apps("corinet")
```

The CLI accepts a genes x samples counts matrix as TSV or CSV, with gene
identifiers as row names. A small example counts matrix is available at
`tests/testdata/counts.tsv` in the package repository for testing.

All subcommands accept the same core arguments:

``` bash
# Full network analysis pipeline — exports 6 TSV files
corinet network --counts counts.tsv --output results/ --n-top 500

# Gene correlation heatmap
corinet heatmap --counts counts.tsv --output results/ --n-top 500

# Force-directed network graph
corinet plot-network --counts counts.tsv --output results/ --n-top 500 --n-plot 80
```

See `corinet --help` and `corinet <subcommand> --help` for full
parameter documentation including correlation method, threshold, and
community detection options.

And that’s that!

## Session Info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] corinet_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.42.0 gtable_0.3.6               
#>  [3] circlize_0.4.18             shape_1.4.6.1              
#>  [5] rjson_0.2.23                xfun_0.57                  
#>  [7] bslib_0.10.0                ggplot2_4.0.3              
#>  [9] GlobalOptions_0.1.4         Biobase_2.72.0             
#> [11] lattice_0.22-9              vctrs_0.7.3                
#> [13] tools_4.6.0                 generics_0.1.4             
#> [15] stats4_4.6.0                parallel_4.6.0             
#> [17] tibble_3.3.1                cluster_2.1.8.2            
#> [19] pkgconfig_2.0.3             Matrix_1.7-5               
#> [21] RColorBrewer_1.1-3          S7_0.2.2                   
#> [23] desc_1.4.3                  S4Vectors_0.50.0           
#> [25] lifecycle_1.0.5             compiler_4.6.0             
#> [27] farver_2.1.2                textshaping_1.0.5          
#> [29] Seqinfo_1.2.0               codetools_0.2-20           
#> [31] ComplexHeatmap_2.28.0       clue_0.3-68                
#> [33] htmltools_0.5.9             sass_0.4.10                
#> [35] yaml_2.3.12                 pkgdown_2.2.0              
#> [37] pillar_1.11.1               crayon_1.5.3               
#> [39] jquerylib_0.1.4             cachem_1.1.0               
#> [41] DelayedArray_0.38.1         iterators_1.0.14           
#> [43] abind_1.4-8                 foreach_1.5.2              
#> [45] tidyselect_1.2.1            digest_0.6.39              
#> [47] dplyr_1.2.1                 labeling_0.4.3             
#> [49] fastmap_1.2.0               grid_4.6.0                 
#> [51] colorspace_2.1-2            cli_3.6.6                  
#> [53] SparseArray_1.12.2          magrittr_2.0.5             
#> [55] S4Arrays_1.12.0             withr_3.0.2                
#> [57] scales_1.4.0                rmarkdown_2.31             
#> [59] XVector_0.52.0              matrixStats_1.5.0          
#> [61] igraph_2.3.0                ragg_1.5.2                 
#> [63] png_0.1-9                   GetoptLong_1.1.1           
#> [65] evaluate_1.0.5              knitr_1.51                 
#> [67] GenomicRanges_1.64.0        IRanges_2.46.0             
#> [69] doParallel_1.0.17           rlang_1.2.0                
#> [71] glue_1.8.1                  BiocGenerics_0.58.0        
#> [73] jsonlite_2.0.0              R6_2.6.1                   
#> [75] MatrixGenerics_1.24.0       systemfonts_1.3.2          
#> [77] fs_2.1.0
```
