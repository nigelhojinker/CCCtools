
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scpeakeR

<!-- badges: start -->

<!-- badges: end -->

The goal of scpeakeR is to provide functions to run [CellChat
(v2)](https://github.com/jinworks/CellChat) and [CellPhoneDB
(v5)](https://github.com/ventolab/CellphoneDB/tree/master) as well as
compare their outputs.

## Installation

### Prerequisites

Before installing **scpeakeR**, be sure to install the below packages
via Bioconductor, as they are dependencies of CellChat:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "BiocNeighbors", "Biobase", "BiocGenerics", "SingleCellExperiment", 
  "S4Vectors", "IRanges", "MatrixGenerics", "ComplexHeatmap"
)

missing_pkgs <- setdiff(bioc_pkgs, rownames(installed.packages()))

if (length(missing_pkgs) > 0) {
  BiocManager::install(missing_pkgs, ask = FALSE, update = TRUE)
}
```

### Download package

Now, you are ready to install the development version of scpeakeR:

``` r
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load_gh("nigelhojinker/scpeakeR")
```

## Example data

`seu.NL` is the Seurat object containing 2,233 cells from non-lesional
human skin from four samples. This dataset was used originally in the
tutorials for CellChat. See this [script for how the data was
generated](data-raw/demo_data.md).

``` r
pacman::p_load_gh("nigelhojinker/scpeakeR")

data(seu.NL)

seu.NL
#> An object of class Seurat 
#> 10353 features across 2233 samples within 1 assay 
#> Active assay: RNA (10353 features, 0 variable features)
#>  1 layer present: data
```

## Functions

In **scpeakeR**, we provide the `scpeaker()` function to perform
CellChat and CellPhoneDB (method 2) analysis directly on their processed
Seurat object.

Please refer to the following links to the full documentations on
running these tools on your dataset:

### Perform CellPhoneDB on Seurat object

- [Running CellChat on Seurat
  object](https://github.com/nigelhojinker/scpeakeR/blob/main/data-raw/Run_CellChat.md)

``` r
# Run this ONCE to create conda environment
config <- file.path(pacman::p_path("scpeakeR"), "data/scpeaker.yaml")
conda_create(envname = "scpeakeR", environment = config)
```

``` r
# Run this each R session when performing CellPhoneDB analysis
use_condaenv("scpeakeR")

# CellPhoneDB analysis -- minimal example
cellphonedb <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")
```

### Perform CellChat on Seurat object

- [Running CellPhoneDB on Seurat
  object](https://github.com/nigelhojinker/scpeakeR/blob/main/data-raw/Run_CellPhoneDB.md)

``` r
# CellChat analysis -- minimal example
cellchat <- scpeaker(seu.NL, labels = "labels", method = "cellchat")
```

### Cross-method analysis

After running BOTH CellChat and CellPhoneDB analysis, we provide the
`crosscheck()` function to compare the results of both CCC tools and
classify the interactions based on the significance testing conducted.
Please refer to the link below on how to compare CellChat and
CellPhoneDB results on your dataset:

- [Comparing CellChat and CellPhoneDB
  results](https://github.com/nigelhojinker/scpeakeR/blob/main/data-raw/Crosscheck.md)

``` r
combine <- crosscheck(cellchat, cellphonedb)
```

### Downstream visualisation and network analysis

The `crosscheck()` function allows users to extract interactions
significant in both CellChat and CellPhoneDB analyses. We provide the
`filter_cellchat()` function that takes in the original CellChat object
and the combined results (`combine`), and returns a new CellChat object
containing only the filtered interactions from cross-method analysis. We
leverage the abundant resource for visualising and conducting network
analysis using CellChat.

Users may refer to the full CellChat vignette
[here](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html).

``` r
cellchat_new <- filter_cellchat(cellchat, combine)
```

## Database harmonisation

In order to identify interactions that map to CellChat and/or
CellPhoneDB, the ligand-receptor pair involved have to share a common
notation. We convert all gene and complex names to uniprot IDs for
unbiased and clear identification of the interacting partners, and used
that uniprot-based identifier to compare interactions from both CCC
tools.

After we harmonized both databases, we compiled all unique interactions
from both CellChat and CellPhoneDB which we coined **scpeakeRDB**, which
contains 3,901 interactions. By default, scpeakeRDB is the database used
for all cell-cell communication inference in this package. However, we
provide users with the flexible to run their analysis with
CellChat/CellPhoneDB-only interactions.

For the full code implementation of from the mapping of interactions in
CellChat and CellPhoneDB, to the final creation of scpeakeRDB, you may
visit the links below:

- [Aligning CellChat and CellPhoneDB interaction databases - Part
  1](https://nigelhojinker.github.io/scpeakeR/data-raw/Database_harmonisation_part1.html)
- [Implementing changes reconstructing CellChatDB.human and CellPhoneDB
  database - Part
  2](https://nigelhojinker.github.io/scpeakeR/data-raw/Database_harmonisation_part2.html)
- [Creating scpeakeRDB after database harmonisation - Part
  3](https://nigelhojinker.github.io/scpeakeR/data-raw/Database_harmonisation_part3.html)

# Final notes

If you have used scpeakeR and any of its related functionalities, please
consider citing our paper:

- Link to upcoming paper

Additionally, please remember to cite the original authors of
[CellChat](https://github.com/jinworks/CellChat) and
[CellPhoneDB](https://github.com/ventolab/CellphoneDB/tree/master) for
their work in building these CCC tools for the community!
