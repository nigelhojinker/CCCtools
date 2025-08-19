
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CCCtools

<!-- badges: start -->

<!-- badges: end -->

The goal of CCCtools is to provide functions to run [CellChat
(v2)](https://github.com/jinworks/CellChat) and [CellPhoneDB
(v5)](https://github.com/ventolab/CellphoneDB/tree/master) as well as
compare their outputs. Gokce was here to check username chnage!

## Installation

### Prerequisites

Before install **CCCtools**, be sure to install the below packages via
Bioconductor, as they are dependencies of CellChat:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "BiocNeighbors", "Biobase", "BiocGenerics", "SingleCellExperiment", 
  "S4Vectors", "IRanges", "MatrixGenerics", "ComplexHeatmap"
)

missing_pkgs <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  BiocManager::install(missing_pkgs, ask = FALSE, update = TRUE)
}
```

### Download package

Now, you are ready to install the development version of CCCtools:

``` r
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load_gh("nigelhojinker/CCCtools")
```

## Example data

`seu.NL` is the Seurat object containing 2,233 cells from non-lesional
human skin from four different patients. This dataset was created from
the data used originally in CellChat tutorial. See this [script for how
the data was generated](data-raw/demo_data.md).

``` r
pacman::p_load_gh("nigelhojinker/CCCtools")

data(seu.NL)

seu.NL
#> An object of class Seurat 
#> 10353 features across 2233 samples within 1 assay 
#> Active assay: RNA (10353 features, 0 variable features)
#>  1 layer present: data
```

## Functions

In **CCCtools**, we provide the functions `run_cellchat()` and
`run_cellphonedb()` for users to perform CellChat and CellPhoneDB
(method 2) analysis directly on their processed Seurat object.

Please refer to the following links to the documentations on running
these tools on your dataset:

- [Running CellChat on Seurat
  object](https://github.com/nigelhojinker/CCCtools/blob/main/data-raw/Run_CellChat.md)
- [Running CellPhoneDB on Seurat
  object](https://github.com/nigelhojinker/CCCtools/blob/main/data-raw/Run_CellPhoneDB.md)

After running BOTH CellChat and CellPhoneDB analysis, we provide the
`crosscheck()` function to compare the results of both CCC tools and
rank the interactions based on the significance testing conducted.
Please refer to the link below on how to compare CellChat and
CellPhoneDB results on your dataset:

- [Comparing CellChat and CellPhoneDB
  results](https://github.com/nigelhojinker/CCCtools/blob/main/data-raw/Crosscheck.md)

## CellPhoneDB-to-CellChatDB Mapping

In order to identify interactions that map to CellChat and/or
CellPhoneDB, the ligand-receptor pair involved have to share a common
notation. We convert all gene and complex names to uniprot IDs for
unbiased and clear identification of the interacting partners, and used
that uniprot-based identifier to compare interactions from both CCC
tools.

Users may view this
[webpage](https://cellphonedb-cellchatdb-mapping.vercel.app/) for full
details on the mapping process used to map interactions across both
databases.

The mapped databases of CellPhoneDB (v5) and CellChat (v2) and are
available by running `data(CPDB)` and `data(CCDB)` respectively. These
databases will be used in the `crosscheck()` function after users have
ran both CellPhoneDB and CellChat on their dataset, and would like to
find out which interactions have been identified by both tools.
