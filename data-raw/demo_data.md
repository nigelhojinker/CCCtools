Demo data for scpeakeR package
================
2025-01-05

# Download the raw data

We will use the same dataset used in the CellChat tutorial. The data is
described [here](https://doi.org/10.6084/m9.figshare.13520015.v1).

``` r
pacman::p_load(tidyverse, janitor, Seurat)

load( url("https://figshare.com/ndownloader/files/25950872"), verbose = TRUE )
# Loading objects:
#   data_humanSkin
```

``` r
names(data_humanSkin)
```

    ## [1] "data" "meta"

The metadata shows that there are 4 patients with non-lesional and 4
patients with lesional samples.

``` r
dim(data_humanSkin$meta)
```

    ## [1] 7563    3

``` r
head(data_humanSkin$meta)
```

    ##                     patient.id condition      labels
    ## S1_AACTCCCAGAGCTGCA   Patient1        LS Inflam. FIB
    ## S1_CAACCAATCCTCATTA   Patient1        LS   FBN1+ FIB
    ## S1_CGCTATCTCCTAGTGA   Patient1        LS Inflam. FIB
    ## S1_ATTTCTGCAGGACGTA   Patient1        LS Inflam. FIB
    ## S1_TGAGCCGAGCTGGAAC   Patient1        LS Inflam. FIB
    ## S1_CAGGTGCAGCCCAACC   Patient1        LS Inflam. FIB

``` r
data_humanSkin$meta %>% 
  count(patient.id, condition) %>% 
  arrange(condition, desc(n))
```

    ##   patient.id condition    n
    ## 1   Patient1        LS 2446
    ## 2   Patient5        LS 1102
    ## 3   Patient2        LS  928
    ## 4   Patient7        LS  535
    ## 5   Patient3        NL 1184
    ## 6  Patient11        NL  655
    ## 7  Patient14        NL  392
    ## 8  Patient15        NL  321

``` r
data_humanSkin$meta %>% 
  tabyl(labels, condition)
```

    ##        labels   LS   NL
    ##     APOE+ FIB 1228 1215
    ##     FBN1+ FIB  813  548
    ##  COL11A1+ FIB  181  196
    ##   Inflam. FIB  484   69
    ##          cDC1  121    7
    ##          cDC2  294   38
    ##            LC   67   20
    ##    Inflam. DC   81    1
    ##            TC  765  212
    ##    Inflam. TC  266   44
    ##    CD40LG+ TC  630  166
    ##           NKT   81   36

The expression data is normalized counts and consists of 17,328 genes
and 7,563 cells.

``` r
dim(data_humanSkin$data)
```

    ## [1] 17328  7563

``` r
data_humanSkin$data[1:5, 1:3]
```

    ## 5 x 3 sparse Matrix of class "dgCMatrix"
    ##          S1_AACTCCCAGAGCTGCA S1_CAACCAATCCTCATTA S1_CGCTATCTCCTAGTGA
    ## A1BG               0.5774867            1.121317                   .
    ## A1BG-AS1           0.5774867            .                          .
    ## A2M                .                    .                          .
    ## A2M-AS1            .                    .                          .
    ## A2ML1              .                    .                          .

# Convert to Seurat object

As most single-cell data is analyzed using Seurat, we built our
functions (i.e. run_cellchat, run_cellphonedb) to use Seurat data as
input. We use the `CreateSeuratObject()` function by the `Seurat`
package to convert our data to a Seurat object. In most cases, this step
is skipped for most users who already have a working Seurat object.

``` r
seu <- CreateSeuratObject(
  counts    = data_humanSkin$data, 
  meta.data = data_humanSkin$meta,
  min.cells    = 75,       # genes expressed in at least 1% of the samples
  min.features = 200
  )

# Set data matrix to data layer correctly
seu[["RNA"]]$data   <- seu[["RNA"]]$counts 
seu[["RNA"]]$counts <-  NULL
```

    ## Warning: Removing default layer, setting default to data

``` r
# Cell composition after QC
seu@meta.data %>% tabyl(labels, condition)
```

    ##        labels   LS   NL
    ##     APOE+ FIB 1103 1038
    ##     FBN1+ FIB  795  509
    ##  COL11A1+ FIB  167  166
    ##   Inflam. FIB  481   62
    ##          cDC1  120    7
    ##          cDC2  294   37
    ##            LC   67   20
    ##    Inflam. DC   80    1
    ##            TC  749  164
    ##    Inflam. TC  262   39
    ##    CD40LG+ TC  610  155
    ##           NKT   80   35

Let us subset the non-lesional and lesional samples separately:

``` r
seu.NL <- subset(seu, condition == "NL")
seu.LS <- subset(seu, condition == "LS")

seu.NL@meta.data <- seu.NL@meta.data %>% select(patient.id, labels)

seu.LS@meta.data <- seu.LS@meta.data %>% select(patient.id, labels)

seu.NL
```

    ## An object of class Seurat 
    ## 10353 features across 2233 samples within 1 assay 
    ## Active assay: RNA (10353 features, 0 variable features)
    ##  1 layer present: data

``` r
seu.LS
```

    ## An object of class Seurat 
    ## 10353 features across 4808 samples within 1 assay 
    ## Active assay: RNA (10353 features, 0 variable features)
    ##  1 layer present: data

Finally, we write this into the data directory:

``` r
usethis::use_data(seu.NL)
# ✔ Saving "seu.NL" to "data/seu.NL.rda".
# ☐ Document your data (see <https://r-pkgs.org/data.html>).
usethis::use_data(seu.LS)
# ✔ Saving "seu.LS" to "data/seu.LS.rda".
# ☐ Document your data (see <https://r-pkgs.org/data.html>
```
