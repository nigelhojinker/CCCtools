# Demo dataset prepared from CellChat tutorial data

The dataset used in the main CellChat tutorial after quality control
filtering (min.cells = 75, min.features = 200)

## Usage

``` r
seu.NL
```

## Format

A Seurat object with 10,353 genes and 2,233 cells. The meta data slot
includes:

- patient.id:

  Patient identifier

- labels:

  Cell type annotation provided by the CellChat authors

## Source

Data was downloaded from https://doi.org/10.6084/m9.figshare.13520015.v1
and processed further (see
https://github.com/Nigel-Ho23/CCCtools/tree/main/data-raw/demo_data.md
for details).

## Examples

``` r

data(seu.NL)

seu.NL@meta.data %>%
  tabyl(labels)
#>        labels    n     percent
#>     APOE+ FIB 1038 0.464845499
#>     FBN1+ FIB  509 0.227944469
#>  COL11A1+ FIB  166 0.074339454
#>   Inflam. FIB   62 0.027765338
#>          cDC1    7 0.003134796
#>          cDC2   37 0.016569637
#>            LC   20 0.008956561
#>    Inflam. DC    1 0.000447828
#>            TC  164 0.073443798
#>    Inflam. TC   39 0.017465293
#>    CD40LG+ TC  155 0.069413345
#>           NKT   35 0.015673981
```
