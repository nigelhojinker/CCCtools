# Extract all relevant interactions from CellChat analysis

Adapted from CellChat's subsetCommunication() function, but retains all
interactions (significant or insignificant) less those that do not pass
min.cells criteria set by filterCommunication()

## Usage

``` r
pull_netslot(object = NULL, thresh = 0.05)
```

## Arguments

- object:

  CellChat object after run_cellchat()

- thresh:

  Threshold of p-value for defining interactions that are statistical
  significant

## Value

A dataframe consisting of inferred interactions for comparison with
CellPhoneDB

## Examples

``` r
if (FALSE) { # \dontrun{
cellchat.res <- pull_netslot(cellchat)
cellchat.res %>% view()
} # }
```
