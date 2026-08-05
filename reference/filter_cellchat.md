# Subsetting CellChat object to contain interactions significant in both CCC tools

Subsetting CellChat object to contain interactions significant in both
CCC tools

## Usage

``` r
filter_cellchat(
  cellchat,
  crosscheck_res,
  category = c("Sig_Both", "Sig_CellChat_Not_Found_CellPhoneDB",
    "Sig_CellPhoneDB_Not_Found_CellChat", "Sig_CellChat_Not_Sig_CellPhoneDB",
    "Sig_CellPhoneDB_Not_Sig_CellChat")
)
```

## Arguments

- cellchat:

  CellChat object obtained from run_cellchat() on a Seurat object.

- crosscheck_res:

  Output data from crosscheck() function.

- category:

  Character value of the type of interaction to filter the cellchat
  object (default: significant in both CellChat and CellPhoneDB).

## Value

A new CellChat object that interactions filtered to your interest.

## Examples

``` r
if (FALSE) { # \dontrun{
## After setting up conda env
cpdb     <- scpeaker(seurat_obj, method = "cellphonedb", ...)
cellchat <- scpeaker(seurat_obj, method = "cellchat",    ...)

combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)

cellchat.new <- filter_cellchat(cellchat = cellchat, crosscheck_res = crosscheck)
} # }
```
