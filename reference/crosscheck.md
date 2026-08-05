# Comparing CellChat and CellPhoneDB results

Comparing CellChat and CellPhoneDB results

## Usage

``` r
crosscheck(
  obj1 = NULL,
  obj2 = NULL,
  comp.type = c("CC_CP", "CP_CC", "CC_CC", "CP_CP"),
  p.cutoff = 0.05,
  return.all = FALSE,
  expand = FALSE
)
```

## Arguments

- obj1:

  CellChat object obtained from run_cellchat() on a Seurat object OR
  output list from run_cellphonedb() on a Seurat object.

- obj2:

  CellChat object obtained from run_cellchat() on a Seurat object OR
  output list from run_cellphonedb() on a Seurat object.

- comp.type:

  Character value stating the type of comparative analysis to run;
  assumes a CellChat-to-CellPhoneDB comparison by default ("CC_CP").

- p.cutoff:

  p-value (default: 0.05)

- return.all:

  Logical value to determine whether to return the full comparative
  analysis results (including non-significant interactions).

- expand:

  Logical value to highlight if interaction is found from cross-database
  analyses. This is used when comparing outputs using different
  databases (default = FALSE).

## Value

A list of five elements: (1) Comparative analysis result (2) 2x2
contingency table of interaction types in both CellChat and CellPhoneDB
(3) Number of significant interactions from both analyses (4) Jaccard
Index for significant interactions inferred (5) Matthew's Correlation
Coefficient to quantify agreement

## Examples

``` r
if (FALSE) { # \dontrun{
## After setting up conda env
cpdb <- run_cellphonedb(seurat.obj, ...)
cellchat <- run_cellchat(seurat.obj, ...)

combine <- crosscheck(obj1 = cellchat, obj2 = cpdb)
} # }
```
