Comparing CellChat and CellPhoneDB results
================
2025-08-19

# Crosscheck function

To proceed with the crosscheck function, users should have ran both
`run_cellchat()` and `run_cellphonedb()` on their Seurat object. If
these steps have not been performed, refer to tutorials to [run
CellChat](https://github.com/nigelhojinker/CCCtools/blob/main/data-raw/Run_CellChat.md)
and [run
CellPhoneDB](https://github.com/nigelhojinker/CCCtools/blob/main/data-raw/Run_CellPhoneDB.md)
at the respective links.

We provide a function `crosscheck()` that takes as input both the
CellChat object from `run_cellchat()` and the list from
`run_cellphonedb()`, as well as an optional `threshold` argument should
users wish to adjust the p-value (default: 0.05). `crosscheck()` returns
a list of six elements: (1) CellChat result (2) CellPhoneDB result (3)
Full combined result including non-significant interactions (4) Filtered
combined result with only significant interactions in at least one tool
(5) 3x3 contingency table of interaction types in both CellChat and
CellPhoneDB (6) CellChat object filtered to contain only significant
interactions in both CellChat and CellPhoneDB. We categorise
interactions in the filtered result in terms of the agreement (and lack
thereof) between CellChat and CellPhoneDB as such:

- High: Interaction is significant in both CellChat and CellPhoneDB
- Mid CellChat: Interaction is found only in CellChat, and p-value \<
  threshold
- Mid CellPhoneDB: Interaction is found only in CellPhoneDB, and p-value
  \< threshold
- Low CellChat: Interaction is significant in CellChat but not in
  CellPhoneDB
- Low CellPhoneDB: Interaction is significant in CellPhoneDB but not in
  CellChat

``` r
data(seu.NL)

cellchat <- run_cellchat(seu.NL, group.by = "labels")

## After setting up python env
cpdb <- run_cellphonedb(seu.NL, labels = "labels")

combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)
# Joining with `by = join_by(interaction_name)`
# Joining with `by = join_by(id_cp_interaction, classification)`
# Removed interactions insignificant in both tools, and interactions only found in one tool and has a p-value >= 0.05 
#  CCCtools summary: 
#            CellPhoneDB
# CellChat      Sig Not Sig Not Found   Sum
#   Sig         221      22       470   713
#   Not Sig     710   31743     29763 62216
#   Not Found   586   21294         0 21880
#   Sum        1517   53059     30233 84809
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
# Created CellChat object with only High confidence interactions.

combine %>% names()
# [1] "cellchat"     "cellphonedb"  "combine.all"  "combine.sig"  "summary"      "cellchat.new"
```
