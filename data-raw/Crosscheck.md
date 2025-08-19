Comparing CellChat and CellPhoneDB results
================
2025-08-19

# Crosscheck function

To proceed with the crosscheck function, users should have ran both
`run_cellchat()` and `run_cellphonedb()` on their Seurat object. If
these steps have not been performed, refer to tutorials to [run
CellChat](../Run_CellChat.md) and [run
CellPhoneDB](../Run_CellPhoneDB.md) at the respective links.

We provide a function `crosscheck()` that takes as input both CellChat
and CellPhoneDB results as dataframes and an optional `threshold`
argument should users wish to adjust the p-value (default: 0.05). The
CellChat result is obtained by running `pull_netslot(cellchat)` and its
CellPhoneDB counterpart is taken from the “pvalues” element of its
output list. `crosscheck()` returns a list of four elements: (1)
Combined result (2) CellChat result (3) CellPhoneDB result (4) Summary
table of containing number of interactions at each confidence level.

We classify the confidence level as such:

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
cellchat.res <- pull_netslot(cellchat)

cpdb.res <- cpdb$pvalues

combine <- crosscheck(cellchat_res = cellchat.res, cellphonedb_res = cpdb.res)
# Joining with `by = join_by(interaction_name)`
# Joining with `by = join_by(id_cp_interaction)`
# Removed interactions insignificant in both tools, and interactions only found in one tool and has a p-value >= 0.05 
#  CCCtools result: 
#       confidence   n percent
#             High 221  0.1105
#     Mid CellChat 470  0.2350
#  Mid CellPhoneDB 586  0.2930
#     Low CellChat  22  0.0110
#  Low CellPhoneDB 701  0.3505

combine %>% names()
# [1] "combine"     "cellchat"    "cellphonedb" "summary"
```
