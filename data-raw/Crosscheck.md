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
a list of three elements: (1) Comparative analysis result (2) 3x3
contingency table of interaction types in both CellChat and CellPhoneDB
(3) CellChat object filtered to contain only significant interactions in
both CellChat and CellPhoneDB. We categorise interactions in the
filtered result in terms of the agreement (and lack thereof) between
CellChat and CellPhoneDB as such:

- Sig_Both: Interaction is significant in both CellChat and CellPhoneDB
- Sig_CellChat_Not_Found_CellPhoneDB: Interaction is found only in
  CellChat, and p-value \< threshold
- Sig_CellPhoneDB_Not_Found_CellChat: Interaction is found only in
  CellPhoneDB, and p-value \< threshold
- Sig_CellChat_Not_Sig_CellPhoneDB: Interaction is significant in
  CellChat but not significant in CellPhoneDB
- Sig_CellPhoneDB_Not_Sig_CellChat: Interaction is significant in
  CellPhoneDB but not significant in CellChat

``` r
data(seu.NL)

cellchat <- run_cellchat(seu.NL, group.by = "labels")

## After setting up python env
cpdb <- run_cellphonedb(seu.NL, labels = "labels")

combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)
# Joining with `by = join_by(interaction_name)`
# Joining with `by = join_by(id_cp_interaction, classification)`
# Comparative analysis done successfully for p-value 0.05 
#  CCCtools summary: 
#            CellPhoneDB
# CellChat      Sig Not_Sig Not_Found   Sum
#   Sig         220      23       470   713
#   Not_Sig     705   31748     29763 62216
#   Not_Found   585   21295         0 21880
#   Sum        1510   53066     30233 84809
# Filtering cellchat object to only interactions significant in both CellChat and CellPhoneDB.
# Filtered object can be obtained from the 'cellchat.new' element of the output.
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
# Created CellChat object with only interactions significant in both CellChat and CellPhoneDB.

combine %>% names()
# [1] "result"       "summary"      "cellchat.new"
```
