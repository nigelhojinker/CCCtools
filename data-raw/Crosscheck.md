Comparing CellChat and CellPhoneDB results
================
2026-07-16

# Crosscheck function

## Comparative analysis across tools

To proceed with the crosscheck function, users should have ran both
`scpeaker(obj, labels, method = "cellchat")` and
`scpeaker(obj, labels, method = "cellphonedb")` on their Seurat object.
If these steps have not been performed, refer to tutorials to [run
CellChat](https://github.com/nigelhojinker/scpeaker/blob/main/data-raw/Run_CellChat.md)
and [run
CellPhoneDB](https://github.com/nigelhojinker/scpeaker/blob/main/data-raw/Run_CellPhoneDB.md)
at the respective links.

We provide a function `crosscheck()` that takes as input both the
CellChat object from `scpeaker(obj, labels, method = "cellchat")` and
the list from `scpeaker(obj, labels, method = "cellphonedb")`, as well
as an optional `threshold` argument should users wish to adjust the
p-value (default: 0.05). `crosscheck()` returns a list of two elements:
(1) Comparative analysis result (2) 2x2 contingency table of interaction
types in both CellChat and CellPhoneDB. We categorise interactions in
the filtered result in terms of the agreement (and lack thereof) between
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
library(scpeaker)
use_condaenv("scpeaker")

data(seu.NL)

cellchat <- scpeaker(seu.NL, labels = "labels", method = "cellchat")
cpdb     <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")

combine <- crosscheck(obj1 = cellchat, obj2 = cpdb)
# Joining with `by = join_by(interaction_name)`
# Joining with `by = join_by(id_cp_interaction)`
# Comparative analysis done successfully for p-value 0.05 
#  scpeaker summary: 
#           CellPhoneDB
# CellChat     Sig   Not_Sig     Sum
#   Sig        841         4     845
#   Not_Sig    554     83993   84547
#   Sum       1395     83997   85392
# Number of significant interactions in both CellChat and CellPhoneDB: 841
# Jaccard Index: 0.6011437
# Mathew’s Correlation coefficient: 0.7720208


combine %>% names()
# [1] "result"    "summary"   "nsig_both" "jaccard"   "mcc"
```

## Comparative analysis within tools

This section is for users who are only interested in investigating the
difference between two Seurat objects (e.g. different conditions,
timepoints etc) that have been processed using the same cell-cell
communication tools. In `crosscheck()`, we provide an argument
`comp.type` where users can specify the tools used for comparative
analysis:

- “CC_CP”: obj1 is a CellChat object and obj2 is a CellPhoneDB object
  (default)
- “CC_CC”: both obj1 and obj2 are CellChat objects
- “CP_CP”: both obj1 and obj2 are CellPhoneDB objects

We provide users with the lesional skin dataset from the CellChat
tutorial, `seu.LS`, as an example here in this tutorial:

``` r
data(seu.NL, seu.LS)
```

### Comparing two CellChat objects

``` r
cellchat.NL <- scpeaker(seu.NL, labels = "labels", method = "cellchat")

cellchat.LS <- scpeaker(seu.LS, labels = "labels", method = "cellchat")

combine <- crosscheck(obj1 = cellchat.LS, obj2 = cellchat.NL, comp.type = "CC_CC")
# Joining with `by = join_by(interaction_name)`
# Joining with `by = join_by(interaction_name)`
# Comparative analysis done successfully for p-value 0.05 
#  scpeaker summary: 
#            cellchat.NL
# cellchat.LS   Sig Not_Sig   Sum
#     Sig       623     641  1264
#     Not_Sig   222   83906 84128
#     Sum       845   84547 85392
# Number of significant interactions in both CellChat analyses: 623 
# Jaccard Index: 0.4192463 
# Matthew's Correlation Coefficient: 0.5981029 

combine %>% names()
# [1] "result"    "summary"   "nsig_both" "jaccard"   "mcc"    
```

### Comparing two CellPhoneDB objects

``` r
cpdb.NL <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")

cpdb.LS <- scpeaker(seu.LS, labels = "labels", method = "cellphonedb")

combine <- crosscheck(obj1 = cpdb.LS, obj2 = cpdb.NL, comp.type = "CP_CP")
# Joining with `by = join_by(id_cp_interaction)`
# Joining with `by = join_by(id_cp_interaction)`
# Comparative analysis done successfully for p-value 0.05 
#  scpeaker summary: 
#          cpdb.NL
# cpdb.LS     Sig Not_Sig   Sum
#   Sig       898     393  1291
#   Not_Sig   503   83598 84101
#   Sum      1401   83991 85392
# Number of significant interactions in both CellPhoneDB analyses: 898 
# Jaccard Index: 0.5005574 
# Matthew's Correlation Coefficient: 0.6624117

combine %>% names()
# [1] "result"    "summary"   "nsig_both" "jaccard"   "mcc"    
```
