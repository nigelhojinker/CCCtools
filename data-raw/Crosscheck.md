Comparing CellChat and CellPhoneDB results
================
2025-01-05

# Crosscheck function

## Comparative analysis across tools

To proceed with the crosscheck function, users should have ran both
`scpeaker(obj, labels, method = "cellchat")` and
`scpeaker(obj, labels, method = "cellphonedb")` on their Seurat object.
If these steps have not been performed, refer to tutorials to [run
CellChat](https://github.com/nigelhojinker/scpeakeR/blob/main/data-raw/Run_CellChat.md)
and [run
CellPhoneDB](https://github.com/nigelhojinker/scpeakeR/blob/main/data-raw/Run_CellPhoneDB.md)
at the respective links.

We provide a function `crosscheck()` that takes as input both the
CellChat object from `scpeaker(obj, labels, method = "cellchat")` and
the list from `scpeaker(obj, labels, method = "cellphonedb")`, as well
as an optional `threshold` argument should users wish to adjust the
p-value (default: 0.05). `crosscheck()` returns a list of two elements:
(1) Comparative analysis result (2) 3x3 contingency table of interaction
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
data(seu.NL)

cellchat <- scpeaker(seu.NL, labels = "labels", method = "cellchat")

# After setting up python env
cpdb <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")

combine <- crosscheck(obj1 = cellchat, obj2 = cpdb)
# Joining with `by = join_by(interaction_name)`
# Joining with `by = join_by(id_cp_interaction)`
# Comparative analysis done successfully for p-value 0.05 
#  scpeakeR summary: 
#            CellPhoneDB
# CellChat      Sig Not_Sig Not_Found   Sum
#   Sig         843       2         0   845
#   Not_Sig      39   83260         0 83299
#   Not_Found   513     735         0  1248
#   Sum        1395   83997         0 85392

combine %>% names()
# [1] "result"       "summary"
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
# Joining with `by = join_by(source, target, interaction_name, interaction_name_2, pathway_name, ligand, receptor, annotation,
# evidence, LR, unique_int)`
# Joining with `by = join_by(source, target, interaction_name, interaction_name_2, pathway_name, ligand, receptor, annotation,
# evidence, LR, unique_int)`
# Comparative analysis done successfully for p-value 0.05 
#  scpeakeR summary: 
#            cellchat.NL
# cellchat.LS   Sig Not_Sig Not_Found   Sum
#   Sig         623     409       232  1264
#   Not_Sig     222   82602      1016 83840
#   Not_Found     0     288         0   288
#   Sum         845   83299      1248 85392

combine %>% names()
# [1] "result"       "summary"
```

### Comparing two CellPhoneDB objects

``` r
cpdb.NL <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")

cpdb.LS <- scpeaker(seu.LS, labels = "labels", method = "cellphonedb")

combine <- crosscheck(obj1 = cpdb.LS, obj2 = cpdb.NL, comp.type = "CP_CP")
# Joining with `by = join_by(id_cp_interaction)`
# Joining with `by = join_by(id_cp_interaction)`
# Joining with `by = join_by(id_cp_interaction, interacting_pair, partner_a, partner_b, gene_a, gene_b, secreted, receptor_a,
# receptor_b, annotation_strategy, is_integrin, directionality, classification, source, target, LR, unique_int)`
# Joining with `by = join_by(id_cp_interaction, interacting_pair, partner_a, partner_b, gene_a, gene_b, secreted, receptor_a,
# receptor_b, annotation_strategy, is_integrin, directionality, classification, source, target, LR, unique_int)`
# Comparative analysis done successfully for p-value 0.05 
#  scpeakeR summary: 
#            cpdb.NL
# cpdb.LS       Sig Not_Sig Not_Found   Sum
#   Sig         893     400         0  1293
#   Not_Sig     502   83597         0 84099
#   Not_Found     0       0         0     0
#   Sum        1395   83997         0 85392

combine %>% names()
# [1] "result"       "summary"
```
