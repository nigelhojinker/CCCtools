Wrapper for running CellChat
================
2025-01-05

## Introduction

CellChat is a popular R package for inferring cell-cell communication
from single-cell and/or spatial transcriptomics data. It offers
functions to visualize the analysis, and the github repository can be
found [**here**](https://github.com/jinworks/CellChat). CellChat has
interaction databases for three different species, of which we will use
HUMAN to find consensus with CellPhoneDB. Like CellPhoneDB, its database
contains only manually curated ligand-receptor interactions and the
subunit architecture of heteromeric complexes are considered. In
additon, CellChat considers co-factors that are involved in the
activation or inhibition of ligand-receptor interactions.

## Running CellChat directly from Seurat object

`scpeaker()` takes as input a Seurat object with the normalized counts
in the data layer, the annotated cell types stored in its meta.data and
the cell-cell communication method to be used (cellchat in this case).
For more details on the cellchat-run in `scpeaker()`, users may refer
[**here**](../R/run_cellchat.R). CellChat offers different methods of
calculating average gene expression per cell type (default: triMean). As
described by the developers of CellChat, triMean is a robust mean method
given by $[(Q1 + 2M + Q3) / 4]$\* that returns fewer but stronger
interactions. If users would like to change the method for calculating
average expression, they may use the `type` parameter to manually
overwrite this.

`scpeaker()` takes in the same arguments as CellChat’s
`createCellChat()`, `subsetDB()`, `computeCommunProb()` and
`filterCommunication()`. Details on the respective functions can be
found by running help(<package_name>) in R.

Lastly, `scpeaker()` runs CellChat analysis on a single dataset up to
and including `filterCommunication()`. For further downstream analyses,
please see the [CellChat vignette for a single
dataset](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html)
for the full tutorial. Note that CellChat is already installed when
users install scpeakeR, therefore, users may immediately run CellChat
functions used in the tutorial(s).

\*: **Q1** First quartile; **M** Median; **Q3** Third quartile

## Tutorial

R set-up and data input ([example data](../demo_data.md))

``` r
pacman::p_load_gh("nigelhojinker/scpeakeR")
setwd(this.path::here())
rm(list = ls())

data("seu.NL")

seu.NL@meta.data %>% head()
```

    ##                     patient.id      labels
    ## S3_ATGAGGGAGTCTTGCA   Patient3   FBN1+ FIB
    ## S3_GTACGTACAAATTGCC   Patient3   FBN1+ FIB
    ## S3_CTCGTCAGTGTTGAGG   Patient3 Inflam. FIB
    ## S3_ATTCTACGTAATCGTC   Patient3   FBN1+ FIB
    ## S3_CTGCCTATCAATCACG   Patient3   FBN1+ FIB
    ## S3_TAGTGGTAGGATGCGT   Patient3   FBN1+ FIB

In `seu.NL`, the metadata column for the annotated cell type labels is
“labels”. As such, we take “labels” as the input for the `labels`
argument in `scpeaker()`.

### Running CellChat with all default parameters

This step typically takes a few minutes to run:

``` r
cellchat <- scpeaker(seu.NL, labels = "labels", method = "cellchat")

# [1] "Create a CellChat object from a Seurat object"
# The `meta.data` slot in the Seurat object is used as cell meta information 
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  APOE+ FIB, FBN1+ FIB, COL11A1+ FIB, Inflam. FIB, cDC1, cDC2, LC, Inflam. DC, TC, Inflam. TC, CD40LG+ TC, NKT 
# triMean is used for calculating the average gene expression per cell group. 
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
# The number of highly variable ligand-receptor pairs used for signaling inference is 588 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-01-05 15:08:19.014103]"
#   |====================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-01-05 15:10:58.284517]"
# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  cDC1, Inflam. DC !   38.4% interactions are removed!
# CellChat analysis completed succesfully.
# Warning message:
# In createCellChat(obj, group.by = labels, assay = assay) :
#   The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`! 
```

### Running CellChat with filtered database

For users interested in only using a subset of CellChat’s database, they
may filter the database in accordance with the `subsetDB()` function by
CellChat. We provide an argument `subsetDB` that is set to FALSE by
default. To filter, simply set `subsetDB = TRUE` in `scpeaker()` and add
in the corresponding arguments as per CellChat’s `subsetDB()` function.

``` r
cellchat <- scpeaker(seu.NL, 
                     labels = "labels",
                     method = "cellchat",
                     subsetDB = TRUE, 
                     search = c("Cell-Cell Contact"), 
                     key = c("annotation"))
```

For more than one level of filtering, `key` takes in a character vector
of the columns in CellChatDB and `search` has to take in a list
corresponding to the item(s) to filter in the keys.

``` r
cellchat <- scpeaker(seu.NL, 
                     labels = "labels", 
                     method = "cellchat",
                     subsetDB = TRUE, 
                     search = list(c("Cell-Cell Contact", "Secreted Signaling"), c("CellChatDB v1")),
                     key = c("annotation", "version"))
```

To remove only non-protein signaling interactions, since a majority of
non-protein signaling interactions consists of metabolic and synaptic
signaling:

``` r
cellchat <- scpeaker(seu.NL, 
                     labels = "labels", 
                     method = "cellchat",
                     subsetDB = TRUE)
```

### Extracting all relevant interactions from CellChat analysis

We define relevant interactions as all inferred interactions from
CellChat - both significant and insignificant - except those that do not
pass the criteria set by min.cells in CellChat’s `filterCommunication()`
function.

We provide a function `pull_netslot()` to execute this extraction and
return a dataframe of all relevant interactions from the analysis. This
function is an edited function based on CellChat’s
`subsetCommunication()`. The rationale for `pull_netslot()` is to obtain
full output comparison from both CellChat and CellPhoneDB with the
reasoning that some interactions in CellChat may be insignificant in
CellChat, but significant in CellPhoneDB, vice versa. This will be
relevant to the `crosscheck()` function after users have ran both
CellChat and CellPhoneDB.

Refer [**here**](../R/pull_netslot.R) for more details.

``` r
cellchat.res <- pull_netslot(cellchat)

cellchat.res %>% glimpse()
# Rows: 84,144
# Columns: 11
# $ source             <fct> APOE+ FIB, FBN1+ FIB, COL11A1+ FIB, Inflam. FIB, cDC1, cDC2, LC, Inflam. DC, TC, Inflam. TC, CD40…
# $ target             <fct> APOE+ FIB, APOE+ FIB, APOE+ FIB, APOE+ FIB, APOE+ FIB, APOE+ FIB, APOE+ FIB, APOE+ FIB, APOE+ FIB…
# $ interaction_name   <fct> TGFB1_TGFBR1_TGFBR2, TGFB1_TGFBR1_TGFBR2, TGFB1_TGFBR1_TGFBR2, TGFB1_TGFBR1_TGFBR2, TGFB1_TGFBR1_…
# $ prob               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
# $ pval               <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
# $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB1 - (TGFBR1+TGFBR2)", "TGFB1 - (TGFBR1+TGFBR2)", "TGFB1 - (TGFBR1…
# $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "…
# $ ligand             <chr> "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1", "TGFB1"…
# $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "TGFbR…
# $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secreted Signaling", "Secreted Signaling", "Secreted…
# $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: …
```
