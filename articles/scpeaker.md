# Getting started with scpeaker

The goal of scpeaker is to provide a harmonised and reproducible
framework for cell-cell communication (CCC) inference involving
single-cell RNA sequencing (scRNA-seq) data. We employ the two most
established tools CellChat (R-based) and CellPhoneDB (Python-based) due
to their innate similarity in predicting CCC, which includes an in-built
ligand-receptor prior knowledge, the handling of heteromeric protein
complexes and statistical analysis via random shuffling of cell type
labels.

scpeaker takes as input a **processed Seurat object with cell type
annotations** for the respective analyses.

## Load package & example data

``` r

pacman::p_load(scpeaker)
```

We use the `seu.NL` example dataset as described in the homepage. The
full script on how the data was processed and generated may be viewed
[here](https://nigelhojinker.github.io/scpeaker/data-raw/demo_data.md).

``` r

data(seu.NL)

seu.NL@meta.data %>% 
  head() %>% 
  knitr::kable()
```

|                     | patient.id | labels      |
|:--------------------|:-----------|:------------|
| S3_ATGAGGGAGTCTTGCA | Patient3   | FBN1+ FIB   |
| S3_GTACGTACAAATTGCC | Patient3   | FBN1+ FIB   |
| S3_CTCGTCAGTGTTGAGG | Patient3   | Inflam. FIB |
| S3_ATTCTACGTAATCGTC | Patient3   | FBN1+ FIB   |
| S3_CTGCCTATCAATCACG | Patient3   | FBN1+ FIB   |
| S3_TAGTGGTAGGATGCGT | Patient3   | FBN1+ FIB   |

``` r

seu.NL@meta.data %>% 
  tabyl(labels) %>% 
  arrange(desc(n)) %>% 
  knitr::kable()
```

| labels       |    n |   percent |
|:-------------|-----:|----------:|
| APOE+ FIB    | 1038 | 0.4648455 |
| FBN1+ FIB    |  509 | 0.2279445 |
| COL11A1+ FIB |  166 | 0.0743395 |
| TC           |  164 | 0.0734438 |
| CD40LG+ TC   |  155 | 0.0694133 |
| Inflam. FIB  |   62 | 0.0277653 |
| Inflam. TC   |   39 | 0.0174653 |
| cDC2         |   37 | 0.0165696 |
| NKT          |   35 | 0.0156740 |
| LC           |   20 | 0.0089566 |
| cDC1         |    7 | 0.0031348 |
| Inflam. DC   |    1 | 0.0004478 |

## CCC analysis

Within the scpeaker package, we provide the general
[`scpeaker()`](https://nigelhojinker.github.io/scpeaker/reference/scpeaker.md)
function to perform **BOTH** CellChat and CellPhoneDB analyses on your
Seurat object. Importantly, scpeaker does not modify the innate
framework that these tools provide, and returns the original predictions
that they infer. By harmonising the ligand-receptor databases, gene
expression averaging strategy and output format that both tools
generate, scpeaker allows you to extract the interactions significantly
tested in CellChat and CellPhoneDB. Additionally, scpeaker accomodates
all arguments used in the statistical analysis of CellChat and
CellPhoneDB within the
[`scpeaker()`](https://nigelhojinker.github.io/scpeaker/reference/scpeaker.md)
function.

We harmonised both the human ligand-receptor databases of CellChat and
CellPhoneDB, providing scpeakerDB which contains 3,901 non-redundant
interactions that have been manually curated by the original authors.
Further, we provide the triMean gene expression averaging strategy for
CellPhoneDB to match that of CellChat, allowing both methods to use the
same robust mean method.

Finally, to optimise the computational efficiency of CellChat’s
analysis, we sped up the `computeCommunProb()` function used in
CellChat, by reimplementing the repeated calculating of gene expression
averages in C++. You may still opt to execute the original CellChat
implemetation by using `fast.mode = FALSE`.

### CellChat analysis

#### scpeaker-optimised

``` r

start_time <- Sys.time()
cellchat <- scpeaker(seu.NL,                  # Seurat object
                     labels   = "labels",     # Cell type label in object's meta data
                     method   = "cellchat",   # CCC method to use
                     type     = "triMean",    # Gene expression averaging strategy (default = triMean)
                     database = "scpeakerDB") # Ligand-receptor prior knowledge to use (default = scpeakerDB)
```

    ## Running cell-cell communication analysis with database: scpeakerDB

    ## [1] "Create a CellChat object from a Seurat object"
    ## The `meta.data` slot in the Seurat object is used as cell meta information

    ## Warning in createCellChat(obj, group.by = labels, assay = assay): The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`!

    ## Set cell identities for the new CellChat object 
    ## The cell groups used for CellChat analysis are  APOE+ FIB, FBN1+ FIB, COL11A1+ FIB, Inflam. FIB, cDC1, cDC2, LC, Inflam. DC, TC, Inflam. TC, CD40LG+ TC, NKT 
    ## triMean is used for calculating the average gene expression per cell group. 
    ## The number of highly variable ligand-receptor pairs used for signaling inference is 588 
    ## [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-08-05 10:03:18.581782]"
    ## [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-08-05 10:03:23.741633]"
    ## The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  cDC1, Inflam. DC !  38.4% interactions are removed!

    ## CellChat analysis completed succesfully.

``` r

end_time <- Sys.time()
print(end_time - start_time)
```

    ## Time difference of 7.049674 secs

#### Original CellChat implementation

``` r

start_time <- Sys.time()
cellchat_ori <- scpeaker(seu.NL,               # Seurat object
                     labels    = "labels",     # Cell type label in object's meta data
                     method    = "cellchat",   #   CCC method to use
                     type      = "triMean",    # Gene expression averaging strategy (default = triMean)
                     database  = "scpeakerDB", # Ligand-receptor prior knowledge to use (default = scpeakerDB)
                     fast.mode = FALSE)        # Original CellChat implementation
```

    ## Running cell-cell communication analysis with database: scpeakerDB

    ## [1] "Create a CellChat object from a Seurat object"
    ## The `meta.data` slot in the Seurat object is used as cell meta information

    ## Warning in createCellChat(obj, group.by = labels, assay = assay): The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`!

    ## Set cell identities for the new CellChat object 
    ## The cell groups used for CellChat analysis are  APOE+ FIB, FBN1+ FIB, COL11A1+ FIB, Inflam. FIB, cDC1, cDC2, LC, Inflam. DC, TC, Inflam. TC, CD40LG+ TC, NKT 
    ## triMean is used for calculating the average gene expression per cell group. 
    ## The number of highly variable ligand-receptor pairs used for signaling inference is 588 
    ## [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-08-05 10:03:24.028109]"
    ## [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-08-05 10:03:37.524316]"
    ## The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  cDC1, Inflam. DC !  38.4% interactions are removed!

    ## CellChat analysis completed succesfully.

``` r

end_time <- Sys.time()
print(end_time - start_time)
```

    ## Time difference of 13.75904 secs

Now that we observe that the runtime in the scpeaker-implemented
CellChat is shorter, how do their results compare? We can utilise the
`subsetCommunication()` function in CellChat to extract all interactions
predicted, alongside their communicaiton probability scores and
p-values.

``` r

subsetCommunication(cellchat) %>% 
  head() %>% 
  knitr::kable()
```

| source | target | ligand | receptor | prob | pval | interaction_name | interaction_name_2 | pathway_name | annotation | evidence |
|:---|:---|:---|:---|---:|---:|:---|:---|:---|:---|:---|
| APOE+ FIB | APOE+ FIB | FGF7 | FGFR1 | 0.0096029 | 0.00 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| Inflam. FIB | APOE+ FIB | FGF7 | FGFR1 | 0.0086075 | 0.01 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| APOE+ FIB | FBN1+ FIB | FGF7 | FGFR1 | 0.0106527 | 0.00 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| Inflam. FIB | FBN1+ FIB | FGF7 | FGFR1 | 0.0095496 | 0.01 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| APOE+ FIB | COL11A1+ FIB | FGF7 | FGFR1 | 0.0102132 | 0.00 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| Inflam. FIB | COL11A1+ FIB | FGF7 | FGFR1 | 0.0091552 | 0.01 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |

``` r

subsetCommunication(cellchat_ori) %>% 
  head() %>% 
  knitr::kable()
```

| source | target | ligand | receptor | prob | pval | interaction_name | interaction_name_2 | pathway_name | annotation | evidence |
|:---|:---|:---|:---|---:|---:|:---|:---|:---|:---|:---|
| APOE+ FIB | APOE+ FIB | FGF7 | FGFR1 | 0.0096029 | 0.00 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| Inflam. FIB | APOE+ FIB | FGF7 | FGFR1 | 0.0086075 | 0.01 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| APOE+ FIB | FBN1+ FIB | FGF7 | FGFR1 | 0.0106527 | 0.00 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| Inflam. FIB | FBN1+ FIB | FGF7 | FGFR1 | 0.0095496 | 0.01 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| APOE+ FIB | COL11A1+ FIB | FGF7 | FGFR1 | 0.0102132 | 0.00 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |
| Inflam. FIB | COL11A1+ FIB | FGF7 | FGFR1 | 0.0091552 | 0.01 | FGF7_FGFR1 | FGF7 - FGFR1 | FGF | Secreted Signaling | PMC: 4393358 |

``` r

identical(subsetCommunication(cellchat), subsetCommunication(cellchat_ori))
```

    ## [1] FALSE

We see that [`identical()`](https://rdrr.io/r/base/identical.html)
returned `FALSE` but the dataframes look almost identical. Rounding the
prob values off to the 15th decimal point, the results obtained from the
scpeaker-implemented CellChat and the original CellChat are the
identical using the same random seed.

``` r

res <- subsetCommunication(cellchat) %>% 
  mutate(prob = round(prob, 15),
         ID   = paste(source, target, ligand, receptor, prob, pval, sep = "|")) %>% 
  pull(ID)

ori <- subsetCommunication(cellchat_ori) %>% 
  mutate(prob = round(prob, 15),
         ID   = paste(source, target, ligand, receptor, prob, pval, sep = "|")) %>% 
  pull(ID)

res %>% head()
```

    ## [1] "APOE+ FIB|APOE+ FIB|FGF7|FGFR1|0.009602894517275|0"        
    ## [2] "Inflam. FIB|APOE+ FIB|FGF7|FGFR1|0.008607528426375|0.01"   
    ## [3] "APOE+ FIB|FBN1+ FIB|FGF7|FGFR1|0.010652686527943|0"        
    ## [4] "Inflam. FIB|FBN1+ FIB|FGF7|FGFR1|0.00954955584343|0.01"    
    ## [5] "APOE+ FIB|COL11A1+ FIB|FGF7|FGFR1|0.01021320237008|0"      
    ## [6] "Inflam. FIB|COL11A1+ FIB|FGF7|FGFR1|0.009155160984692|0.01"

``` r

ori %>% head()
```

    ## [1] "APOE+ FIB|APOE+ FIB|FGF7|FGFR1|0.009602894517275|0"        
    ## [2] "Inflam. FIB|APOE+ FIB|FGF7|FGFR1|0.008607528426375|0.01"   
    ## [3] "APOE+ FIB|FBN1+ FIB|FGF7|FGFR1|0.010652686527943|0"        
    ## [4] "Inflam. FIB|FBN1+ FIB|FGF7|FGFR1|0.00954955584343|0.01"    
    ## [5] "APOE+ FIB|COL11A1+ FIB|FGF7|FGFR1|0.01021320237008|0"      
    ## [6] "Inflam. FIB|COL11A1+ FIB|FGF7|FGFR1|0.009155160984692|0.01"

``` r

identical(res, ori)
```

    ## [1] FALSE

``` r

rm(res, cellchat_ori, ori, start_time, end_time)
```

### CellPhoneDB analysis

As CellPhoneDB is originally a Python package, we leverage the
`reticulate` package to interface with Python within R. We provide a
scpeaker.yaml file with contains the necessary dependencies to build the
required conda environment.

``` r

# Run this once for first conda environment set up
config <- file.path(pacman::p_path("scpeaker"), "data/scpeaker.yaml")
conda_create(envname = "scpeaker", environment = config)
```

``` r

# Run this each R session when conducting CellPhoneDB analysis 
use_condaenv("scpeaker")
```

We are now ready to perform CellPhoneDB statistical analysis:

``` r

cpdb <- scpeaker(seu.NL,                   # Seurat object
                 labels   = "labels",      # Cell type label in object's meta data
                 method   = "cellphonedb", # CCC method to use
                 type     = "triMean",     # Gene expression averaging strategy (default = triMean)
                 database = "scpeakerDB")  # Ligand-receptor prior knowledge to use (default = scpeakerDB)
```

    ## Running cell-cell communication analysis with database: scpeakerDB

    ## use_dir is NULL. Creating temp directory for file creation and storage.

    ## Warning in normalizePath(., winslash = "/"):
    ## path[1]="/var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj": No
    ## such file or directory

    ## Creating input and output files to directory: /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj 
    ##  Note: Make a copy of the temp directory for future use if necessary. 
    ## input_meta.tsv file created and saved to: /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj 
    ## input.h5ad file created and saved to: /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj

    ## Running CellPhoneDB Statistical Analysis with - triMean

    ## Reading user files...
    ## The following user files were loaded successfully:
    ## /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/input.h5ad
    ## /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/input_meta.tsv
    ## [ ][CORE][05/08/26-10:03:42][INFO] type = triMean: Setting threshold to 0.
    ## [ ][CORE][05/08/26-10:03:42][INFO] [Cluster Statistical Analysis] Threshold:0 Iterations:1000 Debug-seed:-1 Threads:4 Precision:3
    ## [ ][CORE][05/08/26-10:03:42][INFO] Running Real Analysis
    ## [ ][CORE][05/08/26-10:03:42][INFO] Running Statistical Analysis
    ## [ ][CORE][05/08/26-10:03:50][INFO] Building Pvalues result
    ## [ ][CORE][05/08/26-10:03:50][INFO] Building results
    ## Saved deconvoluted to /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/statistical_analysis_deconvoluted_08_05_2026_100350.txt
    ## Saved deconvoluted_percents to /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/statistical_analysis_deconvoluted_percents_08_05_2026_100350.txt
    ## Saved means to /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/statistical_analysis_means_08_05_2026_100350.txt
    ## Saved pvalues to /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/statistical_analysis_pvalues_08_05_2026_100350.txt
    ## Saved significant_means to /var/folders/67/jqzxz5vn6c15btlm5m1v3rxm0000gn/T//Rtmph6SRl4/obj/statistical_analysis_significant_means_08_05_2026_100350.txt

    ## Warning in py_to_r.pandas.core.frame.DataFrame(<environment>): index contains
    ## duplicated values: row names not set

    ## Warning in py_to_r.pandas.core.frame.DataFrame(<environment>): index contains
    ## duplicated values: row names not set

    ## CellPhoneDB analysis completed successfully.

Note: The Python numpy random.seed implementation is NOT thread-safe. To
set a random seed for reproducibility purposes, run:

``` r

cpdb <- scpeaker(seu.NL,                     # Seurat object
                 labels     = "labels",      # Cell type label in object's meta data
                 method     = "cellphonedb", # CCC method to use
                 type       = "triMean",     # Gene expression averaging strategy (default = triMean)
                 database   = "scpeakerDB",  # Ligand-receptor prior knowledge to use (default = scpeakerDB)
                 threads    = 1,             # Set number of threads to 1
                 debug_seed = 42)            # Your random seed
```

## Crosscheck analysis

After conducting both CellChat and CellPhoneDB analyses, scpeaker can
determine which interactions are significantly inferred by both tools.
Additionally, we provide functions to filter only these concordant
interactions to a new CellChat object to leverage the extensive
downstream visualisation and network analysis functions already
available in CellChat.

**NOTE:** The filter CellChat object after
[`filter_cellchat()`](https://nigelhojinker.github.io/scpeaker/reference/filter_cellchat.md)
is a processed CellChat object with only interactions significant in
both CellChat and CellPhoneDB. We then ran `computeCommunProbPathway()`,
`aggregateNet()` and `netAnalysis_computeCentrality(slot.name = "netP")`
on the filter object.

``` r

# Compare CellChat and CellPhoneDB results
crosscheck <- crosscheck(cellchat, cpdb)
```

    ## Joining with `by = join_by(interaction_name)`
    ## Joining with `by = join_by(id_cp_interaction)`

    ## Comparative analysis done successfully for p-value 0.05 
    ##  scpeaker summary: 
    ##          CellPhoneDB
    ## CellChat    Sig Not_Sig   Sum
    ##   Sig       843       2   845
    ##   Not_Sig   552   83995 84547
    ##   Sum      1395   83997 85392
    ## Number of significant interactions in both CellChat and CellPhoneDB: 843 
    ## Jaccard Index: 0.6034359 
    ## Matthew's Correlation Coefficient: 0.7738874

``` r

# Filter interactions significant in both analyses and return new CellChat object
cellchat_new <- filter_cellchat(cellchat, crosscheck)
```

    ## Filtering interactions in category: Sig_Both
