Running CellPhoneDB in R
================
2025-01-05

# Introduction

CellPhoneDB is a popular **Python** package developed by [Teichmann
Lab](https://www.teichlab.org/), and currently managed and further
developed by the [Vento-Tormo Lab](https://ventolab.org/) from
CellPhoneDB version 3. It is a useful cell-cell communication tool to
directly look up interactions of interest from single-cell multiomics
data, with their database of only manually curated and reviewed HUMAN
interactions.

A key feature of CellPhoneDB is its inherent function of accounting for
heteromeric complexes by decomplexifying to their subunits and ensuring
that expression of all subunits have to be non-zero and above a
user-defined threshold to be considered present in the interaction. This
accurately recapitulates how interactions can occur in vivo. As such, we
thought to incorporate CellPhoneDB with another popular cell-cell
communication tool, CellChat, that shares the same feature.

We forked the original CellPhoneDB repository and implement an
additional change to allow triMean as a method for computing the average
gene expression per cell type cluster like CellChat. We direct the
installation of the [fork
repository](https://github.com/nigelhojinker/scpeaker-CellphoneDB) with
the changes in the .yaml file mentioned in the section below. triMean is
the default method used in scpeaker for calculating the average gene
expression per cell cluster. Users may switch to the CellPhoneDB-default
method by setting `type = "thresholdedMean`.

# Setting up required environment

As mentioned in the **Introduction**, CellPhoneDB is a Python package
and currently, does not have an R equivalent. Here, we provide a wrapper
function `scpeaker()` for running CellPhoneDB completely within R. To do
this, we require the use of the `reticulate` package (a dependency in
scpeaker) and provide users with the `scpeaker.yaml` file for creating
the CellPhoneDB python environment in R.

## CellPhoneDB conda environment in R

### R set-up

``` r
pacman::p_load_gh("nigelhojinker/scpeaker")
setwd(this.path::here())
rm(list = ls())
```

### Set up conda & python via reticulate

For users **WITHOUT** a prior conda installation, you may run the line
below for installation of miniconda:

``` r
if(!file.exists(conda_binary())) install_miniconda()
```

For users **WITH** conda installed previously, we can check that
`reticulate` automatically detects the installed conda by running:

``` r
conda_binary()
```

#### What if my conda installation is not automatically detected?

If the Conda executable is not automatically detected, located the path
to conda.exe (e.g. using Windows File Explorer) and set it manually in
the .Rprofile.

Open the user profile:

``` r
usethis::edit_r_profile(scope = "user")
```

Add in the following lines (with correct path) and restart Rstudio.

``` r
# In this example, we show the case when conda is installed via miniforge and how users may manually set conda path
if (.Platform$OS.type == "windows") {
  Sys.setenv(
    RETICULATE_CONDA = "C:/users/user/miniforge3/condabin/conda.exe"
  )
}
```

Once installation and conda selection is complete, we can confirm this
by running:

``` r
conda_binary()

## [1] "C:/path/to/miniforge3/condabin/conda.bat"
```

Next, we have to create the CellPhoneDB conda environment with the
necessary Python modules for running the analysis. This is done **only
once**, and typically takes a few minutes to complete:

``` r
config <- file.path( pacman::p_path("scpeaker"), "data/scpeaker.yaml" )

conda_create(envname = "scpeaker", environment = config)
```

Once the cpdb environment has been created, we have to select its python
interpreter:

``` r
use_condaenv("scpeaker")
```

# Running CellPhoneDB

Now, we are ready to run CellPhoneDB in R.

We provide the `run_cellphonedb()` function that can give us output from
CellPhoneDB analysis ([method
2](https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T1_Method2.ipynb)).
`run_cellphonedb()` can take in the same arguments as its Python-based
function, and users may refer to the linked “method 2” for details on
what each argument does.

This function takes in a Seurat object with normalized counts in the
data layer of the RNA assay and annotated cell type labels in its
meta.data.

You may see [here](../R/run_cellphonedb.R) for full details of
`run_cellphonedb()`.

``` r
# Load your Seurat object
data(seu.NL)

seu.NL
```

    ## An object of class Seurat 
    ## 10353 features across 2233 samples within 1 assay 
    ## Active assay: RNA (10353 features, 0 variable features)
    ##  1 layer present: data

``` r
seu.NL@meta.data %>% head()
```

    ##                     patient.id      labels
    ## S3_ATGAGGGAGTCTTGCA   Patient3   FBN1+ FIB
    ## S3_GTACGTACAAATTGCC   Patient3   FBN1+ FIB
    ## S3_CTCGTCAGTGTTGAGG   Patient3 Inflam. FIB
    ## S3_ATTCTACGTAATCGTC   Patient3   FBN1+ FIB
    ## S3_CTGCCTATCAATCACG   Patient3   FBN1+ FIB
    ## S3_TAGTGGTAGGATGCGT   Patient3   FBN1+ FIB

## Running CellPhoneDB with all default parameters

Minimally, `scpeaker()` takes in a Seurat object, the metadata column
name corresponding to the annotated cell types and the cell-cell
communication method (in this case cellphonedb). It creates a temporary
directory (by default) that stores input files required to run
CellPhoneDB, as well as the output .txt files from the analysis. Users
may make a copy of the temporary directory for future use of the data if
necessary. We implemented the triMean gene expression averaging strategy
within CellPhoneDB, which is the default averaging method with
`scpeaker()`.

``` r
cpdb <- scpeaker(obj = seu.NL, labels = "labels", method = "cellphonedb")
# Running cell-cell communication analysis with database: scpeakerDB
# use_dir is NULL. Creating temp directory for file creation and storage.
# Creating input and output files to directory: C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj 
#  Note: Make a copy of the temp directory for future use if necessary. 
# input_meta.tsv file created and saved to: C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj 
# input.h5ad file created and saved to: C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj 
# Running CellPhoneDB Statistical Analysis with - triMean
# Reading user files...
# The following user files were loaded successfully:
# C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj/input.h5ad
# C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj/input_meta.tsv
# [ ][CORE][02/07/26-17:33:05][INFO] type = triMean: Setting threshold to 0.
# [ ][CORE][02/07/26-17:33:05][INFO] [Cluster Statistical Analysis] Threshold:0 Iterations:1000 Debug-seed:-1 Threads:4 Precision:3
# [ ][CORE][02/07/26-17:33:06][INFO] Running Real Analysis
# [ ][CORE][02/07/26-17:33:06][INFO] Running Statistical Analysis
# 100%|██████████| 1000/1000 [00:32<00:00, 30.69it/s][ ][CORE][02/07/26-17:33:39][INFO] Building Pvalues result
# [ ][CORE][02/07/26-17:33:39][INFO] Building results
# 
# Saved deconvoluted to C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj\statistical_analysis_deconvoluted_07_02_2026_173339.txt
# Saved deconvoluted_percents to C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj\statistical_analysis_deconvoluted_percents_07_02_2026_173339.txt
# Saved means to C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj\statistical_analysis_means_07_02_2026_173339.txt
# Saved pvalues to C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj\statistical_analysis_pvalues_07_02_2026_173339.txt
# Saved significant_means to C:\Users\Admin\AppData\Local\Temp\Rtmpknehon/obj\statistical_analysis_significant_means_07_02_2026_173339.txt
# CellPhoneDB analysis completed successfully.
# Warning messages:
# 1: In dir.create(prefix %>% normalizePath(winslash = "/"), recursive = TRUE) :
#   'C:\Users\Admin\AppData\Local\Temp\Rtmpknehon\obj' already exists
# 2: In py_to_r.pandas.core.frame.DataFrame(<environment>) :
#   index contains duplicated values: row names not set
# 3: In py_to_r.pandas.core.frame.DataFrame(<environment>) :
#   index contains duplicated values: row names not set
```

## Running CellPhoneDB with user-defined directory

Should users desire for the files to be created to their directory of
choice, simply add the path to the directory of interest in the argument
`use_dir`.

``` r
cpdb <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb", use_dir = "cpdb/results")
```

## Running CellPhoneDB with adjustable parameters

``` r
cpdb <- scpeaker(seu.NL,
                 labels = "labels",
                 method = "cellphonedb",
                 type = "thresholdedMean",
                 iterations = 123,
                 threshold = 0.2,
                 threads = 1,
                 debug_seed = 42,
                 result_precision = 5,
                 score_interactions = TRUE)
```
