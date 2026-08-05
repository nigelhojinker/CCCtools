# Run CellPhoneDB version 5 in R

Runs CellPhoneDB (v5.0.0) statistical analysis on a Seurat object.

## Usage

``` r
run_cellphonedb(
  obj = NULL,
  labels = NULL,
  type = NULL,
  database = NULL,
  use_dir = NULL,
  ...,
  counts_data = "hgnc_symbol",
  active_tfs_file_path = NULL,
  microenvs_file_path = NULL,
  score_interactions = FALSE,
  threshold = 0.1,
  pvalue = 0.05,
  subsampling = FALSE,
  subsampling_log = FALSE,
  separator = "|",
  debug = FALSE,
  output_suffix = NULL
)
```

## Arguments

- obj:

  Seurat object with normalized counts in the data layer of the RNA
  assay

- labels:

  Metadata column name of cell type annotations

- type:

  Method to compute average gene expression per cell type cluster
  (default: triMean)

- database:

  Choice of database to use for CellChat analysis

- use_dir:

  Input files required to run CellPhoneDB and the output files from the
  analysis are created and saved to a temp file (default). Should users
  desire to have these files be stored in their personal directories,
  they may input their filepath to this argument.

- ...:

  Takes in any arguments from c("iterations", "threads", "debug_seed",
  "result_precision", "subsampling_num_pc", "subsampling_num_cells") for
  user customization, else running on the default values for these
  paramters. See Value below for link to description by CellPhoneDB

- counts_data:

  See Value below for link to description by CellPhoneDB

- active_tfs_file_path:

  See Value below for link to description by CellPhoneDB

- microenvs_file_path:

  See Value below for link to description by CellPhoneDB

- score_interactions:

  See Value below for link to description by CellPhoneDB

- threshold:

  See Value below for link to description by CellPhoneDB

- pvalue:

  See Value below for link to description by CellPhoneDB

- subsampling:

  See Value below for link to description by CellPhoneDB

- subsampling_log:

  See Value below for link to description by CellPhoneDB

- separator:

  See Value below for link to description by CellPhoneDB

- debug:

  See Value below for link to description by CellPhoneDB

- output_suffix:

  See Value below for link to description by CellPhoneDB

## Value

Folder path to CellPhoneDB analysis outputs and result as a list of the
outputs in the R session

For arguments from counts_data to output_suffix, please refer to
https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T1_Method2.ipynb

## Examples

``` r
if (FALSE) { # \dontrun{
## Runs CellPhoneDB statistical analysis with all default parameters
run_cellphonedb(
  obj = seu.NL,
  labels = "labels")

## For users interested in customizing their CellPhoneDB run
run_cellphonedb(seu.NL,
  labels = "labels",
  iterations = 100,
  threshold = 0.2,
  threads = 5,
  debug_seed = 42,
  result_precision = 5,
  score_interactions = TRUE)
} # }
```
