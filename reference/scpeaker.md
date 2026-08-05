# Run CellChat (v2) and/or CellPhoneDB (v5) for cell-cell communication analysis

This function takes as input a Seurat object that has been processed and
contains annotated cell type labels in its meta.data.

It uses scpeakerDB which contains interactions from both CellChat and
CellPhoneDB as the default database, but users may customize the
resource to CellChat/CellPhoneDB only interactions. Note that the
databases used are post-harmonized and slightly differ in the format of
their respective original versions.

## Usage

``` r
scpeaker(
  obj,
  labels,
  method = c("cellchat", "cellphonedb"),
  database = c("scpeakerDB", "CCDB", "CPDB"),
  pvalue = 0.05,
  type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
  assay = "RNA",
  subsetDB = FALSE,
  search = c(),
  key = "annotation",
  non_protein = FALSE,
  threshold = 0.1,
  LR.use = NULL,
  raw.use = TRUE,
  population.size = FALSE,
  distance.use = TRUE,
  interaction.range = 250,
  scale.distance = 0.01,
  k.min = 10,
  contact.dependent = TRUE,
  contact.range = NULL,
  contact.knn.k = NULL,
  contact.dependent.forced = FALSE,
  do.symmetric = TRUE,
  nboot = 100,
  seed.use = 1L,
  Kh = 0.5,
  n = 1,
  min.cells = 10,
  fast.mode = TRUE,
  use_dir = NULL,
  ...,
  counts_data = "hgnc_symbol",
  active_tfs_file_path = NULL,
  microenvs_file_path = NULL,
  score_interactions = FALSE,
  subsampling = FALSE,
  subsampling_log = FALSE,
  separator = "|",
  debug = FALSE,
  output_suffix = NULL
)
```

## Arguments

- obj:

  Seurat object

- labels:

  Metadata column name for the cell type labels

- method:

  Cell-cell communication method for analysis (CellChat or CellPhoneDB)

- database:

  Database to use for cell-cell communication analysis (default:
  scpeakerDB)

- pvalue:

  p-value threshold for significance testing (default: 0.05)

- type:

  method for computing average gene expression per cell type cluster
  (default: triMean)

- assay:

  RNA assay by default

- subsetDB:

  set to TRUE if filtering scpeakerDB (default = FALSE). If TRUE, use
  search, key and non_protein arguments as in help(subsetDB)

- search:

  Taken from CellChat package, run help(subsetDB) for details

- key:

  Taken from CellChat package, run help(subsetDB) for details

- non_protein:

  Taken from CellChat package, run help(subsetDB) for details

- threshold:

  Taken from CellChat package and replaces "trim" argument, run
  help(computeCommunProb) for details

- LR.use:

  Taken from CellChat package, run help(computeCommunProb) for details

- raw.use:

  Taken from CellChat package, run help(computeCommunProb) for details

- population.size:

  Taken from CellChat package, run help(computeCommunProb) for details

- distance.use:

  Taken from CellChat package, run help(computeCommunProb) for details

- interaction.range:

  Taken from CellChat package, run help(computeCommunProb) for details

- scale.distance:

  Taken from CellChat package, run help(computeCommunProb) for details

- k.min:

  Taken from CellChat package, run help(computeCommunProb) for details

- contact.dependent:

  Taken from CellChat package, run help(computeCommunProb) for details

- contact.range:

  Taken from CellChat package, run help(computeCommunProb) for details

- contact.knn.k:

  Taken from CellChat package, run help(computeCommunProb) for details

- contact.dependent.forced:

  Taken from CellChat package, run help(computeCommunProb) for details

- do.symmetric:

  Taken from CellChat package, run help(computeCommunProb) for details

- nboot:

  Taken from CellChat package, run help(computeCommunProb) for details

- seed.use:

  Taken from CellChat package, run help(computeCommunProb) for details

- Kh:

  Taken from CellChat package, run help(computeCommunProb) for details

- n:

  Taken from CellChat package, run help(computeCommunProb) for details

- min.cells:

  Taken from CellChat package, run help(filterCommunication) for details

- fast.mode:

  This argument only applies to CellChat analyses in the current
  version. It runs C++ version of computing average gene expression to
  increase computation speed without changes to prob and pval scores
  (default: TRUE)

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

Cell-cell communication inference results from the chosen method
(CellChat object for CellChat analysis and CellPhoneDB object list).

## CellChat arguments

## CellPhoneDB arguments

## Examples

``` r
if (FALSE) { # \dontrun{
cellchat <- scpeaker(seu.NL, labels = "labels", method = "cellchat")

cpdb     <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")
} # }
```
