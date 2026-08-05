# Run CellChat (v2.2.0) in R

This function takes as input a Seurat object that has been processed and
contains annotated cell type labels in its meta.data.

It runs functions from the CellChat package, namely: createCellChat(),
subsetDB(), subsetData(), identifyOverExpressedGenes(),
identifyOverExpressedInteractions(), computeCommunProb() and
filterCommunication(). run_cellchat() allows users to input arguments
listed in createCellChat(), subsetDB(), computeCommunProb() and
filterCommunication().

## Usage

``` r
run_cellchat(
  obj,
  labels = "ident",
  assay = "RNA",
  database = NULL,
  subsetDB = FALSE,
  search = c(),
  key = "annotation",
  non_protein = FALSE,
  type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
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
  min.cells = 10
)
```

## Arguments

- obj:

  Seurat object

- labels:

  Metadata column name for the cell type labels

- assay:

  RNA assay by default

- database:

  Choice of database to use for CellChat analysis

- subsetDB:

  set to TRUE if filtering scpeakerDB (default = FALSE). If TRUE, use
  search, key and non_protein arguments as in help(subsetDB)

- search:

  taken from CellChat package, run help(subsetDB) for details

- key:

  taken from CellChat package, run help(subsetDB) for details

- non_protein:

  taken from CellChat package, run help(subsetDB) for details

- type:

  taken from CellChat package, run help(computeCommunProb) for details

- threshold:

  taken from CellChat package, run help(computeCommunProb) for details

- LR.use:

  taken from CellChat package, run help(computeCommunProb) for details

- raw.use:

  taken from CellChat package, run help(computeCommunProb) for details

- population.size:

  taken from CellChat package, run help(computeCommunProb) for details

- distance.use:

  taken from CellChat package, run help(computeCommunProb) for details

- interaction.range:

  taken from CellChat package, run help(computeCommunProb) for details

- scale.distance:

  taken from CellChat package, run help(computeCommunProb) for details

- k.min:

  taken from CellChat package, run help(computeCommunProb) for details

- contact.dependent:

  taken from CellChat package, run help(computeCommunProb) for details

- contact.range:

  taken from CellChat package, run help(computeCommunProb) for details

- contact.knn.k:

  taken from CellChat package, run help(computeCommunProb) for details

- contact.dependent.forced:

  taken from CellChat package, run help(computeCommunProb) for details

- do.symmetric:

  taken from CellChat package, run help(computeCommunProb) for details

- nboot:

  taken from CellChat package, run help(computeCommunProb) for details

- seed.use:

  taken from CellChat package, run help(computeCommunProb) for details

- Kh:

  taken from CellChat package, run help(computeCommunProb) for details

- n:

  taken from CellChat package, run help(computeCommunProb) for details

- min.cells:

  taken from CellChat package, run help(filterCommunication) for details

## Value

A CellChat object with computed communication probabilities on the
ligand-receptor level (net slot).

cellchat@net\$prob is a 3-dimensional array consisting of the
source(sender cell), target(receiver cell) and ligand-receptor
interaction as the first, second and third dimension respectively.

cellchat@net\$pval consists of the corresponding p-values computed for
each interaction.

## Examples

``` r
if (FALSE) { # \dontrun{
cellchat <- run_cellchat(seu.NL, group.by = "labels", assay = "RNA")
} # }
```
