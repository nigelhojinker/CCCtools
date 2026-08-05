# Mapped CellPhoneDB (v5) interaction database

This database is adapted from the original CellPhoneDB (v5) database; it
contains information from the original package, but includes uniprot
identifiers for comparison with CellChat results.

## Usage

``` r
CPDB
```

## Format

A dataframe with 2911 interactions and 25 columns:

- id_cp_interaction:

  CellPhoneDB unique interaction ID

- classification:

  Interaction classification

- ligand_name:

  Gene symbol/name of ligand

- ligand:

  Uniprot ID format of ligand

- receptor_name:

  Gene symbol/name of receptor

- receptor:

  Uniprot ID format of receptor

- LR:

  Ligand-receptor pair of interaction in Uniprot format;
  ligand\|receptor separated by "\|"

- complex_ligand:

  Logical value indicating whether ligand is a complex or not

- complex_receptor:

  Logical value indicating whether receptor is a complex or not

- partner_a:

  Uniprot ID format of ligand in original CellPhoneDB

- partner_b:

  Uniprot ID format of receptor in original CellPhoneDB

- protein_name_a:

  Gene symbol/name of ligand in original CellPhoneDB

- protein_name_b:

  Gene symbol/name of receptor in original CellPhoneDB

- annotation_strategy:

  Method of collecting and annotating interaction

- source:

  Source of which interaction was taken from

- is_ppi:

  TRUE if protein-protein interaction

- curator:

  Curator

- reactome_complex:

  Reactome stable identifier for complex

- reactome_reaction:

  Reactome stable identifier for reaction

- reactome_pathway:

  Reactome stable identifier for pathway

- comments:

  Additional information or descriptions on the interaction

- version:

  Version history of interaction in the original CellPhoneDB

- interactors:

  Ligand-receptor pair in gene symbol format

- directionality:

  Interaction type

- modulatory_effect:

  Indicates the type of effect (activatory or inhibitory) interaction
  has

## Source

Database was adapted from the CellPhoneDB (v5) Python package
interaction_table from the zip file at:
https://github.com/ventolab/cellphonedb-data/blob/master/cellphonedb.zip.

## Examples

``` r

data(CPDB)

CPDB %>% view()

```
