# Mapped CellChat (v2) interaction database

This database is adapted from the original CellChat (v2) database; it
contains information from the original package, but includes uniprot
identifiers for comparison with CellPhoneDB results.

## Usage

``` r
CCDB
```

## Format

A dataframe with 3233 interactions and 33 columns:

- interaction_name:

  CellChat interaction name

- pathway_name:

  CellChat pathway name

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

- agonist:

  CellChat agonist - facilitates interaction

- antagonist:

  CellChat antagonist - inhibits interaction

- co_A_receptor:

  CellChat co-activation receptor

- co_I_receptor:

  CellChat co-inhibition receptor

- annotation:

  CellChat interaction category

- interaction_name_2:

  CellChat interaction name 2

- evidence:

  CellChat interaction source

- is_neurotransmitter:

  TRUE if interaction involved in neurotransmission

- ligand.symbol:

  CellChat gene symbol(s) of ligand

- ligand.family:

  CellChat annotation of ligand protein information

- ligand.location:

  Where ligand is found biologically

- ligand.keyword:

  CellChat annotation of ligand properties

- ligand.secreted_type:

  Ligand class if secreted

- ligand.transmembrane:

  TRUE if transmembrane

- receptor.symbol:

  CellChat gene symbol(s) of receptor

- receptor.family:

  CellChat annotation of receptor protein information

- receptor.location:

  Where receptor is found biologically

- receptor.keyword:

  CellChat annotation of receptor properties

- receptor.surfaceome_main:

  Receptor type

- receptor.surfaceome_sub:

  Receptor type breakdown

- receptor.adhesome:

  Patient identifier

- receptor.secreted_type:

  Receptor class if secreted

- receptor.transmembrane:

  TRUE if transmembrane

- version:

  Version of CellChat interaction was uploaded in

## Source

Database was adapted from the CellChat (v2) R package HUMAN database
(CellChatDB.human\$interaction).

## Examples

``` r

data(CCDB)

CCDB %>% view()

```
