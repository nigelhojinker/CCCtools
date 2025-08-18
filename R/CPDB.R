#' Mapped CellPhoneDB (v5) interaction database
#'
#' This database is adapted from the original CellPhoneDB (v5) database;
#' it contains information from the original package, but includes uniprot identifiers for
#' comparison with CellChat results.
#'
#' @format A dataframe with 2911 interactions and 13 columns:
#' \describe{
#'   \item{id_interaction}{Interaction number from 0 to 2910; total 2911 interactions}
#'   \item{id_cp_interaction}{CellPhoneDB unique interaction ID}
#'   \item{source}{Source of which interaction was taken from}
#'   \item{annotation_strategy}{Method of collecting and annotating interaction}
#'   \item{is_ppi}{TRUE if protein-protein interaction}
#'   \item{curator}{Curator}
#'   \item{directionality}{Interaction type}
#'   \item{classification}{Interaction classification}
#'   \item{unique_int}{Uniprot-based interaction identifier}
#'   \item{ligand}{Uniprot ID of ligand}
#'   \item{from}{Gene symbol of ligand}
#'   \item{receptor}{Uniprot ID of receptor}
#'   \item{to}{Gene symbol of receptor}
#' }
#'
#' @examples
#'
#' data(CPDB)
#'
#' CPDB %>% view()
#'
#'
#' @source
#' Database was adapted from the CellPhoneDB (v5) Python package [interaction_table](https://github.com/ventolab/cellphonedb-data/blob/master/cellphonedb.zip) and worked on as described [here](https://cellphonedb-cellchatdb-mapping.vercel.app/)
"CPDB"
