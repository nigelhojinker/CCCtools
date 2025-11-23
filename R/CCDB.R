#' Mapped CellChat (v2) interaction database
#'
#' This database is adapted from the original CellChat (v2) database;
#' it contains information from the original package, but includes uniprot identifiers for comparison with CellPhoneDB results.
#'
#' @format A dataframe with 3233 interactions and 33 columns:
#' \describe{
#'   \item{interaction_name}{CellChat interaction name}
#'   \item{pathway_name}{CellChat pathway name}
#'   \item{ligand_name}{Gene symbol/name of ligand}
#'   \item{ligand}{Uniprot ID format of ligand}
#'   \item{receptor_name}{Gene symbol/name of receptor}
#'   \item{receptor}{Uniprot ID format of receptor}
#'   \item{LR}{Ligand-receptor pair of interaction in Uniprot format; ligand|receptor separated by "|"}
#'   \item{complex_ligand}{Logical value indicating whether ligand is a complex or not}
#'   \item{complex_receptor}{Logical value indicating whether receptor is a complex or not}
#'   \item{agonist}{CellChat agonist - facilitates interaction}
#'   \item{antagonist}{CellChat antagonist - inhibits interaction}
#'   \item{co_A_receptor}{CellChat co-activation receptor}
#'   \item{co_I_receptor}{CellChat co-inhibition receptor}
#'   \item{annotation}{CellChat interaction category}
#'   \item{interaction_name_2}{CellChat interaction name 2}
#'   \item{evidence}{CellChat interaction source}
#'   \item{is_neurotransmitter}{TRUE if interaction involved in neurotransmission}
#'   \item{ligand.symbol}{CellChat gene symbol(s) of ligand}
#'   \item{ligand.family}{CellChat annotation of ligand protein information}
#'   \item{ligand.location}{Where ligand is found biologically}
#'   \item{ligand.keyword}{CellChat annotation of ligand properties}
#'   \item{ligand.secreted_type}{Ligand class if secreted}
#'   \item{ligand.transmembrane}{TRUE if transmembrane}
#'   \item{receptor.symbol}{CellChat gene symbol(s) of receptor}
#'   \item{receptor.family}{CellChat annotation of receptor protein information}
#'   \item{receptor.location}{Where receptor is found biologically}
#'   \item{receptor.keyword}{CellChat annotation of receptor properties}
#'   \item{receptor.surfaceome_main}{Receptor type}
#'   \item{receptor.surfaceome_sub}{Receptor type breakdown}
#'   \item{receptor.adhesome}{Patient identifier}
#'   \item{receptor.secreted_type}{Receptor class if secreted}
#'   \item{receptor.transmembrane}{TRUE if transmembrane}
#'   \item{version}{Version of CellChat interaction was uploaded in}
#' }
#'
#' @examples
#'
#' data(CCDB)
#'
#' CCDB %>% view()
#'
#'
#' @source
#' Database was adapted from the CellChat (v2) R package HUMAN database (CellChatDB.human$interaction).
"CCDB"
