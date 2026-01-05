#' Run CellChat (v2) and/or CellPhoneDB (v5) for cell-cell communication analysis
#'
#' @description
#' This function takes as input a Seurat object that has been processed and contains annotated cell type labels in its meta.data.
#' @description
#' It uses scpeakerDB which contains interactions from both CellChat and CellPhoneDB as the default database, but users may customize the resource to CellChat/CellPhoneDB only interactions.
#'
#' @param obj Seurat object
#' @param labels Metadata column name for the cell type labels
#' @param method Cell-cell communication method for analysis (CellChat or CellPhoneDB)
#' @param database Database to use for cell-cell communication analysis (default: scpeakerDB)
#' @param pvalue p-value threshold for significance testing (default: 0.05)
#' @param type method for computing average gene expression per cell type cluster (default: triMean)
#'
#' @section CellChat arguments:
#'
#' @param assay RNA assay by default
#' @param subsetDB set to TRUE if filtering scpeakerDB (default = FALSE). If TRUE, use search, key and non_protein arguments as in help(subsetDB)
#' @param search taken from CellChat package, run help(subsetDB) for details
#' @param key taken from CellChat package, run help(subsetDB) for details
#' @param non_protein taken from CellChat package, run help(subsetDB) for details
#' @param threshold taken from CellChat package and replaces "trim" argument, run help(computeCommunProb) for details
#' @param LR.use taken from CellChat package, run help(computeCommunProb) for details
#' @param raw.use taken from CellChat package, run help(computeCommunProb) for details
#' @param population.size taken from CellChat package, run help(computeCommunProb) for details
#' @param distance.use taken from CellChat package, run help(computeCommunProb) for details
#' @param interaction.range taken from CellChat package, run help(computeCommunProb) for details
#' @param scale.distance taken from CellChat package, run help(computeCommunProb) for details
#' @param k.min taken from CellChat package, run help(computeCommunProb) for details
#' @param contact.dependent taken from CellChat package, run help(computeCommunProb) for details
#' @param contact.range taken from CellChat package, run help(computeCommunProb) for details
#' @param contact.knn.k taken from CellChat package, run help(computeCommunProb) for details
#' @param contact.dependent.forced taken from CellChat package, run help(computeCommunProb) for details
#' @param do.symmetric taken from CellChat package, run help(computeCommunProb) for details
#' @param nboot taken from CellChat package, run help(computeCommunProb) for details
#' @param seed.use taken from CellChat package, run help(computeCommunProb) for details
#' @param Kh taken from CellChat package, run help(computeCommunProb) for details
#' @param n taken from CellChat package, run help(computeCommunProb) for details
#' @param min.cells taken from CellChat package, run help(filterCommunication) for details
#'
#' @section CellPhoneDB arguments:
#'
#' @param use_dir Input files required to run CellPhoneDB and the output files from the analysis are created and saved to a temp file (default). Should users desire to have these files be stored in their personal directories, they may input their filepath to this argument.
#' @param ... Takes in any arguments from c("iterations", "threads", "debug_seed", "result_precision", "subsampling_num_pc", "subsampling_num_cells") for user customization, else running on the default values for these paramters. See Value below for link to description by CellPhoneDB
#' @param counts_data See Value below for link to description by CellPhoneDB
#' @param active_tfs_file_path See Value below for link to description by CellPhoneDB
#' @param microenvs_file_path See Value below for link to description by CellPhoneDB
#' @param score_interactions See Value below for link to description by CellPhoneDB
#' @param subsampling See Value below for link to description by CellPhoneDB
#' @param subsampling_log See Value below for link to description by CellPhoneDB
#' @param separator See Value below for link to description by CellPhoneDB
#' @param debug See Value below for link to description by CellPhoneDB
#' @param output_suffix See Value below for link to description by CellPhoneDB
#'
#' @returns Cell-cell communication inference results from the chosen method (CellChat object for CellChat analysis and CellPhoneDB object list).
#' @export
#'
#' @examples
#' \dontrun{
#' cellchat <- scpeaker(seu.NL, labels = "labels", method = "cellchat")
#'
#' cpdb     <- scpeaker(seu.NL, labels = "labels", method = "cellphonedb")
#' }
#'

scpeaker <- function(obj, labels, method = c("cellchat", "cellphonedb"), database = c("scpeakerDB", "CCDB", "CPDB"), pvalue = 0.05,
                     assay = "RNA", subsetDB = FALSE, search = c(), key = "annotation", non_protein = FALSE,
                     type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
                     threshold = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE, distance.use = TRUE,
                     interaction.range = 250, scale.distance = 0.01, k.min = 10, contact.dependent = TRUE,
                     contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE,
                     do.symmetric = TRUE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1, min.cells = 10,
                     use_dir = NULL, ..., counts_data = "hgnc_symbol", active_tfs_file_path = NULL, microenvs_file_path = NULL,
                     score_interactions = FALSE, subsampling = FALSE,
                     subsampling_log = FALSE, separator = "|", debug = FALSE, output_suffix = NULL) {

  if (length(method) > 1) stop("Length of method is more than 1.\n Please indicate a method for cell-cell communication analysis - cellchat or cellphonedb")
  if (! method %in% c("cellchat", "cellphonedb")) stop("Invalid input for method.\n Please indicate a method for cell-cell communication analysis - cellchat or cellphonedb")

  database <- match.arg(database)

  if (! database %in% c("scpeakerDB", "CCDB", "CPDB")) stop("Invalid input for database.\n Please indicate a database - scpeakerDB, CCDB, or CPDB")

  type <- match.arg(type)

  if (method == "cellchat"){
    return(run_cellchat(obj, labels = labels, assay = assay,
                        subsetDB = subsetDB, search = search, key = key, non_protein = non_protein, type = type, threshold = threshold,
                        LR.use = LR.use, raw.use = raw.use, population.size = population.size, distance.use = distance.use,
                        interaction.range = interaction.range, scale.distance = scale.distance, k.min = k.min,
                        contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k,
                        contact.dependent.forced = contact.dependent.forced, do.symmetric = do.symmetric, nboot = nboot,
                        seed.use = seed.use, Kh = Kh, n = n, min.cells = min.cells))
  } else if (method == "cellphonedb") {
    return(run_cellphonedb(obj = obj, labels = labels, type = type, use_dir = use_dir, ..., counts_data = counts_data,
                           active_tfs_file_path = active_tfs_file_path, microenvs_file_path = microenvs_file_path,
                           score_interactions = score_interactions, threshold = threshold, pvalue = pvalue,
                           subsampling = subsampling, subsampling_log = subsampling_log, separator = separator, debug = debug,
                           output_suffix = output_suffix))
  }
}
