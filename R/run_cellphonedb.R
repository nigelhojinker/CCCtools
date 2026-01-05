#' Run CellPhoneDB version 5 in R
#'
#' @description
#' Runs CellPhoneDBv5 on a Seurat object using the customized cellphonedbv5a database containing additional interactions from CellChat.
#'
#' @param obj Seurat object with normalized counts in the data layer of the RNA assay
#' @param labels Metadata column name of cell type annotations
#' @param type Method to compute average gene expression per cell type cluster (default: triMean)
#' @param use_dir Input files required to run CellPhoneDB and the output files from the analysis are created and saved to a temp file (default). Should users desire to have these files be stored in their personal directories, they may input their filepath to this argument.
#' @param ... Takes in any arguments from c("iterations", "threads", "debug_seed", "result_precision", "subsampling_num_pc", "subsampling_num_cells") for user customization, else running on the default values for these paramters. See Value below for link to description by CellPhoneDB
#' @param counts_data See Value below for link to description by CellPhoneDB
#' @param active_tfs_file_path See Value below for link to description by CellPhoneDB
#' @param microenvs_file_path See Value below for link to description by CellPhoneDB
#' @param score_interactions See Value below for link to description by CellPhoneDB
#' @param threshold See Value below for link to description by CellPhoneDB
#' @param pvalue See Value below for link to description by CellPhoneDB
#' @param subsampling See Value below for link to description by CellPhoneDB
#' @param subsampling_log See Value below for link to description by CellPhoneDB
#' @param separator See Value below for link to description by CellPhoneDB
#' @param debug See Value below for link to description by CellPhoneDB
#' @param output_suffix See Value below for link to description by CellPhoneDB
#'
#' @returns Folder path to CellPhoneDB analysis outputs and result as a list of the outputs in the R session
#'
#' @returns For arguments from counts_data to output_suffix, please refer to https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T1_Method2.ipynb
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Runs CellPhoneDB method 2 with all default parameters
#' run_cellphonedb(
#'   obj = seu.NL,
#'   labels = "labels")
#'
#' ## For users interested in customizing their CellPhoneDB run
#' run_cellphonedb(seu.NL,
#'   labels = "labels",
#'   iterations = 100,
#'   threshold = 0.2,
#'   threads = 5,
#'   debug_seed = 42,
#'   result_precision = 5,
#'   score_interactions = TRUE)
#' }

run_cellphonedb <- function(obj = NULL, labels = NULL, type = NULL, use_dir = NULL, ...,
                                counts_data = "hgnc_symbol", active_tfs_file_path = NULL, microenvs_file_path = NULL,
                                score_interactions = FALSE, threshold = 0.1, pvalue = 0.05, subsampling = FALSE,
                                subsampling_log = FALSE, separator = "|", debug = FALSE, output_suffix = NULL) {

  if (is.null(obj) || is.null(labels)) {
      stop("obj and labels cannot be NULL. Please input a Seurat object (obj) and a metadata column name (labels).\n")
  }

  if (type == "truncatedMean" | type == "median") {
    stop("truncatedMean and median are currently not available as a mean method. Use either thresholdedMean or triMean.")
  } else if (! type %in% c("thresholdedMean", "triMean")) {
    stop("Unknown type method, use either thresholdedMean or triMean.")
  }

  cpdbPath <- file.path(pacman::p_path("scpeakeR"), "data/scpeakerDB_CP.zip")

  # Use use_dir if provided, else use tempdir with deparse
  if (!is.null(use_dir)) {
    message("use_dir is not NULL. Using user-input directory for file creation and storage.")
    prefix <- use_dir
  } else {
    message("use_dir is NULL. Creating temp directory for file creation and storage.")
    prefix <- file.path(tempdir(), deparse(substitute(obj)))
  }
  dir.create(prefix %>% normalizePath(winslash = "/"), recursive = TRUE)

  cat("Creating input and output files to directory:", prefix, "\n",
      "Note: Make a copy of the temp directory for future use if necessary. \n")

  ## Convert arguments to integers for CellPhoneDB python function
  dots <- list(...)
  to_integer <- c("iterations", "threads", "debug_seed", "result_precision", "subsampling_num_pc", "subsampling_num_cells")

  args <- intersect(names(dots), to_integer)

  dots[args] <- lapply(dots[args], as.integer)

  ## Load required modules
  ad <- import("anndata")
  cpdb.analysis <- import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")

  ## Input file creation
  meta <- obj@meta.data %>%
    rownames_to_column(var = "barcode") %>%
    select(barcode, !!rlang::sym(labels)) %>%
    write_delim(file = file.path(prefix, "input_meta.tsv"), delim = "\t")

  cat("input_meta.tsv file created and saved to:", prefix, "\n")

  normcounts <- obj[["RNA"]]$data

  stopifnot(identical(meta$barcode, colnames(normcounts)))

  genes <- data.frame(gene = rownames(normcounts)) %>%
    column_to_rownames("gene")

  barcode <- data.frame(barcode = colnames(normcounts)) %>%
    column_to_rownames("barcode")

  ## Create AnnData object
  adata <- ad$AnnData(X = normcounts, obs = genes, var = barcode)
  adata <- adata$transpose()
  adata$write_h5ad(file.path(prefix, "input.h5ad"))

  cat("input.h5ad file created and saved to:", prefix, "\n")

  message("Running CellPhoneDB Statistical Analysis with - ", type, "\n")

  ## Run analysis
  cpdb <- do.call(cpdb.analysis$call, c(dots, list(
    cpdb_file_path = cpdbPath,
    meta_file_path = file.path(prefix, "input_meta.tsv"),
    counts_file_path = file.path(prefix, "input.h5ad"),
    counts_data = counts_data,
    active_tfs_file_path = active_tfs_file_path,
    microenvs_file_path = microenvs_file_path,
    score_interactions = score_interactions,
    threshold = threshold,
    pvalue = pvalue,
    subsampling = subsampling,
    subsampling_log = subsampling_log,
    separator = separator,
    debug = debug,
    output_path = prefix,
    output_suffix = output_suffix,
    avg_method = type
  )))

  message("CellPhoneDB analysis completed successfully.")

  # Assign CellPhoneDB class to cpdb out
  class(cpdb) <- "CellPhoneDB"

  return(cpdb)
}
