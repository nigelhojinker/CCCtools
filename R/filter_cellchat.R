#' Subsetting CellChat object to contain interactions significant in both CCC tools
#'
#' @param cellchat CellChat object obtained from run_cellchat() on a Seurat object.
#' @param crosscheck_res Output data from crosscheck() function.
#' @param category Character value of the type of interaction to filter the cellchat object (default: significant in both CellChat and CellPhoneDB).
#'
#' @returns A new CellChat object that interactions filtered to your interest.
#' @export
#'
#' @examples
#' \dontrun{
#' ## After setting up conda env
#' cpdb <- run_cellphonedb(seurat.obj, ...)
#' cellchat <- run_cellchat(seurat.obj, ...)
#'
#' combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)
#'
#' cellchat.new <- filter_cellchat(cellchat = cellchat, crosscheck_res = crosscheck)
#' }


filter_cellchat <- function(cellchat, crosscheck_res, category = c("Sig_Both", "Sig_CellChat_Not_Found_CellPhoneDB", "Sig_CellPhoneDB_Not_Found_CellChat",
                                                                   "Sig_CellChat_Not_Sig_CellPhoneDB", "Sig_CellPhoneDB_Not_Sig_CellChat")) {
  if (! is.data.frame(pull_netslot(cellchat))) {
    stop("Input CellChat result is not a dataframe. Perform run_cellchat() on your Seurat object and set cellchat as your output CellChat object.")
  } else if (! is.data.frame(crosscheck_res$result)) {
    stop("Input crosscheck_res is not a list. Run crosscheck() on your CellChat object and CellPhoneDB output list and extract the result element.")
  }

  category <- match.arg(category)

  to_keep <- crosscheck_res$result %>%
    filter(category == category) %>%
    pull(interaction_name) %>%
    unique()

  cellchat.new <- cellchat
  cellchat.new@net$prob <- cellchat@net$prob[,, to_keep]
  cellchat.new@net$pval <- cellchat@net$pval[,, to_keep]
  cellchat.new@LR$LRsig <- cellchat@LR$LRsig %>%
    filter(interaction_name %in% to_keep)

  cellchat.new <- computeCommunProbPathway(cellchat.new)
  cellchat.new <- aggregateNet(cellchat.new)
  cellchat.new <- netAnalysis_computeCentrality(cellchat.new, slot.name = "netP")

  cat("Filtered CellChat object to only interactions in category:", category, "\n")

  return(cellchat.new)

}



