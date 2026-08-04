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
#' cpdb     <- scpeaker(seurat_obj, method = "cellphonedb", ...)
#' cellchat <- scpeaker(seurat_obj, method = "cellchat",    ...)
#'
#' combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)
#'
#' cellchat.new <- filter_cellchat(cellchat = cellchat, crosscheck_res = crosscheck)
#' }


filter_cellchat <- function(cellchat, crosscheck_res, category = c("Sig_Both", "Sig_CellChat_Not_Found_CellPhoneDB", "Sig_CellPhoneDB_Not_Found_CellChat",
                                                                   "Sig_CellChat_Not_Sig_CellPhoneDB", "Sig_CellPhoneDB_Not_Sig_CellChat")) {
  if (! inherits(cellchat, "CellChat")) {
    stop("Input CellChat result is not a CellChat object. Perform run_cellchat() on your Seurat object and set cellchat as your output CellChat object.")
  } else if (! is.data.frame(crosscheck_res$result)) {
    stop("Input crosscheck_res is not a list. Run crosscheck() on your CellChat object and CellPhoneDB output list and extract the result element.")
  }

  var <- match.arg(category)

  cat("Filtering interactions in category:", var, "\n")

  prob <- cellchat@net$prob
  pval <- cellchat@net$pval

  to_keep <- crosscheck_res$result %>%
    filter(category == var) %>%
    select(source = source_cc, target = target_cc, interaction = interaction_name)

  keep_mask <- array(
    FALSE,
    dim = dim(prob),
    dimnames = dimnames(prob)
  )

  for (i in seq_len(nrow(to_keep))) {
    source      <- to_keep$source[i]
    target      <- to_keep$target[i]
    interaction <- to_keep$interaction[i]

    keep_mask[source, target, interaction] <- TRUE
  }

  cellchat.new <- cellchat

  cellchat.new@net$prob[!(keep_mask)] <- 0
  cellchat.new@net$pval[!(keep_mask)] <- 1

  cellchat.new <- computeCommunProbPathway(cellchat.new)
  cellchat.new <- aggregateNet(cellchat.new)
  cellchat.new <- netAnalysis_computeCentrality(cellchat.new, slot.name = "netP")

  return(cellchat.new)

}



