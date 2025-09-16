#' Subsetting CellChat object to contain interactions significant in both CCC tools
#'
#' @param cellchat CellChat object obtained from run_cellchat() on a Seurat object.
#' @param combine.sig Output data from crosscheck() result element.
#'
#' @returns A new CellChat object that contains only interactions significant in both CellChat and CellPhoneDB.
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
#' cellchat.new <- filter_cellchat(cellchat = cellchat, combine.sig = combine$result)
#' }


filter_cellchat <- function(cellchat, combine.sig) {
  if (! is.data.frame(pull_netslot(cellchat))) {
    stop("Input CellChat result is not a dataframe. Perform run_cellchat() on your Seurat object and set cellchat as your output CellChat object.")
  } else if (! is.data.frame(combine.sig)) {
    stop("Input combine.sig is not a dataframe. Run crosscheck() on your CellChat object and CellPhoneDB output list and extract the result element.")
  }
  
  rel.int <- combine.sig %>%
    filter(confidence == "Sig_Both") %>%
    pull(interaction_name) %>%
    unique()
  
  cellchat.new <- cellchat
  cellchat.new@net$prob <- cellchat@net$prob[,, rel.int]
  cellchat.new@net$pval <- cellchat@net$pval[,, rel.int]
  cellchat.new@LR$LRsig <- cellchat@LR$LRsig %>%
    filter(interaction_name %in% rel.int)
  
  cellchat.new <- computeCommunProbPathway(cellchat.new)
  cellchat.new <- aggregateNet(cellchat.new)
  cellchat.new <- netAnalysis_computeCentrality(cellchat.new, slot.name = "netP")
  
  cat("Created CellChat object with only interactions significant in both CellChat and CellPhoneDB.")
  
  return(cellchat.new)
  
}



