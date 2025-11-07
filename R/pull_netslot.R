#' Extract all relevant interactions from CellChat analysis
#'
#' @description
#' Adapted from CellChat's subsetCommunication() function, but retains all interactions (significant or insignificant) less those that do not pass min.cells criteria set by filterCommunication()
#'
#' @param object CellChat object after run_cellchat()
#' @param thresh Threshold of p-value for defining interactions that are statistical significant
#'
#' @returns A dataframe consisting of inferred interactions for comparison with CellPhoneDB
#' @export
#'
#' @examples
#' \dontrun{
#' cellchat.res <- pull_netslot(cellchat)
#' cellchat.res %>% view()
#' }

pull_netslot <- function (object = NULL, thresh = 0.05) {
  net <- slot(object, "net")

  LR <- object@LR$LRsig

  if (!is.data.frame(net)) {
    prob <- net$prob
    pval <- net$pval
    prob[pval >= thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source", "target", "interaction_name")
    net.pval <- reshape2::melt(pval, value.name = "pval")
    net$pval <- net.pval$pval
    # net <- subset(net, prob > 0) # NOT RUN this line
  }

  if (!("ligand" %in% colnames(net))) {
    col.use <- intersect(c("interaction_name_2", "pathway_name",
                           "ligand", "receptor", "annotation", "evidence"),
                         colnames(LR))
    pairLR <- dplyr::select(LR, col.use)
    idx <- match(net$interaction_name, rownames(pairLR))
    net <- cbind(net, pairLR[idx, ])
  }

  # remove interactions that do not pass criteria (min.cells) in filterCommunication()
  return(net %>% filter( !(prob == 0 & pval < 0.05)))
}
