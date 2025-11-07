#' Run CellChat v2 in R
#'
#' @description
#' This function takes as input a Seurat object that has been processed and contains annotated cell type labels in its meta.data.
#' @description
#' It runs functions from the CellChat package, namely: createCellChat(), subsetDB(), subsetData(), identifyOverExpressedGenes(),
#' identifyOverExpressedInteractions(), computeCommunProb() and filterCommunication(). run_cellchat() allows users to input
#' arguments listed in createCellChat(), subsetDB(), computeCommunProb() and filterCommunication().
#'
#'
#' @param obj Seurat object
#' @param group.by Metadata column name for the cell type labels
#' @param assay RNA assay by default
#' @param subsetDB set to TRUE if filtering CellChatDB.human.new, our customized CellChatDB.human (default = FALSE). If TRUE, use search, key and non_protein arguments as in help(subsetDB)
#' @param search taken from CellChat package, run help(subsetDB) for details
#' @param key taken from CellChat package, run help(subsetDB) for details
#' @param non_protein taken from CellChat package, run help(subsetDB) for details
#' @param type taken from CellChat package, run help(computeCommunProb) for details
#' @param trim taken from CellChat package, run help(computeCommunProb) for details
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
#' @returns A CellChat object with computed communication probabilities on the ligand-receptor level (net slot).
#' @returns cellchat@net$prob is a 3-dimensional array consisting of the source(sender cell), target(receiver cell) and ligand-receptor interaction as the first, second and third dimension respectively.
#' @returns cellchat@net$pval consists of the corresponding p-values computed for each interaction.
#' @export
#'
#' @examples
#' \dontrun{
#' cellchat <- run_cellchat(seu.NL, group.by = "labels", assay = "RNA")
#' }

run_cellchat <- function(obj, group.by = "ident", assay = "RNA",
                         subsetDB = FALSE, search = c(), key = "annotation", non_protein = FALSE,
                         type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
                         trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE, distance.use = TRUE,
                         interaction.range = 250, scale.distance = 0.01, k.min = 10, contact.dependent = TRUE,
                         contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE,
                         do.symmetric = TRUE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1, min.cells = 10) {

  ## Create CellChat object from Seurat object
  cellchat <- createCellChat(obj, group.by = group.by, assay = assay)

  ## Database selection

  # Subset CellChat database if needed
  if (subsetDB){
    cellchat@DB <- subsetDB(CellChatDB.human.new, search = search, key = key, non_protein = non_protein)
  } else {
    cellchat@DB <- CellChatDB.human.new
  }

  ## Run analysis
  cellchat <- subsetData(cellchat) %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = type, trim = trim, LR.use = LR.use, raw.use = raw.use,
                      population.size = population.size, distance.use = distance.use,
                      interaction.range = interaction.range, scale.distance = scale.distance,
                      k.min = k.min, contact.dependent = contact.dependent,
                      contact.range = contact.range, contact.knn.k = contact.knn.k,
                      contact.dependent.forced = contact.dependent.forced,
                      do.symmetric = do.symmetric, nboot = nboot, seed.use = seed.use,
                      Kh = Kh, n = n) %>%
    filterCommunication(min.cells = min.cells)

  message("CellChat analysis completed succesfully.")

  return (cellchat)
}
