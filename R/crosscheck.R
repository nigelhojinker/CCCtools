#' Comparing CellChat and CellPhoneDB results
#'
#' @param obj1 CellChat object obtained from run_cellchat() on a Seurat object OR output list from run_cellphonedb() on a Seurat object.
#' @param obj2 CellChat object obtained from run_cellchat() on a Seurat object OR output list from run_cellphonedb() on a Seurat object.
#' @param comp.type Character value stating the type of comparative analysis to run; assumes a CellChat-to-CellPhoneDB comparison by default ("CC_CP").
#' @param threshold p-value (default: 0.05)
#' @param return.all logical value to determine whether to return the full comparative analysis results (including non-significant interactions).
#'
#' @returns A list of two elements: (1) Comparative analysis result (2) 3x3 contingency table of interaction types in both CellChat and CellPhoneDB.
#' @export
#'
#' @examples
#' \dontrun{
#' ## After setting up conda env
#' cpdb <- run_cellphonedb(seurat.obj, ...)
#' cellchat <- run_cellchat(seurat.obj, ...)
#'
#' combine <- crosscheck(obj1 = cellchat, obj2 = cpdb, comp.type = "CC_CP")
#' }

crosscheck <- function(obj1 = NULL, obj2 = NULL,
                       comp.type = c("CC_CP", "CP_CC", "CC_CC", "CP_CP"),
                       threshold = 0.05, return.all = FALSE) {

  if(is.null(obj1)) stop("obj1 must be supplied.")
  if(is.null(obj2)) stop("obj2 must be supplied.")

  comp.type <- match.arg(comp.type)

  obj1.name <- deparse(substitute(obj1))
  obj2.name <- deparse(substitute(obj2))

  if(comp.type == "CC_CP") {
    if(! inherits(obj1, "CellChat"))    stop(paste(obj1.name, "is not a CellChat object."))
    if(! inherits(obj2, "CellPhoneDB")) stop(paste(obj2.name, "is not a CellPhoneDB object."))

    return(crosscheck_CC_CPDB(obj1, obj2, threshold = threshold, return.all = return.all))
  }

  if(comp.type == "CP_CC") {
    if(! inherits(obj1, "CellPhoneDB")) stop(paste(obj1.name, "is not a CellPhoneDB object."))
    if(! inherits(obj2, "CellChat"))    stop(paste(obj2.name, "is not a CellChat object."))

    return(crosscheck_CC_CPDB(obj2, obj1, threshold = threshold, return.all = return.all))
  }

  if(comp.type == "CC_CC") {
    if(! inherits(obj1, "CellChat"))    stop(paste(obj1.name, "is not a CellChat object."))
    if(! inherits(obj2, "CellChat"))    stop(paste(obj2.name, "is not a CellChat object."))

    return(crosscheck_CC_CC(obj1, obj2, name1 = obj1.name, name2 = obj2.name, threshold = threshold, return.all = return.all))
  }

  if(comp.type == "CP_CP") {
    if(! inherits(obj1, "CellPhoneDB")) stop(paste(obj1.name, "is not a CellPhoneDB object."))
    if(! inherits(obj2, "CellPhoneDB")) stop(paste(obj2.name, "is not a CellPhoneDB object."))

    return(crosscheck_CPDB_CPDB(obj1, obj2, name1 = obj1.name, name2 = obj2.name, threshold = threshold, return.all = return.all))
  }
}
