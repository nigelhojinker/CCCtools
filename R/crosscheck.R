#' Comparing CellChat and CellPhoneDB results
#'
#' @param cellchat_res Dataframe obtained from running pull_netslot() on a CellChat object.
#' @param cellphonedb_res Dataframe of the p-value results from running run_cellphonedb on a Seurat object.
#' @param threshold p-value (default: 0.05)
#'
#' @returns A list of four elements: (1) Combined result (2) CellChat result (3) CellPhoneDB result (4) Summary table of containing number of interactions at each confidence level.
#' @export
#'
#' @examples
#' \dontrun{
#' ## After setting up conda env
#' cpdb <- run_cellphonedb(seurat.obj, ...)
#' cpdb.res <- cpdb$pvalues
#'
#' cellchat <- run_cellchat(seurat.obj, ...)
#' cellchat.res <- pull_netslot(cellchat)
#'
#' combine <- crosscheck(cellchat_res = cellchat.res, cellphonedb_res = cpdb.res)
#' }


crosscheck <- function(cellchat_res, cellphonedb_res, threshold = 0.05) {
  if (! is.data.frame(cellchat_res)) {
    stop("Input CellChat result is not a dataframe. Run pull_netslot() on your CellChat object and set cellchat_res as the output of pull_netslot().")
  } else if (! is.data.frame(cellphonedb_res)) {
    stop("Input CellPhoneDB result is not a dataframe. Run run_cellphonedb() on your Seurat object and set cellphonedb_res as the pvalues dataframe from the output list of run_cellphonedb().")
  }

  data(CCDB, CPDB)

  cellchat_res <- cellchat_res %>%
    left_join(CCDB %>%
                select(interaction_name, unique_int)) %>%
    mutate(unique_int = paste(source, target, unique_int, sep = "_"))

  cellphonedb_res <- cellphonedb_res %>%
    pivot_longer(cols = contains("|"),
                 names_to = "cell_type_pair",
                 values_to = "pval") %>%
    separate(cell_type_pair, into = c("cell_a", "cell_b"), sep = "\\|") %>%
    ## gene-to-uniprot conversion
    mutate(partner_a = gsub(".*:", "", partner_a),
           partner_b = gsub(".*:", "", partner_b)) %>%
    left_join(CPDB %>%
                select(id_cp_interaction, unique_int, ligand, from, receptor, to)) %>%
    ## complex-to-uniprot conversion
    mutate(partner_a = case_when(
      partner_a == from ~ ligand,
      partner_a == to ~ receptor,
      .default = partner_a),
      partner_b = case_when(
        partner_b == from ~ ligand,
        partner_b == to ~ receptor,
        .default = partner_b),
      cpdb_int = paste(partner_a, partner_b, sep = "_")) %>%
    ## creating full unique interaction (with cell type directionality)
    mutate(source = ifelse(unique_int == cpdb_int, cell_a, cell_b),
           target = ifelse(unique_int == cpdb_int, cell_b, cell_a),
           unique_int = paste(source, target, ligand, receptor, sep = "_")) %>%
    relocate(id_cp_interaction, interacting_pair, source, target, ligand, from, receptor, to, unique_int, pval, everything())

  combine <- cellchat_res %>%
    full_join(cellphonedb_res %>% select(id_cp_interaction, interacting_pair, source, target, ligand, from, receptor, to, unique_int, pval),
              by = "unique_int",
              suffix = c(".cc", ".cpdb")) %>%
    filter(! (pval.cc >= threshold & pval.cpdb >= threshold | is.na(pval.cc) & pval.cpdb >= threshold | pval.cc >= threshold & is.na(pval.cpdb)) ) %>%
    mutate(confidence = case_when(
      pval.cc < threshold  & pval.cpdb < threshold  ~ "High",
      pval.cc < threshold  & is.na(pval.cpdb)  ~ "Mid CellChat",
      is.na(pval.cc)  & pval.cpdb < threshold  ~ "Mid CellPhoneDB",
      pval.cc < threshold  & pval.cpdb >= threshold ~ "Low CellChat",
      pval.cc >= threshold & pval.cpdb < threshold  ~ "Low CellPhoneDB"))

  cat("Removed interactions insignificant in both tools, and interactions only found in one tool and has a p-value >=", threshold, "\n",
      "CCCtools result: \n")

  summary <- combine %>%
    mutate(confidence = factor(confidence, levels = c("High", "Mid CellChat", "Mid CellPhoneDB", "Low CellChat", "Low CellPhoneDB"))) %>%
    tabyl(confidence)

  print(summary)

  return(list(combine = combine, cellchat = cellchat_res, cellphonedb = cellphonedb_res, summary = summary))

}




