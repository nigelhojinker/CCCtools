#' Comparing CellChat and CellPhoneDB results
#'
#' @param cellchat CellChat object obtained from run_cellchat() on a Seurat object.
#' @param cellphonedb Output list from run_cellphonedb() on a Seurat object.
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
#' combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)
#' }


crosscheck <- function(cellchat, cellphonedb, threshold = 0.05, return.all = FALSE) {
  if (! is.data.frame(pull_netslot(cellchat))) {
    stop("Input CellChat result is not a dataframe. Perform run_cellchat() on your Seurat object and set cellchat as your output CellChat object.")
  } else if (! is.data.frame(cellphonedb$pvalues)) {
    stop("Input CellPhoneDB result is not a dataframe. Run run_cellphonedb() on your Seurat object and set cellphonedb as the output list of run_cellphonedb().")
  }

  data(CCDB, CPDB)

  cellchat_res <- pull_netslot(cellchat) %>%
    left_join(CCDB %>%
                select(interaction_name, unique_int)) %>%
    mutate(unique_int = paste(source, target, unique_int, sep = "_"))

  cellphonedb_res <- cellphonedb$pvalues %>%
    pivot_longer(cols = contains("|"),
                 names_to = "cell_type_pair",
                 values_to = "pval") %>%
    separate(cell_type_pair, into = c("cell_a", "cell_b"), sep = "\\|") %>%
    ## gene-to-uniprot conversion
    mutate(partner_a = gsub(".*:", "", partner_a),
           partner_b = gsub(".*:", "", partner_b)) %>%
    left_join(CPDB %>%
                select(id_cp_interaction, classification, unique_int, ligand, from, receptor, to)) %>%
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

  levels <- c("Sig", "Not_Sig", "Not_Found")

  combine.all <- cellchat_res %>%
    full_join(cellphonedb_res %>% select(id_cp_interaction, interacting_pair, classification, source, target, ligand, from, receptor, to, unique_int, pval),
              by = "unique_int",
              suffix = c(".cc", ".cpdb")) %>%
    mutate(classification = ifelse(is.na(classification), from, classification),
           pathway_name   = ifelse(is.na(pathway_name), classification, pathway_name),
           confidence = case_when(pval.cc < threshold  & pval.cpdb < threshold  ~ "Sig_Both",
                                  pval.cc < threshold  & is.na(pval.cpdb)  ~ "Sig_CellChat_Not_Found_CellPhoneDB",
                                  is.na(pval.cc)  & pval.cpdb < threshold  ~ "Sig_CellPhoneDB_Not_Found_CellChat",
                                  pval.cc < threshold  & pval.cpdb >= threshold ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                  pval.cc >= threshold & pval.cpdb < threshold  ~ "Sig_CellPhoneDB_Not_Sig_CellChat"),
           CellPhoneDB_status = case_when(pval.cpdb < threshold ~ "Sig",
                                          pval.cpdb >= threshold ~ "Not_Sig",
                                          is.na(pval.cpdb) ~ "Not_Found"),
           CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
           CellChat_status = case_when(pval.cc < threshold ~ "Sig",
                                       pval.cc >= threshold ~ "Not_Sig",
                                       is.na(pval.cc) ~ "Not_Found"),
           CellChat_status = factor(CellChat_status, levels = levels))

  summary <- table(CellChat = combine.all$CellChat_status, CellPhoneDB = combine.all$CellPhoneDB_status) %>% addmargins()

  if (!return.all) {
    res <- combine.all %>%
      filter(!is.na(confidence))
  } else {
    res <- combine.all
  }

  cat("Comparative analysis done successfully for p-value", threshold, "\n",
      "CCCtools summary: \n")

  print(summary)
  
  return(list(result = res, summary = summary))

}



