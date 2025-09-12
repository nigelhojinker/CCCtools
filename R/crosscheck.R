#' Comparing CellChat and CellPhoneDB results
#'
#' @param cellchat CellChat object obtained from run_cellchat() on a Seurat object.
#' @param cellphonedb Output list from run_cellphonedb() on a Seurat object.
#' @param threshold p-value (default: 0.05)
#'
#' @returns A list of six elements: (1) CellChat result (2) CellPhoneDB result (3) Full combined result including non-significant interactions (4) Filtered combined result with only significant interactions in at least one tool (5) 3x3 contingency table of interaction types in both CellChat and CellPhoneDB (6) CellChat object filtered to contain only significant interactions in both CellChat and CellPhoneDB.
#' @export
#'
#' @examples
#' \dontrun{
#' ## After setting up conda env
#' cellchat <- run_cellchat(seurat.obj, ...)
#' cpdb <- run_cellphonedb(seurat.obj, ...)
#'
#' combine <- crosscheck(cellchat = cellchat, cellphonedb = cpdb)
#' }


crosscheck <- function(cellchat, cellphonedb, threshold = 0.05) {
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

  levels <- c("Sig", "Not Sig", "Not Found")

  combine.all <- cellchat_res %>%
    full_join(cellphonedb_res %>% select(id_cp_interaction, interacting_pair, classification, source, target, ligand, from, receptor, to, unique_int, pval),
              by = "unique_int",
              suffix = c(".cc", ".cpdb")) %>%
    mutate(classification = ifelse(is.na(classification), from, classification),
           pathway_name   = ifelse(is.na(pathway_name), classification, pathway_name),
           confidence = case_when(pval.cc < threshold  & pval.cpdb < threshold  ~ "High",
                                  pval.cc < threshold  & is.na(pval.cpdb)  ~ "Mid CellChat",
                                  is.na(pval.cc)  & pval.cpdb < threshold  ~ "Mid CellPhoneDB",
                                  pval.cc < threshold  & pval.cpdb >= threshold ~ "Low CellChat",
                                  pval.cc >= threshold & pval.cpdb < threshold  ~ "Low CellPhoneDB"),
           CellPhoneDB_status = case_when(pval.cpdb < threshold ~ "Sig",
                                          pval.cpdb >= threshold ~ "Not Sig",
                                          is.na(pval.cpdb) ~ "Not Found"),
           CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
           CellChat_status = case_when(pval.cc < threshold ~ "Sig",
                                       pval.cc >= threshold ~ "Not Sig",
                                       is.na(pval.cc) ~ "Not Found"),
           CellChat_status = factor(CellChat_status, levels = levels))

  summary <- table(CellChat = combine.all$CellChat_status, CellPhoneDB = combine.all$CellPhoneDB_status) %>% addmargins()

  combine.sig <- combine.all %>%
    filter(!is.na(confidence))

  cat("Removed interactions insignificant in both tools, and interactions only found in one tool and has a p-value >=", threshold, "\n",
      "CCCtools summary: \n")

  print(summary)

  rel.int <- combine.sig %>%
    filter(confidence == "High") %>%
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

  cat("Created CellChat object with only High confidence interactions.")

  return(list(cellchat = cellchat_res, cellphonedb = cellphonedb_res, combine.all = combine.all, combine.sig = combine.sig, summary = summary, cellchat.new  = cellchat.new))

}



