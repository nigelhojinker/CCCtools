#' Internal helpers for comparative analysis in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CC_CPDB <- function(cellchat, cellphonedb, threshold = 0.05, return.all = FALSE) {

  cellchat_res <- pull_netslot(cellchat) %>%
    left_join(db_map %>% select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  cellphonedb_res <- cellphonedb$pvalues %>% pivot_longer(cols = contains("|"),
                                                          names_to = "cell_type_pair",
                                                          values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(db_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  levels <- c("Sig", "Not_Sig", "Not_Found")

  res <- cellchat_res %>% full_join(cellphonedb_res %>%
                                              select(id_cp_interaction, interacting_pair, source, target,  unique_int, pval),
                                            by = "unique_int",
                                            suffix = c("_cc", "_cpdb")) %>%
    mutate(category = case_when(pval_cc < threshold & pval_cpdb < threshold ~ "Sig_Both",
                                pval_cc < threshold & is.na(pval_cpdb) ~ "Sig_CellChat_Not_Found_CellPhoneDB",
                                is.na(pval_cc) & pval_cpdb < threshold ~ "Sig_CellPhoneDB_Not_Found_CellChat",
                                pval_cc < threshold & pval_cpdb >= threshold ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                pval_cc >= threshold & pval_cpdb < threshold ~ "Sig_CellPhoneDB_Not_Sig_CellChat"),
           CellPhoneDB_status = case_when(pval_cpdb < threshold ~ "Sig",
                                          pval_cpdb >= threshold ~ "Not_Sig", is.na(pval_cpdb) ~ "Not_Found"),
           CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
           CellChat_status = case_when(pval_cc < threshold ~ "Sig",
                                       pval_cc >= threshold ~ "Not_Sig",
                                       is.na(pval_cc) ~ "Not_Found"),
           CellChat_status = factor(CellChat_status, levels = levels))

  summary <- table(CellChat = combine.all$CellChat_status,
                   CellPhoneDB = combine.all$CellPhoneDB_status) %>% addmargins()
  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }
  cat("Comparative analysis done successfully for p-value",
      threshold, "\n", "scpeakeR summary: \n")
  print(summary)
  return(list(result = res, summary = summary))

}

#' Internal helpers for comparative analysis between CellChat objects in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CC_CC <- function(cellchat1, cellchat2, name1, name2, threshold = 0.05, return.all = FALSE) {

  prob1 <- paste0("prob_", name1)
  pval1 <- paste0("pval_", name1)
  prob2 <- paste0("prob_", name2)
  pval2 <- paste0("pval_", name2)
  status1 <- paste0("status_", name1)
  status2 <- paste0("status_", name2)

  levels <- c("Sig", "Not_Sig", "Not_Found")

  ## Process first CellChat object
  cellchat_res1 <- pull_netslot(cellchat1) %>%
    left_join(db_map %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob1 := prob, !!pval1 := pval)

  ## Process second CellChat object
  cellchat_res2 <- pull_netslot(cellchat2) %>%
    left_join(db_map %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob2 := prob, !!pval2 := pval)

  ## Combining results
  ## Combining results
  res <- full_join(cellchat_res1, cellchat_res2, by = "unique_int", suffix = c("_cc1", "_cc2")) %>%
    mutate(interaction_name = coalesce(interaction_name_cc1, interaction_name_cc2),
           pathway_name     = coalesce(pathway_name_cc1, pathway_name_cc2),
           source     = coalesce(source_cc1,   source_cc2),
           target     = coalesce(target_cc1,   target_cc2),
           ligand     = coalesce(ligand_cc1,   ligand_cc2),
           receptor   = coalesce(receptor_cc1, receptor_cc2),
           annotation = coalesce(annotation_cc1, annotation_cc2),
           evidence   = coalesce(evidence_cc1, evidence_cc2)) %>%
    select(interaction_name, pathway_name, source, target, ligand, receptor,
           unique_int, starts_with("prob"), starts_with("pval"), annotation, evidence) %>%
    mutate(category = case_when(!!sym(pval1) < threshold  & !!sym(pval2) < threshold  ~ "Sig_Both",
                                !!sym(pval1) < threshold  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                is.na(!!sym(pval1))  & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                !!sym(pval1) < threshold  & !!sym(pval2) >= threshold ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                !!sym(pval1) >= threshold & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
           !!status1 := factor(
             case_when(!!sym(pval1) < threshold ~ "Sig",
                       !!sym(pval1) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval1)) ~ "Not_Found"),
             levels = levels),
           !!status2 := factor(
             case_when(!!sym(pval2) < threshold ~ "Sig",
                       !!sym(pval2) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval2)) ~ "Not_Found"),
             levels = levels))

  summary <- table(res[[status1]], res[[status2]]) %>%
    addmargins()

  dimnames(summary) <- dimnames(summary) %>%
    setNames(c(name1, name2))

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", threshold, "\n",
      "scpeakeR summary: \n")

  print(summary)

  return(list(result = res, summary = summary))
}

#' Internal helpers for comparative analysis between CellPhoneDB objects in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CPDB_CPDB <- function(cellphonedb1, cellphonedb2, name1, name2, threshold = 0.05, return.all = FALSE) {

  pval1 <- paste0("pval_", name1)
  pval2 <- paste0("pval_", name2)
  status1 <- paste0("status_", name1)
  status2 <- paste0("status_", name2)

  levels <- c("Sig", "Not_Sig", "Not_Found")

  ## CellPhoneDB object 1 processing
  cellphonedb_res1 <- cellphonedb1$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(db_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval1 := pval)

  ## CellPhoneDB object 2 processing
  cellphonedb_res2 <- cellphonedb2$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(db_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval2 := pval)

  ## Combining results
  res <- full_join(cellphonedb_res1, cellphonedb_res2, by = "unique_int", suffix = c("_cp1", "_cp2")) %>%
    mutate(id_cp_interaction = coalesce(id_cp_interaction_cp1, id_cp_interaction_cp2),
           interacting_pair  = coalesce(interacting_pair_cp1, interacting_pair_cp2),
           source     = coalesce(source_cp1,   source_cp2),
           target     = coalesce(target_cp1,   target_cp2),
           partner_a  = coalesce(partner_a_cp1, partner_a_cp2),
           partner_b  = coalesce(partner_b_cp1, partner_b_cp2),
           gene_a = coalesce(gene_a_cp1, gene_b_cp2),
           gene_b = coalesce(gene_b_cp1, gene_b_cp2),
           secreted = coalesce(secreted_cp1, secreted_cp2),
           receptor_a = coalesce(receptor_a_cp1, receptor_b_cp2),
           receptor_b = coalesce(receptor_b_cp1, receptor_b_cp2),
           annotation_strategy = coalesce(annotation_strategy_cp1, annotation_strategy_cp2),
           is_integrin = coalesce(is_integrin_cp1, is_integrin_cp2),
           directionality = coalesce(directionality_cp1, directionality_cp2),
           classification = coalesce(classification_cp1, classification_cp2)) %>%
    select(id_cp_interaction, interacting_pair, source, target, partner_a, partner_b,
           gene_a, gene_b, starts_with("pval"), secreted, receptor_a, receptor_b,
           annotation_strategy, is_integrin, directionality, classification, unique_int) %>%
    mutate(category = case_when(!!sym(pval1) < threshold  & !!sym(pval2) < threshold  ~ "Sig_Both",
                                !!sym(pval1) < threshold  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                is.na(!!sym(pval1))  & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                !!sym(pval1) < threshold  & !!sym(pval2) >= threshold ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                !!sym(pval1) >= threshold & !!sym(pval2) < threshold  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
           !!status1 := factor(
             case_when(!!sym(pval1) < threshold ~ "Sig",
                       !!sym(pval1) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval1)) ~ "Not_Found"),
             levels = levels),
           !!status2 := factor(
             case_when(!!sym(pval2) < threshold ~ "Sig",
                       !!sym(pval2) >= threshold ~ "Not_Sig",
                       is.na(!!sym(pval2)) ~ "Not_Found"),
             levels = levels))

  summary <- table(res[[status1]], res[[status2]]) %>%
    addmargins()

  dimnames(summary) <- dimnames(summary) %>%
    setNames(c(name1, name2))

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", threshold, "\n",
      "scpeakeR summary: \n")

  print(summary)

  return(list(result = res, summary = summary))
}
