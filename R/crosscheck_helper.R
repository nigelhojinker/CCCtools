#' Internal helpers for comparative analysis in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CC_CPDB <- function(cellchat, cellphonedb, threshold = 0.05, return.all = FALSE) {

  cellchat_res <- pull_netslot(cellchat) %>%
    left_join(CCDB %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  cellphonedb_res <- cellphonedb$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(CPDB %>% select(1:7)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  levels <- c("Sig", "Not_Sig", "Not_Found")

  combine.all <- cellchat_res %>%
    full_join(cellphonedb_res %>% select(id_cp_interaction, interacting_pair, classification, source, target, ligand_name, ligand, receptor_name, receptor, unique_int, pval),
              by = "unique_int",
              suffix = c("_cc", "_cpdb")) %>%
    mutate(classification = ifelse(is.na(classification), ligand_name, classification),
           pathway_name   = ifelse(is.na(pathway_name), classification, pathway_name),
           category = case_when(pval_cc < threshold  & pval_cpdb < threshold  ~ "Sig_Both",
                                pval_cc < threshold  & is.na(pval_cpdb)  ~ "Sig_CellChat_Not_Found_CellPhoneDB",
                                is.na(pval_cc)  & pval_cpdb < threshold  ~ "Sig_CellPhoneDB_Not_Found_CellChat",
                                pval_cc < threshold  & pval_cpdb >= threshold ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                pval_cc >= threshold & pval_cpdb < threshold  ~ "Sig_CellPhoneDB_Not_Sig_CellChat"),
           CellPhoneDB_status = case_when(pval_cpdb < threshold ~ "Sig",
                                          pval_cpdb >= threshold ~ "Not_Sig",
                                          is.na(pval_cpdb) ~ "Not_Found"),
           CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
           CellChat_status = case_when(pval_cc < threshold ~ "Sig",
                                       pval_cc >= threshold ~ "Not_Sig",
                                       is.na(pval_cc) ~ "Not_Found"),
           CellChat_status = factor(CellChat_status, levels = levels))

  summary <- table(CellChat = combine.all$CellChat_status, CellPhoneDB = combine.all$CellPhoneDB_status) %>% addmargins()

  if (!return.all) {
    res <- combine.all %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", threshold, "\n",
      "CCCtools summary: \n")

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
    left_join(CCDB %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob1 := prob, !!pval1 := pval)

  ## Process second CellChat object
  cellchat_res2 <- pull_netslot(cellchat2) %>%
    left_join(CCDB %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob2 := prob, !!pval2 := pval)

  ## Combining results
  found_both <- left_join(cellchat_res1, cellchat_res2)

  second_only <- anti_join(cellchat_res2, cellchat_res1)

  res <- bind_rows(found_both, second_only) %>%
    relocate(interaction_name, source, target, ligand, receptor, starts_with("prob"), starts_with("pval"), everything()) %>%
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
      "CCCtools summary: \n")

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
    left_join(CPDB %>% select(1:7)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval1 := pval)

  ## CellPhoneDB object 2 processing
  cellphonedb_res2 <- cellphonedb2$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(CPDB %>% select(1:7)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval2 := pval)

  ## Combining results
  found_both <- left_join(cellphonedb_res1, cellphonedb_res2)

  second_only <- anti_join(cellphonedb_res2, cellphonedb_res1)

  res <- bind_rows(found_both, second_only) %>%
    relocate(id_cp_interaction, interacting_pair, source, target, ligand_name, ligand, receptor_name, receptor, unique_int, starts_with("pval")) %>%
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
      "CCCtools summary: \n")

  print(summary)

  return(list(result = res, summary = summary))
}
