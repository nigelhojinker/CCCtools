#' Internal helpers for comparative analysis in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CC_CPDB <- function(cellchat, cellphonedb, p.cutoff = 0.05, return.all = FALSE,
                               expand = FALSE) {

  cellchat_res <- pull_netslot(cellchat) %>%
    left_join(db_map %>% select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  cellphonedb_res <- cellphonedb$pvalues %>% pivot_longer(cols = contains("|"),
                                                          names_to = "cell_type_pair",
                                                          values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(db_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|"))

  if (expand){
    levels <- c("Sig", "Not_Sig", "Not_Found")

    res <- cellchat_res %>% full_join(cellphonedb_res %>%
                                        select(id_cp_interaction, interacting_pair, source, target,  unique_int, pval),
                                      by = "unique_int",
                                      suffix = c("_cc", "_cpdb")) %>%
      mutate(category = case_when(pval_cc < p.cutoff & pval_cpdb < p.cutoff ~ "Sig_Both",
                                  pval_cc < p.cutoff & is.na(pval_cpdb) ~ "Sig_CellChat_Not_Found_CellPhoneDB",
                                  is.na(pval_cc) & pval_cpdb < p.cutoff ~ "Sig_CellPhoneDB_Not_Found_CellChat",
                                  pval_cc < p.cutoff & pval_cpdb >= p.cutoff ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                  pval_cc >= p.cutoff & pval_cpdb < p.cutoff ~ "Sig_CellPhoneDB_Not_Sig_CellChat"),
             CellPhoneDB_status = case_when(pval_cpdb < p.cutoff ~ "Sig",
                                            pval_cpdb >= p.cutoff ~ "Not_Sig", is.na(pval_cpdb) ~ "Not_Found"),
             CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
             CellChat_status = case_when(pval_cc < p.cutoff ~ "Sig",
                                         pval_cc >= p.cutoff ~ "Not_Sig",
                                         is.na(pval_cc) ~ "Not_Found"),
             CellChat_status = factor(CellChat_status, levels = levels))

    summary <- table(CellChat = res$CellChat_status,
                     CellPhoneDB = res$CellPhoneDB_status) %>% addmargins()

    tp <- summary[1]
    fp <- summary[5] + summary[9]
    tn <- summary[6] + summary[7] + summary[10]
    fn <- summary[2] + summary[3]

    nsig_both <- tp
    jaccard   <- tp / (summary[13] + summary[4] - tp)
    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5

  } else {
    levels <- c("Sig", "Not_Sig")

    res <- cellchat_res %>% full_join(cellphonedb_res %>%
                                        select(id_cp_interaction, interacting_pair, source, target,  unique_int, pval),
                                      by = "unique_int",
                                      suffix = c("_cc", "_cpdb")) %>%
      mutate(category = case_when(pval_cc < p.cutoff & pval_cpdb < p.cutoff ~ "Sig_Both",
                                  pval_cc < p.cutoff & is.na(pval_cpdb) ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                  is.na(pval_cc) & pval_cpdb < p.cutoff ~ "Sig_CellPhoneDB_Not_Sig_CellChat",
                                  pval_cc < p.cutoff & pval_cpdb >= p.cutoff ~ "Sig_CellChat_Not_Sig_CellPhoneDB",
                                  pval_cc >= p.cutoff & pval_cpdb < p.cutoff ~ "Sig_CellPhoneDB_Not_Sig_CellChat"),
             CellPhoneDB_status = case_when(pval_cpdb < p.cutoff ~ "Sig",
                                            (pval_cpdb >= p.cutoff | is.na(pval_cpdb)) ~ "Not_Sig"),
             CellPhoneDB_status = factor(CellPhoneDB_status, levels = levels),
             CellChat_status = case_when(pval_cc < p.cutoff ~ "Sig",
                                         (pval_cc >= p.cutoff | is.na(pval_cc)) ~ "Not_Sig"),
             CellChat_status = factor(CellChat_status, levels = levels))

    summary <- table(CellChat = res$CellChat_status,
                     CellPhoneDB = res$CellPhoneDB_status) %>% addmargins()

    tp <- summary[1]
    fp <- summary[4]
    tn <- summary[5]
    fn <- summary[2]

    nsig_both <- tp
    jaccard   <- tp / (summary[7] + summary[3] - tp)
    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5
  }

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value",
      p.cutoff, "\n", "scpeaker summary: \n")
  print(summary)
  cat("Number of significant interactions in both CellChat and CellPhoneDB:", nsig_both, "\n")
  cat("Jaccard Index:", jaccard, "\n")
  cat("Matthew's Correlation Coefficient:", mcc, "\n")

  return(list(result = res, summary = summary,
              nsig_both = nsig_both, jaccard = jaccard, mcc = mcc))

}

#' Internal helpers for comparative analysis between CellChat objects in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CC_CC <- function(cellchat1, cellchat2, name1, name2, p.cutoff = 0.05, return.all = FALSE, expand = FALSE) {

  prob1 <- paste0("prob_", name1)
  pval1 <- paste0("pval_", name1)
  prob2 <- paste0("prob_", name2)
  pval2 <- paste0("pval_", name2)
  status1 <- paste0("status_", name1)
  status2 <- paste0("status_", name2)

  # Process first CellChat object
  cellchat_res1 <- pull_netslot(cellchat1) %>%
    left_join(db_map %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob1 := prob, !!pval1 := pval)

  # Process second CellChat object
  cellchat_res2 <- pull_netslot(cellchat2) %>%
    left_join(db_map %>%
                select(interaction_name, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!prob2 := prob, !!pval2 := pval)

  if (expand){
    levels <- c("Sig", "Not_Sig", "Not_Found")

    # Combining results
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
      mutate(category = case_when(!!sym(pval1) < p.cutoff  & !!sym(pval2) < p.cutoff  ~ "Sig_Both",
                                  !!sym(pval1) < p.cutoff  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                  is.na(!!sym(pval1))  & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                  !!sym(pval1) < p.cutoff  & !!sym(pval2) >= p.cutoff ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                  !!sym(pval1) >= p.cutoff & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
             !!status1 := factor(
               case_when(!!sym(pval1) < p.cutoff ~ "Sig",
                         !!sym(pval1) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval1)) ~ "Not_Found"),
               levels = levels),
             !!status2 := factor(
               case_when(!!sym(pval2) < p.cutoff ~ "Sig",
                         !!sym(pval2) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval2)) ~ "Not_Found"),
               levels = levels))

    summary <- table(res[[status1]], res[[status2]]) %>%
      addmargins()

    dimnames(summary) <- dimnames(summary) %>%
      setNames(c(name1, name2))

    tp <- summary[1]
    fp <- summary[5] + summary[9]
    tn <- summary[6] + summary[7] + summary[10]
    fn <- summary[2] + summary[3]

    nsig_both <- tp
    jaccard   <- tp / (summary[13] + summary[4] - tp)
    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5

  } else {
    levels <- c("Sig", "Not_Sig")

    # Combining results
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
      mutate(category = case_when(!!sym(pval1) < p.cutoff  & !!sym(pval2) < p.cutoff  ~ "Sig_Both",
                                  !!sym(pval1) < p.cutoff  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                  is.na(!!sym(pval1))  & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Sig_", name1),
                                  !!sym(pval1) < p.cutoff  & !!sym(pval2) >= p.cutoff ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                  !!sym(pval1) >= p.cutoff & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
             !!status1 := factor(
               case_when(!!sym(pval1) < p.cutoff ~ "Sig",
                         !!sym(pval1) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval1)) ~ "Not_Sig"),
               levels = levels),
             !!status2 := factor(
               case_when(!!sym(pval2) < p.cutoff ~ "Sig",
                         !!sym(pval2) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval2)) ~ "Not_Sig"),
               levels = levels))

    summary <- table(res[[status1]], res[[status2]]) %>%
      addmargins()

    dimnames(summary) <- dimnames(summary) %>%
      setNames(c(name1, name2))

    tp <- summary[1]
    fp <- summary[4]
    tn <- summary[5]
    fn <- summary[2]

    nsig_both <- tp
    jaccard   <- tp / (summary[7] + summary[3] - tp)
    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5

  }

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", p.cutoff, "\n",
      "scpeaker summary: \n")
  print(summary)
  cat("Number of significant interactions in both CellChat analyses:", nsig_both, "\n")
  cat("Jaccard Index:", jaccard, "\n")
  cat("Matthew's Correlation Coefficient:", mcc, "\n")

  return(list(result = res, summary = summary,
              nsig_both = nsig_both, jaccard = jaccard, mcc = mcc))
}

#' Internal helpers for comparative analysis between CellPhoneDB objects in crosscheck()
#'
#' @keywords internal
#' @noRd

crosscheck_CPDB_CPDB <- function(cellphonedb1, cellphonedb2, name1, name2, p.cutoff = 0.05, return.all = FALSE, expand = FALSE) {

  pval1 <- paste0("pval_", name1)
  pval2 <- paste0("pval_", name2)
  status1 <- paste0("status_", name1)
  status2 <- paste0("status_", name2)

  # Process first CellPhoneDB object
  cellphonedb_res1 <- cellphonedb1$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(db_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval1 := pval)

  # Process second CellPhoneDB object
  cellphonedb_res2 <- cellphonedb2$pvalues %>%
    pivot_longer(cols = contains("|"), names_to = "cell_type_pair", values_to = "pval") %>%
    separate(cell_type_pair, into = c("source", "target"), sep = "\\|") %>%
    left_join(db_map %>% select(id_cp_interaction, LR)) %>%
    mutate(unique_int = paste(source, target, LR, sep = "|")) %>%
    rename(!!pval2 := pval)

  if (expand){
    levels <- c("Sig", "Not_Sig", "Not_Found")

  # Combining results
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
      mutate(category = case_when(!!sym(pval1) < p.cutoff  & !!sym(pval2) < p.cutoff  ~ "Sig_Both",
                                  !!sym(pval1) < p.cutoff  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                  is.na(!!sym(pval1))  & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                  !!sym(pval1) < p.cutoff  & !!sym(pval2) >= p.cutoff ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                  !!sym(pval1) >= p.cutoff & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
             !!status1 := factor(
               case_when(!!sym(pval1) < p.cutoff ~ "Sig",
                         !!sym(pval1) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval1)) ~ "Not_Found"),
               levels = levels),
             !!status2 := factor(
               case_when(!!sym(pval2) < p.cutoff ~ "Sig",
                         !!sym(pval2) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval2)) ~ "Not_Found"),
               levels = levels))

    summary <- table(res[[status1]], res[[status2]]) %>%
      addmargins()

    dimnames(summary) <- dimnames(summary) %>%
      setNames(c(name1, name2))

    tp <- summary[1]
    fp <- summary[5] + summary[9]
    tn <- summary[6] + summary[7] + summary[10]
    fn <- summary[2] + summary[3]

    nsig_both <- tp
    jaccard   <- tp / (summary[13] + summary[4] - tp)
    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5

  } else {
    levels <- c("Sig", "Not_Sig")

    # Combining results
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
      mutate(category = case_when(!!sym(pval1) < p.cutoff  & !!sym(pval2) < p.cutoff  ~ "Sig_Both",
                                  !!sym(pval1) < p.cutoff  & is.na(!!sym(pval2))  ~ paste0("Sig_", name1, "_Not_Found_", name2),
                                  is.na(!!sym(pval1))  & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Found_", name1),
                                  !!sym(pval1) < p.cutoff  & !!sym(pval2) >= p.cutoff ~ paste0("Sig_", name1, "_Not_Sig_", name2),
                                  !!sym(pval1) >= p.cutoff & !!sym(pval2) < p.cutoff  ~ paste0("Sig_", name2, "_Not_Sig_", name1)),
             !!status1 := factor(
               case_when(!!sym(pval1) < p.cutoff ~ "Sig",
                         !!sym(pval1) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval1)) ~ "Not_Sig"),
               levels = levels),
             !!status2 := factor(
               case_when(!!sym(pval2) < p.cutoff ~ "Sig",
                         !!sym(pval2) >= p.cutoff ~ "Not_Sig",
                         is.na(!!sym(pval2)) ~ "Not_Sig"),
               levels = levels))

    summary <- table(res[[status1]], res[[status2]]) %>%
      addmargins()

    dimnames(summary) <- dimnames(summary) %>%
      setNames(c(name1, name2))

    tp <- summary[1]
    fp <- summary[4]
    tn <- summary[5]
    fn <- summary[2]

    nsig_both <- tp
    jaccard   <- tp / (summary[7] + summary[3] - tp)
    mcc <- (tp * tn - fp * fn) / ( (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) )**0.5

  }

  if (!return.all) {
    res <- res %>%
      filter(!is.na(category))
  }

  cat("Comparative analysis done successfully for p-value", p.cutoff, "\n",
      "scpeaker summary: \n")
  print(summary)
  cat("Number of significant interactions in both CellPhoneDB analyses:", nsig_both, "\n")
  cat("Jaccard Index:", jaccard, "\n")
  cat("Matthew's Correlation Coefficient:", mcc, "\n")

  return(list(result = res, summary = summary,
              nsig_both = nsig_both, jaccard = jaccard, mcc = mcc))
}
