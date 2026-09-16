# ======================================================
# Complete version - all datasets + all groups
# Hard-voting network evaluation | No threshold filtering | No random baseline | AUROC filtering removed
# Datasets: Aridity, Alkalinity, Cold
# All results saved to Hard_Majority_Voting
# Fixes: empty-data check, calc_overlap typo, column-assignment error
# Change: chip file changed from chip_data.txt to chip_seq.csv
# ======================================================
library(vroom)
library(dplyr)
library(pROC)
library(PRROC)
library(tidyverse)

# ======================
# Global parameters
# ======================
datasets   <- c("Aridity", "Alkalinity", "Cold")
root_dir   <- "Hard_Majority_Voting"
chip_file  <- "chip_seq.csv"

# Group mapping: folder name -> group name
group_map <- list(
  c("1_Pair_2Methods",    "Pair_2Methods"),
  c("2_Triple_3Methods",  "Triple_3Methods"),
  c("3_Four_4Methods",    "Four_4Methods"),
  c("4_Five_5Methods",    "Five_5Methods"),
  c("5_All_6Methods",     "All_6Methods")
)

# ======================
# Utility function: compute edge-overlap count (typo fixed)
# ======================
calc_overlap <- function(pred_df, chip_df) {
  pred_pairs <- paste(pred_df$TF, pred_df$Target)
  true_pairs <- paste(chip_df$TF, chip_df$Target)
  sum(pred_pairs %in% true_pairs)
}

# ======================
# Single-network evaluation: AUROC / AUPR
# ======================
evaluate_network <- function(pred_df, chip_df) {
  true_overlap <- calc_overlap(pred_df, chip_df)
  pred_pairs <- paste(pred_df$TF, pred_df$Target)
  true_pairs <- paste(chip_df$TF, chip_df$Target)
  pred_df$label <- as.integer(pred_pairs %in% true_pairs)
  
  if (length(unique(pred_df$label)) < 2) {
    return(list(AUROC = NA, AUPR = NA, True_Overlap = true_overlap))
  }
  
  roc_obj <- roc(pred_df$label, pred_df$EdgeWeight, quiet = TRUE)
  auroc <- auc(roc_obj)
  pr_obj <- pr.curve(scores.class0 = pred_df$EdgeWeight[pred_df$label == 1],
                     scores.class1 = pred_df$EdgeWeight[pred_df$label == 0])
  aupr <- pr_obj$auc.integral
  
  return(list(AUROC = round(auroc,3), AUPR = round(aupr,3), True_Overlap = true_overlap))
}

# ======================
# Classification metrics: Accuracy / F1 / Precision / Recall / TP/FP/FN/TN
# ======================
get_metrics <- function(pred_df, valid_tfs, valid_targets, true_pairs, N_total_chip) {
  pred_filtered <- pred_df %>% filter(TF %in% valid_tfs, Target %in% valid_targets)
  pred_pairs <- unique(paste(pred_filtered$TF, pred_filtered$Target))
  
  TP <- sum(pred_pairs %in% true_pairs)
  FP <- length(pred_pairs) - TP
  FN <- length(true_pairs) - TP
  TN <- max(N_total_chip - TP - FP - FN, 0)
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  F1        <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  ACC       <- (TP + TN) / N_total_chip
  
  return(list(
    Accuracy = round(ACC,4),
    F1 = round(F1,4),
    Precision = round(precision,4),
    Recall = round(recall,4),
    TP = TP, FP = FP, FN = FN, TN = TN,
    N_pred_overlap = length(pred_pairs)
  ))
}

# ======================
# 1. Single dataset: evaluate AUROC/AUPR group by group (AUROC threshold filtering removed)
# ======================
run_eval_single_ds <- function(ds) {
  fusion_root <- file.path(root_dir, paste0("Voting_Sig_Results_", ds, "_0-1"))
  out_dir     <- file.path(root_dir, paste0("Voting_Evaluation_", ds))
  
  if (!file.exists(chip_file)) stop("Missing chip_seq.csv")
  chip_seq <- read_csv(chip_file, col_names = c("TF","Target"), show_col_types = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ds_all <- list()
  
  for(g in group_map){
    folder_name <- g[[1]]
    group_name  <- g[[2]]
    
    cat("\n=========================================\n")
    cat("Evaluating ", ds, ": ", group_name, "\n")
    cat("=========================================\n")
    
    input_dir <- file.path(fusion_root, folder_name)
    out_group <- file.path(out_dir, group_name)
    out_eval  <- file.path(out_group, "Evaluation")
    dir.create(out_eval, recursive = TRUE)
    
    files <- list.files(input_dir, pattern = "Voting_Sig_.*\\.csv$", full.names = TRUE)
    group_res <- list()
    
    for(f in files){
      net_name <- gsub("Voting_Sig_|\\.csv", "", basename(f))
      net <- vroom(f, show_col_types = FALSE)
      total_edges <- nrow(net)
      
      res <- evaluate_network(net, chip_seq)
      row_df <- data.frame(
        Network = net_name,
        AUROC = res$AUROC,
        AUPR = res$AUPR,
        True_Overlap = res$True_Overlap,
        Total_Edges = total_edges
      )
      group_res[[length(group_res)+1]] <- row_df
    }
    
    group_df <- bind_rows(group_res)
    # AUROC threshold filtering removed; all results retained
    write.csv(group_df, file.path(out_eval, paste0("Eval_", group_name, ".csv")), row.names = FALSE)
    
    if(nrow(group_df) > 0){
      ds_all[[length(ds_all)+1]] <- group_df %>% mutate(Group = group_name, Folder = folder_name)
    }
  }
  
  ds_total <- bind_rows(ds_all)
  return(list(ds_total = ds_total, chip = chip_seq))
}

# ======================
# 2. Single dataset: compute the full set of classification metrics
# ======================
run_metrics_single_ds <- function(ds, eval_data, chip_seq) {
  if(nrow(eval_data) == 0){
    cat("\n[Warning] ", ds, " has no evaluation data; skipping metric computation\n")
    return()
  }
  
  OUTPUT_FILE  <- file.path(root_dir, paste0(ds, "_Voting_FINAL.csv"))
  SUMMARY_FILE <- file.path(root_dir, paste0(ds, "_Voting_Summary.csv"))
  
  true_pairs    <- unique(paste(chip_seq$TF, chip_seq$Target))
  valid_tfs     <- unique(chip_seq$TF)
  valid_targets <- unique(chip_seq$Target)
  N_total_chip  <- length(valid_tfs) * length(valid_targets)
  
  df <- eval_data %>% arrange(desc(AUROC))
  
  # Pre-initialize columns to avoid assignment errors
  df$Accuracy      <- NA
  df$F1            <- NA
  df$Precision     <- NA
  df$Recall        <- NA
  df$TP            <- NA
  df$FP            <- NA
  df$FN            <- NA
  df$TN            <- NA
  df$N_pred_overlap<- NA
  
  for(i in 1:nrow(df)){
    folder  <- df$Folder[i]
    netname <- df$Network[i]
    fpath <- file.path(
      root_dir,
      paste0("Voting_Sig_Results_", ds, "_0-1"),
      folder,
      paste0("Voting_Sig_", netname, ".csv")
    )
    if(!file.exists(fpath)) next
    
    net <- vroom(fpath, show_col_types = FALSE)
    m_res <- get_metrics(net, valid_tfs, valid_targets, true_pairs, N_total_chip)
    
    df$Accuracy[i]       <- m_res$Accuracy
    df$F1[i]             <- m_res$F1
    df$Precision[i]      <- m_res$Precision
    df$Recall[i]         <- m_res$Recall
    df$TP[i]             <- m_res$TP
    df$FP[i]             <- m_res$FP
    df$FN[i]             <- m_res$FN
    df$TN[i]             <- m_res$TN
    df$N_pred_overlap[i] <- m_res$N_pred_overlap
  }
  
  write.csv(df, OUTPUT_FILE, row.names = FALSE)
  # Compact summary table
  sum_df <- df %>%
    select(Group, Network, AUROC, AUPR, F1, Precision, Recall, Accuracy) %>%
    arrange(desc(AUROC))
  write.csv(sum_df, SUMMARY_FILE, row.names = FALSE)
  
  cat("\n[Done] ", ds, " classification metrics computed\n")
}

# ======================
# Batch execute all datasets
# ======================
for(ds in datasets){
  cat("\n\n==================================================\n")
  cat("Starting to process dataset: ", ds, "\n")
  cat("==================================================\n")
  
  eval_out <- run_eval_single_ds(ds)
  run_metrics_single_ds(ds, eval_out$ds_total, eval_out$chip)
}

cat("\n[Finished] All datasets and all groups evaluated! Results are under the Hard_Majority_Voting directory\n")
