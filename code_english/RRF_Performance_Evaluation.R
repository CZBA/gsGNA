# ======================================================
# Assessment Pipeline for Multi‑threshold Fusion Networks
# No random baseline, retain top‑10% edges by weight
# Function 1: Generate high‑confidence edges (top‑10% weight edges) and compute AUROC / AUPR
# Function 2: Calculate F1 / ACC / Precision / Recall
# Datasets: Aridity, Alkalinity, Cold
# Gold standard: chip_seq.csv
# Modification: Removed AUROC > 0.7 filtering to avoid empty‑data errors
# ======================================================
library(vroom)
library(dplyr)
library(pROC)
library(PRROC)
library(tidyverse)

datasets <- c("Aridity", "Alkalinity", "Cold")

run_evaluation <- function(dataset) {
  fusion_root   <- paste0("Fusion_Results_", dataset, "_0-1")
  out_dir       <- paste0("Fusion_Evaluation_", dataset, "_Top10")
  chip_file     <- "chip_seq.csv"
  top_ratio     <- 0.1
  
  if (!file.exists(chip_file)) stop(paste("Gold standard file not found:", chip_file))
  chip_seq <- read_csv(chip_file, col_names = TRUE) %>% select(TF = 1, Target = 2)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  process_one_group <- function(folder_name, group_name) {
    cat("\n=========================================\n")
    cat("Evaluating", dataset, ":", group_name, "\n")
    cat("=========================================\n")
    
    input_dir  <- file.path(fusion_root, folder_name)
    out_group  <- file.path(out_dir, group_name)
    out_filt   <- file.path(out_group, "HighConf_Edges")
    out_eval   <- file.path(out_group, "Evaluation")
    lapply(c(out_filt, out_eval), dir.create, recursive = TRUE)
    
    files <- list.files(input_dir, pattern = "Fusion_.*\\.csv$", full.names = TRUE)
    all_eval <- list()
    
    calc_overlap <- function(pred_df, chip_df) {
      pred_pairs <- paste(pred_df$TF, pred_df$Target)
      true_pairs <- paste(chip_df$TF, chip_df$Target)
      sum(pred_pairs %in% true_pairs)
    }
    
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
    
    for (f in files) {
      net_name <- gsub("Fusion_|\\.csv", "", basename(f))
      net <- vroom(f, show_col_types = FALSE)
      fused_total <- nrow(net)
      
      top_n <- max(1, ceiling(fused_total * top_ratio))
      net_high <- net %>% arrange(desc(EdgeWeight)) %>% slice_head(n = top_n)
      
      vroom_write(net_high, file.path(out_filt, paste0(net_name, "_top10.csv")), delim = ",")
      
      res <- evaluate_network(net_high, chip_seq)
      
      row <- data.frame(Network = net_name, TopRatio = top_ratio,
                        AUROC = res$AUROC, AUPR = res$AUPR,
                        True_Overlap = res$True_Overlap,
                        Fused_Edges = fused_total, HighConf_Edges = nrow(net_high))
      all_eval[[length(all_eval)+1]] <- row
    }
    
    final <- bind_rows(all_eval)
    write.csv(final, file.path(out_eval, paste0("Eval_", group_name, ".csv")), row.names = FALSE)
  }
  
  process_one_group("1_Pair_2Methods",    "Pair_2Methods")
  process_one_group("2_Triple_3Methods",  "Triple_3Methods")
  process_one_group("3_Four_4Methods",    "Four_4Methods")
  process_one_group("4_Five_5Methods",    "Five_5Methods")
  process_one_group("5_All_6Methods",     "All_6Methods")
}

# ======================
# Compute secondary metrics (AUROC > 0.7 filter removed)
# ======================
run_metrics <- function(dataset) {
  eval_root    <- paste0("Fusion_Evaluation_", dataset, "_Top10")
  OUTPUT_FILE  <- paste0(dataset, "_FINAL_TOP10.csv")
  SUMMARY_FILE <- paste0(dataset, "_Summary_Rankings_TOP10.csv")
  BASE_HC      <- eval_root
  
  groups <- c("Pair_2Methods","Triple_3Methods","Four_4Methods","Five_5Methods","All_6Methods")
  all_results <- list()
  for (g in groups) {
    f <- file.path(eval_root, g, "Evaluation", paste0("Eval_", g, ".csv"))
    if (file.exists(f)) all_results[[g]] <- read.csv(f) %>% mutate(Group = g)
  }
  all_df <- bind_rows(all_results)
  
  # Keep rows with non‑NA AUROC only; retain all rows by: filtered_df <- all_df
  filtered_df <- all_df %>% filter(!is.na(AUROC)) %>% arrange(desc(AUROC))
  
  # Load gold‑standard ChIP‑seq data
  chip <- read_csv("chip_seq.csv", col_names = c("TF","Target"), show_col_types = FALSE)
  true_pairs <- unique(paste(chip$TF, chip$Target))
  N_true <- length(true_pairs)
  valid_tfs <- unique(chip$TF)
  valid_targets <- unique(chip$Target)
  N_total_chip <- length(valid_tfs) * length(valid_targets)
  
  evaluate_metrics <- function(pred_df, valid_tfs, valid_targets, true_pairs, N_total_chip) {
    pred_filtered <- pred_df %>% filter(TF %in% valid_tfs, Target %in% valid_targets)
    pred_pairs <- unique(paste(pred_filtered$TF, pred_filtered$Target))
    TP <- sum(pred_pairs %in% true_pairs)
    FP <- length(pred_pairs) - TP
    FN <- N_true - TP
    TN <- max(N_total_chip - TP - FP - FN, 0)
    
    precision <- ifelse(TP+FP==0, 0, TP/(TP+FP))
    recall    <- ifelse(TP+FN==0, 0, TP/(TP+FN))
    F1        <- ifelse(precision+recall==0, 0, 2*precision*recall/(precision+recall))
    ACC       <- (TP + TN) / N_total_chip
    
    return(list(Accuracy=round(ACC,4), F1=round(F1,4), Precision=round(precision,4),
                Recall=round(recall,4), TP=TP, FP=FP, FN=FN, TN=TN,
                N_pred_overlap=length(pred_pairs)))
  }
  
  df <- filtered_df
  
  if (nrow(df) == 0) {
    warning(paste("No valid AUROC records found, skip metric computation for", dataset))
    empty_df <- data.frame(Group = character(), Network = character(), TopRatio = numeric(),
                           AUROC = numeric(), AUPR = numeric(), True_Overlap = integer(),
                           Fused_Edges = integer(), HighConf_Edges = integer(),
                           Accuracy = numeric(), F1 = numeric(), Precision = numeric(),
                           Recall = numeric(), TP = integer(), FP = integer(), FN = integer(),
                           TN = integer(), N_pred_overlap = integer())
    write.csv(empty_df, OUTPUT_FILE, row.names = FALSE)
    write.csv(empty_df %>% select(Group, Network, TopRatio, AUROC, F1, Precision, Recall, Accuracy),
              SUMMARY_FILE, row.names = FALSE)
    return()
  }
  
  df$Accuracy <- df$F1 <- df$Precision <- df$Recall <- NA
  df$TP <- df$FP <- df$FN <- df$TN <- df$N_pred_overlap <- NA
  
  for (i in 1:nrow(df)) {
    g <- df$Group[i]
    nw <- df$Network[i]
    fpath <- file.path(BASE_HC, g, "HighConf_Edges", paste0(nw, "_top10.csv"))
    if (!file.exists(fpath)) {
      next
    }
    net_high <- vroom(fpath, show_col_types = FALSE)
    res <- evaluate_metrics(net_high, valid_tfs, valid_targets, true_pairs, N_total_chip)
    
    df$Accuracy[i] <- res$Accuracy
    df$F1[i] <- res$F1
    df$Precision[i] <- res$Precision
    df$Recall[i] <- res$Recall
    df$TP[i] <- res$TP
    df$FP[i] <- res$FP
    df$FN[i] <- res$FN
    df$TN[i] <- res$TN
    df$N_pred_overlap[i] <- res$N_pred_overlap
  }
  
  write.csv(df, OUTPUT_FILE, row.names = FALSE)
  
  s <- df %>% select(Group, Network, TopRatio, AUROC, F1,
                     Precision, Recall, Accuracy) %>%
    arrange(desc(F1))
  write.csv(s, SUMMARY_FILE, row.names = FALSE)
}

# ======================
# Execute full assessment workflow
# ======================
for (ds in datasets) {
  cat("\n\n==================================================\n")
  cat("Processing dataset:", ds, "\n")
  cat("==================================================\n")
  run_evaluation(ds)
  run_metrics(ds)
}

cat("\nWorkflow completed. Assessment outputs generated for Aridity / Alkalinity / Cold using top‑10% edges.\n")
