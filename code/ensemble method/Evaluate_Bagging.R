library(vroom)
library(dplyr)
library(pROC)
library(PRROC)
library(tidyr)

# ====================== Global parameters ======================
datasets   <- c("Aridity", "Alkalinity", "Cold")
fusion_root_main   <- "Bagging_Unsup_AlgoBootstrap"
eval_out_root      <- "Bagging_Unsup_Evaluation_AllEdges"
chip_file  <- "chip_seq.csv"

# Group mapping: ensemble folder name <-> group alias
group_map <- data.frame(
  folder = c("Pair2", "Triple3", "Four4", "Five5", "All6"),
  group_name = c("Pair_2Methods", "Triple_3Methods", "Four_4Methods", "Five_5Methods", "All_6Methods")
)

# ====================== Low-level evaluation utility functions ======================
calc_overlap <- function(pred_df, chip_df) {
  pred_pairs <- paste(pred_df$TF, pred_df$Target)
  true_pairs <- paste(chip_df$TF, chip_df$Target)
  sum(pred_pairs %in% true_pairs)
}

evaluate_network <- function(pred_df, chip_df) {
  # Core fix 1: unify the sample space, deduplicate keeping the max weight per TF-Target pair
  pred_df <- pred_df %>%
    group_by(TF, Target) %>%
    slice_max(EdgeWeight, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  pred_pairs <- paste(pred_df$TF, pred_df$Target)
  true_pairs <- paste(chip_df$TF, chip_df$Target)
  pred_df$label <- as.integer(pred_pairs %in% true_pairs)
  
  if (length(unique(pred_df$label)) < 2) {
    return(list(AUROC = NA, AUPR = NA, True_Overlap = sum(pred_pairs %in% true_pairs)))
  }
  
  roc_obj <- roc(pred_df$label, pred_df$EdgeWeight, quiet = TRUE)
  auroc <- auc(roc_obj)
  
  # Comment distinguishes positive/negative samples to prevent swapped arguments
  pr_obj <- pr.curve(
    scores.class0 = pred_df$EdgeWeight[pred_df$label == 1], # positive-sample scores
    scores.class1 = pred_df$EdgeWeight[pred_df$label == 0]  # negative-sample scores
  )
  aupr <- pr_obj$auc.integral
  
  return(list(AUROC = round(auroc,3), AUPR = round(aupr,3), True_Overlap = sum(pred_pairs %in% true_pairs)))
}

# Core fix 2: reconstruct the TN calculation logic, rename the ambiguous variable
get_metrics <- function(pred_df, valid_tfs, valid_targets, true_pairs, total_candidate_pairs) {
  pred_filtered <- pred_df %>% filter(TF %in% valid_tfs, Target %in% valid_targets)
  pred_pairs <- unique(paste(pred_filtered$TF, pred_filtered$Target))
  
  P <- length(true_pairs)                  # total number of true positive interaction pairs
  all_neg_candidates <- total_candidate_pairs - P # all negative candidate interaction pairs
  
  TP <- sum(pred_pairs %in% true_pairs)
  FP <- length(pred_pairs) - TP
  FN <- P - TP
  TN <- all_neg_candidates - FP            # standardized TN formula
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  F1        <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  ACC       <- (TP + TN) / total_candidate_pairs
  
  return(list(
    Accuracy = round(ACC,4),
    F1 = round(F1,4),
    Precision = round(precision,4),
    Recall = round(recall,4),
    TP = TP, FP = FP, FN = FN, TN = TN,
    N_pred_overlap = length(pred_pairs)
  ))
}

# ====================== No-filter function, returns all original edges ======================
filter_top10pct <- function(df) {
  return(df)
}

# ====================== Single dataset: batch AUROC/AUPR/overlap computation ======================
run_eval_single_ds <- function(ds, chip_seq) {
  ds_fusion_dir <- file.path(fusion_root_main, ds)
  ds_eval_root  <- file.path(eval_out_root, paste0("Eval_", ds))
  dir.create(ds_eval_root, recursive = TRUE, showWarnings = FALSE)
  
  ds_all_metrics <- list()
  
  for (g_idx in 1:nrow(group_map)) {
    folder_nm  <- group_map$folder[g_idx]
    group_nm   <- group_map$group_name[g_idx]
    cat("\n=========================================\n")
    cat("Dataset", ds, " | Group: ", group_nm, "\n")
    cat("=========================================\n")
    
    input_dir <- file.path(ds_fusion_dir, folder_nm)
    if (!dir.exists(input_dir)) {
      cat("  [Warning] Directory does not exist; skipping: ", input_dir, "\n")
      next
    }
    
    group_eval_dir <- file.path(ds_eval_root, group_nm)
    dir.create(file.path(group_eval_dir, "Evaluation"), recursive = TRUE, showWarnings = FALSE)
    
    csv_files <- list.files(
      path = input_dir,
      pattern = "^Bagging_.*\\.csv$",
      full.names = TRUE
    )
    if (length(csv_files) == 0) {
      cat("  [Warning] No result files in the current group\n")
      next
    }
    
    group_metric_list <- list()
    for (f in csv_files) {
      net_name <- gsub("^Bagging_|\\.csv$", "", basename(f))
      net_df <- vroom(f, show_col_types = FALSE)
      
      # New: clean weight anomalies NA / Inf / -Inf
      net_df <- net_df %>% filter(!is.na(EdgeWeight), is.finite(EdgeWeight))
      
      # Use the raw data directly without any reduction
      net_all <- filter_top10pct(net_df)
      edge_num <- nrow(net_all)
      
      if (edge_num == 0) {
        cat("  [Warning] ", net_name, " has no valid edges; skipping\n")
        next
      }
      
      eval_res <- evaluate_network(net_all, chip_seq)
      line_df <- data.frame(
        Network = net_name,
        GroupFolder = folder_nm,
        GroupName = group_nm,
        AUROC = eval_res$AUROC,
        AUPR = eval_res$AUPR,
        True_ChIP_Overlap = eval_res$True_Overlap,
        Total_All_Edges = edge_num
      )
      group_metric_list[[length(group_metric_list)+1]] <- line_df
    }
    
    group_df <- bind_rows(group_metric_list)
    write.csv(
      group_df,
      file.path(group_eval_dir, "Evaluation", paste0("Eval_", group_nm, ".csv")),
      row.names = FALSE
    )
    if (nrow(group_df) > 0) {
      ds_all_metrics[[length(ds_all_metrics)+1]] <- group_df
    }
  }
  
  ds_total_df <- bind_rows(ds_all_metrics)
  return(ds_total_df)
}

# ====================== Single dataset: supplementary F1/Precision/Recall/TP/TN metrics ======================
run_extra_metrics <- function(ds, full_eval_df, chip_seq) {
  if (nrow(full_eval_df) == 0) {
    cat("\n[Warning] ", ds, " has no valid evaluation data; skipping classification-metric computation\n")
    return(invisible(NULL))
  }
  
  full_out_csv  <- file.path(eval_out_root, paste0(ds, "_Bagging_FINAL_AllEdges.csv"))
  summary_csv   <- file.path(eval_out_root, paste0(ds, "_Bagging_Summary_AllEdges.csv"))
  
  true_pairs    <- unique(paste(chip_seq$TF, chip_seq$Target))
  valid_tfs     <- unique(chip_seq$TF)
  valid_targets <- unique(chip_seq$Target)
  total_candidate_pairs  <- length(valid_tfs) * length(valid_targets) # renamed ambiguous variable
  
  df <- full_eval_df %>% arrange(desc(AUROC))
  
  # Initialize metric columns
  df$Accuracy      <- NA
  df$F1            <- NA
  df$Precision     <- NA
  df$Recall        <- NA
  df$TP            <- NA
  df$FP            <- NA
  df$FN            <- NA
  df$TN            <- NA
  df$N_pred_overlap<- NA
  df$Dataset       <- ds # New dataset identifier column
  
  for(i in 1:nrow(df)){
    folder  <- df$GroupFolder[i]
    netname <- df$Network[i]
    
    fpath <- file.path(
      fusion_root_main,
      ds,
      folder,
      paste0("Bagging_", netname, ".csv")
    )
    if(!file.exists(fpath)) {
      cat("[Warning] Missing network file: ", fpath, "\n") # Print missing log
      next
    }
    
    net_raw <- vroom(fpath, show_col_types = FALSE)
    # Clean weight anomalies
    net_raw <- net_raw %>% filter(!is.na(EdgeWeight), is.finite(EdgeWeight))
    net_all <- filter_top10pct(net_raw)
    
    if(nrow(net_all) == 0) next
    
    # Pass the corrected variable name total_candidate_pairs
    m_res <- get_metrics(net_all, valid_tfs, valid_targets, true_pairs, total_candidate_pairs)
    
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
  
  write.csv(df, full_out_csv, row.names = FALSE)
  
  # Compact summary table, Dataset field placed first for easy merging
  sum_df <- df %>%
    select(Dataset, GroupName, Network, AUROC, AUPR, F1, Precision, Recall, Accuracy, Total_All_Edges) %>%
    arrange(desc(AUROC))
  write.csv(sum_df, summary_csv, row.names = FALSE)
  
  cat("\n[Done] ", ds, " all-edge evaluation metrics computed\n")
}

# ====================== Main batch workflow ======================
if (!file.exists(chip_file)) stop("Missing chip_seq.csv")
chip_global <- vroom(chip_file, show_col_types = FALSE) %>%
  mutate(
    TF = toupper(gene_id),
    Target = toupper(TARGET)
  )

for(ds in datasets){
  cat("\n\n==================================================\n")
  cat("Starting all-edge evaluation of Bagging | Dataset: ", ds, "\n")
  cat("==================================================\n")
  
  eval_all_df <- run_eval_single_ds(ds, chip_global)
  run_extra_metrics(ds, eval_all_df, chip_global)
}

cat("\n[Finished] All datasets completed: no filtering applied, evaluated with all edges of the networks\nOutput directory: ", eval_out_root, "\n")
