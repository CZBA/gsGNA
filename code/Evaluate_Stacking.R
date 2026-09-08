library(vroom)
library(dplyr)
library(pROC)
library(PRROC)
library(tidyr)

# ====================== Global tunable parameters ======================
datasets   <- c("Aridity", "Alkalinity", "Cold")
fusion_root_main   <- "Stacking_Supervised_ChIP"
eval_out_root      <- "Stacking_Sup_Evaluation_Top10Pct"
chip_file  <- "chip_seq.csv"
top_pct    <- 0.10    # Filter the top X%, unified entry
seed_fix   <- 123     # Fixed random seed for reproducibility

# Group mapping
group_map <- data.frame(
  folder = c("1_Pair_2Methods", "2_Triple_3Methods", "3_Four_4Methods", "4_Five_5Methods", "5_All_6Methods"),
  group_name = c("Pair_2Methods", "Triple_3Methods", "Four_4Methods", "Five_5Methods", "All_6Methods")
)

# ====================== Low-level utility functions ======================
# Compute the number of predicted edges overlapping with ChIP true edges
calc_overlap <- function(pred_df, chip_df) {
  pred_pairs <- paste(pred_df$TF, pred_df$Target)
  true_pairs <- paste(chip_df$TF, chip_df$Target)
  sum(pred_pairs %in% true_pairs)
}

# Compute AUROC, AUPR, and the number of ChIP-overlapping edges
evaluate_network <- function(pred_df, chip_df) {
  true_overlap <- calc_overlap(pred_df, chip_df)
  pred_pairs <- paste(pred_df$TF, pred_df$Target)
  true_pairs <- paste(chip_df$TF, chip_df$Target)
  pred_df$label <- as.integer(pred_pairs %in% true_pairs)
  
  # Set metrics to NA if positive/negative samples are incomplete
  if (length(unique(pred_df$label)) < 2) {
    return(list(AUROC = NA, AUPR = NA, True_Overlap = true_overlap))
  }
  
  roc_obj <- roc(pred_df$label, pred_df$Stacking_Score, quiet = TRUE)
  auroc <- auc(roc_obj)
  
  # pr.curve: scores.class0 = positive samples, scores.class1 = negative samples
  pr_obj <- pr.curve(
    scores.class0 = pred_df$Stacking_Score[pred_df$label == 1],
    scores.class1 = pred_df$Stacking_Score[pred_df$label == 0]
  )
  aupr <- pr_obj$auc.integral
  
  return(list(AUROC = round(auroc,3), AUPR = round(aupr,3), True_Overlap = true_overlap))
}

# Compute Precision/Recall/F1/ACC/confusion matrix
get_metrics <- function(pred_filtered, valid_tfs, valid_targets, true_pairs) {
  pred_pairs <- unique(paste(pred_filtered$TF, pred_filtered$Target))
  
  TP <- sum(pred_pairs %in% true_pairs)
  FP <- length(pred_pairs) - TP
  true_positive_total <- length(true_pairs)
  FN <- true_positive_total - TP
  
  # Full Cartesian product of TF-Target
  all_tf_target_pairs <- length(valid_tfs) * length(valid_targets)
  all_neg <- all_tf_target_pairs - true_positive_total
  TN <- max(all_neg - FP, 0)
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  F1        <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  ACC       <- (TP + TN) / all_tf_target_pairs
  
  return(list(
    Accuracy = round(ACC,4),
    F1 = round(F1,4),
    Precision = round(precision,4),
    Recall = round(recall,4),
    TP = TP, FP = FP, FN = FN, TN = TN,
    N_pred_overlap = length(pred_pairs)
  ))
}

# ====================== Unified filter function: remove zero scores + keep the top top_pct high-scoring edges ======================
filter_top_pct <- function(df, pct = top_pct, seed = seed_fix) {
  set.seed(seed)
  df <- df %>% filter(Stacking_Score > 0)
  if (nrow(df) == 0) return(df)
  
  df <- df %>% arrange(desc(Stacking_Score))
  total_n <- nrow(df)
  keep_n  <- ceiling(total_n * pct)
  df_top <- df[1:keep_n, ]
  return(df_top)
}

# ====================== Single-dataset batch AUROC/AUPR computation, cache filtered subsets ======================
run_eval_single_ds <- function(ds, chip_seq) {
  ds_fusion_dir <- file.path(fusion_root_main, ds)
  ds_eval_root  <- file.path(eval_out_root, paste0("Eval_", ds))
  ds_cache_dir  <- file.path(ds_eval_root, "Cache_FilteredSubset")
  
  # Fix: create directories twice, cannot pass a vector
  dir.create(ds_eval_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_cache_dir, recursive = TRUE, showWarnings = FALSE)
  
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
      pattern = "^SupStacking_.*\\.csv$",
      full.names = TRUE
    )
    if (length(csv_files) == 0) {
      cat("  [Warning] No Stacking result files in the current group\n")
      next
    }
    
    group_metric_list <- list()
    for (f in csv_files) {
      tryCatch({
        net_name <- gsub("^SupStacking_|\\.csv$", "", basename(f))
        net_df <- vroom(f, show_col_types = FALSE)
        
        # Filter the high-scoring subset
        net_top <- filter_top_pct(net_df)
        edge_num <- nrow(net_top)
        
        if (edge_num == 0) {
          cat("  [Warning] ", net_name, " has no valid edges after filtering; skipping\n")
          next
        }
        
        # Cache the filtered subset; F1 computation later avoids re-reading the large original file
        cache_rds <- file.path(ds_cache_dir, paste0(ds, "_", group_nm, "_", net_name, ".rds"))
        saveRDS(net_top, file = cache_rds)
        
        # AUROC/AUPR evaluation
        eval_res <- evaluate_network(net_top, chip_seq)
        line_df <- data.frame(
          Network = net_name,
          GroupFolder = folder_nm,
          GroupName = group_nm,
          CacheRDS = cache_rds,
          AUROC = eval_res$AUROC,
          AUPR = eval_res$AUPR,
          True_ChIP_Overlap = eval_res$True_Overlap,
          Total_TopPct_Edges = edge_num
        )
        group_metric_list <- append(group_metric_list, list(line_df))
      }, error = function(e) {
        cat("  [Error] Failed to process file: ", f, " error: ", e$message, "\n")
      })
    }
    
    group_df <- bind_rows(group_metric_list)
    if (nrow(group_df) > 0) {
      write.csv(
        group_df,
        file.path(group_eval_dir, "Evaluation", paste0("Eval_", group_nm, ".csv")),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
      ds_all_metrics <- append(ds_all_metrics, list(group_df))
    }
  }
  
  ds_total_df <- bind_rows(ds_all_metrics)
  return(ds_total_df)
}

# ====================== Read cached subsets, supplement F1/Precision/Recall/ACC ======================
run_extra_metrics <- function(ds, full_eval_df, chip_seq) {
  if (nrow(full_eval_df) == 0) {
    cat("\n[Warning] ", ds, " has no valid evaluation data; skipping classification-metric computation\n")
    return(invisible(NULL))
  }
  
  full_out_csv  <- file.path(eval_out_root, paste0(ds, "_Stacking_FINAL_TopPct.csv"))
  summary_csv   <- file.path(eval_out_root, paste0(ds, "_Stacking_Summary_TopPct.csv"))
  
  true_pairs    <- unique(paste(chip_seq$TF, chip_seq$Target))
  valid_tfs     <- unique(chip_seq$TF)
  valid_targets <- unique(chip_seq$Target)
  
  df <- full_eval_df %>% arrange(desc(AUROC))
  
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
    cache_path <- df$CacheRDS[i]
    if(!file.exists(cache_path)) next
    
    tryCatch({
      net_top <- readRDS(cache_path)
      m_res <- get_metrics(net_top, valid_tfs, valid_targets, true_pairs)
      
      df$Accuracy[i]       <- m_res$Accuracy
      df$F1[i]             <- m_res$F1
      df$Precision[i]      <- m_res$Precision
      df$Recall[i]         <- m_res$Recall
      df$TP[i]             <- m_res$TP
      df$FP[i]             <- m_res$FP
      df$FN[i]             <- m_res$FN
      df$TN[i]             <- m_res$TN
      df$N_pred_overlap[i] <- m_res$N_pred_overlap
    }, error = function(e) {
      cat("  [Error] Cache read failed: ", cache_path, e$message, "\n")
    })
  }
  
  # Full metrics table
  write.csv(df, full_out_csv, row.names = FALSE, fileEncoding = "UTF-8")
  
  # Compact summary table
  sum_df <- df %>%
    select(GroupName, Network, AUROC, AUPR, F1, Precision, Recall, Accuracy, Total_TopPct_Edges) %>%
    arrange(desc(AUPR))
  write.csv(sum_df, summary_csv, row.names = FALSE, fileEncoding = "UTF-8")
  
  cat("\n[Done] ", ds, " Stacking top", round(top_pct*100), "% high-confidence edge evaluation completed\n")
}

# ====================== Main entry ======================
# 1. Validate the existence of the ChIP file
if (!file.exists(chip_file)) stop("Missing file: ", chip_file)

# 2. Validate the required ChIP columns
chip_raw <- vroom(chip_file, show_col_types = FALSE)
required_chip_cols <- c("gene_id", "TARGET")
missing_cols <- setdiff(required_chip_cols, colnames(chip_raw))
if (length(missing_cols) > 0) {
  stop(paste0("ChIP file is missing required columns: ", paste(missing_cols, collapse = ",")))
}

# 3. Unify TF/Target to uppercase for matching
chip_global <- chip_raw %>%
  mutate(
    TF = toupper(gene_id),
    Target = toupper(TARGET)
  )

# 4. Iterate over all datasets for evaluation
for(ds in datasets){
  cat("\n\n==================================================\n")
  cat("Starting top", round(top_pct*100), "% high-confidence edge evaluation of Stacking (supervised ChIP) | Dataset: ", ds, "\n")
  cat("==================================================\n")
  
  eval_all_df <- run_eval_single_ds(ds, chip_global)
  run_extra_metrics(ds, eval_all_df, chip_global)
}

cat("\n[Finished] All datasets completed: filter out Stacking_Score=0 edges + only the top", round(top_pct*100), "% high-confidence regulatory edges evaluated\nOutput root directory: ", eval_out_root, "\n")
