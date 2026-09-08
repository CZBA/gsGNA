library(vroom)
library(dplyr)
library(pROC)
library(PRROC)
library(tidyr)

# ====================== Centralized global constants ======================
datasets           <- c("Aridity", "Alkalinity", "Cold")
fusion_root_main   <- "FeatureFusion_XGB_GRN_Opt"
eval_out_root      <- "XGB_Fusion_Evaluation_AllNonZeroEdges"
chip_file          <- "chip_seq.csv"
SCORE_COL          <- "XGB_Score"
FILE_PREFIX        <- "XGBFusion_"
group_map <- data.frame(
  folder = c("1_Pair_2Methods", "2_Triple_3Methods", "3_Four_4Methods", "4_Five_5Methods", "5_All_6Methods"),
  group_name = c("Pair_2Methods", "Triple_3Methods", "Four_4Methods", "Five_5Methods", "All_6Methods"),
  stringsAsFactors = FALSE
)
REQUIRED_COLS      <- c("TF", "Target", SCORE_COL)

# ====================== Utility function 1: preprocessing: filter zero scores + uppercase + deduplicate duplicate edges ======================
preprocess_net <- function(df, score_col = SCORE_COL) {
  df <- df %>%
    filter(.data[[score_col]] > 0) %>%
    mutate(
      TF = toupper(TF),
      Target = toupper(Target)
    ) %>%
    distinct(TF, Target, .keep_all = TRUE)
  return(df)
}

# ====================== Utility function 2: ChIP overlap count ======================
calc_overlap <- function(pred_df, chip_dt) {
  inner_join(pred_df, chip_dt, by = c("TF", "Target")) %>% nrow()
}

# ====================== Utility function 3: AUROC / AUPR calculation, returns the network with labels ======================
evaluate_network <- function(pred_df, chip_dt) {
  pred_label <- pred_df %>%
    left_join(chip_dt %>% mutate(label = 1), by = c("TF", "Target")) %>%
    mutate(label = replace_na(label, 0))
  
  true_overlap <- sum(pred_label$label == 1)
  if (length(unique(pred_label$label)) < 2) {
    return(list(AUROC = NA, AUPR = NA, True_Overlap = true_overlap, net_df = pred_label))
  }
  
  roc_obj <- roc(pred_label$label, pred_label[[SCORE_COL]], quiet = TRUE)
  auroc <- auc(roc_obj)
  pr_pos <- pred_label[[SCORE_COL]][pred_label$label == 1]
  pr_neg <- pred_label[[SCORE_COL]][pred_label$label == 0]
  pr_obj <- pr.curve(scores.class0 = pr_pos, scores.class1 = pr_neg)
  aupr <- pr_obj$auc.integral
  
  return(list(
    AUROC = round(as.numeric(auroc), 3),
    AUPR = round(aupr, 3),
    True_Overlap = true_overlap,
    net_df = pred_label
  ))
}

# ====================== Utility function 4: correct confusion matrix (with TN, Accuracy) ======================
# Candidate space = all TFs in ChIP x all Targets in ChIP; TP/FP/FN/TN are counted only within this space
get_metrics_correct <- function(pred_label_df, chip_dt) {
  # Extract the full TF and Target sets of the gold-standard space
  chip_tfs <- unique(chip_dt$TF)
  chip_tgts <- unique(chip_dt$Target)
  # All theoretical pairs within this space
  all_candidate_pairs <- expand.grid(
    TF = chip_tfs,
    Target = chip_tgts,
    stringsAsFactors = FALSE
  ) %>% distinct(TF, Target)
  total_candidate_space <- nrow(all_candidate_pairs)
  
  # True positive set
  true_pos_pairs <- chip_dt %>% distinct(TF, Target)
  total_true_pos <- nrow(true_pos_pairs)
  
  # Restrict predicted edges to the candidate space
  pred_in_space <- inner_join(pred_label_df, all_candidate_pairs, by = c("TF", "Target"))
  TP <- sum(pred_in_space$label == 1)
  FP <- sum(pred_in_space$label == 0)
  FN <- total_true_pos - TP
  # TN = total candidate space - predicted edges (TP+FP) - missed positives (FN)
  TN <- total_candidate_space - TP - FP - FN
  TN <- max(TN, 0)
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  F1        <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  ACC       <- ifelse(total_candidate_space == 0, 0, (TP + TN) / total_candidate_space)
  
  return(list(
    Accuracy = round(ACC, 4),
    F1 = round(F1, 4),
    Precision = round(precision, 4),
    Recall = round(recall, 4),
    TP = TP, FP = FP, FN = FN, TN = TN,
    N_pred_overlap = nrow(pred_in_space),
    Total_Candidate_Space = total_candidate_space
  ))
}

# ====================== Single-dataset batch AUROC/AUPR computation + cache networks ======================
run_eval_single_ds <- function(ds, chip_dt) {
  ds_fusion_dir <- file.path(fusion_root_main, ds)
  ds_eval_root  <- file.path(eval_out_root, paste0("Eval_", ds))
  dir.create(ds_eval_root, recursive = TRUE, showWarnings = FALSE)
  
  ds_all_metrics <- list()
  ds_net_cache   <- list() # Cache labeled networks to avoid second file reads
  
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
    
    # ========== Fixed reading logic: match ds_XGBFusion_xxx.csv ==========
    csv_files <- list.files(
      path = input_dir,
      pattern = "XGBFusion.*\\.csv$",
      ignore.case = TRUE,
      full.names = TRUE
    )
    if (length(csv_files) == 0) {
      cat("  [Warning] No XGB fusion result files in the current group\n")
      next
    }
    
    group_metric_list <- list()
    for (f in csv_files) {
      base_fn <- tools::file_path_sans_ext(basename(f))
      # Strip the preceding dataset prefix, keep only the combination name as the Network name
      net_name <- sub("^.*_XGBFusion_", "", base_fn)
      
      net_df_raw <- vroom(f, show_col_types = FALSE, progress = FALSE)
      
      # Column-name validation
      if (!all(REQUIRED_COLS %in% colnames(net_df_raw))) {
        cat("  [Warning] ", net_name, " lacks TF/Target/XGB_Score; skipping\n")
        next
      }
      
      net_valid <- preprocess_net(net_df_raw)
      edge_num <- nrow(net_valid)
      if (edge_num == 0) {
        cat("  [Warning] ", net_name, " has no edges with score > 0; skipping\n")
        next
      }
      
      eval_res <- evaluate_network(net_valid, chip_dt)
      cache_key <- paste(ds, group_nm, net_name, sep = "|")
      ds_net_cache[[cache_key]] <- eval_res$net_df
      
      line_df <- data.frame(
        CacheKey = cache_key,
        Network = net_name,
        GroupFolder = folder_nm,
        GroupName = group_nm,
        AUROC = eval_res$AUROC,
        AUPR = eval_res$AUPR,
        True_ChIP_Overlap = eval_res$True_Overlap,
        Total_Valid_Edges = edge_num,
        stringsAsFactors = FALSE
      )
      group_metric_list[[length(group_metric_list)+1]] <- line_df
    }
    
    group_df <- bind_rows(group_metric_list)
    if (nrow(group_df) > 0) {
      vroom_write(
        group_df,
        file.path(group_eval_dir, "Evaluation", paste0("Eval_", group_nm, ".csv")),
        delim = ",", progress = FALSE
      )
      ds_all_metrics[[length(ds_all_metrics)+1]] <- group_df
    }
  }
  
  ds_total_df <- bind_rows(ds_all_metrics)
  return(list(metric_df = ds_total_df, net_cache = ds_net_cache))
}

# ====================== Supplementary F1/Precision/Recall/TN/Acc classification metrics ======================
run_extra_metrics <- function(ds, eval_list, chip_dt) {
  full_eval_df <- eval_list$metric_df
  net_cache    <- eval_list$net_cache
  
  if (nrow(full_eval_df) == 0) {
    cat("\n[Warning] ", ds, " has no valid evaluation data; skipping classification-metric computation\n")
    return(invisible(NULL))
  }
  
  full_out_csv  <- file.path(eval_out_root, paste0(ds, "_XGB_FINAL_AllNonZero.csv"))
  summary_csv   <- file.path(eval_out_root, paste0(ds, "_XGB_Summary_AllNonZero.csv"))
  
  df <- full_eval_df %>% arrange(desc(AUPR))
  
  # New metric columns
  df$Accuracy      <- NA
  df$F1            <- NA
  df$Precision     <- NA
  df$Recall        <- NA
  df$TP            <- NA
  df$FP            <- NA
  df$FN            <- NA
  df$TN            <- NA
  df$N_pred_overlap<- NA
  df$Total_Candidate_Space <- NA
  
  for(i in 1:nrow(df)){
    ck <- df$CacheKey[i]
    if (!ck %in% names(net_cache)) next
    net_label_df <- net_cache[[ck]]
    
    m_res <- get_metrics_correct(net_label_df, chip_dt)
    
    df$Accuracy[i]             <- m_res$Accuracy
    df$F1[i]                   <- m_res$F1
    df$Precision[i]            <- m_res$Precision
    df$Recall[i]               <- m_res$Recall
    df$TP[i]                   <- m_res$TP
    df$FP[i]                   <- m_res$FP
    df$FN[i]                   <- m_res$FN
    df$TN[i]                   <- m_res$TN
    df$N_pred_overlap[i]       <- m_res$N_pred_overlap
    df$Total_Candidate_Space[i]<- m_res$Total_Candidate_Space
  }
  
  vroom_write(df, full_out_csv, delim = ",", progress = FALSE)
  
  sum_df <- df %>%
    select(
      GroupName, Network, AUROC, AUPR, F1, Precision, Recall, Accuracy,
      TP, FP, FN, TN, Total_Valid_Edges, Total_Candidate_Space
    ) %>%
    arrange(desc(AUPR))
  vroom_write(sum_df, summary_csv, delim = ",", progress = FALSE)
  
  cat("\n[Done] ", ds, " XGB-fusion all-valid-edge evaluation completed (with TN/Acc)\n")
}

# ====================== Main batch workflow ======================
if (!file.exists(chip_file)) stop("Missing chip_seq.csv")
# Read and preprocess the ChIP gold standard
chip_global <- vroom(chip_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    TF = toupper(gene_id),
    Target = toupper(TARGET)
  ) %>%
  select(TF, Target) %>%
  drop_na(TF, Target) %>%
  filter(TF != Target) %>%
  distinct(TF, Target)

for(ds in datasets){
  cat("\n\n==================================================\n")
  cat("Starting all-valid-edge evaluation of XGB feature fusion | Dataset: ", ds, "\n")
  cat("==================================================\n")
  
  eval_res_list <- run_eval_single_ds(ds, chip_global)
  run_extra_metrics(ds, eval_res_list, chip_global)
}

cat("\n[Finished] All datasets completed: only XGB_Score=0 filtered, all non-zero scoring regulatory edges evaluated\nOutput root directory: ", eval_out_root, "\n")
