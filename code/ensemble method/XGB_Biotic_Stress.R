library(dplyr)
library(tidyr)
library(vroom)
library(xgboost)
library(caret)
library(pROC)
library(PRROC)
library(matrixStats)  # New: for colMins/colMaxs

# ====================== Global parameters ======================
set.seed(123)
sup_neg_ratio    <- 1
only_all6        <- FALSE
normalize_weights <- TRUE
test_ratio       <- 0.2
top_percent      <- 0.1        # Keep the top 10% of high-scoring edges
write_full_all_edges <- FALSE
compress_output <- FALSE

# XGB hyperparameters
xgb_nrounds      <- 500
xgb_early_stop   <- 30
xgb_max_depth    <- 3
xgb_eta          <- 0.15
xgb_subsample    <- 0.8
xgb_colsample    <- 0.8

chip_file_path   <- "chip_seq.csv"
root_dir         <- "FeatureFusion_XGB_GRN_Opt"
method_folders   <- c("GENIE3", "Kboost", "GRNBoost", "3DCEMA", "DeepRIG", "IGEGRNs")
# =====================================================

# New: fill_na_zero function definition
fill_na_zero <- function(df, cols) {
  # Replace NA with 0 in the specified columns
  for (col in cols) {
    if (col %in% colnames(df)) {
      df[[col]][is.na(df[[col]])] <- 0
    }
  }
  return(df)
}

# Write function
safe_write_csv <- function(df, filepath) {
  if (nrow(df) == 0) return(invisible(NULL))
  vroom_write(df, filepath, delim = ",", progress = FALSE)
}

# Read a single-model GRN file
clean_read <- function(path, delim = ",") {
  df <- vroom(path, delim = delim, show_col_types = FALSE, progress = FALSE)
  df <- df[, 1:3]
  colnames(df) <- c("TF", "Target", "Weight")
  df <- df %>% drop_na(TF, Target, Weight) %>% filter(is.finite(Weight))
  df <- df %>% filter(TF != Target)
  df <- df %>% group_by(TF, Target) %>% summarise(Weight = max(Weight), .groups = "drop")
  gc()
  return(df)
}

# Normalization function
normalize_by_train_stats <- function(train_mat, full_mat) {
  col_min <- colMins(train_mat, na.rm = TRUE)
  col_max <- colMaxs(train_mat, na.rm = TRUE)
  range_val <- col_max - col_min
  range_val[range_val == 0] <- 1
  full_scaled <- sweep(full_mat, 2, col_min, "-")
  full_scaled <- sweep(full_scaled, 2, range_val, "/")
  full_scaled[is.na(full_scaled)] <- 0
  return(full_scaled)
}

# Read the ChIP gold standard
read_chip_pos <- function(chip_path) {
  chip_df <- vroom(chip_path, show_col_types = FALSE, progress = FALSE)
  chip_df <- chip_df[, 1:2]
  colnames(chip_df) <- c("TF", "Target")
  chip_df <- chip_df %>% drop_na(TF, Target) %>% filter(TF != Target) %>% distinct(TF, Target)
  gc()
  return(chip_df)
}

# AUC/PR calculation
calc_metric <- function(test_df, score_col = "XGB_Score") {
  test_df <- test_df %>% filter(!is.na(!!sym(score_col)))
  lab <- test_df$label
  scr <- test_df[[score_col]]
  if (length(unique(lab)) < 2 || sum(lab == 1) == 0 || sum(lab == 0) == 0) {
    return(c(AUC_ROC = NA, AUC_PR = NA))
  }
  roc_obj <- roc(lab, scr, quiet = TRUE)
  auc_roc <- round(as.numeric(auc(roc_obj)), 4)
  pr_obj <- pr.curve(scores.class0 = scr[lab == 1], scores.class1 = scr[lab == 0], curve = FALSE)
  auc_pr <- round(pr_obj$auc.integral, 4)
  return(c(AUC_ROC = auc_roc, AUC_PR = auc_pr))
}

# Main XGB training function
pure_xgb_grn_optim <- function(method_list, all_base_wide, chip_pos, neg_ratio = 1, normalize_w = TRUE, test_frac = 0.2) {
  method_names <- names(method_list)
  feat_cols <- method_names
  wide_subset <- all_base_wide %>% select(TF, Target, all_of(feat_cols)) %>% distinct(TF, Target, .keep_all = TRUE)
  wide_subset <- fill_na_zero(wide_subset, feat_cols)  # Now defined above
  
  all_edge_key <- wide_subset %>% select(TF, Target)
  all_pos_edge <- inner_join(chip_pos, all_edge_key, by = c("TF", "Target"))
  n_pos_total <- nrow(all_pos_edge)
  if (n_pos_total == 0) {
    warning("No matching ChIP positive edges; skipping")
    return(list(result_full = data.frame(), result_save = data.frame(), metrics = c(AUC_ROC = NA, AUC_PR = NA), sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, best_nrounds=0)))
  }
  
  all_neg_candidate <- anti_join(all_edge_key, chip_pos, by = c("TF", "Target"))
  if (nrow(all_neg_candidate) == 0) {
    warning("No negative-sample pool; skipping")
    return(list(result_full = data.frame(), result_save = data.frame(), metrics = c(AUC_ROC = NA, AUC_PR = NA), sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, best_nrounds=0)))
  }
  
  # Stratified train/test split
  n_pos_test <- max(1, floor(n_pos_total * test_frac))
  n_pos_train <- n_pos_total - n_pos_test
  if (n_pos_train < 2) {
    warning("Too few positive samples to train")
    return(list(result_full = data.frame(), result_save = data.frame(), metrics = c(AUC_ROC = NA, AUC_PR = NA), sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, best_nrounds=0)))
  }
  
  if (n_distinct(all_pos_edge$TF) >= 2) {
    train_idx <- createDataPartition(all_pos_edge$TF, p = 1 - test_frac, list = FALSE)
    pos_train <- all_pos_edge[train_idx, ]
    pos_test <- all_pos_edge[-train_idx, ]
  } else {
    pos_idx <- sample(seq_len(n_pos_total), n_pos_test)
    pos_test <- all_pos_edge[pos_idx, ]
    pos_train <- all_pos_edge[-pos_idx, ]
  }
  
  n_neg_total <- nrow(all_neg_candidate)
  n_neg_test <- min(n_pos_test, n_neg_total)
  if (n_neg_test > 0 && n_distinct(all_neg_candidate$TF) >= 2) {
    neg_test_idx <- createDataPartition(all_neg_candidate$TF, p = n_neg_test / n_neg_total, list = FALSE)
    neg_test <- all_neg_candidate[neg_test_idx, ]
    neg_train_pool <- all_neg_candidate[-neg_test_idx, ]
  } else {
    neg_idx <- sample(seq_len(n_neg_total), n_neg_test)
    neg_test <- all_neg_candidate[neg_idx, ]
    neg_train_pool <- all_neg_candidate[-neg_idx, ]
  }
  
  n_neg_train_need <- floor(n_pos_train * neg_ratio)
  if (nrow(neg_train_pool) < n_neg_train_need) {
    neg_train <- neg_train_pool
  } else {
    set.seed(123)
    neg_train <- slice_sample(neg_train_pool, n = n_neg_train_need)
  }
  
  tr_pos_n <- nrow(pos_train)
  tr_neg_n <- nrow(neg_train)
  te_pos_n <- nrow(pos_test)
  te_neg_n <- nrow(neg_test)
  
  train_set <- bind_rows(pos_train %>% mutate(label = 1), neg_train %>% mutate(label = 0))
  if (length(unique(train_set$label)) < 2) {
    warning("The training set lacks one class label")
    return(list(result_full = data.frame(), result_save = data.frame(), metrics = c(AUC_ROC = NA, AUC_PR = NA), sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, te_neg=te_neg_n, best_nrounds=0)))
  }
  
  pred_set <- wide_subset %>% select(TF, Target) %>% distinct()
  train_feat_raw <- left_join(train_set, wide_subset, by = c("TF", "Target"))
  pred_feat_raw  <- left_join(pred_set, wide_subset, by = c("TF", "Target"))
  
  x_train_raw <- as.matrix(train_feat_raw[, feat_cols, drop = FALSE])
  x_pred_raw  <- as.matrix(pred_feat_raw[, feat_cols, drop = FALSE])
  
  # Filter zero-variance features
  zero_var <- apply(x_train_raw, 2, var, na.rm = TRUE) == 0
  if (any(zero_var)) {
    feat_cols <- feat_cols[!zero_var]
    x_train_raw <- x_train_raw[, !zero_var, drop = FALSE]
    x_pred_raw  <- x_pred_raw[, !zero_var, drop = FALSE]
  }
  if (ncol(x_train_raw) == 0) {
    warning("No usable features after filtering")
    return(list(result_full = data.frame(), result_save = data.frame(), metrics = c(AUC_ROC = NA, AUC_PR = NA), sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, te_neg=te_neg_n, best_nrounds=0)))
  }
  
  if (normalize_w) {
    x_train <- normalize_by_train_stats(x_train_raw, x_train_raw)
    x_pred  <- normalize_by_train_stats(x_train_raw, x_pred_raw)
  } else {
    x_train <- x_train_raw
    x_pred  <- x_pred_raw
  }
  y_train <- train_feat_raw$label
  
  dtrain <- xgb.DMatrix(data = x_train, label = y_train)
  xgb_params <- list(
    objective = "binary:logistic", 
    eval_metric = "logloss", 
    max_depth = xgb_max_depth, 
    eta = xgb_eta, 
    subsample = xgb_subsample, 
    colsample_bytree = xgb_colsample, 
    seed = 123
  )
  
  cv_folds <- min(5, max(2, floor(nrow(train_set) / 2)))
  cv_result <- xgb.cv(
    params = xgb_params, 
    data = dtrain, 
    nrounds = xgb_nrounds, 
    nfold = cv_folds, 
    early_stopping_rounds = xgb_early_stop, 
    verbose = 0, 
    prediction = FALSE, 
    stratified = TRUE
  )
  best_nrounds <- ifelse(is.null(cv_result$best_iteration) || cv_result$best_iteration == 0, 
                         xgb_nrounds, 
                         cv_result$best_iteration)
  
  fit <- xgb.train(params = xgb_params, data = dtrain, nrounds = best_nrounds, verbose = 0)
  
  dpred <- xgb.DMatrix(x_pred)
  pred_set$XGB_Score <- round(predict(fit, newdata = dpred), 4)
  result_full <- pred_set %>% select(TF, Target, XGB_Score) %>% arrange(desc(XGB_Score))
  
  # Test-set evaluation
  test_set <- bind_rows(pos_test %>% mutate(label=1), neg_test %>% mutate(label=0))
  test_feat <- inner_join(test_set, result_full, by = c("TF", "Target"))
  metrics <- calc_metric(test_feat, "XGB_Score")
  
  # Keep the top 10%
  total_all <- nrow(result_full)
  keep_num <- floor(total_all * top_percent)
  if (keep_num < 1) keep_num <- 1
  result_save <- slice_head(result_full, n = keep_num)
  
  cat(sprintf("    Training set: %d pos, %d neg | Test set: %d pos, %d neg | Best nrounds: %d\n", 
              tr_pos_n, tr_neg_n, te_pos_n, te_neg_n, best_nrounds))
  cat(sprintf("    Test AUC-ROC = %.4f, AUC-PR = %.4f | Output top%.0f%% edges: %d / %d\n", 
              metrics["AUC_ROC"], metrics["AUC_PR"], top_percent*100, nrow(result_save), total_all))
  
  gc()
  return(list(result_full = result_full, result_save = result_save, metrics = metrics, 
              sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, 
                              te_neg=te_neg_n, best_nrounds = best_nrounds)))
}

# Exception-capture wrapper
safe_run_xgb <- function(method_subset, all_base_wide, chip_pos, comb_name, neg_ratio, normalize_w, test_frac, dataset_name) {
  res <- tryCatch({
    pure_xgb_grn_optim(method_subset, all_base_wide, chip_pos, neg_ratio = neg_ratio, 
                       normalize_w = normalize_w, test_frac = test_frac)
  }, error = function(e) {
    cat("  >> Combination", comb_name, "error:", conditionMessage(e), "\n")
    return(list(result_full = data.frame(), result_save = data.frame(), 
                metrics = c(AUC_ROC = NA, AUC_PR = NA), 
                sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, best_nrounds=0)))
  })
  
  info <- res$sample_info
  m_val <- res$metrics
  out_df <- data.frame(
    Dataset = dataset_name, 
    Combination = comb_name, 
    NumMethods = length(method_subset),
    Methods = paste(sort(names(method_subset)), collapse = ";"),
    Train_Pos = info["tr_pos"], 
    Train_Neg = info["tr_neg"],
    Test_Pos = info["te_pos"], 
    Test_Neg = info["te_neg"],
    Best_Nrounds = info["best_nrounds"], 
    AUC_ROC = m_val["AUC_ROC"], 
    AUC_PR = m_val["AUC_PR"],
    stringsAsFactors = FALSE
  )
  return(list(out_df = out_df, pred_save = res$result_save, pred_full = res$result_full))
}

# Build the wide table
build_wide_table <- function(all_methods, method_names) {
  all_edges <- lapply(all_methods, function(df) df %>% select(TF, Target, Weight)) %>% 
    bind_rows(.id = "Method") %>%
    pivot_wider(id_cols = c(TF, Target), names_from = Method, values_from = Weight, values_fill = 0)
  missing_cols <- setdiff(method_names, colnames(all_edges))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) all_edges[[col]] <- 0
  }
  return(all_edges)
}

# ====================== Main program ======================
# Update: use the new dataset names
datasets <- c("MO", "NC", "RSV")
chip_gt_pos <- read_chip_pos(chip_file_path)
eval_summary <- data.frame()

for (ds in datasets) {
  cat("\n=============================================\n")
  cat("Optimized single-layer XGB feature fusion: ", ds, "\n")
  cat("=============================================\n")
  set.seed(123)
  
  # Update: dataset names in the file paths
  base_files <- list(
    GENIE3 = file.path("GENIE3", paste0(ds, "_significant_edges.txt")),
    Kboost = file.path("Kboost", paste0(ds, "_significant_edges.csv")),
    GRNBoost = file.path("GRNBoost", paste0(ds, "_significant_edges.csv")),
    `3DCEMA` = file.path("3DCEMA", paste0(ds, "_significant_edges.csv")),
    DeepRIG = file.path("DeepRIG", paste0(ds, "_significant_edges.csv")),
    IGEGRNs = file.path("IGEGRNs", paste0(ds, "_significant_edges.csv"))
  )
  missing <- base_files[!sapply(base_files, file.exists)]
  if (length(missing) > 0) {
    warning(paste0("Dataset ", ds, " missing files: ", paste(names(missing), collapse = ", ")))
    next
  }
  
  genie3     <- clean_read(base_files$GENIE3, delim = "\t")
  kboost     <- clean_read(base_files$Kboost, delim = ",")
  grnboost   <- clean_read(base_files$GRNBoost, delim = ",")
  dcema3     <- clean_read(base_files$`3DCEMA`, delim = ",")
  deeprig    <- clean_read(base_files$DeepRIG, delim = ",")
  igegrn     <- clean_read(base_files$IGEGRNs, delim = ",")
  
  all_methods <- list(GENIE3 = genie3, Kboost = kboost, GRNBoost = grnboost, 
                      `3DCEMA` = dcema3, DeepRIG = deeprig, IGEGRNs = igegrn)
  method_names <- names(all_methods)
  
  all_base_wide <- build_wide_table(all_methods, method_names)
  rm(genie3, kboost, grnboost, dcema3, deeprig, igegrn)
  gc()
  
  total_pos <- nrow(inner_join(chip_gt_pos, all_base_wide[,c("TF","Target")], by=c("TF","Target")))
  if(total_pos < 20) {
    warning(sprintf("Dataset %s has only %d positive edges; skipping", ds, total_pos))
    next
  }
  
  # Directories
  out_root <- file.path(root_dir, ds)
  dir_pair    <- file.path(out_root, "1_Pair_2Methods")
  dir_triple  <- file.path(out_root, "2_Triple_3Methods")
  dir_4way    <- file.path(out_root, "3_Four_4Methods")
  dir_5way    <- file.path(out_root, "4_Five_5Methods")
  dir_all6    <- file.path(out_root, "5_All_6Methods")
  all_dirs <- c(out_root, dir_pair, dir_triple, dir_4way, dir_5way, dir_all6)
  for (d in all_dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  
  process_combinations <- function(comb_list, output_dir) {
    batch_metrics <- list()
    idx <- 1
    for (nm in comb_list) {
      nm_sorted <- sort(nm)
      comb_name <- paste(nm_sorted, collapse = "_")
      cat("  Processing combination: ", comb_name, "\n")
      run_res <- safe_run_xgb(all_methods[nm], all_base_wide, chip_gt_pos, comb_name, 
                              sup_neg_ratio, normalize_weights, test_ratio, ds)
      if (nrow(run_res$pred_save) > 0) {
        out_file <- file.path(output_dir, paste0(ds, "_XGBFusion_", comb_name, ".csv"))
        safe_write_csv(run_res$pred_save, out_file)
        cat("    [OK] Saved top-10% high-scoring edges (", nrow(run_res$pred_save), ")\n")
      }
      if (write_full_all_edges && nrow(run_res$pred_full) > 0) {
        full_file <- file.path(output_dir, paste0(ds, "_XGBFusion_FULL_", comb_name, ".csv"))
        safe_write_csv(run_res$pred_full, full_file)
      }
      batch_metrics[[idx]] <- run_res$out_df
      idx <- idx + 1
    }
    bind_rows(batch_metrics)
  }
  
  ds_metrics <- data.frame()
  if (only_all6) {
    ds_metrics <- process_combinations(list(method_names), dir_all6)
  } else {
    m2 <- process_combinations(combn(method_names, 2, simplify = FALSE), dir_pair)
    m3 <- process_combinations(combn(method_names, 3, simplify = FALSE), dir_triple)
    m4 <- process_combinations(combn(method_names, 4, simplify = FALSE), dir_4way)
    m5 <- process_combinations(combn(method_names, 5, simplify = FALSE), dir_5way)
    m6 <- process_combinations(list(method_names), dir_all6)
    ds_metrics <- bind_rows(m2, m3, m4, m5, m6)
  }
  
  eval_summary <- bind_rows(eval_summary, ds_metrics)
  cat("\n[Done]", ds, "dataset computed; output path: ", out_root, "\n")
  gc()
}

# Output the summary metrics
if (nrow(eval_summary) > 0) {
  eval_summary <- eval_summary %>%
    select(Dataset, Combination, NumMethods, Methods, Train_Pos, Train_Neg, 
           Test_Pos, Test_Neg, Best_Nrounds, AUC_ROC, AUC_PR) %>%
    arrange(Dataset, desc(AUC_PR), NumMethods)
  sum_file <- file.path(root_dir, "XGB_Fusion_Evaluation_Summary.csv")
  safe_write_csv(eval_summary, sum_file)
  cat("\n[Finished] All datasets completed! Summary metrics saved\n")
  cat("Total combinations:", nrow(eval_summary), "\n")
  print(table(eval_summary$Dataset))
  print(table(eval_summary$NumMethods))
} else {
  cat("\nNo valid run combinations\n")
}
