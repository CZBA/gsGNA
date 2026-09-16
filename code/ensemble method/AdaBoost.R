# ============================================================================
# AdaBoost-based GRN Feature Fusion Pipeline for Rice Stress Response
# ============================================================================
# Load required packages
library(dplyr)
library(tidyr)
library(vroom)
library(ada)
library(caret)
library(pROC)
library(PRROC)
library(data.table)
library(rpart)

# Install missing dependency digest automatically
if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest")
}
library(digest)

# ====================== Global Parameters ======================
set.seed(123)
sup_neg_ratio    <- 1
normalize_weights <- TRUE
test_ratio       <- 0.2
save_all_edges   <- TRUE   # Keep all edges without quantile filtering

# AdaBoost hyperparameters
ada_maxiter      <- 150
ada_nu           <- 0.5
ada_type         <- "real"
ada_maxdepth     <- 2
ada_minsplit     <- 100
ada_cp           <- 0.05

# Performance & memory settings
batch_size       <- 20
cv_folds         <- 3
cv_iter_step     <- 20
min_samples_cv   <- 1000
memory_limit_gb  <- 8
chip_file_path   <- "chip_seq.csv"
root_dir         <- "FeatureFusion_AdaBoost_Optimized"
all_method_names <- c("GENIE3", "Kboost", "GRNBoost", "3DCEMA", "DeepRIG", "IGEGRNs")

# ====================== Memory Management ======================
#' Check memory usage and trigger garbage collection if exceeding threshold
check_memory <- function() {
  if (requireNamespace("pryr", quietly = TRUE)) {
    mem_used <- pryr::mem_used()
    mem_gb <- as.numeric(mem_used) / 1e9
    if (mem_gb > memory_limit_gb * 0.8) {
      cat("Memory usage reaches", round(mem_gb, 2), "GB, running garbage collection\n")
      gc()
    }
  } else {
    gc()
  }
}

# ====================== Utility Functions ======================
#' Safe CSV writing with try-catch error handling
#' @param df Input data frame
#' @param filepath Output file path
safe_write_csv <- function(df, filepath) {
  if (nrow(df) == 0) return(invisible(NULL))
  tryCatch({
    vroom_write(df, filepath, delim = ",", progress = FALSE)
  }, error = function(e) {
    cat("File writing failed: ", filepath, "Error: ", conditionMessage(e), "\n")
  })
}

#' Read GRN edge table, retain TF-Target-Weight and deduplicate by max weight
#' @param path File path of GRN edge file
#' @param delim Delimiter character
clean_read <- function(path, delim = ",") {
  df <- fread(path, select = 1:3, data.table = FALSE)
  colnames(df) <- c("TF", "Target", "Weight")
  df <- df %>%
    filter(!is.na(TF), !is.na(Target), is.finite(Weight), TF != Target) %>%
    group_by(TF, Target) %>%
    summarise(Weight = max(Weight), .groups = "drop")
  gc()
  return(df)
}

#' Min-Max normalization, use statistics only from training set to avoid data leakage
#' @param train_mat Training feature matrix
#' @param full_mat Full feature matrix to scale
normalize_by_train_stats <- function(train_mat, full_mat) {
  col_min <- apply(train_mat, 2, min, na.rm = TRUE)
  col_max <- apply(train_mat, 2, max, na.rm = TRUE)
  range_val <- col_max - col_min
  range_val[range_val == 0] <- 1
  
  full_scaled <- sweep(full_mat, 2, col_min, "-")
  full_scaled <- sweep(full_scaled, 2, range_val, "/")
  full_scaled[is.na(full_scaled)] <- 0
  return(full_scaled)
}

#' Read ChIP-seq gold-standard positive TF-target interactions
#' @param chip_path Path to ChIP-seq ground truth file
read_chip_pos <- function(chip_path) {
  chip_df <- fread(chip_path, select = 1:2, data.table = FALSE)
  colnames(chip_df) <- c("TF", "Target")
  chip_df <- chip_df %>%
    filter(!is.na(TF), !is.na(Target), TF != Target) %>%
    distinct(TF, Target)
  gc()
  return(chip_df)
}

#' Calculate AUC-ROC and AUC-PR for prediction evaluation
#' @param test_df Test dataset with true labels and prediction scores
#' @param score_col Column name of prediction score
calc_metric <- function(test_df, score_col = "AdaBoost_Score") {
  test_df <- test_df %>% filter(!is.na(!!sym(score_col)))
  lab <- test_df$label
  scr <- test_df[[score_col]]
  
  if (length(unique(lab)) < 2 || sum(lab == 1) == 0 || sum(lab == 0) == 0) {
    return(c(AUC_ROC = NA, AUC_PR = NA))
  }
  
  roc_obj <- roc(lab, scr, quiet = TRUE, direction = "<")
  auc_roc <- round(as.numeric(auc(roc_obj)), 4)
  
  pr_obj <- tryCatch({
    pr.curve(scores.class0 = scr[lab == 1], scores.class1 = scr[lab == 0], curve = FALSE)
  }, error = function(e) NULL)
  
  auc_pr <- ifelse(is.null(pr_obj), NA, round(pr_obj$auc.integral, 4))
  return(c(AUC_ROC = auc_roc, AUC_PR = auc_pr))
}

# ====================== Feature Construction ======================
#' Build wide-format feature table for selected model combination
#' @param method_subset List of selected base GRN models
#' @param all_methods_list Full list of base GRN results
build_subset_wide <- function(method_subset, all_methods_list) {
  sub_methods <- names(method_subset)
  all_edges <- NULL
  for (m in sub_methods) {
    df_temp <- all_methods_list[[m]] %>%
      select(TF, Target, Weight) %>%
      rename(!!paste0(m, "_W") := Weight) %>%
      mutate(!!paste0(m, "_Exist") := 1)
    
    if (is.null(all_edges)) {
      all_edges <- df_temp
    } else {
      all_edges <- full_join(all_edges, df_temp, by = c("TF", "Target"))
    }
  }
  
  for (m in sub_methods) {
    w_col <- paste0(m, "_W")
    e_col <- paste0(m, "_Exist")
    all_edges[[w_col]][is.na(all_edges[[w_col]])] <- 0
    all_edges[[e_col]][is.na(all_edges[[e_col]])] <- 0
  }
  
  feat_names <- c(paste0(sub_methods, "_W"), paste0(sub_methods, "_Exist"))
  return(list(wide = all_edges, feat_cols = feat_names))
}

#' Safe stratified sampling with fixed seed
#' @param df Input data frame
#' @param n Sample size
#' @param seed Random seed
efficient_sample <- function(df, n, seed) {
  if (nrow(df) <= n) return(df)
  set.seed(seed)
  df[sample(nrow(df), n), ]
}

#' Split train/test by both TF space and Target space for strict out-of-distribution evaluation
#' @param pos_df Positive interaction table
#' @param test_frac Fraction for test set
split_tf_target_two_space <- function(pos_df, test_frac = 0.2) {
  all_tfs <- unique(pos_df$TF)
  all_tars <- unique(pos_df$Target)
  
  test_tfs <- sample(all_tfs, size = max(1, floor(length(all_tfs) * test_frac)))
  test_tars <- sample(all_tars, size = max(1, floor(length(all_tars) * test_frac)))
  
  train_tfs <- setdiff(all_tfs, test_tfs)
  train_tars <- setdiff(all_tars, test_tars)
  
  pos_train <- pos_df %>% filter(TF %in% train_tfs & Target %in% train_tars)
  pos_test  <- pos_df %>% filter(TF %in% test_tfs & Target %in% test_tars)
  
  # Fallback: split only by TF if double-space split fails
  if(nrow(pos_train) < 1 || nrow(pos_test) < 1){
    pos_split_old <- split_by_gene_single_tf(pos_df, test_frac)
    return(pos_split_old)
  }
  return(list(train_pos = pos_train, test_pos = pos_test,
              train_tfs = train_tfs, train_tars = train_tars,
              test_tfs = test_tfs, test_tars = test_tars))
}

#' Fallback split: split dataset only on TF IDs
#' @param pos_df Positive interaction table
#' @param test_frac Fraction for test set
split_by_gene_single_tf <- function(pos_df, test_frac=0.2){
  all_tfs <- unique(pos_df$TF)
  test_tfs <- sample(all_tfs, max(1, floor(length(all_tfs)*test_frac)))
  train_tfs <- setdiff(all_tfs, test_tfs)
  list(train_pos = pos_df %>% filter(TF %in% train_tfs),
       test_pos = pos_df %>% filter(TF %in% test_tfs))
}

# ====================== CV Hyperparameter Optimization (Maximize AUC-PR) ======================
#' Cross-validation to select optimal AdaBoost iteration number
#' @param train_df Training dataset
#' @param max_iter Maximum candidate iterations
#' @param cv_folds Number of CV folds
#' @param iter_step Step size for iteration grid search
cv_optimize_ada <- function(train_df, max_iter, cv_folds = 3, iter_step = 20) {
  set.seed(123)
  if (nrow(train_df) < min_samples_cv) return(max_iter)
  
  folds <- createFolds(train_df$label, k = cv_folds, list = TRUE)
  iter_grid <- seq(10, min(max_iter, 100), by = iter_step)
  if (length(iter_grid) == 0) return(max_iter)
  
  best_iter <- max_iter
  best_score <- -Inf
  
  for (iter in iter_grid) {
    cv_scores <- numeric(cv_folds)
    for (fold in 1:cv_folds) {
      val_idx <- folds[[fold]]
      tr_cv <- train_df[-val_idx, ]
      va_cv <- train_df[val_idx, ]
      
      if (sum(tr_cv$label == 1) < 2 || sum(tr_cv$label == 0) < 2) {
        cv_scores[fold] <- NA
        next
      }
      tryCatch({
        ada_cv <- ada(label ~ ., data = tr_cv, iter = iter, nu = ada_nu, type = ada_type,
                      control = rpart.control(maxdepth = ada_maxdepth, cp = ada_cp, minsplit = ada_minsplit))
        prob_cv <- predict(ada_cv, va_cv[, -ncol(va_cv)], type = "prob")[, 2]
        pr_cv <- pr.curve(scores.class0 = prob_cv[va_cv$label == 1], scores.class1 = prob_cv[va_cv$label == 0], curve = FALSE)
        cv_scores[fold] <- pr_cv$auc.integral
      }, error = function(e) cv_scores[fold] <- NA)
    }
    mean_score <- mean(cv_scores, na.rm = TRUE)
    if (!is.na(mean_score) && mean_score > best_score) {
      best_score <- mean_score
      best_iter <- iter
    }
  }
  return(best_iter)
}

# ====================== Core Training Function ======================
#' Train AdaBoost fusion model for one combination of base GRN methods
#' @param method_subset Selected subset of base GRN methods
#' @param all_methods_list List containing all base GRN outputs
#' @param chip_pos ChIP-seq gold-standard positive edges
#' @param ds_name Dataset name
#' @param comb_name Name for current method combination
#' @param neg_ratio Negative / positive sampling ratio
#' @param normalize_w Whether to perform min-max normalization
#' @param test_frac Fraction of test set
pure_ada_grn_optimized <- function(method_subset, all_methods_list, chip_pos, ds_name, comb_name,
                                   neg_ratio = 1, normalize_w = TRUE, test_frac = 0.2) {
  sub_res <- build_subset_wide(method_subset, all_methods_list)
  wide_sub <- sub_res$wide
  feat_cols <- sub_res$feat_cols
  
  # Generate unique seed for each combination via hash digest
  hash_raw <- digest::digest(paste0(ds_name, comb_name), serialize = FALSE)
  hash_num <- sum(as.integer(charToRaw(hash_raw)))
  comb_seed <- abs(hash_num) %% 999999
  comb_seed <- max(1, comb_seed)
  
  all_edge_key <- wide_sub %>% select(TF, Target)
  all_pos_edge <- inner_join(chip_pos, all_edge_key, by = c("TF", "Target"))
  n_pos_total <- nrow(all_pos_edge)
  
  if (n_pos_total == 0) {
    return(list(result_full = data.frame(), result_save = data.frame(),
                metrics = c(AUC_ROC = NA, AUC_PR = NA),
                sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, 
                                best_iter=0, real_neg_ratio=NA, 
                                total_edges=0, keep_edges=0)))
  }
  
  all_neg_candidate <- anti_join(all_edge_key, chip_pos, by = c("TF", "Target"))
  if (nrow(all_neg_candidate) == 0) {
    return(list(result_full = data.frame(), result_save = data.frame(),
                metrics = c(AUC_ROC = NA, AUC_PR = NA),
                sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, 
                                best_iter=0, real_neg_ratio=NA, 
                                total_edges=0, keep_edges=0)))
  }
  
  # Two-space split for TF and Target
  pos_split <- split_tf_target_two_space(all_pos_edge, test_frac)
  pos_train <- pos_split$train_pos
  pos_test  <- pos_split$test_pos
  n_pos_train <- nrow(pos_train)
  n_pos_test <- nrow(pos_test)
  
  if (n_pos_train < 2) {
    return(list(result_full = data.frame(), result_save = data.frame(),
                metrics = c(AUC_ROC = NA, AUC_PR = NA),
                sample_info = c(tr_pos=n_pos_train, tr_neg=0, te_pos=n_pos_test, 
                                te_neg=0, best_iter=0, real_neg_ratio=NA, 
                                total_edges=0, keep_edges=0)))
  }
  
  # Negative sampling
  n_neg_train <- min(floor(n_pos_train * neg_ratio), nrow(all_neg_candidate))
  neg_train <- efficient_sample(all_neg_candidate, n_neg_train, seed = comb_seed)
  real_neg_ratio <- nrow(neg_train) / n_pos_train
  
  train_set <- bind_rows(
    pos_train %>% mutate(label = 1),
    neg_train %>% mutate(label = 0)
  )
  
  if (length(unique(train_set$label)) < 2) {
    return(list(result_full = data.frame(), result_save = data.frame(),
                metrics = c(AUC_ROC = NA, AUC_PR = NA),
                sample_info = c(tr_pos=n_pos_train, tr_neg=nrow(neg_train), 
                                te_pos=n_pos_test, te_neg=nrow(all_neg_candidate),
                                best_iter=0, real_neg_ratio=real_neg_ratio, 
                                total_edges=0, keep_edges=0)))
  }
  
  pred_set <- wide_sub %>% select(TF, Target) %>% distinct()
  train_feat <- left_join(train_set, wide_sub, by = c("TF", "Target"))
  pred_feat <- left_join(pred_set, wide_sub, by = c("TF", "Target"))
  
  x_train_raw <- as.matrix(train_feat[, feat_cols, drop = FALSE])
  x_pred_raw <- as.matrix(pred_feat[, feat_cols, drop = FALSE])
  
  # Filter zero-variance features
  feat_var <- apply(x_train_raw, 2, var, na.rm = TRUE)
  keep_feat <- feat_var > 1e-8
  if (sum(keep_feat) < length(feat_cols)) {
    feat_cols <- feat_cols[keep_feat]
    x_train_raw <- x_train_raw[, keep_feat, drop = FALSE]
    x_pred_raw <- x_pred_raw[, keep_feat, drop = FALSE]
  }
  
  if (ncol(x_train_raw) == 0) {
    return(list(result_full = data.frame(), result_save = data.frame(),
                metrics = c(AUC_ROC = NA, AUC_PR = NA),
                sample_info = c(tr_pos=n_pos_train, tr_neg=nrow(neg_train), 
                                te_pos=n_pos_test, te_neg=nrow(all_neg_candidate),
                                best_iter=0, real_neg_ratio=real_neg_ratio, 
                                total_edges=0, keep_edges=0)))
  }
  
  # Feature normalization
  if (normalize_w) {
    x_train <- normalize_by_train_stats(x_train_raw, x_train_raw)
    x_pred <- normalize_by_train_stats(x_train_raw, x_pred_raw)
  } else {
    x_train <- x_train_raw
    x_pred <- x_pred_raw
  }
  
  train_df <- as.data.frame(x_train)
  train_df$label <- as.factor(train_feat$label)
  best_iter <- cv_optimize_ada(train_df, ada_maxiter, cv_folds, cv_iter_step)
  
  # Train final AdaBoost model
  set.seed(comb_seed)
  fit <- ada(label ~ ., data = train_df, iter = best_iter, 
             nu = ada_nu, type = ada_type,
             control = rpart.control(maxdepth = ada_maxdepth, cp = ada_cp, minsplit = ada_minsplit))
  
  # Predict scores for all candidate edges
  pred_df <- as.data.frame(x_pred)
  pred_probs <- predict(fit, newdata = pred_df, type = "prob")[, 2]
  pred_set$AdaBoost_Score <- round(pred_probs, 4)
  result_full <- pred_set %>% select(TF, Target, AdaBoost_Score) %>% arrange(desc(AdaBoost_Score))
  
  # Model evaluation on test set
  test_set <- bind_rows(pos_test %>% mutate(label = 1), all_neg_candidate %>% mutate(label = 0))
  test_feat <- inner_join(test_set, result_full, by = c("TF", "Target"))
  metrics <- calc_metric(test_feat, "AdaBoost_Score")
  
  # Save all edges (quantile filtering disabled by save_all_edges flag)
  total_all <- nrow(result_full)
  if(save_all_edges){
    result_save <- result_full
  }else{
    threshold <- quantile(result_full$AdaBoost_Score, probs = 0.9)
    result_save <- result_full %>% filter(AdaBoost_Score >= threshold)
  }
  keep_num <- nrow(result_save)
  
  cat(sprintf("    [%s] Train: %d pos / %.1f neg | Test: %d pos / %d neg | Iter: %d | Seed: %d\n",
              paste(names(method_subset), collapse=","), 
              n_pos_train, real_neg_ratio, n_pos_test, nrow(all_neg_candidate), best_iter, comb_seed))
  cat(sprintf("    AUC-ROC=%.4f, AUC-PR=%.4f | SaveAllEdges=%s | Keep: %d/%d\n",
              metrics["AUC_ROC"], metrics["AUC_PR"], save_all_edges, keep_num, total_all))
  
  gc()
  return(list(result_full = result_full, result_save = result_save, metrics = metrics,
              sample_info = c(tr_pos=n_pos_train, tr_neg=nrow(neg_train), 
                              te_pos=n_pos_test, te_neg=nrow(all_neg_candidate),
                              best_iter = best_iter, 
                              real_neg_ratio = round(real_neg_ratio, 3),
                              total_edges = total_all, keep_edges = keep_num)))
}

# ====================== Wrapper with Exception Handling ======================
#' Wrapper function to catch errors for each model combination
#' @param method_subset Subset of base GRN methods
#' @param all_methods_list Full list of base GRN results
#' @param comb_name Combination identifier
#' @param ds_name Dataset name
#' @param neg_ratio Negative sampling ratio
#' @param normalize_w Whether normalize features
#' @param test_frac Test proportion
safe_run_ada <- function(method_subset, all_methods_list, chip_pos, comb_name, ds_name,
                         neg_ratio, normalize_w, test_frac) {
  res <- tryCatch({
    pure_ada_grn_optimized(method_subset, all_methods_list, chip_pos, comb_name, ds_name,
                           neg_ratio = neg_ratio, normalize_w = normalize_w, 
                           test_frac = test_frac)
  }, error = function(e) {
    cat("  Combination", comb_name, "error:", conditionMessage(e), "\n")
    return(list(result_full = data.frame(), result_save = data.frame(),
                metrics = c(AUC_ROC = NA, AUC_PR = NA),
                sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0, 
                                best_iter=0, real_neg_ratio=NA, 
                                total_edges=0, keep_edges=0)))
  })
  
  info <- res$sample_info
  m_val <- res$metrics
  
  out_df <- data.frame(
    Dataset = ds_name,
    Combination = comb_name,
    NumMethods = length(method_subset),
    Methods = paste(sort(names(method_subset)), collapse = ";"),
    Train_Pos = info["tr_pos"],
    Train_Neg = info["tr_neg"],
    Test_Pos = info["te_pos"],
    Test_All_Neg = info["te_neg"],
    Real_Neg_Ratio = info["real_neg_ratio"],
    Best_Iter = info["best_iter"],
    Total_All_Edges = info["total_edges"],
    Filter_Keep_Edges = info["keep_edges"],
    AUC_ROC = m_val["AUC_ROC"],
    AUC_PR = m_val["AUC_PR"],
    stringsAsFactors = FALSE
  )
  return(list(out_df = out_df, pred_save = res$pred_save, pred_full = res$result_full))
}

# ====================== Serial Batch Processing ======================
#' Serial loop over all method combinations
#' @param comb_list List of model combinations
#' @param output_dir Output directory
#' @param all_methods_list All base GRN results
#' @param chip_gt_pos ChIP-seq gold standard
#' @param ds Dataset name
process_combinations_serial <- function(comb_list, output_dir, all_methods_list, chip_gt_pos, ds) {
  batch_results <- data.frame()
  for(nm in comb_list){
    gc()
    nm_sorted <- sort(nm)
    comb_name <- paste(nm_sorted, collapse = "_")
    out_file <- file.path(output_dir, paste0(ds, "_AdaBoostFusion_", comb_name, ".csv"))
    # Skip if output file already exists
    if(file.exists(out_file)){
      cat("  Skip existing combination:", comb_name, "\n")
      next
    }
    cat("  Processing:", comb_name, "\n")
    run_res <- safe_run_ada(all_methods_list[nm], all_methods_list, chip_gt_pos, comb_name, ds,
                            sup_neg_ratio, normalize_weights, test_ratio)
    if (nrow(run_res$pred_save) > 0) {
      safe_write_csv(run_res$pred_save, out_file)
      cat("  Saved:", comb_name, "(", nrow(run_res$pred_save), "edges)\n")
    }
    batch_results <- bind_rows(batch_results, run_res$out_df)
  }
  return(batch_results)
}

# ====================== Main Program ======================
cat("\n=============================================\n")
cat("AdaBoost GRN Fusion Serial Pipeline (No checkpoint, retain all edges) Start\n")
cat("=============================================\n")
datasets <- c("Aridity", "Alkalinity", "Cold")
chip_gt_pos <- read_chip_pos(chip_file_path)
eval_summary <- data.frame()
total_start_time <- Sys.time()

for (ds in datasets) {
  cat("\n=============================================\n")
  cat("Processing Dataset:", ds, "\n")
  cat("=============================================\n")
  set.seed(123)
  
  # File paths for base GRN outputs
  base_files <- list(
    GENIE3 = file.path("GENIE3", paste0(ds, "_significant_edges.txt")),
    Kboost = file.path("Kboost", paste0(ds, "_significant_edges.csv")),
    GRNBoost = file.path("GRNBoost", paste0(ds, "_significant_edges.csv")),
    `3DCEMA` = file.path("3DCEMA", paste0(ds, "_significant_edges.csv")),
    DeepRIG = file.path("DeepRIG", paste0(ds, "_significant_edges.csv")),
    IGEGRNs = file.path("IGEGRNs", paste0(ds, "_significant_edges.csv"))
  )
  
  # Check missing input files
  missing <- base_files[!sapply(base_files, file.exists)]
  if (length(missing) > 0) {
    warning(paste0("Dataset ", ds, " missing files: ", paste(names(missing), collapse = ", ")))
    next
  }
  
  cat("Reading base GRN datasets...\n")
  all_methods <- list(
    GENIE3 = clean_read(base_files$GENIE3, delim = "\t"),
    Kboost = clean_read(base_files$Kboost, delim = ","),
    GRNBoost = clean_read(base_files$GRNBoost, delim = ","),
    `3DCEMA` = clean_read(base_files$`3DCEMA`, delim = ","),
    DeepRIG = clean_read(base_files$DeepRIG, delim = ","),
    IGEGRNs = clean_read(base_files$IGEGRNs, delim = ",")
  )
  
  # Create hierarchical output directories for 2~6 model combinations
  out_root <- file.path(root_dir, ds)
  dir_pair <- file.path(out_root, "1_Pair_2Methods")
  dir_triple <- file.path(out_root, "2_Triple_3Methods")
  dir_4way <- file.path(out_root, "3_Four_4Methods")
  dir_5way <- file.path(out_root, "4_Five_5Methods")
  dir_all6 <- file.path(out_root, "5_All_6Methods")
  dir_map <- list("2" = dir_pair, "3" = dir_triple, "4" = dir_4way, "5" = dir_5way, "6" = dir_all6)
  all_dirs <- c(out_root, unlist(dir_map))
  for (d in all_dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  
  # Generate all combinations (2/3/4/5/6 base models)
  comb_2 <- combn(all_method_names, 2, simplify = FALSE)
  comb_3 <- combn(all_method_names, 3, simplify = FALSE)
  comb_4 <- combn(all_method_names, 4, simplify = FALSE)
  comb_5 <- combn(all_method_names, 5, simplify = FALSE)
  comb_6 <- list(all_method_names)
  comb_list <- c(comb_2, comb_3, comb_4, comb_5, comb_6)
  
  cat("Total combinations for current dataset:", length(comb_list), "\n")
  cat("No checkpoint RDS, iterate all combinations, skip existing CSV only\n")
  cat("Total combinations to run:", length(comb_list), "\n")
  
  # Process groups by combination size
  for (size_str in names(dir_map)) {
    size <- as.integer(size_str)
    size_comb <- comb_list[sapply(comb_list, length) == size]
    if (length(size_comb) == 0) next
    cat("\nProcessing", size, "-model fusion combinations, total", length(size_comb), "groups\n")
    n_batches <- ceiling(length(size_comb) / batch_size)
    
    for (batch in 1:n_batches) {
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, length(size_comb))
      batch_comb <- size_comb[start_idx:end_idx]
      cat("  Batch", batch, "/", n_batches, "running", length(batch_comb), "combinations\n")
      batch_start <- Sys.time()
      
      batch_results <- process_combinations_serial(batch_comb, dir_map[[size_str]], all_methods, chip_gt_pos, ds)
      
      if (nrow(batch_results) > 0) {
        eval_summary <- bind_rows(eval_summary, batch_results)
      }
      
      batch_end <- Sys.time()
      cat("  Batch finished, time cost:", round(difftime(batch_end, batch_start, units="mins"), 2), "minutes\n")
      gc()
      check_memory()
    }
  }
  
  # Release large objects
  rm(all_methods, comb_list)
  gc()
  cat("\n✅ Dataset [", ds, "] all combinations completed\n")
}
total_end_time <- Sys.time()
cat("\n=============================================\n")
cat("All stress datasets finished! Total runtime:", round(difftime(total_end_time, total_start_time, units="hours"), 2), "hours\n")
cat("=============================================\n")

# Export global evaluation summary table
if (nrow(eval_summary) > 0) {
  eval_summary <- eval_summary %>% arrange(Dataset, desc(AUC_PR), NumMethods)
  sum_file <- file.path(root_dir, "AdaBoost_Fusion_Summary.csv")
  safe_write_csv(eval_summary, sum_file)
  
  cat("\nSummary statistics:\n")
  cat("Total fusion combinations: ", nrow(eval_summary), "\n")
  cat("Combination count per stress dataset:\n")
  print(table(eval_summary$Dataset))
  cat("Combination count by number of base models:\n")
  print(table(eval_summary$NumMethods))
  
  # Best combination per dataset (ranked by AUC-PR)
  best_comb <- eval_summary %>%
    group_by(Dataset) %>%
    slice_max(AUC_PR, n = 1, with_ties = FALSE) %>%
    select(Dataset, Combination, NumMethods, AUC_ROC, AUC_PR)
  cat("\n🏆 Best fusion combination for each dataset:\n")
  print(best_comb, row.names = FALSE)
  
  # Global best combination
  overall_best <- eval_summary %>%
    slice_max(AUC_PR, n = 1) %>%
    select(Dataset, Combination, NumMethods, AUC_ROC, AUC_PR)
  cat("\n🌟 Global optimal combination:\n")
  print(overall_best, row.names = FALSE)
} else {
  cat("\n⚠️ No valid evaluation results. Please check input GRN files and ChIP file\n")
}
cat("\n✅ Pipeline completed normally\n")
