library(dplyr)
library(tidyr)
library(vroom)
library(glmnet)
library(readr)
library(matrixStats)
library(pROC)
library(PRROC)
library(progress)  # progress bar

# ============================== Global fixed parameters & path configuration ==============================
set.seed(123)                      # Global base seed
sup_neg_ratio    <- 1              # Positive:negative sample ratio in training set (1:1)
only_all6        <- FALSE
normalize_weights <- TRUE
test_ratio       <- 0.2            # Test set split ratio (based on positive samples)
filter_zero_neg  <- TRUE           # Whether to filter negative samples with all-model weights equal to 0
chip_file_path   <- "chip_seq.csv" # Path to the gold-standard ChIP file
root_dir         <- "Stacking_Supervised_ChIP"
method_folders   <- c("GENIE3", "Kboost", "GRNBoost", "3DCEMA", "DeepRIG", "IGEGRNs")
# =====================================================================================

# Read and normalize a single-model GRN file
clean_read <- function(path, delim = ",") {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  df <- vroom(path, delim = delim, show_col_types = FALSE)
  df <- df[, 1:3]
  colnames(df) <- c("TF", "Target", "Weight")
  df <- df %>% drop_na(TF, Target, Weight)
  df <- df %>% filter(is.finite(Weight))
  df$TF <- toupper(df$TF)
  df$Target <- toupper(df$Target)
  df <- df %>% filter(TF != Target)
  df <- df %>% 
    group_by(TF, Target) %>% 
    summarise(Weight = max(Weight), .groups = "drop")
  return(df)
}

# Fill NA weights in the wide table with 0
fill_na_zero <- function(df, cols) {
  df <- df %>%
    mutate(across(all_of(cols), ~ ifelse(is.na(.), 0, .)))
  return(df)
}

# 0-1 min-max normalization (optimized)
normalize_by_train_stats <- function(train_mat, full_mat) {
  col_min <- colMins(train_mat, na.rm = TRUE)
  col_max <- colMaxs(train_mat, na.rm = TRUE)
  range_val <- col_max - col_min
  range_val[range_val == 0] <- 1
  
  # Optimization: use matrix operations instead of sweep
  full_scaled <- t((t(full_mat) - col_min) / range_val)
  full_scaled[is.na(full_scaled)] <- 0
  full_scaled[!is.finite(full_scaled)] <- 0
  return(full_scaled)
}

# Read the ChIP-seq gold standard (positive edges only)
read_chip_pos <- function(chip_path) {
  if (!file.exists(chip_path)) stop(paste("ChIP file missing:", chip_path))
  chip_df <- vroom(chip_path, show_col_types = FALSE)
  chip_df <- chip_df[, 1:2]
  colnames(chip_df) <- c("TF", "Target")
  chip_df <- chip_df %>% drop_na(TF, Target)
  chip_df$TF <- toupper(chip_df$TF)
  chip_df$Target <- toupper(chip_df$Target)
  chip_df <- chip_df %>% filter(TF != Target) %>% distinct(TF, Target)
  return(chip_df)
}

# Generate stratified CV fold ids (fixed version)
get_stratified_foldid <- function(y, nfolds) {
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)
  
  # Ensure each class has at least as many samples as folds
  if (length(pos_idx) < nfolds || length(neg_idx) < nfolds) {
    stop("The number of samples in one class is less than the number of folds")
  }
  
  # Randomly shuffle and assign folds
  fold_pos <- sample(rep(1:nfolds, length.out = length(pos_idx)))
  fold_neg <- sample(rep(1:nfolds, length.out = length(neg_idx)))
  
  foldid <- integer(length(y))
  foldid[pos_idx] <- fold_pos
  foldid[neg_idx] <- fold_neg
  return(foldid)
}

# ==================== Core function: fixed version ====================
stacking_supervised_chip <- function(method_list, all_base_wide, chip_pos,
                                     neg_ratio = 1, normalize_w = TRUE,
                                     test_frac = 0.2, local_seed = 999, 
                                     filter_zero_neg = TRUE, 
                                     zero_var_threshold = 1e-12) {
  set.seed(local_seed)
  method_names <- names(method_list)
  if (length(method_names) == 0) stop("The model combination is empty")
  
  feat_cols <- method_names
  wide_subset <- all_base_wide %>%
    select(TF, Target, all_of(feat_cols)) %>%
    distinct(TF, Target, .keep_all = TRUE)
  wide_subset <- fill_na_zero(wide_subset, feat_cols)
  
  all_edge_key <- wide_subset %>% select(TF, Target)
  all_pos_edge <- inner_join(chip_pos, all_edge_key, by = c("TF", "Target"))
  n_pos_total <- nrow(all_pos_edge)
  
  if (n_pos_total == 0) {
    warning("No matching ChIP positive edges for the current model combination; skipping")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=NA)
    ))
  }
  
  # Skip if there are too few positive samples
  if (n_pos_total < 3) {
    warning(sprintf("Too few positive samples (%d), cannot train effectively", n_pos_total))
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=NA)
    ))
  }
  
  all_neg_candidate <- anti_join(all_edge_key, chip_pos, by = c("TF", "Target"))
  
  if (filter_zero_neg) {
    neg_with_feat <- all_neg_candidate %>% 
      left_join(wide_subset, by = c("TF", "Target"))
    
    # Compute the sum of weights using a more reasonable threshold
    sum_w <- rowSums(neg_with_feat[, feat_cols, drop = FALSE])
    # Adjust the threshold dynamically based on the data distribution
    weight_quantile <- quantile(sum_w, probs = 0.1, na.rm = TRUE)
    threshold <- max(1e-8, weight_quantile * 0.01)  # Use 1% of the 10% quantile as threshold
    
    keep_idx <- sum_w > threshold
    cat(sprintf("    Filtering all-zero negative samples: %d / %d kept (threshold=%.2e)\n", 
                sum(keep_idx), nrow(all_neg_candidate), threshold))
    all_neg_candidate <- neg_with_feat[keep_idx, c("TF","Target")]
  }
  
  if (nrow(all_neg_candidate) == 0) {
    warning("No usable negative samples after filtering; skipping")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=NA)
    ))
  }
  
  # ========== 1. Fixed split of positive/negative samples, permanent train/test isolation ==========
  n_pos_test <- max(1, floor(n_pos_total * test_frac))
  # Ensure the training set has at least 3 positive samples
  if (n_pos_total - n_pos_test < 3) {
    n_pos_test <- max(1, n_pos_total - 3)
  }
  n_pos_train <- n_pos_total - n_pos_test
  
  if (n_pos_train < 1) {
    warning("Too few positive samples to split a test set; skipping")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=NA)
    ))
  }
  
  pos_idx <- sample(1:n_pos_total, n_pos_test)
  pos_test  <- all_pos_edge[pos_idx, ]
  pos_train <- all_pos_edge[-pos_idx, ]
  
  n_neg_total <- nrow(all_neg_candidate)
  # Test-set negative samples: use all available negatives or the same count as positives (whichever is smaller)
  n_neg_test <- min(n_pos_test * 2, n_neg_total)  # At least 2x positive samples
  n_neg_test <- max(1, n_neg_test)
  
  if (n_neg_test > n_neg_total) {
    warning(sprintf("Test negative-sample demand %d exceeds available negatives %d; using all negatives as the test set", 
                    n_neg_test, n_neg_total))
    n_neg_test <- n_neg_total
  }
  
  neg_idx <- sample(1:n_neg_total, n_neg_test)
  neg_test  <- all_neg_candidate[neg_idx, ]
  neg_train_pool <- all_neg_candidate[-neg_idx, ]
  
  n_neg_train_need <- floor(n_pos_train * neg_ratio)
  if (nrow(neg_train_pool) < n_neg_train_need) {
    warning(sprintf("Insufficient training negatives (need %d, have %d); using all available negatives", 
                    n_neg_train_need, nrow(neg_train_pool)))
    neg_train <- neg_train_pool
  } else {
    neg_train <- slice_sample(neg_train_pool, n = n_neg_train_need)
  }
  
  tr_pos_n <- nrow(pos_train)
  tr_neg_n <- nrow(neg_train)
  te_pos_n <- nrow(pos_test)
  te_neg_n <- nrow(neg_test)
  actual_neg_ratio <- round(tr_neg_n / tr_pos_n, 3)
  
  pos_train_labeled <- pos_train %>% mutate(label = 1)
  neg_train_labeled <- neg_train %>% mutate(label = 0)
  train_set <- bind_rows(pos_train_labeled, neg_train_labeled)
  
  if (length(unique(train_set$label)) < 2) {
    warning("The training set lacks one class of samples; cannot train")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, te_neg=te_neg_n),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=actual_neg_ratio)
    ))
  }
  
  pred_set <- wide_subset %>% select(TF, Target) %>% distinct()
  train_feat_raw <- left_join(train_set, wide_subset, by = c("TF", "Target"))
  pred_feat_raw  <- left_join(pred_set, wide_subset, by = c("TF", "Target"))
  
  x_train_raw <- as.matrix(train_feat_raw[, feat_cols, drop = FALSE])
  x_pred_raw  <- as.matrix(pred_feat_raw[, feat_cols, drop = FALSE])
  
  # Check whether features exist
  if (ncol(x_train_raw) == 0) {
    warning("Training features are empty; skipping")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, te_neg=te_neg_n),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=actual_neg_ratio)
    ))
  }
  
  # Remove zero-variance features
  zero_var <- apply(x_train_raw, 2, function(v) var(v, na.rm = TRUE) < zero_var_threshold)
  zero_feat_names <- paste(feat_cols[zero_var], collapse = ";")
  if (any(zero_var)) {
    cat("    Removing zero-variance features: ", zero_feat_names, "\n")
    x_train_raw <- x_train_raw[, !zero_var, drop = FALSE]
    x_pred_raw  <- x_pred_raw[, !zero_var, drop = FALSE]
  }
  
  # Check feature count again
  if (ncol(x_train_raw) == 0) {
    warning("No features remain after removing zero-variance features; skipping")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, te_neg=te_neg_n),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats=zero_feat_names, Actual_Neg_Ratio=actual_neg_ratio)
    ))
  }
  
  # Normalization
  if (normalize_w) {
    x_train <- normalize_by_train_stats(x_train_raw, x_train_raw)
    x_pred  <- normalize_by_train_stats(x_train_raw, x_pred_raw)
  } else {
    x_train <- x_train_raw
    x_pred  <- x_pred_raw
  }
  y_train <- train_feat_raw$label
  
  # Stratified CV to select lambda
  best_lambda <- 1e-4
  n_train <- length(y_train)  # Fix: define n_train
  pos_num <- sum(y_train == 1)
  neg_num <- sum(y_train == 0)
  
  # Check whether cross-validation conditions are met
  if (n_train >= 6 && pos_num >= 3 && neg_num >= 3) {
    nfolds <- min(5, pos_num, neg_num)
    nfolds <- max(3, nfolds)
    
    # Ensure each fold has enough samples
    if (pos_num >= nfolds && neg_num >= nfolds) {
      tryCatch({
        foldid <- get_stratified_foldid(y_train, nfolds)
        cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0,
                            nfolds = nfolds, foldid = foldid, 
                            type.measure = "auc",
                            lambda = 10^seq(-3, 1, length = 20))  # More reasonable lambda range
        if (!is.null(cv_fit$lambda.min)) {
          best_lambda <- cv_fit$lambda.min
          cat("    CV optimal lambda =", best_lambda, "\n")
        }
      }, error = function(e) {
        warning("Stratified cross-validation failed; using fallback lambda=1e-4:", conditionMessage(e))
      })
    } else {
      warning(sprintf("Insufficient samples for %d-fold CV (pos:%d neg:%d); using fallback lambda=1e-4", 
                      nfolds, pos_num, neg_num))
    }
  } else {
    warning(sprintf("Insufficient samples (total:%d pos:%d neg:%d); using fallback lambda=1e-4", 
                    n_train, pos_num, neg_num))
  }
  
  # Train the final model
  final_model <- glmnet(x_train, y_train, family = "binomial", alpha = 0, 
                        lambda = best_lambda)
  
  # Predict and clamp to 0~1
  raw_pred <- as.vector(predict(final_model, newx = x_pred, type = "response"))
  pred_set$Stacking_Score <- round(pmax(0, pmin(1, raw_pred)), 4)
  result <- pred_set %>% select(TF, Target, Stacking_Score) %>% arrange(desc(Stacking_Score))
  
  # ========== Test-set evaluation ==========
  test_set <- bind_rows(pos_test %>% mutate(label = 1), neg_test %>% mutate(label = 0))
  test_feat <- left_join(test_set, pred_set, by = c("TF", "Target"))
  auc_roc <- NA
  auc_pr  <- NA
  
  if (nrow(test_feat) > 0 && length(unique(test_feat$label)) >= 2) {
    lab_vec <- test_feat$label
    scr_vec <- test_feat$Stacking_Score
    
    # Check whether valid predictions exist
    if (sum(!is.na(scr_vec)) > 0 && sum(lab_vec == 1) > 0 && sum(lab_vec == 0) > 0) {
      tryCatch({
        roc_obj <- roc(lab_vec, scr_vec, quiet = TRUE)
        auc_roc <- round(as.numeric(auc(roc_obj)), 4)
        pr_obj <- pr.curve(
          scores.class0 = scr_vec[lab_vec == 1],
          scores.class1 = scr_vec[lab_vec == 0],
          curve = FALSE
        )
        auc_pr <- round(pr_obj$auc.integral, 4)
      }, error = function(e) {
        warning("AUC evaluation failed:", conditionMessage(e))
      })
    }
  }
  
  cat(sprintf("    Training set: %d pos, %d neg | Test set: %d pos, %d neg\n", 
              tr_pos_n, tr_neg_n, te_pos_n, te_neg_n))
  cat(sprintf("    Test AUC-ROC = %.4f, AUC-PR = %.4f\n", auc_roc, auc_pr))
  
  return(list(
    result = result,
    metrics = c(AUC_ROC = auc_roc, AUC_PR = auc_pr),
    sample_info = c(tr_pos=tr_pos_n, tr_neg=tr_neg_n, te_pos=te_pos_n, te_neg=te_neg_n),
    extra_info = c(Used_Lambda=best_lambda, Zero_Var_Feats=zero_feat_names, 
                   Actual_Neg_Ratio=actual_neg_ratio)
  ))
}

# Exception-capture wrapper function
safe_run_stacking <- function(method_subset, all_base_wide, chip_pos, comb_name,
                              neg_ratio, normalize_w, test_frac, seed, 
                              filter_zero_neg) {
  res <- tryCatch({
    stacking_supervised_chip(method_subset, all_base_wide, chip_pos,
                             neg_ratio = neg_ratio, normalize_w = normalize_w,
                             test_frac = test_frac, local_seed = seed, 
                             filter_zero_neg = filter_zero_neg)
  }, error = function(e) {
    cat("  >> Combination", comb_name, "error:", conditionMessage(e), "\n")
    return(list(
      result = data.frame(TF = character(), Target = character(), Stacking_Score = numeric()),
      metrics = c(AUC_ROC = NA, AUC_PR = NA),
      sample_info = c(tr_pos=0, tr_neg=0, te_pos=0, te_neg=0),
      extra_info = c(Used_Lambda=NA, Zero_Var_Feats="", Actual_Neg_Ratio=NA)
    ))
  })
  
  info <- res$sample_info
  m_val <- res$metrics
  extra <- res$extra_info
  out_df <- data.frame(
    Combination = comb_name,
    NumMethods = length(method_subset),
    Methods = paste(sort(names(method_subset)), collapse = ";"),
    Train_Pos = info["tr_pos"],
    Train_Neg = info["tr_neg"],
    Test_Pos = info["te_pos"],
    Test_Neg = info["te_neg"],
    Actual_Neg_Ratio = extra["Actual_Neg_Ratio"],
    Used_Lambda = extra["Used_Lambda"],
    Zero_Var_Feats = extra["Zero_Var_Feats"],
    AUC_ROC = m_val["AUC_ROC"],
    AUC_PR = m_val["AUC_PR"],
    Test_Balanced = TRUE,
    stringsAsFactors = FALSE
  )
  return(list(out_df = out_df, pred_df = res$result))
}

# ======================================================
# Main program: batch process multiple stress datasets
# ======================================================
datasets <- c("Aridity", "Alkalinity", "Cold")
chip_gt_pos <- read_chip_pos(chip_file_path)
eval_summary <- data.frame()

for (ds_idx in seq_along(datasets)) {
  ds <- datasets[ds_idx]
  cat("\n=============================================\n")
  cat("Supervised Stacking (ChIP gold standard): ", ds, "\n")
  cat("=============================================\n")
  set.seed(123)
  
  base_files <- list(
    GENIE3    = file.path("GENIE3", paste0(ds, "_significant_edges.txt")),
    Kboost    = file.path("Kboost", paste0(ds, "_significant_edges.csv")),
    GRNBoost  = file.path("GRNBoost", paste0(ds, "_significant_edges.csv")),
    `3DCEMA`  = file.path("3DCEMA", paste0(ds, "_significant_edges.csv")),
    DeepRIG   = file.path("DeepRIG", paste0(ds, "_significant_edges.csv")),
    IGEGRNs   = file.path("IGEGRNs", paste0(ds, "_significant_edges.csv"))
  )
  missing <- base_files[!sapply(base_files, file.exists)]
  if (length(missing) > 0) {
    stop(paste0("Dataset ", ds, " missing files: ", paste(names(missing), collapse = ", ")))
  }
  
  genie3     <- clean_read(base_files$GENIE3, delim = "\t")
  kboost     <- clean_read(base_files$Kboost, delim = ",")
  grnboost   <- clean_read(base_files$GRNBoost, delim = ",")
  dcema3     <- clean_read(base_files$`3DCEMA`, delim = ",")
  deeprig    <- clean_read(base_files$DeepRIG, delim = ",")
  igegrn     <- clean_read(base_files$IGEGRNs, delim = ",")
  
  all_methods <- list(
    GENIE3    = genie3,
    Kboost    = kboost,
    GRNBoost  = grnboost,
    `3DCEMA`  = dcema3,
    DeepRIG   = deeprig,
    IGEGRNs   = igegrn
  )
  method_names <- names(all_methods)
  
  # ========== Fix: correctly build the wide table ==========
  cat("  Building wide table...\n")
  all_long <- bind_rows(lapply(names(all_methods), function(name) {
    df <- all_methods[[name]]
    df %>%
      select(TF, Target, Weight) %>%
      mutate(Model = name)
  }))
  
  # Use pivot_wider to build the wide table
  all_base_wide <- all_long %>%
    pivot_wider(
      names_from = Model, 
      values_from = Weight, 
      values_fill = 0
    ) %>%
    filter(TF != Target) %>%
    distinct(TF, Target, .keep_all = TRUE)
  
  # Check whether the wide table was built successfully
  if (nrow(all_base_wide) == 0) {
    warning(sprintf("The wide table for dataset %s is empty; skipping", ds))
    next
  }
  
  cat(sprintf("  Wide table built: %d edges, %d features\n", 
              nrow(all_base_wide), ncol(all_base_wide) - 2))
  
  total_pos <- nrow(inner_join(chip_gt_pos, all_base_wide[,c("TF","Target")], 
                               by=c("TF","Target")))
  if(total_pos < 5){
    warning(sprintf("Dataset %s matches only %d ChIP positive edges; too few samples; skipping all combinations", ds, total_pos))
    next
  }
  cat(sprintf("  Matched %d ChIP positive edges\n", total_pos))
  
  # Create output directories
  out_root <- file.path(root_dir, ds)
  dir_pair    <- file.path(out_root, "1_Pair_2Methods")
  dir_triple  <- file.path(out_root, "2_Triple_3Methods")
  dir_4way    <- file.path(out_root, "3_Four_4Methods")
  dir_5way    <- file.path(out_root, "4_Five_5Methods")
  dir_all6    <- file.path(out_root, "5_All_6Methods")
  all_dirs <- c(out_root, dir_pair, dir_triple, dir_4way, dir_5way, dir_all6)
  for (d in all_dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  
  # Generic function to process combinations
  process_combinations <- function(comb_list, output_dir, seed_base) {
    batch_metrics <- list()
    idx <- 1
    total_comb <- length(comb_list)
    cat("  Total combinations in this batch: ", total_comb, "\n")
    
    # Add a progress bar
    pb <- progress_bar$new(
      format = "    processing [:bar] :percent remaining: :eta",
      total = total_comb, 
      clear = FALSE,
      width = 60
    )
    
    for (nm in comb_list) {
      pb$tick()
      nm_sorted <- sort(nm)
      comb_name <- paste(nm_sorted, collapse = "_")
      
      # Print detailed info only occasionally to avoid excessive output
      if (idx %% 10 == 1 || idx == total_comb) {
        cat("\n  Processing combination [", idx, "/", total_comb, "]: ", comb_name, "\n", sep="")
      }
      
      comb_seed <- seed_base + idx * 100
      run_res <- safe_run_stacking(all_methods[nm], all_base_wide, chip_gt_pos, 
                                   comb_name, sup_neg_ratio, normalize_weights, 
                                   test_ratio, seed = comb_seed, 
                                   filter_zero_neg = filter_zero_neg)
      
      if (nrow(run_res$pred_df) > 0) {
        out_file <- file.path(output_dir, paste0("SupStacking_", comb_name, ".csv"))
        write.csv(run_res$pred_df, out_file, row.names = FALSE)
      }
      batch_metrics[[idx]] <- run_res$out_df
      idx <- idx + 1
    }
    cat("\n")  # newline
    bind_rows(batch_metrics)
  }
  
  ds_metrics <- data.frame()
  base_seed <- ds_idx * 10000
  
  if (only_all6) {
    cat("Running only the full 6-model combination\n")
    ds_metrics <- process_combinations(list(method_names), dir_all6, 
                                       seed_base = base_seed)
  } else {
    cat("Iterating over all 2~6-model combinations\n")
    # Generate all combinations
    cat("  Generating 2-method combinations...\n")
    m2_list <- combn(method_names, 2, simplify = FALSE)
    cat("  Generating 3-method combinations...\n")
    m3_list <- combn(method_names, 3, simplify = FALSE)
    cat("  Generating 4-method combinations...\n")
    m4_list <- combn(method_names, 4, simplify = FALSE)
    cat("  Generating 5-method combinations...\n")
    m5_list <- combn(method_names, 5, simplify = FALSE)
    
    m2 <- process_combinations(m2_list, dir_pair, base_seed + 1000)
    m3 <- process_combinations(m3_list, dir_triple, base_seed + 2000)
    m4 <- process_combinations(m4_list, dir_4way, base_seed + 3000)
    m5 <- process_combinations(m5_list, dir_5way, base_seed + 4000)
    m6 <- process_combinations(list(method_names), dir_all6, base_seed + 5000)
    ds_metrics <- bind_rows(m2, m3, m4, m5, m6)
  }
  
  ds_metrics$Dataset <- ds
  eval_summary <- bind_rows(eval_summary, ds_metrics)
  cat("\n[Done]", ds, "dataset: all combinations computed, output path: ", out_root, "\n")
  
  # Clear memory
  rm(all_methods, all_base_wide, all_long)
  gc()
}

# Save and sort the evaluation summary
eval_summary <- eval_summary %>%
  select(Dataset, Combination, NumMethods, Methods,
         Train_Pos, Train_Neg, Test_Pos, Test_Neg, Actual_Neg_Ratio,
         Used_Lambda, Zero_Var_Feats, Test_Balanced, AUC_ROC, AUC_PR) %>%
  arrange(Dataset, NumMethods, Combination)

write.csv(eval_summary, file.path(root_dir, "Evaluation_Results.csv"), row.names = FALSE)

cat("\n[Finished] Supervised Stacking for all stress datasets completed!\n")
cat("Evaluation summary table: ", file.path(root_dir, "Evaluation_Results.csv"), "\n")
cat("Total valid combinations: ", nrow(eval_summary), "\n")

valid_roc <- eval_summary[!is.na(eval_summary$AUC_ROC), ]
if (nrow(valid_roc) > 0) {
  cat("  Valid AUC-ROC range: ", paste(range(valid_roc$AUC_ROC, na.rm = TRUE), collapse = " - "), "\n")
  cat("  Valid AUC-PR range: ", paste(range(valid_roc$AUC_PR, na.rm = TRUE), collapse = " - "), "\n")
  
  # Output the best combination
  best_roc <- valid_roc[which.max(valid_roc$AUC_ROC), ]
  cat("Best AUC-ROC combination: ", best_roc$Combination, 
      " (AUC-ROC:", best_roc$AUC_ROC, ", AUC-PR:", best_roc$AUC_PR, ")\n")
} else {
  cat("No valid evaluation results\n")
}
