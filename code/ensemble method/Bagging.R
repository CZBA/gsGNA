library(dplyr)
library(tidyr)
library(vroom)
library(matrixStats)

# ============================== Global configuration (centralized path/parameter management) ==============================
set.seed(123)
# Run mode switch: debug=25 / mid=200 / final=500
run_mode <- "final"
# Automatically select the matching Bootstrap iteration count for each mode
B <- dplyr::case_when(
  run_mode == "debug" ~ 25,    # Breiman L. Bagging Predictors, Machine Learning, 1996
  run_mode == "mid"  ~ 200,    # Lau M et al. Scientific Reports, 2023
  run_mode == "final"~ 500     # Ghosh Roy et al. Bioinformatics, 2020
)
normalize_weights <- TRUE   # Global normalization switch
top_pct <- 0.1               # Keep the top 10% of high-weight edges
# Input root directory
input_dirs <- list(
  GENIE3   = "GENIE3",
  Kboost   = "Kboost",
  GRNBoost2= "GRNBoost",
  DCEMA    = "3DCEMA",
  DeepRIG  = "DeepRIG",
  IGE      = "IGEGRNs"
)
# Output root directory
out_root_main <- "Bagging_Unsup_AlgoBootstrap"
# Dataset list
datasets <- c("Aridity", "Alkalinity", "Cold")
# Print current parameters for reproducibility
cat("==== Current Run Configuration ====\n")
cat("Run mode: ", run_mode, "\n")
cat("Bootstrap iterations B =", B, "\n")
cat("Min-Max normalization: ", normalize_weights, "\n")
cat("Keeping top", top_pct*100, "% high-weight edges, removing edges with weight = 0\n")
cat("==================================\n")
# ==========================================================================================

# ------------------------------------------------
# Utility 1: read + clean a single algorithm's results
# ------------------------------------------------
clean_read <- function(path, delim = ",") {
  df <- vroom(path, delim = delim, show_col_types = FALSE)
  df <- df[, 1:3]
  colnames(df) <- c("TF", "Target", "Weight")
  df <- df %>% drop_na(TF, Target, Weight)
  df$TF <- toupper(df$TF)
  df$Target <- toupper(df$Target)
  df <- df %>% filter(TF != Target)
  # For duplicate edges, keep the maximum weight
  df <- df %>%
    group_by(TF, Target) %>%
    summarise(Weight = max(Weight), .groups = "drop")
  return(df)
}

# ------------------------------------------------
# Utility 2: fill missing values with 0
# ------------------------------------------------
fill_na_zero <- function(df, feat_cols) {
  df <- df %>%
    mutate(across(all_of(feat_cols), ~ suppressWarnings(as.numeric(.)))) %>%
    mutate(across(all_of(feat_cols), ~ ifelse(is.na(.), 0, .)))
  return(df)
}

# ------------------------------------------------
# Utility 3: vectorized min-max normalization (optimized)
# ------------------------------------------------
normalize_cols <- function(df, cols) {
  if (length(cols) == 0) return(df)
  mat <- as.matrix(df[, cols])
  mat <- apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  })
  df[, cols] <- mat
  return(df)
}

# ------------------------------------------------
# Core: unsupervised Bagging ensemble (based on algorithm-column Bootstrap)
# New post-processing: remove zero-weight edges, keep only the top 10% high-weight regulatory edges
# ------------------------------------------------
bagging_fusion_unsup <- function(wide_df, feat_cols) {
  M <- length(feat_cols)
  if (M == 0) stop("No algorithm feature columns provided")
  N <- nrow(wide_df)
  if (N == 0) {
    warning("No candidate edges")
    return(data.frame(TF = character(), Target = character(), EdgeWeight = numeric()))
  }
  
  weight_mat <- as.matrix(wide_df[, feat_cols, drop = FALSE])
  final_scores <- numeric(N)
  
  for (b in 1:B) {
    # Sample algorithm columns with replacement
    sampled_cols <- sample(1:M, size = M, replace = TRUE)
    local_scores <- rowMeans(weight_mat[, sampled_cols, drop = FALSE])
    final_scores <- final_scores + local_scores
  }
  
  final_scores <- final_scores / B
  res <- wide_df %>%
    select(TF, Target) %>%
    mutate(EdgeWeight = final_scores) %>%
    arrange(desc(EdgeWeight))
  
  # ========== New filtering logic ==========
  # 1. Remove edges with weight equal to 0
  res <- res %>% filter(EdgeWeight > 0)
  if(nrow(res) == 0){
    warning("No edges remain after zero-weight filtering")
    return(res)
  }
  # 2. Keep only the top 10% of high-weight edges
  keep_num <- ceiling(nrow(res) * top_pct)
  res <- res %>% slice_head(n = keep_num)
  # ==================================
  
  return(res)
}

# ------------------------------------------------
# Exception-capture wrapper
# ------------------------------------------------
safe_bagging <- function(wide_df, feat_cols) {
  tryCatch(
    bagging_fusion_unsup(wide_df = wide_df, feat_cols = feat_cols),
    error = function(e) {
      cat("  Error: ", e$message, "\n")
      return(data.frame(TF = character(), Target = character(), EdgeWeight = numeric()))
    }
  )
}

# ------------------------------------------------
# Run + save wrapper: file name includes the B-round marker
# ------------------------------------------------
run_and_save <- function(wide_df, feat_cols, comb_name, out_dir) {
  cat("  Processing: ", comb_name, "\n")
  start_t <- Sys.time()
  out_res <- safe_bagging(wide_df, feat_cols)
  end_t <- Sys.time()
  cat("  Elapsed: ", round(end_t - start_t, 2), "s | Output edges: ", nrow(out_res), "\n")
  
  if (nrow(out_res) > 0) {
    fpath <- file.path(out_dir, paste0("Bagging_", comb_name, "_B", B, "_Top10pct.csv"))
    write.csv(out_res, fpath, row.names = FALSE)
  } else {
    cat("  Warning: combination", comb_name, "has no valid edges after filtering, not saved\n")
  }
}

# ======================================================
# Main workflow: batch process datasets one by one
# ======================================================
for(ds in datasets) {
  cat("\n=============================================\n")
  cat("Unsupervised Algorithm-Column Bootstrap Bagging | Dataset: ", ds, " | B =", B, " | Output Top10% non-zero-weight edges\n")
  cat("=============================================\n")
  
  # 1. Read all 6 algorithm results for the current dataset at once
  cat("  Reading raw data of each algorithm...\n")
  genie3     <- clean_read(file.path(input_dirs$GENIE3,    paste0(ds, "_significant_edges.txt")), delim = "\t")
  kboost     <- clean_read(file.path(input_dirs$Kboost,    paste0(ds, "_significant_edges.csv")))
  grnboost2  <- clean_read(file.path(input_dirs$GRNBoost2,paste0(ds, "_significant_edges.csv")))
  dcema      <- clean_read(file.path(input_dirs$DCEMA,     paste0(ds, "_significant_edges.csv")))
  deeprig    <- clean_read(file.path(input_dirs$DeepRIG,   paste0(ds, "_significant_edges.csv")))
  ige        <- clean_read(file.path(input_dirs$IGE,       paste0(ds, "_significant_edges.csv")))
  
  all_methods <- list(
    GENIE3=genie3, Kboost=kboost, GRNBoost2=grnboost2,
    DCEMA=dcema, DeepRIG=deeprig, IGE=ige
  )
  method_names <- names(all_methods)
  
  # 2. Build the global edge set + global wide table (done once, reused by all combinations)
  cat("  Building global feature wide table...\n")
  all_edges <- bind_rows(lapply(all_methods, function(df) df[, c("TF", "Target")])) %>%
    distinct() %>%
    filter(TF != Target)
  
  # Append weight columns of all algorithms
  for (m in method_names) {
    col_name <- paste0("feat_", m)
    all_edges <- all_edges %>%
      left_join(all_methods[[m]] %>% select(TF, Target, Weight),
                by = c("TF", "Target")) %>%
      rename(!!col_name := Weight)
  }
  feat_all <- paste0("feat_", method_names)
  
  # 3. Global zero-filling + one-shot global normalization (to keep all combinations on the same scale)
  all_edges <- fill_na_zero(all_edges, feat_all)
  if (normalize_weights) {
    all_edges <- normalize_cols(all_edges, feat_all)
  }
  
  # 4. Create hierarchical output directories
  out_root <- file.path(out_root_main, ds)
  dir_pair   <- file.path(out_root, "Pair2")
  dir_triple <- file.path(out_root, "Triple3")
  dir_4way   <- file.path(out_root, "Four4")
  dir_5way   <- file.path(out_root, "Five5")
  dir_all6   <- file.path(out_root, "All6")
  
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  lapply(c(dir_pair, dir_triple, dir_4way, dir_5way, dir_all6),
         dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # 5. Iterate over all combinations
  cat("  Starting to iterate over algorithm combinations...\n")
  
  # 2-algorithm combinations
  comb2 <- combn(method_names, 2, simplify = FALSE)
  for (nm in comb2) {
    cols <- paste0("feat_", nm)
    run_and_save(all_edges, cols, paste(nm, collapse = "+"), dir_pair)
  }
  
  # 3-algorithm combinations
  comb3 <- combn(method_names, 3, simplify = FALSE)
  for (nm in comb3) {
    cols <- paste0("feat_", nm)
    run_and_save(all_edges, cols, paste(nm, collapse = "+"), dir_triple)
  }
  
  # 4-algorithm combinations
  comb4 <- combn(method_names, 4, simplify = FALSE)
  for (nm in comb4) {
    cols <- paste0("feat_", nm)
    run_and_save(all_edges, cols, paste(nm, collapse = "+"), dir_4way)
  }
  
  # 5-algorithm combinations
  comb5 <- combn(method_names, 5, simplify = FALSE)
  for (nm in comb5) {
    cols <- paste0("feat_", nm)
    run_and_save(all_edges, cols, paste(nm, collapse = "+"), dir_5way)
  }
  
  # All-6-algorithm combination
  run_and_save(all_edges, feat_all, "All6", dir_all6)
  
  cat("\n[Done] ", ds, " dataset: all combinations processed, B =", B, ", only top-10% non-zero-weight edges output\n")
}

cat("\n[Finished] All datasets: unsupervised algorithm-column Bootstrap Bagging completed! Current iterations B =", B, "\n")
