# Load required packages (keep only the essential ones)
library(tidyverse)
library(pROC)      # Compute AUROC
library(PRROC)     # Compute AUPR
library(vroom)     # Fast data reading
library(ROSE) 

# ------------------------------------------------------------------------------
# 1. Data-loading function (compatible with csv/txt suffixes)
# ------------------------------------------------------------------------------
load_grn_data <- function(file_path) {
  vroom(file_path, show_col_types = FALSE)[, 1:3] %>%
    rename(TF = 1, Target = 2, EdgeWeight = 3) %>%  # Unify column names
    drop_na() %>%                                    # Remove rows containing NA
    distinct(TF, Target, .keep_all = TRUE)           # Deduplicate: keep unique TF-target pairs
}

# ChIP-seq gold standard (shared globally)
chip_seq <- read_csv("chip_seq.csv", col_names = TRUE) %>% 
  dplyr::select(TF = 1, Target = 2) %>%
  distinct()

# ------------------------------------------------------------------------------
# 2. Overlap statistics, null-model permutation, and comprehensive evaluation functions (original logic fully retained)
# ------------------------------------------------------------------------------
calc_overlap <- function(pred_df, chip_df, top_n = Inf) { # Inf = all
  pred_top <- pred_df %>%
    arrange(desc(EdgeWeight)) %>%
    slice_head(n = top_n) %>%
    mutate(pair = paste(TF, Target, sep = "_")) %>% pull(pair)
  
  true_pairs <- paste(chip_df$TF, chip_df$Target, sep = "_")
  
  # Debug printing
  cat("  Total predicted edges:", length(pred_top), "\n")
  cat("  Total ChIP edges:", length(true_pairs), "\n")
  cat("  Number of shared TFs:", length(intersect(pred_df$TF, chip_df$TF)), "\n")
  cat("  True overlap count:", sum(pred_top %in% true_pairs), "\n")
  
  overlap <- sum(pred_top %in% true_pairs)
  return(overlap)
}

# Permutation null model
null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  return(list(true = true_overlap, null = null_dist, pvalue = p_val))
}

# AUROC+AUPR+permutation test comprehensive evaluation
evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  pred_df <- pred_df %>%
    mutate(label = ifelse(paste(TF, Target) %in% paste(chip_df$TF, chip_df$Target), 1, 0))
  
  roc_obj <- roc(pred_df$label, pred_df$EdgeWeight, quiet = TRUE)
  auroc <- auc(roc_obj)
  
  pr_obj <- pr.curve(
    scores.class0 = pred_df$EdgeWeight[pred_df$label == 1],
    scores.class1 = pred_df$EdgeWeight[pred_df$label == 0]
  )
  aupr <- pr_obj$auc.integral
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    auroc = auroc,
    aupr = aupr,
    true_overlap = null_res$true,
    null_mean = mean(null_res$null),
    p_value = null_res$pvalue
  ))
}

# ------------------------------------------------------------------------------
# 3. Batch configuration: 6 algorithms + 3 stresses
# ------------------------------------------------------------------------------
# Algorithm folder names & file suffixes
method_config <- list(
  "IGEGRNs" = list(suffix = ".csv"),
  "3DCEMA" = list(suffix = ".csv"),
  "DeepRIG" = list(suffix = ".csv"),
  "Kboost" = list(suffix = ".csv"),
  "GRNBoost" = list(suffix = ".csv"),
  "GENIE3" = list(suffix = ".txt")
)

# Stress-name list
stress_names <- c("Alkalinity", "Aridity", "Cold")

# ========== Fix point 1: pre-initialize an empty data frame with complete column names ==========
all_results <- tibble(
  Method = character(),
  Stress = character(),
  AUROC = numeric(),
  AUPR = numeric(),
  True_Overlap = integer(),
  Null_Mean = numeric(),
  P_Value = character()
)

# Double loop: outer loop over algorithms, inner loop over stresses
for (method_name in names(method_config)) {
  cat("\n==================== Starting evaluation of algorithm: ", method_name, "====================\n")
  file_suffix <- method_config[[method_name]]$suffix
  
  for (stress in stress_names) {
    # Concatenate the file path
    file_path <- paste0(method_name, "/", stress, "_significant_edges", file_suffix)
    cat("\n>> Current: ", method_name, " | ", stress, " | file path: ", file_path, "\n")
    
    # Read the GRN network for the current algorithm + stress
    grn_data <- load_grn_data(file_path)
    
    # Perform evaluation (1000 permutations)
    eval_res <- evaluate_with_null(
      pred_df = grn_data,
      chip_df = chip_seq,
      top_n = Inf,
      n_perm = 1000
    )
    
    # Convert the P value to scientific notation with 2 significant digits
    p_sci <- format(eval_res$p_value, scientific = TRUE, digits = 2)
    
    # ========== Fix point 2: use bind_rows for concatenation, compatible with the initialized tibble ==========
    temp_row <- tibble(
      Method = method_name,
      Stress = stress,
      AUROC = round(eval_res$auroc, 3),
      AUPR = round(eval_res$aupr, 3),
      True_Overlap = eval_res$true_overlap,
      Null_Mean = round(eval_res$null_mean, 2),
      P_Value = p_sci
    )
    all_results <- bind_rows(all_results, temp_row)
    
    cat(stress, "evaluation completed, P =", p_sci, "\n")
  }
}

# ------------------------------------------------------------------------------
# 4. Output and save all results
# ------------------------------------------------------------------------------
# Print the full results table
print(all_results, row.names = FALSE)

# Export to csv (containing all metrics for 6 algorithms x 3 stresses)
write.csv(all_results, "six_methods_stress_evaluation.csv", row.names = FALSE)
cat("\nAll algorithm + stress evaluations completed; results saved to six_methods_stress_evaluation.csv\n")
