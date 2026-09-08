# ==============================
# Quantile‑Rank Fusion (QRF) Batch Pipeline: Alkalinity, Aridity, Cold
# Equal‑weight fusion, sum of weights within each ensemble = 1, final score range 0~1
# Output directory: Fusion_Results_{dataset}_0‑1
# Subdirectories: 1_Pair_2Methods / 2_Triple_3Methods / 3_Four_4Methods / 4_Five_5Methods / 5_All_6Methods
# Directly compatible with downstream run_evaluation() assessment script
# ==============================
library(dplyr)
library(tidyr)
library(vroom)

# ------------------------------------------------
# 1. Core fusion function (Quantile‑Rank Fusion, QRF)
# ------------------------------------------------
quantile_rank_fusion <- function(method_list, weights = NULL) {
  
  process_one <- function(df) {
    df <- df[, 1:3]
    colnames(df) <- c("TF", "Target", "Weight")
    df <- df[!is.na(df$Weight) & df$Weight != "", ]
    
    df <- df %>%
      arrange(desc(as.numeric(Weight))) %>%
      mutate(quant_score = 1 - (row_number() / n())) %>%
      select(TF, Target, quant_score)
    
    return(df)
  }
  
  processed_list <- lapply(method_list, process_one)
  
  all_edges <- bind_rows(lapply(processed_list, function(x) x[, 1:2])) %>% distinct()
  
  for (i in seq_along(processed_list)) {
    current <- processed_list[[i]]
    col_name <- paste0("s", i)
    all_edges[[col_name]] <- current$quant_score[match(
      interaction(all_edges$TF, all_edges$Target),
      interaction(current$TF, current$Target)
    )]
    all_edges[[col_name]][is.na(all_edges[[col_name]])] <- 0
  }
  
  # Equal‑weight setting: sum of weights within ensemble equals 1
  if (is.null(weights)) {
    weights <- rep(1 / length(processed_list), length(processed_list))
  }
  
  score_cols <- paste0("s", seq_along(processed_list))
  all_edges$final_score <- as.matrix(all_edges[, score_cols]) %*% matrix(weights)
  
  result <- all_edges %>%
    arrange(desc(final_score)) %>%
    select(TF, Target, EdgeWeight = final_score)
  
  return(result)
}

# ------------------------------------------------
# 2. Input reading and cleaning function
# ------------------------------------------------
clean_read <- function(path, delim = ",") {
  df <- vroom(path, delim = delim, show_col_types = FALSE)
  df <- df[, 1:3]
  colnames(df) <- c("TF", "Target", "Weight")
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$Weight), ]
  return(df)
}

# ------------------------------------------------
# 3. Main function for processing single dataset
# ------------------------------------------------
run_qrf_fusion <- function(dataset) {
  cat("\n=============================================\n")
  cat("Processing dataset: ", dataset, "\n")
  cat("=============================================\n")
  
  # Load six GRN inference outputs for current stress condition
  genie3     <- clean_read(paste0("GENIE3/", dataset, "_significant_edges.txt"), delim = "\t")
  kboost     <- clean_read(paste0("Kboost/", dataset, "_significant_edges.csv"))
  grnboost2  <- clean_read(paste0("GRNBoost/", dataset, "_significant_edges.csv"))
  dcema      <- clean_read(paste0("3DCEMA/", dataset, "_significant_edges.csv"))
  deeprig    <- clean_read(paste0("DeepRIG/", dataset, "_significant_edges.csv"))
  ige        <- clean_read(paste0("IGEGRNs/", dataset, "_significant_edges.csv"))
  
  all_methods <- list(
    GENIE3    = genie3,
    Kboost    = kboost,
    GRNBoost2 = grnboost2,
    DCEMA     = dcema,
    DeepRIG   = deeprig,
    IGE       = ige
  )
  
  # Create output directory tree
  out_root <- paste0("Fusion_Results_", dataset, "_0-1")
  dir_pair    <- file.path(out_root, "1_Pair_2Methods")
  dir_triple  <- file.path(out_root, "2_Triple_3Methods")
  dir_4way    <- file.path(out_root, "3_Four_4Methods")
  dir_5way    <- file.path(out_root, "4_Five_5Methods")
  dir_all6    <- file.path(out_root, "5_All_6Methods")
  
  lapply(c(out_root, dir_pair, dir_triple, dir_4way, dir_5way, dir_all6), function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  })
  
  # All 2‑model combinations
  cat("Computing all 2‑model combinations\n")
  comb2 <- combn(names(all_methods), 2, simplify = FALSE)
  for (nm in comb2) {
    res <- quantile_rank_fusion(all_methods[nm])
    write.csv(res, file.path(dir_pair, paste0("Fusion_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # All 3‑model combinations
  cat("Computing all 3‑model combinations\n")
  comb3 <- combn(names(all_methods), 3, simplify = FALSE)
  for (nm in comb3) {
    res <- quantile_rank_fusion(all_methods[nm])
    write.csv(res, file.path(dir_triple, paste0("Fusion_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # All 4‑model combinations
  cat("Computing all 4‑model combinations\n")
  comb4 <- combn(names(all_methods), 4, simplify = FALSE)
  for (nm in comb4) {
    res <- quantile_rank_fusion(all_methods[nm])
    write.csv(res, file.path(dir_4way, paste0("Fusion_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # All 5‑model combinations
  cat("Computing all 5‑model combinations\n")
  comb5 <- combn(names(all_methods), 5, simplify = FALSE)
  for (nm in comb5) {
    res <- quantile_rank_fusion(all_methods[nm])
    write.csv(res, file.path(dir_5way, paste0("Fusion_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # Full 6‑model ensemble
  cat("Computing full 6‑model fusion\n")
  all6 <- quantile_rank_fusion(all_methods)
  write.csv(all6, file.path(dir_all6, "Fusion_All6Methods.csv"), row.names = FALSE)
  
  cat("Completed QRF fusion for ", dataset, ". Output directory: ", out_root, "\n")
}

# ------------------------------------------------
# Run pipeline for all three stress datasets
# ------------------------------------------------
datasets <- c("Alkalinity", "Aridity", "Cold")
for (ds in datasets) {
  run_qrf_fusion(ds)
}

cat("\nAll‑dataset QRF fusion finished. Proceed to run_evaluation() and run_metrics() for downstream assessment.\n")
