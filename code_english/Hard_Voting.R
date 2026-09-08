# ==============================
# Weighted hard-voting fusion + significant-edge filtering (Scheme 2: retain original weights, continuous score)
# Datasets: Aridity, Alkalinity, Cold
# Rules: 1. Hard-voting filter vote_score > 0.5  2. Final score uses the multi-algorithm weight mean (continuous value)
# Pair + triple + quadruple + quintuple + all-6 combinations, automatically organized into subfolders
# Results uniformly saved to the Hard_Majority_Voting root directory
# ==============================
library(dplyr)
library(tidyr)
library(vroom)

# ------------------------------------------------
# 1. Weighted hard-voting fusion function (Scheme 2 modified version)
# ------------------------------------------------
hard_voting_fusion <- function(method_list) {
  
  process_one <- function(df) {
    df <- df[, 1:3]
    colnames(df) <- c("TF", "Target", "Weight")
    df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$Weight), ]
    df$vote <- 1          # Detection marker
    return(df)
  }
  
  processed_list <- lapply(method_list, process_one)
  all_edges <- bind_rows(lapply(processed_list, function(x) x[, 1:2])) %>% distinct()
  n_methods <- length(processed_list)
  
  # Append voting columns + original weight columns
  for (i in seq_along(processed_list)) {
    current <- processed_list[[i]]
    key <- interaction(all_edges$TF, all_edges$Target)
    curr_key <- interaction(current$TF, current$Target)
    
    # 0/1 voting column
    vote_col <- paste0("v", i)
    all_edges[[vote_col]] <- current$vote[match(key, curr_key)]
    all_edges[[vote_col]][is.na(all_edges[[vote_col]])] <- 0
    
    # Original weight column, filled with 0 when not detected
    weight_col <- paste0("w", i)
    all_edges[[weight_col]] <- current$Weight[match(key, curr_key)]
    all_edges[[weight_col]][is.na(all_edges[[weight_col]])] <- 0
  }
  
  # Compute the voting score (original filtering rule unchanged)
  vote_cols <- paste0("v", seq_along(processed_list))
  all_edges$total_vote <- rowSums(all_edges[, vote_cols])
  all_edges$vote_score <- all_edges$total_vote / n_methods
  
  # Compute the weight mean as a continuous confidence score (core optimization)
  weight_cols <- paste0("w", seq_along(processed_list))
  all_edges$mean_weight <- rowMeans(all_edges[, weight_cols])
  
  # Majority-voting filter: vote_score > 0.5, sorted by fused weight in descending order
  result <- all_edges %>%
    filter(vote_score > 0.5) %>%
    arrange(desc(mean_weight)) %>%
    select(TF, Target, EdgeWeight = mean_weight)
  
  return(result)
}

# ------------------------------------------------
# 2. Read and clean function
# ------------------------------------------------
clean_read <- function(path, delim = ",") {
  df <- vroom(path, delim = delim, show_col_types = FALSE)
  df <- df[, 1:3]
  colnames(df) <- c("TF", "Target", "Weight")
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$Weight), ]
  return(df)
}

# ======================================================
# Batch processing: the three abiotic-stress datasets Aridity / Alkalinity / Cold
# ======================================================
datasets <- c("Aridity", "Alkalinity", "Cold")

for(ds in datasets){
  
  cat("\n=============================================\n")
  cat("Processing dataset: ", ds, " (weighted hard voting + significant-edge filtering)\n")
  cat("=============================================\n")
  
  # Read the results of the 6 algorithms
  genie3     <- clean_read(paste0("GENIE3/", ds, "_significant_edges.txt"), delim = "\t")
  kboost     <- clean_read(paste0("Kboost/", ds, "_significant_edges.csv"))
  grnboost2  <- clean_read(paste0("GRNBoost/", ds, "_significant_edges.csv"))
  dcema      <- clean_read(paste0("3DCEMA/", ds, "_significant_edges.csv"))
  deeprig    <- clean_read(paste0("DeepRIG/", ds, "_significant_edges.csv"))
  ige        <- clean_read(paste0("IGEGRNs/", ds, "_significant_edges.csv"))
  
  all_methods <- list(
    GENIE3    = genie3,
    Kboost    = kboost,
    GRNBoost2 = grnboost2,
    DCEMA     = dcema,
    DeepRIG   = deeprig,
    IGE       = ige
  )
  
  # Top-level root folder: Hard_Majority_Voting
  top_root <- "Hard_Majority_Voting"
  # Dataset subdirectory
  ds_root <- paste0("Voting_Sig_Results_", ds, "_0-1")
  # Concatenate the full path
  out_root <- file.path(top_root, ds_root)
  
  dir_pair    <- file.path(out_root, "1_Pair_2Methods")
  dir_triple  <- file.path(out_root, "2_Triple_3Methods")
  dir_4way    <- file.path(out_root, "3_Four_4Methods")
  dir_5way    <- file.path(out_root, "4_Five_5Methods")
  dir_all6    <- file.path(out_root, "5_All_6Methods")
  
  # Batch-create multi-level folders
  lapply(c(out_root, dir_pair, dir_triple, dir_4way, dir_5way, dir_all6), function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  })
  
  cat("Starting weighted hard-voting fusion and significant-edge filtering...\n")
  
  # 1. Pair combinations
  comb2 <- combn(names(all_methods), 2, simplify = FALSE)
  for (nm in comb2) {
    res <- hard_voting_fusion(all_methods[nm])
    write.csv(res, file.path(dir_pair, paste0("Voting_Sig_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # 2. Triple combinations
  comb3 <- combn(names(all_methods), 3, simplify = FALSE)
  for (nm in comb3) {
    res <- hard_voting_fusion(all_methods[nm])
    write.csv(res, file.path(dir_triple, paste0("Voting_Sig_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # 3. Quadruple combinations
  comb4 <- combn(names(all_methods), 4, simplify = FALSE)
  for (nm in comb4) {
    res <- hard_voting_fusion(all_methods[nm])
    write.csv(res, file.path(dir_4way, paste0("Voting_Sig_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # 4. Quintuple combinations
  comb5 <- combn(names(all_methods), 5, simplify = FALSE)
  for (nm in comb5) {
    res <- hard_voting_fusion(all_methods[nm])
    write.csv(res, file.path(dir_5way, paste0("Voting_Sig_", paste(nm, collapse = "+"), ".csv")), row.names = FALSE)
  }
  
  # 5. All 6 algorithms
  all6 <- hard_voting_fusion(all_methods)
  write.csv(all6, file.path(dir_all6, "Voting_Sig_All6Methods.csv"), row.names = FALSE)
  
  cat("\n[Done]", ds, " processed; weighted-voting significant edges output\n")
}

cat("\n[Finished] Fusion processing for Aridity / Alkalinity / Cold completed!\n")
