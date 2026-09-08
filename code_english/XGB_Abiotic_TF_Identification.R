library(dplyr)
library(tidyr)
library(readr)
library(vroom)
library(igraph)
library(tidyverse)
library(gridExtra)
library(VennDiagram)

# ========== 2. Rewrite the data-reading function: adapted to XGB-fusion CSVs with header (TF,Target,XGB_Score) ==========
# New data source: your three XGB-fusion GRN CSVs, which carry the header TF/Target/XGB_Score
load_grn_data <- function(file_path) {
  # Read the fusion regulatory edge table with header
  raw_data <- vroom(
    file_path, 
    show_col_types = FALSE
  )
  
  # Rename to unified column names, compatible with all downstream analysis functions
  result <- raw_data %>%
    rename(TF = TF, Target = Target, EdgeWeight = XGB_Score) %>%
    drop_na() %>%                                    # Remove rows containing NA
    distinct(TF, Target, .keep_all = TRUE)           # Deduplicate TF-Target pairs, keep weights
  
  # Validate data
  if(nrow(result) == 0) {
    warning(paste("File", file_path, "has no valid data after reading"))
  }
  
  return(result)
}

# ========== 3. Read the new XGB-fusion GRN files (replacing the original significant_edges.csv) ==========
Z_ll_Alkalinity <- load_grn_data("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv")
Z_ll_Aridity <- load_grn_data("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv")
Z_ll_Cold <- load_grn_data("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv")

# Verify the data-reading results
cat("=== XGB-fusion GRN data reading verification ===\n")
cat("Alkalinity fusion network rows: ", nrow(Z_ll_Alkalinity), "\n")
cat("Alkalinity column names: ", paste(colnames(Z_ll_Alkalinity), collapse = ", "), "\n\n")

# Count the number of TFs under each stress
alkalinity_tf_count <- Z_ll_Alkalinity %>% distinct(TF) %>% nrow()
aridity_tf_count <- Z_ll_Aridity %>% distinct(TF) %>% nrow()
cold_tf_count <- Z_ll_Cold %>% distinct(TF) %>% nrow()

cat("=== TF count statistics across stress fusion networks ===\n")
cat("Total TFs under Alkalinity stress: ", alkalinity_tf_count, "\n")
cat("Total TFs under Aridity stress: ", aridity_tf_count, "\n")
cat("Total TFs under Cold stress: ", cold_tf_count, "\n\n")

# ========== 4. Degree centrality calculation & plotting functions (unchanged, fully compatible with new data) ==========
cal_deg <- function(link_list) {
  link_list %>% 
    group_by(TF) %>% 
    summarise(degree = n(), .groups = "drop") %>% 
    arrange(desc(degree))
}

# Plot degree-centrality bar chart
plot_deg_centrality <- function(deg_data, stress_name) {
  ggplot(data = deg_data) +
    geom_col(mapping = aes(x = reorder(TF, -degree), y = degree), fill = "#2E86AB") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(
      x = paste("TFs in", stress_name), 
      y = "Degree (Number of target genes)",
      title = paste("TF Degree Centrality |", stress_name, "XGB Fusion GRN")
    ) +
    coord_cartesian(ylim = c(0, 20000))
}

# Compute degree centrality under each stress
deg_Aridity  <- cal_deg(Z_ll_Aridity)
deg_Alkalinity <- cal_deg(Z_ll_Alkalinity)
deg_Cold <- cal_deg(Z_ll_Cold)

# Plot
p_Aridity <- plot_deg_centrality(deg_Aridity, "Aridity")
p_Alkalinity <- plot_deg_centrality(deg_Alkalinity, "Alkalinity")
p_Cold <- plot_deg_centrality(deg_Cold, "Cold")

# Combine and save (file name includes the XGBFusion marker to distinguish old results)
combined_plot <- grid.arrange(p_Aridity, p_Alkalinity, p_Cold, nrow = 1)
ggsave("XGBFusion_degree_centrality_plots.pdf",
       combined_plot, 
       width = 15, 
       height = 5,
       device = "pdf",
       dpi = 300,
       units = "in")

# ========== 5. Poisson-distribution screening of hub TFs (function unchanged, compatible with new data) ==========
find_key_TFs <- function(deg_data) {
  if(nrow(deg_data) == 0) {
    warning("Input data is empty; returning an empty data frame")
    return(data.frame(TF = character(), degree = integer(), p_value = numeric(), adj_p_value = numeric()))
  }
  
  mean_degree <- mean(deg_data$degree, na.rm = TRUE)
  deg_data$p_value <- ppois(deg_data$degree, lambda = mean_degree, lower.tail = FALSE)
  deg_data$adj_p_value <- p.adjust(deg_data$p_value, method = "bonferroni")
  key_TFs <- deg_data %>% filter(adj_p_value < 0.05)
  return(key_TFs)
}

# Screen key hub TFs
key_TFs_Aridity <- find_key_TFs(deg_Aridity)
key_TFs_Alkalinity <- find_key_TFs(deg_Alkalinity)
key_TFs_Cold <- find_key_TFs(deg_Cold)

# Extract TF-name vectors
P_key_TFs_Aridity <- key_TFs_Aridity$TF
P_key_TFs_Alkalinity <- key_TFs_Alkalinity$TF
P_key_TFs_Cold <- key_TFs_Cold$TF

# Compute shared TFs between conditions
P_Aridity_Alkalinity <- intersect(P_key_TFs_Aridity, P_key_TFs_Alkalinity)
P_Aridity_Cold <- intersect(P_key_TFs_Aridity, P_key_TFs_Cold)
P_Alkalinity_Cold <- intersect(P_key_TFs_Alkalinity, P_key_TFs_Cold)
P_all <- intersect(P_Aridity_Alkalinity, P_key_TFs_Cold)

# ========== 6. Venn diagram (file names replaced with fusion versions) ==========
venn_plot <- venn.diagram(
  x = list(P_key_TFs_Aridity, P_key_TFs_Alkalinity, P_key_TFs_Cold),
  category.names = c("Aridity", "Alkalinity", "Cold"),
  filename = NULL,
  col = "transparent",
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.5,
  cex = 0.8,
  fontfamily = "sans",
  cat.cex = 0.8,
  cat.fontfamily = "sans"
)

pdf("XGBFusion_poisson_venn.pdf", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()

# ========== 7. Output TF list files (all file names carry the XGBFusion prefix to distinguish old files) ==========
write.table(P_key_TFs_Aridity, file = "XGBFusion_P_key_TFs_Aridity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_key_TFs_Alkalinity, file = "XGBFusion_P_key_TFs_Alkalinity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_key_TFs_Cold, file = "XGBFusion_P_key_TFs_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Output shared TFs
write.table(P_Aridity_Alkalinity, file = "XGBFusion_P_Aridity_Alkalinity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_Aridity_Cold, file = "XGBFusion_P_Aridity_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_Alkalinity_Cold, file = "XGBFusion_P_Alkalinity_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_all, file = "XGBFusion_P_all_common_TF.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ========== 8. Extract key TF-Target regulatory edges (with XGB weights) ==========
extract_key_edges <- function(link_list, key_TFs) {
  link_list %>% 
    filter(TF %in% key_TFs)
}

get_key_tf_targets <- function(link_list, key_tfs) {
  link_list %>% 
    filter(TF %in% key_tfs) %>%
    dplyr::select(TF, Target, EdgeWeight) %>%
    arrange(TF, desc(EdgeWeight))
}

# Regulatory relationships of single-stress specific hub TFs
key_tf_targets_Aridity <- get_key_tf_targets(Z_ll_Aridity, P_key_TFs_Aridity)
key_tf_targets_Alkalinity <- get_key_tf_targets(Z_ll_Alkalinity, P_key_TFs_Alkalinity)
key_tf_targets_Cold <- get_key_tf_targets(Z_ll_Cold, P_key_TFs_Cold)

# Regulatory relationships of core TFs shared across all three stresses
common_tf_targets_Aridity <- get_key_tf_targets(Z_ll_Aridity, P_all)
common_tf_targets_Alkalinity <- get_key_tf_targets(Z_ll_Alkalinity, P_all)
common_tf_targets_Cold <- get_key_tf_targets(Z_ll_Cold, P_all)

# Output TSV results (with fusion marker)
write_tsv(key_tf_targets_Aridity, "XGBFusion_key_TF_targets_Aridity.tsv")
write_tsv(key_tf_targets_Alkalinity, "XGBFusion_key_TF_targets_Alkalinity.tsv")
write_tsv(key_tf_targets_Cold, "XGBFusion_key_TF_targets_Cold.tsv")

write_tsv(common_tf_targets_Aridity, "XGBFusion_common_TF_targets_Aridity.tsv")
write_tsv(common_tf_targets_Alkalinity, "XGBFusion_common_TF_targets_Alkalinity.tsv")
write_tsv(common_tf_targets_Cold, "XGBFusion_common_TF_targets_Cold.tsv")

# Print statistics
cat("=== XGB fusion network: number of key TF-target regulatory pairs ===\n")
cat("Aridity: ", nrow(key_tf_targets_Aridity), "\n")
cat("Alkalinity: ", nrow(key_tf_targets_Alkalinity), "\n")
cat("Cold: ", nrow(key_tf_targets_Cold), "\n")

# 1. Extract target genes of the common key TFs under each stress (with enhanced deduplication and empty-value removal)
# Function: get all target genes of a specific TF set under a given stress (deduplicate + remove empty values)
get_tf_targets_unique <- function(link_list, tf_set) {
  link_list %>% 
    filter(TF %in% tf_set) %>%  # Filter common key TFs
    pull(Target) %>%            # Extract target genes
    unique() %>%                # Deduplicate: keep unique target genes
    na.omit() %>%               # Remove NA values (empty values)
    .[. != ""]                  # Remove empty strings (avoid invalid values like "")
}

# Extract target genes of the common key TFs under the three stresses (deduplication and empty-value removal included)
common_tf_targets_Aridity_unique <- get_tf_targets_unique(Z_ll_Aridity, P_all)
common_tf_targets_Alkalinity_unique <- get_tf_targets_unique(Z_ll_Alkalinity, P_all)
common_tf_targets_Cold_unique <- get_tf_targets_unique(Z_ll_Cold, P_all)


# 2. Re-check and clean the target-gene lists (double safeguard)
# Remove possible remaining empty or invalid values
clean_targets <- function(targets) {
  targets %>% 
    unique() %>% 
    na.omit() %>% 
    .[. != ""]
}

# Apply a second cleaning to the three lists
common_tf_targets_Aridity_unique <- clean_targets(common_tf_targets_Aridity_unique)
common_tf_targets_Alkalinity_unique <- clean_targets(common_tf_targets_Alkalinity_unique)
common_tf_targets_Cold_unique <- clean_targets(common_tf_targets_Cold_unique)


# 3. Compute shared and specific target genes (based on cleaned data)
# Target genes shared by all three conditions
targets_all_common <- intersect(
  intersect(common_tf_targets_Aridity_unique, common_tf_targets_Alkalinity_unique),
  common_tf_targets_Cold_unique
) %>% clean_targets()  # Clean the result again

# Target genes shared only by two conditions
targets_Aridity_Alkalinity_only <- setdiff(
  intersect(common_tf_targets_Aridity_unique, common_tf_targets_Alkalinity_unique),
  targets_all_common
) %>% clean_targets()

targets_Aridity_Cold_only <- setdiff(
  intersect(common_tf_targets_Aridity_unique, common_tf_targets_Cold_unique),
  targets_all_common
) %>% clean_targets()

targets_Alkalinity_Cold_only <- setdiff(
  intersect(common_tf_targets_Alkalinity_unique, common_tf_targets_Cold_unique),
  targets_all_common
) %>% clean_targets()

# Target genes specific to each stress
targets_Aridity_unique <- setdiff(
  common_tf_targets_Aridity_unique,
  union(common_tf_targets_Alkalinity_unique, common_tf_targets_Cold_unique)
) %>% clean_targets()

targets_Alkalinity_unique <- setdiff(
  common_tf_targets_Alkalinity_unique,
  union(common_tf_targets_Aridity_unique, common_tf_targets_Cold_unique)
) %>% clean_targets()

targets_Cold_unique <- setdiff(
  common_tf_targets_Cold_unique,
  union(common_tf_targets_Aridity_unique, common_tf_targets_Alkalinity_unique)
) %>% clean_targets()


# 4. Output results to files (ensure no empty files)
if (!dir.exists("common_TF_targets_analysis")) {
  dir.create("common_TF_targets_analysis")
}

# Check for emptiness before output; write a hint message when empty (avoid empty files)
write_clean_targets <- function(targets, file_path) {
  if (length(targets) == 0) {
    writeLines("No valid targets found", file_path)  # Write a hint when the list is empty
  } else {
    writeLines(targets, file_path)
  }
}

# Output target genes shared by all three conditions
write_clean_targets(targets_all_common, 
                    "common_TF_targets_analysis/targets_all_common.txt")

# Output target genes shared by two conditions
write_clean_targets(targets_Aridity_Alkalinity_only, 
                    "common_TF_targets_analysis/targets_Aridity_Alkalinity_only.txt")
write_clean_targets(targets_Aridity_Cold_only, 
                    "common_TF_targets_analysis/targets_Aridity_Cold_only.txt")
write_clean_targets(targets_Alkalinity_Cold_only, 
                    "common_TF_targets_analysis/targets_Alkalinity_Cold_only.txt")

# Output target genes specific to each stress
write_clean_targets(targets_Aridity_unique, 
                    "common_TF_targets_analysis/targets_Aridity_unique.txt")
write_clean_targets(targets_Alkalinity_unique, 
                    "common_TF_targets_analysis/targets_Alkalinity_unique.txt")
write_clean_targets(targets_Cold_unique, 
                    "common_TF_targets_analysis/targets_Cold_unique.txt")


# 5. Count and print the results (with empty-value checks)
cat("Target-gene statistics of common key TFs (after deduplication and empty-value removal):\n")
cat("Target genes shared by all three stresses: ", length(targets_all_common), "\n")
cat("Target genes shared only by Aridity and Alkalinity: ", length(targets_Aridity_Alkalinity_only), "\n")
cat("Target genes shared only by Aridity and Cold: ", length(targets_Aridity_Cold_only), "\n")
cat("Target genes shared only by Alkalinity and Cold: ", length(targets_Alkalinity_Cold_only), "\n")
cat("Target genes specific to Aridity: ", length(targets_Aridity_unique), "\n")
cat("Target genes specific to Alkalinity: ", length(targets_Alkalinity_unique), "\n")
cat("Target genes specific to Cold: ", length(targets_Cold_unique), "\n")

# Hint message update
cat("\nAnalysis completed!\n")
cat("1. Target-gene lists saved to the common_TF_targets_analysis folder\n")
cat("2. Venn diagram saved as common_TF_targets_analysis/common_TF_targets_venn.pdf\n")
