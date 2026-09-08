# Load required packages
library(dplyr)
library(tidyr)
library(vroom)
library(tidyverse)
library(VennDiagram)
library(readxl)

# ------------------------------------------------------------------------------
# New: read the rice reference TF list Osj_TF_list.txt, TF names are in the second column
# ------------------------------------------------------------------------------
# Note: adjust delim according to your file separator; commonly tab "\t", comma ",", or space " "
# First try tab; if it fails, modify the delim parameter yourself
tf_ref <- vroom("Osj_TF_list.txt", delim = "\t", show_col_types = FALSE)
# Extract the second column as the standard TF-name vector
valid_tf_pool <- pull(tf_ref, 2) %>% unique()

cat("Total number of TFs in the reference TF library Osj_TF_list.txt: ", length(valid_tf_pool), "\n")

# ------------------------------------------------------------------------------
# 1. Data reading and preprocessing (unified format, deduplication)
# ------------------------------------------------------------------------------
load_grn_data <- function(file_path) {
  vroom(file_path, show_col_types = FALSE)[, 1:3] %>%
    rename(TF = TF, Target = Target, EdgeWeight = XGB_Score) %>%
    drop_na() %>%
    distinct(TF, Target, .keep_all = TRUE)
}

# Read the three abiotic-stress XGB fusion network files
Z_ll_Alkalinity <- load_grn_data("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv")
Z_ll_Aridity <- load_grn_data("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv")
Z_ll_Cold <- load_grn_data("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv")

# ------------------------------------------------------------------------------
# New function: filter the GRN by the reference TF list, keep only real TFs, and output validation info
# ------------------------------------------------------------------------------
filter_grn_by_tf_ref <- function(grn_df, stress_name){
  all_tf_in_net <- distinct(grn_df, TF) %>% pull(TF)
  # Matching
  true_tf <- intersect(all_tf_in_net, valid_tf_pool)
  fake_tf <- setdiff(all_tf_in_net, valid_tf_pool)
  
  cat("==================== ", stress_name, " ====================\n")
  cat("Number of independent TFs in the raw network: ", length(all_tf_in_net), "\n")
  cat("Number of valid TFs matching Osj_TF_list: ", length(true_tf), "\n")
  cat("Number of TFs NOT in Osj_TF_list: ", length(fake_tf), "\n")
  
  if(length(fake_tf) > 0){
    cat("[Unmatched TF list]\n")
    print(fake_tf)
    # Export the unmatched TFs to a text file for later checking of ID-naming differences
    writeLines(fake_tf, con = paste0(stress_name,"_unmatched_TF.txt"))
  }
  
  # Filter the network to keep only real TFs in the reference library
  grn_filtered <- grn_df %>% filter(TF %in% valid_tf_pool)
  return(grn_filtered)
}

# Perform the filtering (key! all subsequent degree computations and key-TF screening use the filtered network)
Z_ll_Alkalinity_filtered <- filter_grn_by_tf_ref(Z_ll_Alkalinity, "Alkalinity")
Z_ll_Aridity_filtered   <- filter_grn_by_tf_ref(Z_ll_Aridity, "Aridity")
Z_ll_Cold_filtered      <- filter_grn_by_tf_ref(Z_ll_Cold, "Cold")

# ------------------------------------------------------------------------------
# 2. Compute the total TFs, key TFs, and non-key TFs (all using filtered data!)
# ------------------------------------------------------------------------------
# 2.1 Extract all valid TFs under each stress
tf_alk_all <- Z_ll_Alkalinity_filtered %>% distinct(TF) %>% pull(TF)
tf_ari_all <- Z_ll_Aridity_filtered %>% distinct(TF) %>% pull(TF)
tf_cold_all <- Z_ll_Cold_filtered %>% distinct(TF) %>% pull(TF)

# 2.2 Degree-centrality computation
cal_deg <- function(link_list) {
  link_list %>% 
    group_by(TF) %>% 
    summarise(degree = n(), .groups = "drop") %>%
    arrange(desc(degree))
}

deg_alk <- cal_deg(Z_ll_Alkalinity_filtered)
deg_ari <- cal_deg(Z_ll_Aridity_filtered)
deg_cold <- cal_deg(Z_ll_Cold_filtered)

# 2.3 Poisson-distribution screening of key TFs
find_key_TFs <- function(deg_data) {
  if (nrow(deg_data) == 0) {
    warning("Degree-centrality data is empty; returning an empty list")
    return(character(0))
  }
  mean_deg <- mean(deg_data$degree, na.rm = TRUE)
  deg_data$p_val <- ppois(deg_data$degree, lambda = mean_deg, lower.tail = FALSE)
  deg_data$adj_p_val <- p.adjust(deg_data$p_val, method = "bonferroni")
  deg_data %>% filter(adj_p_val < 0.05) %>% pull(TF)
}

# Extract the key TFs under each abiotic stress (based on the real TF subset)
tf_alk_key <- find_key_TFs(deg_alk)
tf_ari_key <- find_key_TFs(deg_ari)
tf_cold_key <- find_key_TFs(deg_cold)

# 2.4 Summarize the plotting data
combined_data <- bind_rows(
  data.frame(
    Stress = "Alkalinity",
    TF_Category = c("Key TFs", "Non-key TFs"),
    Count = c(length(tf_alk_key), length(tf_alk_all) - length(tf_alk_key))
  ),
  data.frame(
    Stress = "Aridity",
    TF_Category = c("Key TFs", "Non-key TFs"),
    Count = c(length(tf_ari_key), length(tf_ari_all) - length(tf_ari_key))
  ),
  data.frame(
    Stress = "Cold",
    TF_Category = c("Key TFs", "Non-key TFs"),
    Count = c(length(tf_cold_key), length(tf_cold_all) - length(tf_cold_key))
  )
)

print("Summary table of TF statistics per stress after filtering:")
print(combined_data)

# ------------------------------------------------------------------------------
# 3. Draw the stacked bar chart
# ------------------------------------------------------------------------------
key_color <- "#DC143C"
non_key_color <- "#2E86AB"

p_combined <- ggplot(combined_data, aes(
  x = Stress,
  y = Count,
  fill = TF_Category,
  label = Count
)) +
  geom_col(width = 0.4, color = "black", alpha = 0.9) +
  geom_text(
    position = position_stack(vjust = 0.5),
    size = 4.2, 
    fontface = "bold", 
    color = "white"
  ) +
  scale_fill_manual(values = c("Key TFs" = key_color, "Non-key TFs" = non_key_color)) +
  labs(
    title = "Transcription Factors Distribution Across Abiotic Stresses",
    x = "Abiotic Stress Type",
    y = "Number of Transcription Factors",
    fill = "TF Category"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  ylim(0, max(aggregate(Count ~ Stress, combined_data, sum)$Count) * 1.1)

ggsave(
  filename = "Abiotic_TFs_Combined_Stacked_FilteredByOsjTF.pdf",
  plot = p_combined,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white",
  device = "pdf"
)

cat("==============================================\n")
cat("Analysis completed!\n")
cat("1. TFs not in Osj_TF_list have been exported separately as *unmatched_TF.txt\n")
cat("2. The plot uses networks filtered/validated by the reference TF list\n")
cat("3. Figure file: Abiotic_TFs_Combined_Stacked_FilteredByOsjTF.pdf\n")

#=========================================================================
# Extract the stress-specific key TFs and export them as txt
#=========================================================================
# The three key-TF sets
key_tf_sets <- list(
  Alkalinity = tf_alk_key,
  Aridity    = tf_ari_key,
  Cold       = tf_cold_key
)

A <- key_tf_sets$Alkalinity
B <- key_tf_sets$Aridity
C <- key_tf_sets$Cold

# Key TFs unique to each stress
Only_Alk_KeyTF <- sort(setdiff(A, union(B, C)))
Only_Ari_KeyTF <- sort(setdiff(B, union(A, C)))
Only_Cold_KeyTF <- sort(setdiff(C, union(A, B)))

# Export as text
writeLines(Only_Alk_KeyTF, con = "Only_Alkalinity_KeyTF.txt")
writeLines(Only_Ari_KeyTF, con = "Only_Aridity_KeyTF.txt")
writeLines(Only_Cold_KeyTF, con = "Only_Cold_KeyTF.txt")

# Console output of the counts
cat("\n==================== Stress-specific key TF counts ====================\n")
cat("Key TFs unique to Alkalinity: ", length(Only_Alk_KeyTF), "\n")
cat("Key TFs unique to Aridity: ", length(Only_Ari_KeyTF), "\n")
cat("Key TFs unique to Cold: ", length(Only_Cold_KeyTF), "\n")
