library(tidyverse)
library(vroom) 

# ------------------------------------------------------------------------------
# 1. Data loading and preprocessing (extract TFs per stress and filter TFBS)
# ------------------------------------------------------------------------------
load_grn_data <- function(file_path) {
  vroom(file_path, show_col_types = FALSE)[, 1:3] %>%
    rename(TF = 1, Target = 2, EdgeWeight = 3) %>%
    drop_na() %>%
    distinct(TF, Target, .keep_all = TRUE)
}

# Load GRN data (note: make sure the file paths are correct)
Alkalinity_data <- load_grn_data("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv")
Aridity_data <- load_grn_data("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv")
Cold_data <- load_grn_data("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv")

# Extract the TFs of each stress separately (deduplicated)
alkalinity_tfs <- unique(Alkalinity_data$TF)
aridity_tfs <- unique(Aridity_data$TF)
cold_tfs <- unique(Cold_data$TF)

message(sprintf("-> Extracted %d unique TFs from the Alkalinity GRN", length(alkalinity_tfs)))
message(sprintf("-> Extracted %d unique TFs from the Aridity GRN", length(aridity_tfs)))
message(sprintf("-> Extracted %d unique TFs from the Cold GRN", length(cold_tfs)))

# First load the complete TFBS data
tfbs_full <- read_tsv("TFTarget.txt", col_names = TRUE) %>% 
  select(TF = 1, Target = 2) %>% 
  distinct()

# Filter the corresponding TFBS pairs by the TFs of each stress
tfbs_alkalinity <- tfbs_full %>% filter(TF %in% alkalinity_tfs)
tfbs_aridity <- tfbs_full %>% filter(TF %in% aridity_tfs)
tfbs_cold <- tfbs_full %>% filter(TF %in% cold_tfs)

# Generate the regulatory pairs of each stress's TFBS and their counts
tfbs_alkalinity_pairs <- paste(tfbs_alkalinity$TF, tfbs_alkalinity$Target, sep = "_")
tfbs_aridity_pairs <- paste(tfbs_aridity$TF, tfbs_aridity$Target, sep = "_")
tfbs_cold_pairs <- paste(tfbs_cold$TF, tfbs_cold$Target, sep = "_")

# Total TFBS validation set per stress & the number validated by the GRN
tfbs_alkalinity_total <- length(tfbs_alkalinity_pairs)
tfbs_aridity_total <- length(tfbs_aridity_pairs)
tfbs_cold_total <- length(tfbs_cold_pairs)

alkalinity_pairs <- paste(Alkalinity_data$TF, Alkalinity_data$Target, sep = "_")
aridity_pairs <- paste(Aridity_data$TF, Aridity_data$Target, sep = "_")
cold_pairs <- paste(Cold_data$TF, Cold_data$Target, sep = "_")

alkalinity_verified <- sum(alkalinity_pairs %in% tfbs_alkalinity_pairs)
aridity_verified <- sum(aridity_pairs %in% tfbs_aridity_pairs)
cold_verified <- sum(cold_pairs %in% tfbs_cold_pairs)

message(sprintf("-> Alkalinity TFBS validation-set total: %d | Number verified by the GRN: %d", tfbs_alkalinity_total, alkalinity_verified))
message(sprintf("-> Aridity TFBS validation-set total: %d | Number verified by the GRN: %d", tfbs_aridity_total, aridity_verified))
message(sprintf("-> Cold TFBS validation-set total: %d | Number verified by the GRN: %d", tfbs_cold_total, cold_verified))

# ------------------------------------------------------------------------------
# 2. Build the comparison dataset (TFBS total vs GRN verified counts)
# ------------------------------------------------------------------------------
comparison_data <- tibble(
  condition = factor(c("Alkalinity", "Aridity", "Cold"), 
                     levels = c("Alkalinity", "Aridity", "Cold")),
  tfbs_total = c(tfbs_alkalinity_total, tfbs_aridity_total, tfbs_cold_total),
  grn_verified = c(alkalinity_verified, aridity_verified, cold_verified)
) %>%
  # Reshape the data into long format, suitable for a stacked plot
  pivot_longer(-condition, names_to = "data_type", values_to = "count") %>%
  mutate(
    data_type = factor(data_type,
                       levels = c("grn_verified", "tfbs_total"),
                       labels = c("GRN Verified by TFBS", "TFBS"))  # Key modification: TFBS Validation Set Total changed to TFBS
  )

# View the comparison data
message("\n-> Comparison table of TFBS total vs GRN verified counts:")
print(comparison_data)

# ------------------------------------------------------------------------------
# 3. Draw the stacked plot (TFBS total vs GRN verified counts, red-blue color scheme + border)
# ------------------------------------------------------------------------------
# Define the red-blue color palette (labels exactly match those in the dataset to avoid mismatches)
color_palette <- c(
  "GRN Verified by TFBS" = "#DC143C",    # GRN verified count in red
  "TFBS" = "#2E86AB"                    # TFBS total in blue (corresponding to the simplified label)
)

ggplot(comparison_data, aes(x = condition, y = count, fill = data_type)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.4) +
  scale_fill_manual(values = color_palette) +  # Key: use the custom red-blue color palette
  labs(
    x = "Stress Condition",
    y = "Number of Regulatory Pairs",
    fill = "Data Category",
    title = "TFBS Validation of Gene Regulatory Networks"
  ) +
  theme_minimal() +
  # Core modification: add border + optimize axis-line styles
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    # Axis-line styles
    axis.line = element_line(color = "black", linewidth = 0.5),  
    axis.ticks = element_line(color = "black", linewidth = 0.5), 
    axis.ticks.length = unit(2, "mm"),                          
    panel.grid = element_blank(),                               
    # Key: add an overall border around the image (all four sides)
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    # Ensure the border displays completely
    plot.margin = margin(10, 10, 10, 10, "pt")
  ) +
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black",
    fontface = "bold"
  )

# Save the image as a PDF
ggsave("TFBS_Total_vs_GRN_Verified_Stacked_Barplot_RedBlue.pdf", 
       width = 6.5, height = 7, dpi = 300, bg = "white",
       device = "pdf",  # Specify PDF output
       useDingbats = FALSE)  # Avoid font-compatibility issues in PDFs
message("\n-> Red-blue stacked plot saved as PDF (with axis lines + overall border): TFBS_Total_vs_GRN_Verified_Stacked_Barplot_RedBlue.pdf")
