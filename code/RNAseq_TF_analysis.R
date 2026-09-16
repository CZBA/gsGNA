# Load required packages
library(tidyverse)
library(dplyr)
library(readxl)
library(GenomicRanges)
library(ggplot2)
library(openxlsx)
library(edgeR)
library(org.Osativa.eg.db)
library(VennDiagram)
library(gridExtra)

#==================== Custom utility functions ====================
#' Read excel file and select specified columns
#' @param file_path Path to input excel file
#' @param cols Vector of column indices or column names to keep
#' @return Tibble with selected columns
read_excel_select <- function(file_path, cols) {
  read_excel(file_path) %>% dplyr::select(all_of(cols))
}

#' Remove rows where all numeric columns equal zero
#' @param data Input expression tibble
#' @return Filtered tibble without all‑zero rows
remove_zero_rows <- function(data) {
  data %>% filter(rowSums(dplyr::select(., where(is.numeric))) != 0)
}

#' Filter input data to retain only transcription factors
#' @param data Input gene expression data
#' @param gene_col Character string of gene ID column name
#' @param tf_list Tibble containing TF Gene_ID column
#' @return Subsetted tibble containing only TF genes
filter_TF <- function(data, gene_col, tf_list) {
  data %>% filter(!!sym(gene_col) %in% tf_list$Gene_ID)
}

#' Convert rice gene ID using org.Osativa.eg.db
#' @param id_data Single‑column tibble containing input gene IDs
#' @param from_type Source keytype in org.Osativa.eg.db, e.g. "RAP"
#' @param to_type Target keytype in org.Osativa.eg.db, e.g. "GID"
#' @return Tibble with converted IDs; keep original ID when conversion returns NA
convert_and_replace_id <- function(id_data, from_type = "RAP", to_type = "GID") {
  gene_ids <- pull(id_data, 1)
  original_colname <- colnames(id_data)[1]
  
  if (length(gene_ids) == 0) {
    stop("Input ID data is empty, please check your input!")
  }
  
  converted_ids <- mapIds(
    org.Osativa.eg.db,
    keys = gene_ids,
    keytype = from_type,
    column = to_type,
    multiVals = "first"
  )
  
  # Retain original ID for NA conversion results
  result <- tibble(
    !!original_colname := ifelse(is.na(converted_ids), gene_ids, converted_ids)
  )
  return(result)
}

#' Draw 3‑set Venn diagram and export to PDF
#' @param id_list Named list of gene ID vectors for Venn comparison
#' @param title Plot main title
#' @param fill_colors Vector of fill colors for each set
#' @param filename Output PDF file name (without .pdf suffix)
#' @return Venn diagram gList object
draw_venn <- function(id_list, title, fill_colors, filename) {
  venn_plot <- venn.diagram(
    x = id_list,
    filename = NULL,
    col = "black",
    fill = fill_colors,
    alpha = 0.6,
    label.col = "black",
    cex = 1.2,
    fontfamily = "serif",
    cat.col = "black",
    cat.cex = 1.2,
    cat.fontfamily = "serif",
    main = title,
    main.cex = 1.5,
    main.fontfamily = "serif"
  )
  
  pdf(paste0(filename, ".pdf"), width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
  
  return(venn_plot)
}

#' Extract unique‑only genes for three‑set comparison
#' @param set_list List containing three gene ID vectors
#' @return List of genes exclusively present in each single set
get_venn_unique_three <- function(set_list){
  A <- set_list[[1]]
  B <- set_list[[2]]
  C <- set_list[[3]]
  
  res <- list(
    Only_Alkalinity = setdiff(A, union(B,C)),
    Only_Aridity    = setdiff(B, union(A,C)),
    Only_Cold       = setdiff(C, union(A,B))
  )
  return(res)
}

#' Export gene subsets into separate plain text files
#' @param subset_list Named list of gene ID vectors
#' @param prefix Output file name prefix
export_unique_txt <- function(subset_list, prefix){
  purrr::walk(names(subset_list), function(nm){
    gene_vec <- sort(subset_list[[nm]])
    writeLines(gene_vec, con = paste0(prefix,"_",nm,".txt"))
  })
}

#==================== 1. Import expression datasets ====================
# Alkalinity stress dataset
cpm_Alkalinity <- read_excel("GSE104928/GSE104928.xlsx") %>% remove_zero_rows()

# Aridity stress dataset
cpm_Aridity <- read_excel_select("GSE121303_Processed_data.xlsx", c(2, 21,22,23)) %>% remove_zero_rows()
cpm_Aridity <- cpm_Aridity %>% distinct(Gene_ID, .keep_all = TRUE)

# Cold stress dataset
cpm_Cold <- read_excel_select("GSE112547.xlsx", 1:6) %>% remove_zero_rows()

# Import rice transcription factor reference list
rice_TF <- read_tsv("Osj_TF_list.txt", col_names = TRUE) %>% dplyr::select(Gene_ID)

#==================== 2. Gene ID conversion (RAP → GID for aridity dataset) ====================
# Extract first ID column and perform ID conversion for aridity data
cpm_Aridity_col1 <- cpm_Aridity %>% dplyr::select(1)
cpm_Aridity_col1_replaced <- convert_and_replace_id(
  id_data = cpm_Aridity_col1,
  from_type = "RAP",
  to_type = "GID"
)

# Alkalinity and cold datasets do not require ID conversion
cpm_Alkalinity_col1 <- cpm_Alkalinity %>% dplyr::select(1)
cpm_Cold_col1       <- cpm_Cold %>% dplyr::select(1)

#==================== 3. Subset transcription factor genes ====================
# Filter TF genes for each stress condition
TF_Alkalinity <- filter_TF(cpm_Alkalinity, "geneID", rice_TF)
TF_Aridity    <- filter_TF(cpm_Aridity_col1_replaced, "Gene_ID", rice_TF)
TF_Cold       <- filter_TF(cpm_Cold, "Gene_ID", rice_TF)

# Extract ID‑only tibble for downstream set operations
TF_Alkalinity_col1 <- TF_Alkalinity %>% dplyr::select(1)
TF_Aridity_col1_replaced <- TF_Aridity %>% dplyr::select(1)
TF_Cold_col1       <- TF_Cold %>% dplyr::select(1)

#==================== 4. Preview ID conversion results ====================
cat("===== ID conversion preview =====\n")
cat("Original aridity IDs (top 5 rows):\n")
print(head(cpm_Aridity_col1, 5))
cat("\nConverted aridity IDs (top 5 rows):\n")
print(head(cpm_Aridity_col1_replaced, 5))

#==================== 5. Construct unique ID sets for Venn diagram ====================
## All expressed genes from RNA‑seq
rna_alkalinity_ids <- pull(cpm_Alkalinity_col1, 1) %>% unique()
rna_aridity_ids    <- pull(cpm_Aridity_col1_replaced, 1) %>% unique()
rna_cold_ids       <- pull(cpm_Cold_col1, 1) %>% unique()

## Transcription factor genes
tf_alkalinity_ids <- pull(TF_Alkalinity_col1, 1) %>% unique()
tf_aridity_ids    <- pull(TF_Aridity_col1_replaced, 1) %>% unique()
tf_cold_ids       <- pull(TF_Cold_col1, 1) %>% unique()

#==================== 6. Plot 3‑set Venn diagrams ====================
# Venn for all RNA‑seq detected genes
rna_colors <- c("#0072B2", "#D55E00", "#009E73")
rna_id_list <- list(
  Alkalinity = rna_alkalinity_ids,
  Aridity    = rna_aridity_ids,
  Cold       = rna_cold_ids
)
rna_venn <- draw_venn(
  id_list = rna_id_list,
  title = "RNA‑seq Data (Converted IDs)",
  fill_colors = rna_colors,
  filename = "RNA_Data_Venn"
)

# Venn for transcription factor genes
tf_colors <- c("#CC79A7", "#F0E442", "#56B4E9")
tf_id_list <- list(
  Alkalinity = tf_alkalinity_ids,
  Aridity    = tf_aridity_ids,
  Cold       = tf_cold_ids
)
tf_venn <- draw_venn(
  id_list = tf_id_list,
  title = "TF Data (Converted IDs)",
  fill_colors = tf_colors,
  filename = "TF_Data_Venn"
)

#==================== 7. Summary statistics of gene counts ====================
gene_counts <- tibble(
  Group = c("Alkalinity", "Aridity", "Cold"),
  RNA_Count = c(
    length(rna_alkalinity_ids),
    length(rna_aridity_ids),
    length(rna_cold_ids)
  ),
  TF_Count = c(
    length(tf_alkalinity_ids),
    length(tf_aridity_ids),
    length(tf_cold_ids)
  )
)

cat("\n===== Gene count summary =====\n")
print(gene_counts)

#==================== 8. ggplot visualization ====================
## Grouped bar plot: RNA‑seq vs TF gene count
gene_counts_long <- gene_counts %>%
  pivot_longer(
    cols = c(RNA_Count, TF_Count),
    names_to = "DataType",
    values_to = "GeneNumber"
  ) %>%
  mutate(
    DataType = factor(DataType, levels = c("RNA_Count", "TF_Count"), labels = c("RNA‑seq", "TF")),
    Group    = factor(Group, levels = c("Alkalinity", "Aridity", "Cold"))
  )

bar_plot <- ggplot(gene_counts_long, aes(x = Group, y = GeneNumber, fill = DataType)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("#0072B2", "#CC79A7")) +
  geom_text(
    aes(label = GeneNumber),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = "Stress Treatment",
    y = "Number of Genes",
    fill = "Data Type",
    title = "Comparison of Gene Numbers (RNA‑seq vs TF)",
    subtitle = "After ID Conversion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

ggsave(
  filename = "RNA_TF_GeneCount_BarPlot.pdf",
  plot = bar_plot,
  width = 10,
  height = 7,
  dpi = 300
)

## Bar plot for TF proportion
gene_counts_ratio <- gene_counts %>%
  mutate(
    TF_Ratio = (TF_Count / RNA_Count) * 100,
    Group = factor(Group, levels = c("Alkalinity", "Aridity", "Cold"))
  )

ratio_plot <- ggplot(gene_counts_ratio, aes(x = Group, y = TF_Ratio, fill = Group)) +
  geom_col(width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  geom_text(
    aes(label = sprintf("%.1f%%", TF_Ratio)),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = "Stress Treatment",
    y = "TF Genes Ratio (%)",
    title = "TF Genes Proportion in RNA‑seq Data",
    subtitle = "After ID Conversion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave(
  filename = "TF_Ratio_BarPlot.pdf",
  plot = ratio_plot,
  width = 8,
  height = 7,
  dpi = 300
)

## Arrange two bar plots into one combined figure
combined_plot <- grid.arrange(bar_plot, ratio_plot, ncol = 2)
ggsave(
  filename = "RNA_TF_Combined_BarPlots.pdf",
  plot = combined_plot,
  width = 18,
  height = 7,
  dpi = 300
)

print(bar_plot)
print(ratio_plot)

#==================== 9. Selective environment cleanup ====================
keep_vars <- c(
  "read_excel_select", "remove_zero_rows", "filter_TF", "convert_and_replace_id", "draw_venn",
  "get_venn_unique_three", "export_unique_txt",
  "rna_venn", "tf_venn", "bar_plot", "ratio_plot", "combined_plot",
  "gene_counts", "gene_counts_ratio"
)
all_vars <- ls(envir = globalenv())
remove_vars <- setdiff(all_vars, keep_vars)
rm(list = remove_vars)

#==================== Export stress‑specific unique TF gene lists ====================
tf_unique_subsets <- get_venn_unique_three(tf_id_list)
export_unique_txt(tf_unique_subsets, prefix = "TF_Genes")

cat("\n===== Count of stress‑exclusive TF genes =====\n")
print(purrr::map_int(tf_unique_subsets, length))
