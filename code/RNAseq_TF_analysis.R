
# Load required packages
library(tidyverse)
library(readxl)
library(GenomicRanges)
library(ggplot2)
library(openxlsx)
library(edgeR)
library(org.Osativa.eg.db)
library(VennDiagram)
library(gridExtra)

#==================== 1. Define utility functions ====================

# Read Excel file and select specified columns
read_excel_select <- function(file_path, cols) {
  read_excel(file_path) %>%
    dplyr::select(all_of(cols))
}

# Remove rows where all numeric expression columns are zero
remove_zero_rows <- function(data) {
  
  numeric_data <- data %>%
    dplyr::select(where(is.numeric))
  
  if (ncol(numeric_data) == 0) {
    stop("No numeric expression columns were detected.")
  }
  
  data %>%
    filter(rowSums(numeric_data, na.rm = TRUE) != 0)
}

# Filter genes using a TF reference list
filter_TF <- function(data, gene_col, tf_list) {
  
  data %>%
    filter(
      !is.na(.data[[gene_col]]),
      .data[[gene_col]] %in% tf_list$Gene_ID
    )
}

# Convert rice gene IDs using org.Osativa.eg.db
convert_and_replace_id <- function(
    id_data,
    from_type = "RAP",
    to_type = "GID") {
  
  gene_ids <- as.character(pull(id_data, 1))
  original_colname <- colnames(id_data)[1]
  
  if (length(gene_ids) == 0) {
    stop("Input gene ID data is empty. Please check your input.")
  }
  
  converted_ids <- mapIds(
    org.Osativa.eg.db,
    keys = unique(gene_ids),
    keytype = from_type,
    column = to_type,
    multiVals = "first"
  )
  
  converted_ids <- converted_ids[
    match(gene_ids, names(converted_ids))
  ]
  
  result <- tibble(
    !!original_colname := ifelse(
      is.na(converted_ids),
      gene_ids,
      unname(converted_ids)
    )
  )
  
  return(result)
}

# Draw three-set Venn diagram
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
  
  pdf(
    paste0(filename, ".pdf"),
    width = 8,
    height = 8
  )
  
  grid::grid.draw(venn_plot)
  dev.off()
  
  return(venn_plot)
}

# Extract genes exclusively present in each stress condition
get_venn_unique_three <- function(set_list) {
  
  A <- set_list[[1]]
  B <- set_list[[2]]
  C <- set_list[[3]]
  
  list(
    Only_Alkaline = setdiff(A, union(B, C)),
    Only_Drought = setdiff(B, union(A, C)),
    Only_Cold = setdiff(C, union(A, B))
  )
}

# Export gene lists
export_unique_txt <- function(subset_list, prefix) {
  
  purrr::walk(names(subset_list), function(nm) {
    
    gene_vec <- sort(subset_list[[nm]])
    
    writeLines(
      gene_vec,
      con = paste0(prefix, "_", nm, ".txt")
    )
  })
}

#==================== 2. Import expression datasets ====================

# Alkaline stress dataset
expr_Alkaline <- read_excel(
  "GSE104928/GSE104928.xlsx"
) %>%
  remove_zero_rows()

# Drought stress dataset
expr_Drought <- read_excel_select(
  "GSE121303_Processed_data.xlsx",
  c(2, 21, 22, 23)
) %>%
  remove_zero_rows() %>%
  distinct(Gene_ID, .keep_all = TRUE)

# Cold stress dataset
expr_Cold <- read_excel_select(
  "GSE112547.xlsx",
  1:6
) %>%
  remove_zero_rows()

#==================== 3. Import TF reference lists ====================

# TF reference list for alkaline and cold stress
rice_TF <- read_tsv(
  "Osj_TF_list.txt",
  col_names = TRUE
) %>%
  dplyr::select(Gene_ID) %>%
  filter(
    !is.na(Gene_ID),
    Gene_ID != "None"
  ) %>%
  distinct(Gene_ID)

# TF reference list for drought stress
rice_TF_1 <- read_excel(
  "Osj_TF_list.xlsx"
) %>%
  dplyr::select(Gene_ID) %>%
  filter(
    !is.na(Gene_ID),
    Gene_ID != "None"
  ) %>%
  distinct(Gene_ID)

#==================== 4. TF screening ====================

# Alkaline stress
TF_Alkaline <- filter_TF(
  expr_Alkaline,
  "geneID",
  rice_TF
)

# Drought stress
# Screen TFs using original Gene_ID before ID conversion
TF_Drought_original <- filter_TF(
  expr_Drought,
  "Gene_ID",
  rice_TF_1
)

# Cold stress
TF_Cold <- filter_TF(
  expr_Cold,
  "Gene_ID",
  rice_TF
)

#==================== 5. Count TF screening results ====================

# Count non-missing TF gene ID records
lib_Alkalinity <- sum(
  !is.na(TF_Alkaline$geneID)
)

lib_Drought <- sum(
  !is.na(TF_Drought_original$Gene_ID)
)

lib_Cold <- sum(
  !is.na(TF_Cold$Gene_ID)
)

lib_tf <- tibble(
  Stress = c("Alkalinity", "Drought", "Cold"),
  Count = c(
    lib_Alkalinity,
    lib_Drought,
    lib_Cold
  )
)

cat("\n===== TF screening count summary =====\n")
print(lib_tf)

#==================== 6. Gene ID conversion ====================

# Convert all drought expression gene IDs
expr_Drought_col1 <- expr_Drought %>%
  dplyr::select(Gene_ID)

expr_Drought_col1_replaced <- convert_and_replace_id(
  id_data = expr_Drought_col1,
  from_type = "RAP",
  to_type = "GID"
)

# Convert drought TF IDs using the same mapping
TF_Drought_col1 <- TF_Drought_original %>%
  dplyr::select(Gene_ID)

TF_Drought_col1_replaced <- convert_and_replace_id(
  id_data = TF_Drought_col1,
  from_type = "RAP",
  to_type = "GID"
)

# Alkaline and cold datasets do not require ID conversion
expr_Alkaline_col1 <- expr_Alkaline %>%
  dplyr::select(1)

expr_Cold_col1 <- expr_Cold %>%
  dplyr::select(1)

# Extract TF IDs
TF_Alkaline_col1 <- TF_Alkaline %>%
  dplyr::select(1)

TF_Cold_col1 <- TF_Cold %>%
  dplyr::select(1)

#==================== 7. Preview ID conversion ====================

cat("\n===== ID conversion preview =====\n")

cat("\nOriginal drought IDs (top 5 rows):\n")
print(head(expr_Drought_col1, 5))

cat("\nConverted drought IDs (top 5 rows):\n")
print(head(expr_Drought_col1_replaced, 5))

cat("\nConverted drought TF IDs (top 5 rows):\n")
print(head(TF_Drought_col1_replaced, 5))

#==================== 8. Construct unique gene ID sets ====================

# All expressed genes
rna_alkaline_ids <- unique(
  pull(expr_Alkaline_col1, 1)
)

rna_drought_ids <- unique(
  pull(expr_Drought_col1_replaced, 1)
)

rna_cold_ids <- unique(
  pull(expr_Cold_col1, 1)
)

# Transcription factor genes
tf_alkaline_ids <- unique(
  pull(TF_Alkaline_col1, 1)
)

tf_drought_ids <- unique(
  pull(TF_Drought_col1_replaced, 1)
)

tf_cold_ids <- unique(
  pull(TF_Cold_col1, 1)
)

# Remove missing IDs
rna_alkaline_ids <- na.omit(rna_alkaline_ids)
rna_drought_ids <- na.omit(rna_drought_ids)
rna_cold_ids <- na.omit(rna_cold_ids)

tf_alkaline_ids <- na.omit(tf_alkaline_ids)
tf_drought_ids <- na.omit(tf_drought_ids)
tf_cold_ids <- na.omit(tf_cold_ids)

#==================== 9. Plot TF screening count ====================

p1 <- ggplot(
  data = lib_tf,
  mapping = aes(x = Stress, y = Count)
) +
  geom_col() +
  geom_label(
    aes(label = Count),
    size = 5
  ) +
  labs(
    x = NULL,
    y = "TF Screening Count"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    panel.grid = element_blank()
  )

ggsave(
  filename = "TF_Screening_count.pdf",
  plot = p1,
  width = 7,
  height = 7,
  units = "in"
)

print(p1)

#==================== 10. Plot RNA-seq and TF Venn diagrams ====================

# RNA-seq gene Venn diagram
rna_colors <- c(
  "#0072B2",
  "#D55E00",
  "#009E73"
)

rna_id_list <- list(
  Alkaline = rna_alkaline_ids,
  Drought = rna_drought_ids,
  Cold = rna_cold_ids
)

rna_venn <- draw_venn(
  id_list = rna_id_list,
  title = "RNA-seq Data (Converted IDs)",
  fill_colors = rna_colors,
  filename = "RNA_Data_Venn"
)

# TF gene Venn diagram
tf_colors <- c(
  "#CC79A7",
  "#F0E442",
  "#56B4E9"
)

tf_id_list <- list(
  Alkaline = tf_alkaline_ids,
  Drought = tf_drought_ids,
  Cold = tf_cold_ids
)

tf_venn <- draw_venn(
  id_list = tf_id_list,
  title = "TF Data (Converted IDs)",
  fill_colors = tf_colors,
  filename = "TF_Data_Venn"
)

#==================== 11. Summary statistics ====================

gene_counts <- tibble(
  Group = c("Alkaline", "Drought", "Cold"),
  RNA_Count = c(
    length(rna_alkaline_ids),
    length(rna_drought_ids),
    length(rna_cold_ids)
  ),
  TF_Count = c(
    length(tf_alkaline_ids),
    length(tf_drought_ids),
    length(tf_cold_ids)
  )
)

cat("\n===== Gene count summary =====\n")
print(gene_counts)

#==================== 12. RNA-seq and TF gene count visualization ====================

gene_counts_long <- gene_counts %>%
  pivot_longer(
    cols = c(RNA_Count, TF_Count),
    names_to = "DataType",
    values_to = "GeneNumber"
  ) %>%
  mutate(
    DataType = factor(
      DataType,
      levels = c("RNA_Count", "TF_Count"),
      labels = c("RNA-seq", "TF")
    ),
    Group = factor(
      Group,
      levels = c("Alkaline", "Drought", "Cold")
    )
  )

bar_plot <- ggplot(
  gene_counts_long,
  aes(
    x = Group,
    y = GeneNumber,
    fill = DataType
  )
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "RNA-seq" = "#0072B2",
      "TF" = "#CC79A7"
    )
  ) +
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
    title = "Comparison of Gene Numbers (RNA-seq vs TF)",
    subtitle = "After ID Conversion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 16,
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(size = 12),
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
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

#==================== 13. TF proportion visualization ====================

gene_counts_ratio <- gene_counts %>%
  mutate(
    TF_Ratio = (TF_Count / RNA_Count) * 100,
    Group = factor(
      Group,
      levels = c("Alkaline", "Drought", "Cold")
    )
  )

ratio_plot <- ggplot(
  gene_counts_ratio,
  aes(
    x = Group,
    y = TF_Ratio,
    fill = Group
  )
) +
  geom_col(
    width = 0.7,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Alkaline" = "#0072B2",
      "Drought" = "#D55E00",
      "Cold" = "#009E73"
    )
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", TF_Ratio)
    ),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = "Stress Treatment",
    y = "TF Genes Ratio (%)",
    title = "TF Genes Proportion in RNA-seq Data",
    subtitle = "After ID Conversion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 16,
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
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

#==================== 14. Combine bar plots ====================

combined_plot <- grid.arrange(
  bar_plot,
  ratio_plot,
  ncol = 2
)

ggsave(
  filename = "RNA_TF_Combined_BarPlots.pdf",
  plot = combined_plot,
  width = 18,
  height = 7,
  dpi = 300
)

print(bar_plot)
print(ratio_plot)

#==================== 15. Export stress-specific unique TF genes ====================

tf_unique_subsets <- get_venn_unique_three(
  tf_id_list
)

export_unique_txt(
  tf_unique_subsets,
  prefix = "TF_Genes"
)

cat("\n===== Count of stress-exclusive TF genes =====\n")

print(
  purrr::map_int(
    tf_unique_subsets,
    length
  )
)

#==================== 16. Export summary tables ====================

write.xlsx(
  gene_counts,
  file = "RNA_TF_Gene_Counts.xlsx",
  rowNames = FALSE
)

write.xlsx(
  lib_tf,
  file = "TF_Screening_Counts.xlsx",
  rowNames = FALSE
)

cat("\n===== Analysis completed successfully =====\n")
