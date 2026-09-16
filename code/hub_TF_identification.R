library(dplyr)
library(tidyr)
library(readr)
library(vroom)
library(igraph)
library(tidyverse)
library(gridExtra)
library(VennDiagram)

# ========== 2. Load GRN CSV with header (TF,Target,XGB_Score) ==========
load_grn_data <- function(file_path) {
  raw_data <- vroom(
    file_path,
    show_col_types = FALSE
  )
  
  result <- raw_data %>%
    rename(TF = TF, Target = Target, EdgeWeight = XGB_Score) %>%
    drop_na() %>%
    distinct(TF, Target, .keep_all = TRUE)
  
  if(nrow(result) == 0) {
    warning(paste("File", file_path, "has zero valid rows after filtering"))
  }
  return(result)
}

# ========== 3. Load XGB‑fusion GRN datasets ==========
Z_ll_Alkalinity <- load_grn_data("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv")
Z_ll_Aridity    <- load_grn_data("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv")
Z_ll_Cold       <- load_grn_data("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv")

cat("=== XGB‑fusion GRN loading summary ===\n")
cat("Alkalinity edge count: ", nrow(Z_ll_Alkalinity), "\n")
cat("Alkalinity columns: ", paste(colnames(Z_ll_Alkalinity), collapse = ", "), "\n\n")

alkalinity_tf_count <- Z_ll_Alkalinity %>% distinct(TF) %>% nrow()
aridity_tf_count    <- Z_ll_Aridity %>% distinct(TF) %>% nrow()
cold_tf_count       <- Z_ll_Cold %>% distinct(TF) %>% nrow()

cat("=== TF count per stress network ===\n")
cat("Alkalinity TF count: ", alkalinity_tf_count, "\n")
cat("Aridity TF count: ", aridity_tf_count, "\n")
cat("Cold TF count: ", cold_tf_count, "\n\n")

# ========== 4. Degree centrality calculation & plotting ==========
cal_deg <- function(link_list) {
  link_list %>%
    group_by(TF) %>%
    summarise(degree = n(), .groups = "drop") %>%
    arrange(desc(degree))
}

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

deg_Aridity    <- cal_deg(Z_ll_Aridity)
deg_Alkalinity <- cal_deg(Z_ll_Alkalinity)
deg_Cold       <- cal_deg(Z_ll_Cold)

p_Aridity    <- plot_deg_centrality(deg_Aridity, "Aridity")
p_Alkalinity <- plot_deg_centrality(deg_Alkalinity, "Alkalinity")
p_Cold       <- plot_deg_centrality(deg_Cold, "Cold")

combined_plot <- grid.arrange(p_Aridity, p_Alkalinity, p_Cold, nrow = 1)
ggsave("XGBFusion_degree_centrality_plots.pdf",
       combined_plot,
       width = 15,
       height = 5,
       device = "pdf",
       dpi = 300,
       units = "in")

# ========== 5. Poisson‑based hub‑TF detection ==========
find_key_TFs <- function(deg_data) {
  if(nrow(deg_data) == 0) {
    warning("Input degree data is empty")
    return(data.frame(TF = character(), degree = integer(), p_value = numeric(), adj_p_value = numeric()))
  }
  mean_degree <- mean(deg_data$degree, na.rm = TRUE)
  deg_data$p_value <- ppois(deg_data$degree, lambda = mean_degree, lower.tail = FALSE)
  deg_data$adj_p_value <- p.adjust(deg_data$p_value, method = "bonferroni")
  key_TFs <- deg_data %>% filter(adj_p_value < 0.05)
  return(key_TFs)
}

key_TFs_Aridity    <- find_key_TFs(deg_Aridity)
key_TFs_Alkalinity <- find_key_TFs(deg_Alkalinity)
key_TFs_Cold       <- find_key_TFs(deg_Cold)

P_key_TFs_Aridity    <- key_TFs_Aridity$TF
P_key_TFs_Alkalinity <- key_TFs_Alkalinity$TF
P_key_TFs_Cold       <- key_TFs_Cold$TF

P_Aridity_Alkalinity <- intersect(P_key_TFs_Aridity, P_key_TFs_Alkalinity)
P_Aridity_Cold       <- intersect(P_key_TFs_Aridity, P_key_TFs_Cold)
P_Alkalinity_Cold    <- intersect(P_key_TFs_Alkalinity, P_key_TFs_Cold)
P_all                <- intersect(P_Aridity_Alkalinity, P_key_TFs_Cold)

# ========== 6. Venn diagram for hub‑TFs ==========
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

# ========== 7. Export hub‑TF ID lists ==========
write.table(P_key_TFs_Aridity, file = "XGBFusion_P_key_TFs_Aridity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_key_TFs_Alkalinity, file = "XGBFusion_P_key_TFs_Alkalinity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_key_TFs_Cold, file = "XGBFusion_P_key_TFs_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(P_Aridity_Alkalinity, file = "XGBFusion_P_Aridity_Alkalinity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_Aridity_Cold, file = "XGBFusion_P_Aridity_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_Alkalinity_Cold, file = "XGBFusion_P_Alkalinity_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_all, file = "XGBFusion_P_all_common_TF.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ========== 8. Extract hub‑TF‑target edges ==========
extract_key_edges <- function(link_list, key_TFs) {
  link_list %>% filter(TF %in% key_TFs)
}

get_key_tf_targets <- function(link_list, key_tfs) {
  link_list %>%
    filter(TF %in% key_tfs) %>%
    dplyr::select(TF, Target, EdgeWeight) %>%
    arrange(TF, desc(EdgeWeight))
}

key_tf_targets_Aridity    <- get_key_tf_targets(Z_ll_Aridity, P_key_TFs_Aridity)
key_tf_targets_Alkalinity <- get_key_tf_targets(Z_ll_Alkalinity, P_key_TFs_Alkalinity)
key_tf_targets_Cold       <- get_key_tf_targets(Z_ll_Cold, P_key_TFs_Cold)

common_tf_targets_Aridity    <- get_key_tf_targets(Z_ll_Aridity, P_all)
common_tf_targets_Alkalinity <- get_key_tf_targets(Z_ll_Alkalinity, P_all)
common_tf_targets_Cold       <- get_key_tf_targets(Z_ll_Cold, P_all)

write_tsv(key_tf_targets_Aridity, "XGBFusion_key_TF_targets_Aridity.tsv")
write_tsv(key_tf_targets_Alkalinity, "XGBFusion_key_TF_targets_Alkalinity.tsv")
write_tsv(key_tf_targets_Cold, "XGBFusion_key_TF_targets_Cold.tsv")
write_tsv(common_tf_targets_Aridity, "XGBFusion_common_TF_targets_Aridity.tsv")
write_tsv(common_tf_targets_Alkalinity, "XGBFusion_common_TF_targets_Alkalinity.tsv")
write_tsv(common_tf_targets_Cold, "XGBFusion_common_TF_targets_Cold.tsv")

cat("=== XGB‑fusion GRN hub‑TF‑target pair counts ===\n")
cat("Aridity：", nrow(key_tf_targets_Aridity), "\n")
cat("Alkalinity：", nrow(key_tf_targets_Alkalinity), "\n")
cat("Cold：", nrow(key_tf_targets_Cold), "\n")

get_tf_targets_unique <- function(link_list, tf_set) {
  link_list %>%
    filter(TF %in% tf_set) %>%
    pull(Target) %>%
    unique() %>%
    na.omit() %>%
    .[. != ""]
}

common_tf_targets_Aridity_unique    <- get_tf_targets_unique(Z_ll_Aridity, P_all)
common_tf_targets_Alkalinity_unique <- get_tf_targets_unique(Z_ll_Alkalinity, P_all)
common_tf_targets_Cold_unique       <- get_tf_targets_unique(Z_ll_Cold, P_all)

clean_targets <- function(targets) {
  targets %>%
    unique() %>%
    na.omit() %>%
    .[. != ""]
}

common_tf_targets_Aridity_unique    <- clean_targets(common_tf_targets_Aridity_unique)
common_tf_targets_Alkalinity_unique <- clean_targets(common_tf_targets_Alkalinity_unique)
common_tf_targets_Cold_unique       <- clean_targets(common_tf_targets_Cold_unique)

targets_all_common <- intersect(
  intersect(common_tf_targets_Aridity_unique, common_tf_targets_Alkalinity_unique),
  common_tf_targets_Cold_unique
) %>% clean_targets()

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

if (!dir.exists("common_TF_targets_analysis")) {
  dir.create("common_TF_targets_analysis")
}

write_clean_targets <- function(targets, file_path) {
  if (length(targets) == 0) {
    writeLines("No valid targets found", file_path)
  } else {
    writeLines(targets, file_path)
  }
}

write_clean_targets(targets_all_common, "common_TF_targets_analysis/targets_all_common.txt")
write_clean_targets(targets_Aridity_Alkalinity_only, "common_TF_targets_analysis/targets_Aridity_Alkalinity_only.txt")
write_clean_targets(targets_Aridity_Cold_only, "common_TF_targets_analysis/targets_Aridity_Cold_only.txt")
write_clean_targets(targets_Alkalinity_Cold_only, "common_TF_targets_analysis/targets_Alkalinity_Cold_only.txt")
write_clean_targets(targets_Aridity_unique, "common_TF_targets_analysis/targets_Aridity_unique.txt")
write_clean_targets(targets_Alkalinity_unique, "common_TF_targets_analysis/targets_Alkalinity_unique.txt")
write_clean_targets(targets_Cold_unique, "common_TF_targets_analysis/targets_Cold_unique.txt")

cat("=== Shared‑hub‑TF target‑gene statistics (deduplicated & cleaned) ===\n")
cat("3‑stress common targets: ", length(targets_all_common), "\n")
cat("Aridity & Alkalinity only: ", length(targets_Aridity_Alkalinity_only), "\n")
cat("Aridity & Cold only: ", length(targets_Aridity_Cold_only), "\n")
cat("Alkalinity & Cold only: ", length(targets_Alkalinity_Cold_only), "\n")
cat("Aridity‑specific targets: ", length(targets_Aridity_unique), "\n")
cat("Alkalinity‑specific targets: ", length(targets_Alkalinity_unique), "\n")
cat("Cold‑specific targets: ", length(targets_Cold_unique), "\n")

cat("\nAnalysis finished.\n")
cat("Target‑gene lists saved in common_TF_targets_analysis/\n")
