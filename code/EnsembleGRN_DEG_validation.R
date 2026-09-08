# ============================================================================
# Load required R packages
# ============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("org.Osativa.eg.db", quietly = TRUE)) BiocManager::install("org.Osativa.eg.db")
library(tidyverse)
library(vroom)
library(readxl)
library(writexl)
library(VennDiagram)
library(org.Osativa.eg.db)
library(grid)

# ============================================================================
# Utility function: Rice gene ID conversion (RAP → GID)
# Retained for dependency compatibility, NOT invoked in this workflow
# ============================================================================
rice_id_convert <- function(id_vec, fromType = "RAP", toType = "GID") {
  id_unique <- unique(id_vec)
  id_map <- mapIds(org.Osativa.eg.db,
                   keys = id_unique,
                   keytype = fromType,
                   column = toType,
                   multiVals = "first")
  # Match to original input order; NA for unmatched entries
  id_converted <- id_map[match(id_vec, names(id_map))]
  return(id_converted)
}

# ============================================================================
# Load GRN network data (CSV without header; first row contains raw data)
# EdgeWeight remains character type to suppress parsing warnings
# ============================================================================
load_grn_data <- function(file_path) {
  raw_data <- tryCatch({
    vroom::vroom(file_path, show_col_types = FALSE, col_names = FALSE,
                 col_types = cols(.default = "c"))
  }, error = function(e) {
    message("vroom read failed, fallback to read.csv: ", e$message)
    read.csv(file_path, stringsAsFactors = FALSE, na.strings = c("", "NA"),
             header = FALSE, colClasses = "character")
  })
  if (ncol(raw_data) < 3) stop("Input file contains fewer than 3 columns!")
  grn_data <- raw_data[, 1:3, drop = FALSE]
  colnames(grn_data) <- c("TF", "Target", "EdgeWeight")
  grn_data <- grn_data %>%
    mutate(TF = str_trim(TF), Target = str_trim(Target)) %>%
    drop_na() %>%
    distinct(TF, Target, .keep_all = TRUE)
  grn_data$pair_id <- paste(grn_data$TF, grn_data$Target, sep = "_")
  return(grn_data)
}

# ============================================================================
# Load DEG table (1st column = TF, 2nd column = Target)
# ============================================================================
load_deg_data <- function(file_path) {
  df <- read_xlsx(file_path) %>% drop_na()
  colnames(df) <- tolower(colnames(df))
  if (all(c("tf", "target") %in% colnames(df))) {
    df <- df %>% dplyr::select(tf, target) %>% dplyr::rename(TF = tf, Target = target)
  } else if ("tf" %in% colnames(df) & any(str_detect(colnames(df), "gene|gid|rap"))) {
    target_col <- str_subset(colnames(df), "gene|gid|rap")[1]
    df <- df %>% dplyr::select(tf, all_of(target_col)) %>% dplyr::rename(TF = tf, Target = target)
  } else {
    df <- df %>% dplyr::select(TF = 1, Target = 2)
    message("⚠️ Standard column names not detected; assign col1=TF, col2=Target by default")
  }
  df <- df %>%
    mutate(TF = str_trim(TF), Target = str_trim(Target)) %>%
    distinct(TF, Target, .keep_all = FALSE)
  return(df)
}

# ============================================================================
# Hypergeometric enrichment test for a single given transcription factor
# @param grn_df GRN edge data frame
# @param deg_df DEG pair data frame
# @param background_targets Background gene set (universe)
# @param tf_name Target TF identifier
# @return Tibble containing test statistics
# ============================================================================
single_tf_hyper_test <- function(grn_df, deg_df, background_targets, tf_name) {
  grn_tf_edges <- grn_df %>% filter(TF == tf_name)
  grn_targets <- unique(grn_tf_edges$Target)
  deg_targets <- unique(deg_df$Target)
  
  # Background universe: union of all target genes from input GRNs
  all_targets <- unique(background_targets)
  N <- 1 * length(all_targets)
  n <- length(grn_targets)
  
  # K: number of DEG targets present within background universe
  deg_in_background <- intersect(deg_targets, all_targets)
  K <- length(deg_in_background)
  
  # m: overlapping genes between GRN targets and DEG targets
  m <- sum(grn_targets %in% deg_targets)
  
  if (m == 0) {
    p_raw <- 1; fold_enrich <- 0; expect <- 0
  } else {
    expect <- n * (K / N)
    p_raw <- phyper(q = m - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
    fold_enrich <- m / expect
  }
  
  return(tibble(
    TF = tf_name,
    N = N,
    K = K,
    n = n,
    m = m,
    Expected = expect,
    Fold_Enrichment = fold_enrich,
    P_raw = p_raw
  ))
}

# ============================================================================
# Main workflow
# ============================================================================
message("=== Loading GRN data (no edge‑level filtering applied) ===")
Alkalinity_data <- load_grn_data("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv")
Aridity_data <- load_grn_data("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv")
Cold_data <- load_grn_data("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv")
message(sprintf("Alkalinity GRN edge count: %d", nrow(Alkalinity_data)))
message(sprintf("Aridity GRN edge count: %d", nrow(Aridity_data)))
message(sprintf("Cold GRN edge count: %d", nrow(Cold_data)))

message("=== Loading DEG datasets ===")
Alkalinity_DEG_raw <- load_deg_data("AlkalinityDEG.xlsx")
Aridity_DEG_raw <- load_deg_data("AridityDEG.xlsx")
Cold_DEG_raw <- load_deg_data("ColdDEG.xlsx")

# ====================== ID handling & unified filtering rules ======================
# Rule: Skip RAP‑GID ID conversion for all three stress datasets
# Filter rule: drop NA, deduplicate pairs; retain only Target genes starting with "LOC"
filter_loc_target <- function(df){
  df %>%
    drop_na(TF, Target) %>%
    distinct(TF, Target, .keep_all = FALSE) %>%
    filter(str_starts(Target, pattern = "^LOC"))
}

# 1. Alkalinity stress
message("✅ Alkalinity DEG: keep raw IDs; retain only LOC‑prefixed targets")
Alkalinity_DEG <- filter_loc_target(Alkalinity_DEG_raw)
message(sprintf("Alkalinity DEG: raw %d pairs, %d pairs retained after LOC filter", nrow(Alkalinity_DEG_raw), nrow(Alkalinity_DEG)))

# 2. Aridity stress
message("✅ Aridity DEG: keep raw IDs; retain only LOC‑prefixed targets")
Aridity_DEG <- filter_loc_target(Aridity_DEG_raw)
message(sprintf("Aridity DEG: raw %d pairs, %d pairs retained after LOC filter", nrow(Aridity_DEG_raw), nrow(Aridity_DEG)))

# 3. Cold stress
message("✅ Cold DEG: keep raw IDs; retain only LOC‑prefixed targets")
Cold_DEG <- filter_loc_target(Cold_DEG_raw)
message(sprintf("Cold DEG: raw %d pairs, %d pairs retained after LOC filter", nrow(Cold_DEG_raw), nrow(Cold_DEG)))

# ==========================================================================
# Extract unique TF identifiers from DEG tables
alk_tf <- unique(Alkalinity_DEG$TF)
ari_tf <- unique(Aridity_DEG$TF)
cold_tf <- unique(Cold_DEG$TF)

message(sprintf("Alkalinity DEG TF: %s (DEG pairs: %d)", alk_tf[1], nrow(Alkalinity_DEG)))
message(sprintf("Aridity DEG TF: %s (DEG pairs: %d)", ari_tf[1], nrow(Aridity_DEG)))
message(sprintf("Cold DEG TF: %s (DEG pairs: %d)", cold_tf[1], nrow(Cold_DEG)))

# Global background universe: union of all target genes across three GRNs
global_background <- unique(c(
  Alkalinity_data$Target, Aridity_data$Target, Cold_data$Target
))
message(sprintf("Global background gene count (GRN‑derived targets only): %d", length(global_background)))

alk_background <- global_background
ari_background <- global_background
cold_background <- global_background

# Hypergeometric enrichment test for each TF
alk_res <- single_tf_hyper_test(Alkalinity_data, Alkalinity_DEG, alk_background, alk_tf[1]) %>%
  mutate(Stress = "Alkaline")
ari_res <- single_tf_hyper_test(Aridity_data, Aridity_DEG, ari_background, ari_tf[1]) %>%
  mutate(Stress = "Aridity")
cold_res <- single_tf_hyper_test(Cold_data, Cold_DEG, cold_background, cold_tf[1]) %>%
  mutate(Stress = "Cold")

tf_enrich <- bind_rows(alk_res, ari_res, cold_res) %>%
  mutate(
    P_adjust_FDR = p.adjust(P_raw, method = "fdr"),
    Significance = case_when(
      P_adjust_FDR < 0.001 ~ "***",
      P_adjust_FDR < 0.01  ~ "**",
      P_adjust_FDR < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  mutate(Gene = TF)

# ============================================================================
# Compile summary statistics table
# ============================================================================
validation_stats <- tibble(
  Gene = c(alk_tf[1], ari_tf[1], cold_tf[1]),
  Stress = c("Alkaline stress", "Aridity stress", "Cold stress"),
  GRN_TF_Edges = c(
    nrow(Alkalinity_data %>% filter(TF == alk_tf[1])),
    nrow(Aridity_data %>% filter(TF == ari_tf[1])),
    nrow(Cold_data %>% filter(TF == cold_tf[1]))
  ),
  DEG_Edges = c(nrow(Alkalinity_DEG), nrow(Aridity_DEG), nrow(Cold_DEG)),
  Overlap = c(
    sum(Alkalinity_data$Target[Alkalinity_data$TF == alk_tf[1]] %in% Alkalinity_DEG$Target),
    sum(Aridity_data$Target[Aridity_data$TF == ari_tf[1]] %in% Aridity_DEG$Target),
    sum(Cold_data$Target[Cold_data$TF == cold_tf[1]] %in% Cold_DEG$Target)
  )
) %>%
  mutate(
    GRN_Match_Rate = ifelse(GRN_TF_Edges > 0, round(Overlap / GRN_TF_Edges * 100, 2), NA),
    DEG_Utilization_Rate = ifelse(DEG_Edges > 0, round(Overlap / DEG_Edges * 100, 2), NA)
  ) %>%
  left_join(tf_enrich %>% dplyr::select(Gene, N, K, n, Expected, Fold_Enrichment,
                                        P_raw, P_adjust_FDR, Significance),
            by = "Gene") %>%
  dplyr::select(Gene, Stress, GRN_TF_Edges, DEG_Edges, Overlap,
                GRN_Match_Rate, DEG_Utilization_Rate,
                N, K, n, Expected, Fold_Enrichment,
                P_raw, P_adjust_FDR, Significance)

message("\n=== Hypergeometric test results (background = GRN‑only targets, no GRN edge filter) ===")
print(validation_stats, width = Inf)
write_tsv(validation_stats, "GRN_DEG_GRNonlyBG_NoFilter_Hypergeometric_Results_NoIDConvert.tsv")
message("✅ Statistics saved: GRN_DEG_GRNonlyBG_NoFilter_Hypergeometric_Results_NoIDConvert.tsv")

# ============================================================================
# Venn‑diagram plotting module
# ============================================================================
message("\n🔍 Generating Venn diagrams with enrichment annotations ...")
grn_colors <- c("#E63946", "#457B9D", "#1D3557")
deg_color <- "#F1FAEE"
plot_w <- 8; plot_h <- 6

enrich_display <- tf_enrich %>%
  mutate(Stress_clean = case_when(
    Stress == "Alkaline" ~ "Alkaline stress",
    Stress == "Aridity"  ~ "Aridity stress",
    Stress == "Cold"     ~ "Cold stress"
  ))

#' Draw pairwise Venn for GRN‑predicted targets vs over‑expression DEG targets
#' @param gene_title TF gene name displayed on plot
#' @param stress_en Stress condition label
#' @param grn_df GRN edge table
#' @param deg_df DEG pair table
#' @param fill_color Fill color for GRN set
#' @param out_name Output file base name
#' @param enrich_df Table containing fold‑enrichment and significance
#' @return Summary tibble of overlap statistics
draw_tf_target_venn <- function(gene_title, stress_en, grn_df, deg_df, fill_color, out_name, enrich_df) {
  common_tf <- intersect(grn_df$TF, deg_df$TF)
  if (length(common_tf) == 0) {
    message(paste0("⚠️  ", gene_title, " (", stress_en, "): no shared TF between GRN and DEG"))
    return(NULL)
  }
  grn_target <- grn_df %>% filter(TF %in% common_tf) %>% pull(Target) %>% unique()
  deg_target <- deg_df %>% filter(TF %in% common_tf) %>% pull(Target) %>% unique()
  n_grn <- length(grn_target)
  n_deg <- length(deg_target)
  n_inter <- length(intersect(grn_target, deg_target))
  
  enrich <- enrich_df %>% filter(Stress_clean == stress_en)
  fe <- round(enrich$Fold_Enrichment, 2)
  sig <- enrich$Significance
  
  pdf(paste0(out_name, "_NoIDConvert.pdf"), width = plot_w, height = plot_h)
  draw.pairwise.venn(
    area1 = n_grn, area2 = n_deg, cross.area = n_inter,
    category = c(paste0("Targets in ", stress_en), "DEGs in OE Lines"),
    ellipse = TRUE, rotation.ratio = 0.5, scaled = FALSE, offset = 0.15,
    fill = c(fill_color, deg_color), col = "black", lwd = 1.2, alpha = 0.75,
    fontface = "bold", fontsize = 12,
    cat.cex = c(1.2, 1.2), cat.fontface = "bold", cat.dist = 0.03,
    cat.pos = c(0, 180), cat.col = "black"
  )
  grid.text(paste0(gene_title, "\nFE = ", fe, "  ", sig),
            x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))
  dev.off()
  
  return(tibble(
    Gene = gene_title, Stress = stress_en, Common_TF = length(common_tf),
    GRN_Targets = n_grn, OE_DEGs = n_deg, Overlap_Genes = n_inter,
    Overlap_Rate = round(n_inter / max(n_grn, 1) * 100, 2)
  ))
}

# Plot pairwise Venn for three stress conditions
OsWRKY1_stats <- draw_tf_target_venn(
  gene_title = alk_tf[1], stress_en = "Alkaline stress",
  grn_df = Alkalinity_data, deg_df = Alkalinity_DEG,
  fill_color = grn_colors[1], out_name = paste0(alk_tf[1], "_Alkaline_Venn"),
  enrich_df = enrich_display
)

OsbZIP62_stats <- draw_tf_target_venn(
  gene_title = ari_tf[1], stress_en = "Aridity stress",
  grn_df = Aridity_data, deg_df = Aridity_DEG,
  fill_color = grn_colors[2], out_name = paste0(ari_tf[1], "_Aridity_Venn"),
  enrich_df = enrich_display
)

OsWRKY70_stats <- draw_tf_target_venn(
  gene_title = cold_tf[1], stress_en = "Cold stress",
  grn_df = Cold_data, deg_df = Cold_DEG,
  fill_color = grn_colors[3], out_name = paste0(cold_tf[1], "_Cold_Venn"),
  enrich_df = enrich_display
)

venn_stats <- bind_rows(OsWRKY1_stats, OsbZIP62_stats, OsWRKY70_stats)
if (nrow(venn_stats) > 0) {
  message("\n=== Venn‑diagram target‑gene overlap statistics ===")
  print(venn_stats, n = 3)
  write_tsv(venn_stats, "Venn_Gene_Overlap_GRNonlyBG_NoFilter_NoIDConvert.tsv")
}

message("\n✅ Pipeline finished")
message("💡 Unified rule for three stress datasets: skip RAP‑GID ID conversion, keep raw identifiers")
message("💡 Filtering: drop NA & duplicate pairs; retain only Target genes starting with LOC")
message("💡 Background universe: aggregated target genes from all input GRNs; no GRN‑edge filtering")
message("💡 Significance markers: *** p<0.001, ** p<0.01, * p<0.05, ns non‑significant")
