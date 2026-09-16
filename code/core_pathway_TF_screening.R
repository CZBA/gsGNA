# ===================== Load required R packages =====================
library(org.Osativa.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(VennDiagram)
library(vroom)     # Fast data reading
library(ggplot2)
library(readr)

# ===================== 0. Create output directory =====================
result_dir <- "Core_pathway_TF_screen_results"
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
  message(paste("Created output directory: ", result_dir))
} else {
  message(paste("Output directory already exists: ", result_dir))
}

# ===================== 1. Read global pathway gene sets (no ID conversion, use 3rd column) =====================
read_pathway_gene <- function(file_path) {
  gene_vec <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)[[3]]
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""] %>% unique()
  return(gene_vec)
}

# Global pathway gene sets (shared across all stress conditions)
MAPK_all_gid <- read_pathway_gene("MAPK_pathway_gene_anno.txt")
Hormone_all_gid <- read_pathway_gene("Hormone_pathway_gene_anno.txt")

# Check intersection between two pathways
cat("===== Intersection check: MAPK and hormone pathway gene sets =====\n")
intersection_gid <- intersect(MAPK_all_gid, Hormone_all_gid)
cat(paste0("Total MAPK pathway genes: ", length(MAPK_all_gid), "\n"))
cat(paste0("Total hormone pathway genes: ", length(Hormone_all_gid), "\n"))
cat(paste0("Intersection gene count: ", length(intersection_gid), "\n"))
if (length(intersection_gid) > 0) {
  cat(paste0("Example intersection genes: ", paste(head(intersection_gid, 5), collapse = ","), "\n"))
}
cat("\n")

# ===================== 2. Read key‑TF lists (for exclusion only, not used in downstream analysis) =====================
read_key_tfs <- function(file_path) {
  key_tfs <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
  key_tfs <- key_tfs[!is.na(key_tfs) & key_tfs != ""] %>% unique() %>% as.character()
  return(key_tfs)
}

key_TFs_Alkalinity <- read_key_tfs("Abiotic_TF_identification/XGBFusion_P_key_TFs_Alkalinity.txt")
key_TFs_Aridity   <- read_key_tfs("Abiotic_TF_identification/XGBFusion_P_key_TFs_Aridity.txt")
key_TFs_Cold      <- read_key_tfs("Abiotic_TF_identification/XGBFusion_P_key_TFs_Cold.txt")

key_tfs_list <- list(
  Alkalinity = key_TFs_Alkalinity,
  Aridity    = key_TFs_Aridity,
  Cold       = key_TFs_Cold
)

# ===================== 3. Compile target gene sets for each stress (share global pathway genes) =====================
target_genes <- list(
  Aridity    = c(MAPK_all_gid, Hormone_all_gid) %>% unique() %>% na.omit(),
  Alkalinity = c(MAPK_all_gid, Hormone_all_gid) %>% unique() %>% na.omit(),
  Cold       = c(MAPK_all_gid, Hormone_all_gid) %>% unique() %>% na.omit()
)

# ===================== 4. Load GRN edge tables for each stress =====================
edges_data <- list(
  Alkalinity = vroom("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv", show_col_types = FALSE)[,1:3],
  Aridity    = vroom("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv", show_col_types = FALSE)[,1:3],
  Cold       = vroom("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv", show_col_types = FALSE)[,1:3]
)

for (stress in names(edges_data)) {
  colnames(edges_data[[stress]]) <- c("source", "target", "weight")
  edges_data[[stress]]$target  <- as.character(edges_data[[stress]]$target)
  edges_data[[stress]]$source <- as.character(edges_data[[stress]]$source)
}

# ===================== 5. Matching function: retain non‑key TFs only =====================
#' Filter GRN edges targeting pathway genes and exclude predefined key‑TFs
#' @param edges_df GRN edge data frame
#' @param target_gids Combined gene set of MAPK and hormone pathways
#' @param key_tfs Vector of TFs to be excluded
#' @param stress Stress condition label
#' @param mapk_gid MAPK pathway gene vector
#' @param hormone_gid Hormone pathway gene vector
#' @return List containing filtered edges, TF list and pathway annotation table
match_only_nonkey_tf <- function(edges_df, target_gids, key_tfs, stress, mapk_gid, hormone_gid) {
  # Filter edges whose targets belong to pathway gene sets
  matched_by_target <- edges_df %>%
    filter(target %in% target_gids) %>%
    distinct(source, target, weight, .keep_all = TRUE)
  
  # Exclude key‑TFs, keep only non‑key TFs
  matched_non_key_tf <- matched_by_target %>%
    filter(!source %in% key_tfs)
  
  # Count distinct targets regulated by each non‑key TF
  non_key_tf_count <- matched_non_key_tf %>%
    group_by(source) %>%
    summarise(target_num = n_distinct(target))
  
  # Assign pathway association for each non‑key TF
  non_key_tf_pathway_association <- matched_non_key_tf %>%
    group_by(source) %>%
    mutate(
      is_mapk_target    = target %in% mapk_gid,
      is_hormone_target = target %in% hormone_gid
    ) %>%
    summarise(
      mapk_target_num     = sum(is_mapk_target),
      hormone_target_num  = sum(is_hormone_target),
      total_target_num    = n_distinct(target),
      pathway = case_when(
        mapk_target_num > 0 & hormone_target_num == 0 ~ "MAPK_signaling_pathway",
        mapk_target_num == 0 & hormone_target_num > 0 ~ "Plant_hormone_signal_transduction",
        mapk_target_num > 0 & hormone_target_num > 0 ~ "Dual_pathways(MAPK+Hormone)",
        TRUE ~ "Unknown"
      ),
      mapk_target_ratio    = ifelse(total_target_num > 0, mapk_target_num / total_target_num, 0),
      hormone_target_ratio = ifelse(total_target_num > 0, hormone_target_num / total_target_num, 0)
    ) %>%
    ungroup()
  
  # Console summary for current stress
  cat(paste0("===== ", stress, " stress: non‑key TF pathway classification =====\n"))
  pathway_dist <- table(non_key_tf_pathway_association$pathway)
  for (path in names(pathway_dist)) {
    cat(paste0("  ", path, "：", pathway_dist[path], " non‑key TFs\n"))
  }
  cat(paste0("  Total non‑key TFs under this stress: ", nrow(non_key_tf_pathway_association), "\n\n"))
  
  matched_non_key_tfs <- matched_non_key_tf$source %>% unique() %>% na.omit()
  non_key_tf_with_pathway <- non_key_tf_pathway_association %>%
    filter(source %in% matched_non_key_tfs)
  
  return(list(
    non_key_tf_matched_pairs = matched_non_key_tf,
    non_key_tfs              = matched_non_key_tfs,
    non_key_tf_pathway_info  = non_key_tf_with_pathway
  ))
}

# ===================== 6. Run matching pipeline for three stress conditions =====================
nonkey_results <- list()

# Aridity
nonkey_results[["Aridity"]] <- match_only_nonkey_tf(
  edges_df = edges_data[["Aridity"]],
  target_gids = target_genes[["Aridity"]],
  key_tfs = key_tfs_list[["Aridity"]],
  stress = "Aridity",
  mapk_gid = MAPK_all_gid,
  hormone_gid = Hormone_all_gid
)

# Alkalinity
nonkey_results[["Alkalinity"]] <- match_only_nonkey_tf(
  edges_df = edges_data[["Alkalinity"]],
  target_gids = target_genes[["Alkalinity"]],
  key_tfs = key_tfs_list[["Alkalinity"]],
  stress = "Alkalinity",
  mapk_gid = MAPK_all_gid,
  hormone_gid = Hormone_all_gid
)

# Cold
nonkey_results[["Cold"]] <- match_only_nonkey_tf(
  edges_df = edges_data[["Cold"]],
  target_gids = target_genes[["Cold"]],
  key_tfs = key_tfs_list[["Cold"]],
  stress = "Cold",
  mapk_gid = MAPK_all_gid,
  hormone_gid = Hormone_all_gid
)

# ===================== 7. Merge non‑key TFs across stresses; remove filter of >=2 stresses =====================
ntf_alk <- nonkey_results[["Alkalinity"]]$non_key_tf_pathway_info %>% mutate(stress = "Alkalinity", TF_type = "NonKey_TF")
ntf_ari <- nonkey_results[["Aridity"]]$non_key_tf_pathway_info %>% mutate(stress = "Aridity", TF_type = "NonKey_TF")
ntf_cold <- nonkey_results[["Cold"]]$non_key_tf_pathway_info %>% mutate(stress = "Cold", TF_type = "NonKey_TF")

merged_nonkey_all <- bind_rows(ntf_alk, ntf_ari, ntf_cold)
unique_all_nonkey_tfs <- merged_nonkey_all$source %>% unique() %>% na.omit()

# Count stress occurrence for each TF (for statistics only, not used as filter)
tf_stress_count <- merged_nonkey_all %>%
  distinct(source, stress) %>%
  count(source, name = "stress_occur")

# Keep all dual‑pathway non‑key TFs, no stress‑occurrence threshold
conserved_dual_tf_df <- merged_nonkey_all %>%
  filter(pathway == "Dual_pathways(MAPK+Hormone)")
unique_conserved_dual_tfs <- unique(conserved_dual_tf_df$source)

# TFs present in all three stresses (for supplementary statistics only)
common_all3_tf <- tf_stress_count %>% filter(stress_occur == 3) %>% pull(source)
common_dual_3stress_df <- merged_nonkey_all %>%
  filter(source %in% common_all3_tf, pathway == "Dual_pathways(MAPK+Hormone)")
unique_common_dual_3stress <- unique(common_dual_3stress_df$source)

# Global console summary
cat("===== Merged statistics for non‑key TFs across three stresses =====\n")
cat(paste0("Non‑key TFs (Alkalinity): ", nrow(ntf_alk), "\n"))
cat(paste0("Non‑key TFs (Aridity): ", nrow(ntf_ari), "\n"))
cat(paste0("Non‑key TFs (Cold): ", nrow(ntf_cold), "\n"))
cat(paste0("Deduplicated total non‑key TFs across stresses: ", length(unique_all_nonkey_tfs), "\n"))
cat(paste0("All dual‑pathway non‑key TFs (no ≥2‑stress filter): ", length(unique_conserved_dual_tfs), "\n"))
cat(paste0("Dual‑pathway TFs shared by all three stresses: ", length(unique_common_dual_3stress), "\n\n"))

conserved_dual_path_stat <- conserved_dual_tf_df %>%
  distinct(source, pathway) %>%
  group_by(pathway) %>%
  summarise(conserved_tf_num = n_distinct(source), .groups = "drop")

cat("【Final screening】Pathway distribution of dual‑pathway non‑key TFs:\n")
for (i in 1:nrow(conserved_dual_path_stat)) {
  cat(paste0("  - ", conserved_dual_path_stat$pathway[i], "：", conserved_dual_path_stat$conserved_tf_num, " TFs\n"))
}
cat("\n")

# ===================== 8. Export files for all dual‑pathway non‑key TFs =====================
# 1. Full annotation table for dual‑pathway TFs
write.csv(conserved_dual_tf_df,
          file.path(result_dir, "All_DualNonKey_TFs_FullInfo.csv"),
          row.names = FALSE)

# 2. ID list of all dual‑pathway non‑key TFs
if (length(unique_conserved_dual_tfs) > 0) {
  write.table(data.frame(TF_ID = unique_conserved_dual_tfs),
              file.path(result_dir, "All_DualNonKey_TFs_List.txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# 3. Supplementary list: dual‑pathway TFs present in all three stresses
if (length(unique_common_dual_3stress) > 0) {
  write.table(data.frame(TF_ID = unique_common_dual_3stress),
              file.path(result_dir, "Common_DualNonKey_TFs_All_3_stresses.txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# 4. Export TF lists split by pathway
conserved_tf_path_unique <- conserved_dual_tf_df %>% distinct(source, pathway)
conserved_split_path <- split(conserved_tf_path_unique, conserved_tf_path_unique$pathway)
for (path in names(conserved_split_path)) {
  fn <- gsub("/|\\(|\\)", "_", path)
  unique_tf_df <- conserved_split_path[[path]] %>% distinct(source)
  write.table(unique_tf_df,
              file.path(result_dir, paste0("All_", fn, "_NonKeyTFs.txt")),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# ===================== Generate annotated TF‑target table with pathway column =====================
# 1. Read BioMart annotation table (1st=RAP ID, 3rd=gene symbol/annotation)
BioMart_data <- read_delim("BioMart_data.txt", delim = "\t", show_col_types = FALSE)
colnames(BioMart_data)[1] <- "RAP"
colnames(BioMart_data)[3] <- "Gene_Annotation"
BioMart_data <- BioMart_data %>%
  select(RAP, Gene_Annotation) %>%
  distinct(RAP, Gene_Annotation) %>%
  drop_na()

# 2. Target TF set: all dual‑pathway non‑key TFs
target_tf_set <- unique_conserved_dual_tfs

# 3. Extract edges and attach pathway labels
edge_ari <- nonkey_results[["Aridity"]]$non_key_tf_matched_pairs %>%
  filter(source %in% target_tf_set) %>%
  left_join(conserved_dual_tf_df %>% filter(stress=="Aridity") %>% select(source, pathway), by=c("source"="source")) %>%
  mutate(Stress = "Aridity")

edge_alk <- nonkey_results[["Alkalinity"]]$non_key_tf_matched_pairs %>%
  filter(source %in% target_tf_set) %>%
  left_join(conserved_dual_tf_df %>% filter(stress=="Alkalinity") %>% select(source, pathway), by=c("source"="source")) %>%
  mutate(Stress = "Alkalinity")

edge_cold <- nonkey_results[["Cold"]]$non_key_tf_matched_pairs %>%
  filter(source %in% target_tf_set) %>%
  left_join(conserved_dual_tf_df %>% filter(stress=="Cold") %>% select(source, pathway), by=c("source"="source")) %>%
  mutate(Stress = "Cold")

all_edges_raw <- bind_rows(edge_ari, edge_alk, edge_cold) %>%
  rename(TF_GID = source, Target_GID = target) %>%
  select(TF_GID, Target_GID, pathway) %>% distinct()

# 4. GID‑to‑RAP ID mapping
all_gid_vec <- unique(c(all_edges_raw$TF_GID, all_edges_raw$Target_GID))
gid2rap_map <- bitr(geneID = all_gid_vec, fromType = "GID", toType = "RAP", OrgDb = org.Osativa.eg.db)
colnames(gid2rap_map) <- c("GID", "RAP")
gid2rap_map <- gid2rap_map %>% distinct(GID, RAP) %>% drop_na()

# 5. Annotate TF
final_table <- all_edges_raw %>%
  left_join(gid2rap_map, by = c("TF_GID" = "GID"), relationship = "many‑to‑many") %>%
  left_join(BioMart_data, by = c("RAP" = "RAP"), relationship = "many‑to‑many") %>%
  rename(TF_Annotation = Gene_Annotation) %>%
  select(-RAP)

# 6. Annotate Target
final_table <- final_table %>%
  left_join(gid2rap_map, by = c("Target_GID" = "GID"), relationship = "many‑to‑many") %>%
  left_join(BioMart_data, by = c("RAP" = "RAP"), relationship = "many‑to‑many") %>%
  rename(Target_Annotation = Gene_Annotation) %>%
  select(-RAP)

# 7. Remove NA entries, re‑order columns
final_table <- final_table %>%
  select(TF_GID, TF_Annotation, Target_GID, Target_Annotation, pathway) %>%
  drop_na(TF_Annotation, Target_Annotation) %>%
  distinct()

# 8. Export annotated table
out_file <- file.path(result_dir, "TF_Target_GID_Annotation_Pathway_Table_AllDualTF.csv")
write.csv(final_table, out_file, row.names = FALSE)

cat("\n========== Annotated TF‑target table (with pathway) exported ==========\n")
cat("Output file: ", out_file, "\n")
cat("Columns: TF_GID | TF_Annotation | Target_GID | Target_Annotation | pathway\n")
cat("Filter rules:\n1. Non‑key TF + dual‑pathway TF (no ≥2‑stress restriction)\n2. Rows without gene annotation are removed\n3. pathway label: Dual_pathways(MAPK+Hormone)\n")

# ===================== Export full target‑gene list regulated by dual‑pathway TFs =====================
all_dual_tf_edges <- bind_rows(
  nonkey_results[["Alkalinity"]]$non_key_tf_matched_pairs,
  nonkey_results[["Aridity"]]$non_key_tf_matched_pairs,
  nonkey_results[["Cold"]]$non_key_tf_matched_pairs
) %>%
  filter(source %in% unique_conserved_dual_tfs) %>%
  distinct(source, target) %>%
  rename(TF_GID = source, Target_GID = target)

# Assign target‑gene pathway membership
all_dual_tf_edges <- all_dual_tf_edges %>%
  mutate(
    Target_Pathway = case_when(
      Target_GID %in% intersect(MAPK_all_gid, Hormone_all_gid) ~ "MAPK & Hormone",
      Target_GID %in% MAPK_all_gid ~ "MAPK_signaling_pathway",
      Target_GID %in% Hormone_all_gid ~ "Plant_hormone_signal_transduction"
    )
  )

# TF annotation lookup
tf_anno <- all_dual_tf_edges %>%
  left_join(gid2rap_map, by = c("TF_GID" = "GID"), relationship = "many‑to‑many") %>%
  left_join(BioMart_data, by = "RAP", relationship = "many‑to‑many") %>%
  rename(TF_Annotation = Gene_Annotation) %>%
  select(TF_GID, TF_Annotation) %>%
  distinct()

# Target annotation lookup
target_anno <- all_dual_tf_edges %>%
  left_join(gid2rap_map, by = c("Target_GID" = "GID"), relationship = "many‑to‑many") %>%
  left_join(BioMart_data, by = "RAP", relationship = "many‑to‑many") %>%
  rename(Target_Annotation = Gene_Annotation) %>%
  select(Target_GID, Target_Annotation) %>%
  distinct()

# Merge annotation tables
tf_target_full_list <- all_dual_tf_edges %>%
  left_join(tf_anno, by = "TF_GID") %>%
  left_join(target_anno, by = "Target_GID") %>%
  select(TF_GID, TF_Annotation, Target_GID, Target_Annotation, Target_Pathway) %>%
  drop_na(TF_Annotation, Target_Annotation) %>%
  distinct() %>%
  arrange(TF_GID, Target_GID)

write.csv(tf_target_full_list,
          file.path(result_dir, "All_DualTF_All_TargetGenes_List.csv"),
          row.names = FALSE)

cat("=============================================\n")
cat(paste0("Unique TF‑target pairs regulated by dual‑pathway non‑key TFs: ", nrow(tf_target_full_list), "\n"))
cat(paste0("Distinct regulated target genes: ", n_distinct(tf_target_full_list$Target_GID), "\n"))
cat("Output file: All_DualTF_All_TargetGenes_List.csv\n")
cat("Columns: TF_GID, TF_Annotation, Target_GID, Target_Annotation, Target_Pathway\n")
cat("=============================================\n")

# ===================== Statistics for annotated valid dual‑pathway TFs =====================
valid_dual_tf_anno <- tf_target_full_list %>%
  select(TF_GID, TF_Annotation) %>%
  distinct()

cat("=============================================\n")
cat(paste0("Raw dual‑pathway non‑key TFs (GID list, before annotation filter): ", length(unique_conserved_dual_tfs), "\n"))
cat(paste0("Valid TFs after BioMart annotation & NA removal: ", nrow(valid_dual_tf_anno), "\n"))
cat("Valid TF list (GID + annotation):\n")
print(valid_dual_tf_anno, row.names = FALSE)

write.csv(valid_dual_tf_anno,
          file.path(result_dir, "All_DualNonKey_TFs_Annotated_Unique.csv"),
          row.names = FALSE)
cat("\nExport annotated deduplicated TF list: All_DualNonKey_TFs_Annotated_Unique.csv\n")
cat("=============================================\n")
