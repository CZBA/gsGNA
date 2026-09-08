# Load required R packages
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(VennDiagram)
library(igraph)
library(tidyverse)
library(miRspongeR)
library(ggplot2)
library(cowplot)
library(vroom)
library(ReactomePA)

# -------------------------- 2. Read network data (retain first three columns only) --------------------------
# Load edge‑list for three stress conditions, keep only first three columns
Z_ll_Alkalinity <- vroom("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv", delim = ",") %>% select(1:3)
Z_ll_Aridity   <- vroom("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv", delim = ",") %>% select(1:3)
Z_ll_Cold      <- vroom("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv", delim = ",") %>% select(1:3)

# -------------------------- 3. Define core functions (process each clustering method individually) --------------------------
#' Detect modules using a single community‑detection algorithm
#' @param network_data input edge table
#' @param method community detection method name
#' @return list of gene modules
get_single_method_modules <- function(network_data, method) {
  netModule(
    spongenetwork = network_data[, 1:2],  # only node pairs for network construction
    method = method,
    directed = FALSE,                     # treat GRN as undirected for module detection
    modulesize = 1                        # retain all raw modules; filter size afterwards
  )
}

#' Iterate over multiple community‑detection methods
#' @param network_data input edge table
#' @param methods vector of algorithm names
#' @return named list, each element stores module results from one algorithm
get_modules <- function(network_data, methods) {
  module_list <- list()
  for (method in methods) {
    cat(sprintf("→ Running module detection using %s ...\n", method))
    module_list[[method]] <- get_single_method_modules(network_data, method)
  }
  return(module_list)
}

#' Filter modules: keep modules with gene count > 10
#' @param modules nested list output from get_modules()
#' @return filtered nested module list
filter_large_modules <- function(modules) {
  lapply(modules, function(alg_modules) {
    alg_modules[sapply(alg_modules, length) > 10]
  })
}

#' Save detected modules with hierarchical folder structure (stress / algorithm)
#' @param valid_modules filtered module list
#' @param stress_name name of stress condition
save_modules_with_folder <- function(valid_modules, stress_name) {
  # Level‑1 folder for given stress condition
  stress_folder <- sprintf("%s_module", stress_name)
  dir.create(stress_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Iterate each clustering algorithm, create sub‑folder
  for (method in names(valid_modules)) {
    method_folder <- file.path(stress_folder, method)
    dir.create(method_folder, recursive = TRUE, showWarnings = FALSE)
    
    method_mods <- valid_modules[[method]]
    if (length(method_mods) == 0) {
      cat(sprintf("→ %s stress: no valid modules (>10 genes) for %s, skip\n", stress_name, method))
      next
    }
    
    # Write each module to individual text file
    for (mod_idx in seq_along(method_mods)) {
      mod_genes <- method_mods[[mod_idx]]
      file_path <- file.path(method_folder, sprintf("module_%d.txt", mod_idx))
      
      write.table(
        x = data.frame(gene = mod_genes),
        file = file_path,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
    }
    cat(sprintf("→ %s stress | %s: %d valid modules saved to %s\n",
                stress_name, method, length(method_mods), method_folder))
  }
  cat(sprintf("=== All modules for %s stress written into: %s\n\n", stress_name, stress_folder))
}

# -------------------------- 4. Batch process three stress datasets --------------------------
# Eight community‑detection algorithms applied
used_methods <- c("MCODE", "infomap", "prop", "eigen", "louvain", "walktrap", "FN", "MCL")

# 4.1 Alkalinity stress
cat("Start processing Alkalinity stress dataset ...\n")
alkalinity_all_modules  <- get_modules(Z_ll_Alkalinity, used_methods)
alkalinity_large_modules<- filter_large_modules(alkalinity_all_modules)
save_modules_with_folder(alkalinity_large_modules, stress_name = "Alkalinity")

# 4.2 Aridity stress
cat("Start processing Aridity stress dataset ...\n")
aridity_all_modules  <- get_modules(Z_ll_Aridity, used_methods)
aridity_large_modules<- filter_large_modules(aridity_all_modules)
save_modules_with_folder(aridity_large_modules, stress_name = "Aridity")

# 4.3 Cold stress
cat("Start processing Cold stress dataset ...\n")
cold_all_modules  <- get_modules(Z_ll_Cold, used_methods)
cold_large_modules<- filter_large_modules(cold_all_modules)
save_modules_with_folder(cold_large_modules, stress_name = "Cold")

# -------------------------- 5. Summarize valid module counts across stresses --------------------------
get_module_count <- function(large_modules) {
  sapply(large_modules, length)
}

stats_summary <- data.frame(
  Algorithm = used_methods,
  Alkalinity_valid_modules = get_module_count(alkalinity_large_modules)[used_methods],
  Aridity_valid_modules    = get_module_count(aridity_large_modules)[used_methods],
  Cold_valid_modules       = get_module_count(cold_large_modules)[used_methods]
)

cat("=== Summary: valid module count (size >10 genes) across three stress conditions ===\n")
print(stats_summary, row.names = FALSE)
