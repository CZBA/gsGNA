# Load packages
library(org.Osativa.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(VennDiagram)
library(ggplot2)  # Explicitly load ggplot2 to ensure plotting functions work

# Read the processed gene data
Aridity_Gene <- read.table("非生物胁迫TF鉴定/XGBFusion_P_key_TFs_Aridity.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
Alkalinity_Gene <- read.table("非生物胁迫TF鉴定/XGBFusion_P_key_TFs_Alkalinity.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
Cold_Gene <- read.table("非生物胁迫TF鉴定/XGBFusion_P_key_TFs_Cold.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
All_Gene <- read.table("非生物胁迫TF鉴定/XGBFusion_P_all_common_TF.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]

# Gene-ID conversion (GID to GO)
Aridity_Gene_list <- bitr(Aridity_Gene, fromType = "GID",
                          toType = "GO", OrgDb = org.Osativa.eg.db)

Alkalinity_Gene_list <- bitr(Alkalinity_Gene, fromType = "GID",
                             toType = "GO", OrgDb = org.Osativa.eg.db)

Cold_Gene_list <- bitr(Cold_Gene, fromType = "GID",
                       toType = "GO", OrgDb = org.Osativa.eg.db)

All_Gene_list <- bitr(All_Gene, fromType = "GID",
                      toType = "GO", OrgDb = org.Osativa.eg.db)

# Define the enrichment-analysis and PDF-visualization function
perform_enrichment_analysis <- function(gene_list, condition_name) {
  # Perform GO enrichment analysis (BP: biological process)
  ego <- enrichGO(gene          = gene_list$GO,
                  OrgDb         = org.Osativa.eg.db,
                  ont           = "BP",  
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  keyType       = "GO"
  )
  
  # Generate a Barplot and save as PDF (high resolution)
  pdf(paste0(condition_name, "_GO_barplot.pdf"), width = 10, height = 8)  # Set width/height adapted to GO-term length
  print(barplot(ego, showCategory = 10) + 
          ggtitle(paste(condition_name, "GO Enrichment (BP)")) +
          theme(plot.title = element_text(hjust = 0.5)))  # Center the title
  dev.off()
  
  # Generate a Dotplot and save as PDF
  pdf(paste0(condition_name, "_GO_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(ego, showCategory = 10) + 
          ggtitle(paste(condition_name, "GO Enrichment (BP)")) +
          theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  return(list(GO = ego))
}

# Run enrichment analysis for each gene list and output PDFs
Aridity_results <- perform_enrichment_analysis(Aridity_Gene_list, "Aridity")
Alkalinity_results <- perform_enrichment_analysis(Alkalinity_Gene_list, "Alkalinity")
Cold_results <- perform_enrichment_analysis(Cold_Gene_list, "Cold")
All_results <- perform_enrichment_analysis(All_Gene_list, "Common_tf")
