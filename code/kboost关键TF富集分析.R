# 加载包
library(org.Osativa.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(VennDiagram)
library(ggplot2)  # 显式加载ggplot2，确保绘图功能正常

# 读取处理后的基因数据
Aridity_Gene <- read.table("TF/P_key_TFs_Aridity.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
Alkalinity_Gene <- read.table("TF/P_key_TFs_Alkalinity.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
Cold_Gene <- read.table("TF/P_key_TFs_Cold.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
All_Gene <- read.table("TF/P_all.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]

# 基因ID转换（GID转GO）
Aridity_Gene_list <- bitr(Aridity_Gene, fromType = "GID",
                          toType = "GO", OrgDb = org.Osativa.eg.db)

Alkalinity_Gene_list <- bitr(Alkalinity_Gene, fromType = "GID",
                             toType = "GO", OrgDb = org.Osativa.eg.db)

Cold_Gene_list <- bitr(Cold_Gene, fromType = "GID",
                       toType = "GO", OrgDb = org.Osativa.eg.db)

All_Gene_list <- bitr(All_Gene, fromType = "GID",
                      toType = "GO", OrgDb = org.Osativa.eg.db)

# 定义富集分析和PDF可视化函数
perform_enrichment_analysis <- function(gene_list, condition_name) {
  # 进行GO富集分析（BP：生物学过程）
  ego <- enrichGO(gene          = gene_list$GO,
                  OrgDb         = org.Osativa.eg.db,
                  ont           = "BP",  
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  keyType       = "GO"
  )
  
  # 生成Barplot并保存为PDF（高分辨率）
  pdf(paste0(condition_name, "_GO_barplot.pdf"), width = 10, height = 8)  # 设置宽高，适配GO术语长度
  print(barplot(ego, showCategory = 10) + 
          ggtitle(paste(condition_name, "GO Enrichment (BP)")) +
          theme(plot.title = element_text(hjust = 0.5)))  # 标题居中
  dev.off()
  
  # 生成Dotplot并保存为PDF
  pdf(paste0(condition_name, "_GO_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(ego, showCategory = 10) + 
          ggtitle(paste(condition_name, "GO Enrichment (BP)")) +
          theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  return(list(GO = ego))
}

# 对每个基因列表执行富集分析并输出PDF
Aridity_results <- perform_enrichment_analysis(Aridity_Gene_list, "Aridity")
Alkalinity_results <- perform_enrichment_analysis(Alkalinity_Gene_list, "Alkalinity")
Cold_results <- perform_enrichment_analysis(Cold_Gene_list, "Cold")
All_results <- perform_enrichment_analysis(All_Gene_list, "Common_tf")