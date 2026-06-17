library(dplyr)
library(tidyr)
library(readr)
library(vroom)
library(igraph)
library(tidyverse)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(org.Osativa.eg.db)
library(clusterProfiler)
library(readxl)  
library(VennDiagram)



# 读取数据函数（彻底修复列名问题，不依赖临时列名）
load_grn_data <- function(file_path) {
  # 读取数据，不使用文件中的表头，直接按列位置处理
  raw_data <- vroom(
    file_path, 
    show_col_types = FALSE,
    col_names = FALSE,  # 忽略文件中的任何表头
    n_max = Inf
  )
  
  # 只保留前3列，并重命名（直接通过位置，不依赖列名）
  result <- raw_data[, 1:3] %>%
    set_names(c("TF", "Target", "EdgeWeight")) %>%  # 直接设置最终列名
    drop_na() %>%                                    # 删除含NA的行
    distinct(TF, Target, .keep_all = TRUE)           # 去重
  
  # 校验数据
  if(nrow(result) == 0) {
    warning(paste("文件", file_path, "读取后无有效数据"))
  }
  
  return(result)
}

# 读取边数据（现在不会再出现V1找不到的问题）
Z_ll_Alkalinity <- load_grn_data("Alkalinity_significant_edges.csv")
Z_ll_Aridity <- load_grn_data("Aridity_significant_edges.csv")
Z_ll_Cold <- load_grn_data("Cold_significant_edges.csv")

# 验证数据读取结果（关键：确认数据正确加载）
cat("=== 数据读取验证 ===\n")
cat("Alkalinity数据行数：", nrow(Z_ll_Alkalinity), "\n")
cat("Alkalinity数据列名：", paste(colnames(Z_ll_Alkalinity), collapse = ", "), "\n\n")

# 计算各胁迫下的TF数量
alkalinity_tf_count <- Z_ll_Alkalinity %>% distinct(TF) %>% nrow()
aridity_tf_count <- Z_ll_Aridity %>% distinct(TF) %>% nrow()
cold_tf_count <- Z_ll_Cold %>% distinct(TF) %>% nrow()

cat("=== TF数量统计 ===\n")
cat("Alkalinity胁迫的TF数量：", alkalinity_tf_count, "\n")
cat("Aridity胁迫的TF数量：", aridity_tf_count, "\n")
cat("Cold胁迫的TF数量：", cold_tf_count, "\n\n")

# 计算度中心性的函数
cal_deg <- function(link_list) {
  link_list %>% 
    group_by(TF) %>% 
    summarise(degree = n(), .groups = "drop") %>% 
    arrange(desc(degree))
}

# 绘制度中心性柱状图函数
plot_deg_centrality <- function(deg_data, stress_name) {
  ggplot(data = deg_data) +
    geom_col(mapping = aes(x = reorder(TF, -degree), y = degree)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      x = paste("TFs in", stress_name), 
      y = "Degree (Number of targets)",
      title = paste("TF Degree Centrality in", stress_name)
    ) +
    coord_cartesian(ylim = c(0, 35000))
}

# 计算每种胁迫下的度中心性
deg_Aridity  <- cal_deg(Z_ll_Aridity)
deg_Alkalinity <- cal_deg(Z_ll_Alkalinity)
deg_Cold <- cal_deg(Z_ll_Cold)

# 绘制度中心性柱状图
p_Aridity <- plot_deg_centrality(deg_Aridity, "Aridity")
p_Alkalinity <- plot_deg_centrality(deg_Alkalinity, "Alkalinity")
p_Cold <- plot_deg_centrality(deg_Cold, "Cold")

# 组合三个图并保存为PDF
combined_plot <- grid.arrange(p_Aridity, p_Alkalinity, p_Cold, nrow = 1)
ggsave("degree_centrality_plots.pdf",
       combined_plot, 
       width = 15, 
       height = 5,
       device = "pdf",
       dpi = 300,
       units = "in")

# 泊松分布计算关键转录因子的函数
find_key_TFs <- function(deg_data) {
  if(nrow(deg_data) == 0) {
    warning("输入数据为空，返回空数据框")
    return(data.frame(TF = character(), degree = integer(), p_value = numeric(), adj_p_value = numeric()))
  }
  
  # 计算平均度中心性
  mean_degree <- mean(deg_data$degree, na.rm = TRUE)
  
  # 计算泊松p值（使用生存函数计算极端值概率）
  deg_data$p_value <- ppois(deg_data$degree, lambda = mean_degree, lower.tail = FALSE)
  
  # 校正p值（Bonferroni校正）
  deg_data$adj_p_value <- p.adjust(deg_data$p_value, method = "bonferroni")
  
  # 筛选显著转录因子（p值小于0.05）
  key_TFs <- deg_data %>% filter(adj_p_value < 0.05)
  return(key_TFs)
}

# 找出每种胁迫下的关键转录因子
key_TFs_Aridity <- find_key_TFs(deg_Aridity)
key_TFs_Alkalinity <- find_key_TFs(deg_Alkalinity)
key_TFs_Cold <- find_key_TFs(deg_Cold)

# 提取泊松分布筛选出的转录因子名称
P_key_TFs_Aridity <- key_TFs_Aridity$TF
P_key_TFs_Alkalinity <- key_TFs_Alkalinity$TF
P_key_TFs_Cold <- key_TFs_Cold$TF

# 计算共有转录因子
P_Aridity_Alkalinity <- intersect(P_key_TFs_Aridity, P_key_TFs_Alkalinity)
P_Aridity_Cold <- intersect(P_key_TFs_Aridity, P_key_TFs_Cold)
P_Alkalinity_Cold <- intersect(P_key_TFs_Alkalinity, P_key_TFs_Cold)
P_all <- intersect(P_Aridity_Alkalinity, P_key_TFs_Cold)

# 提取关键TF边列表的函数
extract_key_edges <- function(link_list, key_TFs) {
  link_list %>% 
    filter(TF %in% key_TFs)
}

# 提取每种胁迫下的关键TF边列表
key_edges_Aridity <- extract_key_edges(Z_ll_Aridity, P_all)
key_edges_Alkalinity <- extract_key_edges(Z_ll_Alkalinity, P_all)
key_edges_Cold <- extract_key_edges(Z_ll_Cold, P_all)

# 绘制泊松分布方法的韦恩图并保存为PDF
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

pdf("poisson_venn.pdf", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()
# 设置工作目录
setwd("F:/20242015110015/水稻基因调控网络/KBoost/TF")
# 输出关键转录因子结果
write.table(P_key_TFs_Aridity, file = "P_key_TFs_Aridity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_key_TFs_Alkalinity, file = "P_key_TFs_Alkalinity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_key_TFs_Cold, file = "P_key_TFs_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 输出共有转录因子结果
write.table(P_Aridity_Alkalinity, file = "P_Aridity_Alkalinity.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_Aridity_Cold, file = "P_Aridity_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_Alkalinity_Cold, file = "P_Alkalinity_Cold.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(P_all, file = "P_all.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 筛选每种胁迫下关键TF及其靶基因（包含边权重）
get_key_tf_targets <- function(link_list, key_tfs) {
  link_list %>% 
    filter(TF %in% key_tfs) %>%
    dplyr::select(TF, Target, EdgeWeight) %>%
    arrange(TF, desc(EdgeWeight))
}

# 提取各胁迫下关键TF的靶基因关系
key_tf_targets_Aridity <- get_key_tf_targets(Z_ll_Aridity, P_key_TFs_Aridity)
key_tf_targets_Alkalinity <- get_key_tf_targets(Z_ll_Alkalinity, P_key_TFs_Alkalinity)
key_tf_targets_Cold <- get_key_tf_targets(Z_ll_Cold, P_key_TFs_Cold)

# 提取共有关键TF在各胁迫下的靶基因关系
common_tf_targets_Aridity <- get_key_tf_targets(Z_ll_Aridity, P_all)
common_tf_targets_Alkalinity <- get_key_tf_targets(Z_ll_Alkalinity, P_all)
common_tf_targets_Cold <- get_key_tf_targets(Z_ll_Cold, P_all)

# 输出结果到文件
write_tsv(key_tf_targets_Aridity, "key_TF_targets_Aridity.tsv")
write_tsv(key_tf_targets_Alkalinity, "key_TF_targets_Alkalinity.tsv")
write_tsv(key_tf_targets_Cold, "key_TF_targets_Cold.tsv")

# 输出共有关键TF的靶基因关系
write_tsv(common_tf_targets_Aridity, "common_TF_targets_Aridity.tsv")
write_tsv(common_tf_targets_Alkalinity, "common_TF_targets_Alkalinity.tsv")
write_tsv(common_tf_targets_Cold, "common_TF_targets_Cold.tsv")

# 打印各胁迫下关键TF-靶基因对数量
cat("=== 关键TF-靶基因对数量 ===\n")
cat("Aridity胁迫下关键TF-靶基因对数量：", nrow(key_tf_targets_Aridity), "\n")
cat("Alkalinity胁迫下关键TF-靶基因对数量：", nrow(key_tf_targets_Alkalinity), "\n")
cat("Cold胁迫下关键TF-靶基因对数量：", nrow(key_tf_targets_Cold), "\n")

