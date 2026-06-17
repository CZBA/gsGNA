# 方法1：使用tibble包的column_to_rownames函数
library(tibble)
load("K_alkalinity_result.RData")
load("K_aridity_result.RData")
load("K_cold_result.RData")
# 读取数据
# 读取数据并设置行名（行名为基因）
alkalinity_GRN <- K_alkalinity$GRN
aridity_GRN <- K_aridity$GRN
cold_GRN <- K_cold$GRN

# 转换为矩阵（行：基因，列：TF）
alkalinity_GRN <- as.matrix(alkalinity_GRN)
aridity_GRN <- as.matrix(aridity_GRN)
cold_GRN <- as.matrix(cold_GRN)

generate_edge_list <- function(grn_matrix) {
  # 找到非零调控关系（TF->基因）
  edges <- which(grn_matrix != 0, arr.ind = TRUE)
  if(nrow(edges) == 0) return(data.frame(source=character(), target=character(), weight=numeric()))
  
  # 源为TF（列名），目标为基因（行名）
  data.frame(
    source = colnames(grn_matrix)[edges[, 2]],  # TF名称（列名）
    target = rownames(grn_matrix)[edges[, 1]],  # 基因名称（行名）
    weight = grn_matrix[edges]
  )
}

rm(K_alkalinity,K_aridity,K_cold)

# 生成边列表
alkalinity_edges <- generate_edge_list(alkalinity_GRN)
aridity_edges <- generate_edge_list(aridity_GRN)
cold_edges <- generate_edge_list(cold_GRN)
# 输出为CSV文件，不包含行索引
write.csv(alkalinity_edges, "alkalinity_edges.csv", row.names = FALSE)
write.csv(aridity_edges, "aridity_edges.csv", row.names = FALSE)
write.csv(cold_edges, "cold_edges.csv", row.names = FALSE)

# 安装并加载必要的库（首次使用需安装）
# install.packages(c("vroom", "ggplot2", "dplyr", "tidyr"))
library(vroom)
library(ggplot2)
library(dplyr)
library(tidyr)

# 使用vroom快速读取CSV文件（假设权重列名为"weight"，请根据实际调整）
alkalinity_edges <- vroom("alkalinity_edges.csv", show_col_types = FALSE)
aridity_edges <- vroom("aridity_edges.csv", show_col_types = FALSE)
cold_edges <- vroom("cold_edges.csv", show_col_types = FALSE)

# 数据合并与预处理
combined_data <- bind_rows(
  alkalinity_edges %>% select(weight) %>% mutate(type = "碱度"),
  aridity_edges %>% select(weight) %>% mutate(type = "干旱度"),
  cold_edges %>% select(weight) %>% mutate(type = "寒冷度")
) %>%
  mutate(
    weight = pmin(pmax(weight, 0), 1)  # 将权重限制在0-1范围内
  )

# 绘制权重区间分布图
ggplot(combined_data, aes(x = weight, fill = type)) +
  geom_histogram(
    position = "dodge",  # 并列显示（避免重叠，也可改为"identity"加alpha=0.7重叠）
    bins = 15,           # 区间数量
    color = "white",     # 边框颜色
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +  # 自定义颜色
  labs(
    title = "不同因子的权重区间分布",
    x = "权重值",
    y = "频数",
    fill = "因子类型"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.key.size = unit(0.8, "cm"),
    axis.text = element_text(color = "black")
  ) +
  xlim(0, 1)
