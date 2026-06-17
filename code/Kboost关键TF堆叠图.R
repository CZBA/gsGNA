# 加载必需包
library(dplyr)
library(tidyr)
library(vroom)
library(tidyverse)
library(VennDiagram)
library(readxl)

# ------------------------------------------------------------------------------
# 1. 数据读取与预处理（统一格式、去重）
# ------------------------------------------------------------------------------
# 读取基因调控网络边数据（TF-靶基因-边权重）
load_grn_data <- function(file_path) {
  vroom(file_path, show_col_types = FALSE)[, 1:3] %>%
    rename(TF = 1, Target = 2, EdgeWeight = 3) %>%  # 统一列名
    drop_na() %>%                                    # 删除NA行
    distinct(TF, Target, .keep_all = TRUE)           # 去重TF-靶基因对
}

# 读取三种胁迫数据（确保文件路径正确）
Z_ll_Alkalinity <- load_grn_data("Alkalinity_significant_edges.csv")
Z_ll_Aridity <- load_grn_data("Aridity_significant_edges.csv")
Z_ll_Cold <- load_grn_data("Cold_significant_edges.csv")

# ------------------------------------------------------------------------------
# 2. 计算总TF、关键TF、非关键TF数量
# ------------------------------------------------------------------------------
# 2.1 提取各胁迫下的所有TF（去重）
tf_alkalinity_all <- Z_ll_Alkalinity %>% distinct(TF) %>% pull(TF)
tf_aridity_all <- Z_ll_Aridity %>% distinct(TF) %>% pull(TF)
tf_cold_all <- Z_ll_Cold %>% distinct(TF) %>% pull(TF)

# 2.2 度中心性计算（为筛选关键TF做准备）
cal_deg <- function(link_list) {
  link_list %>% 
    group_by(TF) %>% 
    summarise(degree = n(), .groups = "drop") %>%  # 统计每个TF的靶基因数（度）
    arrange(desc(degree))
}

deg_Alkalinity <- cal_deg(Z_ll_Alkalinity)
deg_Aridity <- cal_deg(Z_ll_Aridity)
deg_Cold <- cal_deg(Z_ll_Cold)

# 2.3 泊松分布筛选关键TF（右尾概率+Bonferroni校正）
find_key_TFs <- function(deg_data) {
  if (nrow(deg_data) == 0) {
    warning("度中心性数据为空，返回空列表")
    return(character(0))
  }
  mean_deg <- mean(deg_data$degree, na.rm = TRUE)
  deg_data$p_val <- ppois(deg_data$degree, lambda = mean_deg, lower.tail = FALSE)
  deg_data$adj_p_val <- p.adjust(deg_data$p_val, method = "bonferroni")
  deg_data %>% filter(adj_p_val < 0.05) %>% pull(TF)  # 仅返回关键TF名称
}

# 提取各胁迫下的关键TF
tf_alkalinity_key <- find_key_TFs(deg_Alkalinity)
tf_aridity_key <- find_key_TFs(deg_Aridity)
tf_cold_key <- find_key_TFs(deg_Cold)

# 2.4 合并所有胁迫的数据为一个数据框（用于合并绘图）
combined_data <- bind_rows(
  data.frame(
    Stress = "Alkalinity",
    TF_Category = c("Key TFs", "Non-key TFs"),
    Count = c(length(tf_alkalinity_key), length(tf_alkalinity_all) - length(tf_alkalinity_key))
  ),
  data.frame(
    Stress = "Aridity",
    TF_Category = c("Key TFs", "Non-key TFs"),
    Count = c(length(tf_aridity_key), length(tf_aridity_all) - length(tf_aridity_key))
  ),
  data.frame(
    Stress = "Cold",
    TF_Category = c("Key TFs", "Non-key TFs"),
    Count = c(length(tf_cold_key), length(tf_cold_all) - length(tf_cold_key))
  )
)

# ------------------------------------------------------------------------------
# 3. 绘制合并的堆叠柱状图
# ------------------------------------------------------------------------------
# 定义颜色
key_color <- "#DC143C"
non_key_color <- "#2E86AB"

# 绘制合并的堆叠图
p_combined <- ggplot(combined_data, aes(
  x = Stress,
  y = Count,
  fill = TF_Category,
  label = Count
)) +
  geom_col(width = 0.7, color = "black", alpha = 0.9) +
  geom_text(
    position = position_stack(vjust = 0.5),
    size = 4.2, 
    fontface = "bold", 
    color = "white"
  ) +
  scale_fill_manual(values = c("Key TFs" = key_color, "Non-key TFs" = non_key_color)) +
  labs(
    title = "Transcription Factors Distribution Across Different Stresses",
    x = "Stress Type",
    y = "Number of Transcription Factors",
    fill = "TF Category"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,    # 标题居中
      size = 14,      # 标题大小
      face = "bold",  # 标题加粗
      margin = margin(b = 15)  # 标题底部边距
    ),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 11),  # 坐标轴文本大小
    legend.title = element_text(size = 11, face = "bold"),  # 图例标题
    legend.text = element_text(size = 10),  # 图例文本
    legend.position = "top"  # 图例位置在顶部
  ) +
  ylim(0, max(aggregate(Count ~ Stress, combined_data, sum)$Count) * 1.1)  # 统一Y轴范围

# 保存合并后的图为PDF格式（核心修改处）
ggsave(
  filename = "TFs_Combined_Stacked.pdf",  # 文件名改为pdf
  plot = p_combined,
  width = 8,  # 适当加宽画布以容纳三个胁迫
  height = 6,
  dpi = 300,   # PDF为矢量图，dpi不影响清晰度，但保留参数保证兼容性
  bg = "white",
  device = "pdf"  # 显式指定保存为PDF格式
)

# 运行完成提示
cat("分析完成！\n")
cat("1. 合并的堆叠柱状图已保存为 TFs_Combined_Stacked.pdf\n")
