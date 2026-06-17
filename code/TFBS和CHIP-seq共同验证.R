library(tidyverse)
library(vroom) 

# ------------------------------------------------------------------------------
# 1. 数据加载与预处理（分胁迫提取TF+筛选TFBS）
# ------------------------------------------------------------------------------
load_grn_data <- function(file_path) {
  vroom(file_path, show_col_types = FALSE)[, 1:3] %>%
    rename(TF = 1, Target = 2, EdgeWeight = 3) %>%
    drop_na() %>%
    distinct(TF, Target, .keep_all = TRUE)
}

# 加载GRN数据（注意：请确保文件路径正确）
Alkalinity_data <- load_grn_data("Alkalinity_significant_edges.csv")
Aridity_data <- load_grn_data("Aridity_significant_edges.csv")
Cold_data <- load_grn_data("Cold_significant_edges.csv")

# 分别提取三个胁迫各自的TF（去重）
alkalinity_tfs <- unique(Alkalinity_data$TF)
aridity_tfs <- unique(Aridity_data$TF)
cold_tfs <- unique(Cold_data$TF)

message(sprintf("→ 碱胁迫（Alkalinity）GRN中提取到 %d 个独特TF", length(alkalinity_tfs)))
message(sprintf("→ 干旱胁迫（Aridity）GRN中提取到 %d 个独特TF", length(aridity_tfs)))
message(sprintf("→ 冷胁迫（Cold）GRN中提取到 %d 个独特TF", length(cold_tfs)))

# 先加载完整TFBS数据
tfbs_full <- read_tsv("TFTarget.txt", col_names = TRUE) %>% 
  select(TF = 1, Target = 2) %>% 
  distinct()

# 按各胁迫TF筛选对应的TFBS对
tfbs_alkalinity <- tfbs_full %>% filter(TF %in% alkalinity_tfs)
tfbs_aridity <- tfbs_full %>% filter(TF %in% aridity_tfs)
tfbs_cold <- tfbs_full %>% filter(TF %in% cold_tfs)

# 生成各胁迫TFBS对应的调控对及数量统计
tfbs_alkalinity_pairs <- paste(tfbs_alkalinity$TF, tfbs_alkalinity$Target, sep = "_")
tfbs_aridity_pairs <- paste(tfbs_aridity$TF, tfbs_aridity$Target, sep = "_")
tfbs_cold_pairs <- paste(tfbs_cold$TF, tfbs_cold$Target, sep = "_")

# 各胁迫TFBS验证集总量 & GRN被验证的数量
tfbs_alkalinity_total <- length(tfbs_alkalinity_pairs)
tfbs_aridity_total <- length(tfbs_aridity_pairs)
tfbs_cold_total <- length(tfbs_cold_pairs)

alkalinity_pairs <- paste(Alkalinity_data$TF, Alkalinity_data$Target, sep = "_")
aridity_pairs <- paste(Aridity_data$TF, Aridity_data$Target, sep = "_")
cold_pairs <- paste(Cold_data$TF, Cold_data$Target, sep = "_")

alkalinity_verified <- sum(alkalinity_pairs %in% tfbs_alkalinity_pairs)
aridity_verified <- sum(aridity_pairs %in% tfbs_aridity_pairs)
cold_verified <- sum(cold_pairs %in% tfbs_cold_pairs)

message(sprintf("→ 碱胁迫TFBS验证集总量：%d | GRN被验证数量：%d", tfbs_alkalinity_total, alkalinity_verified))
message(sprintf("→ 干旱胁迫TFBS验证集总量：%d | GRN被验证数量：%d", tfbs_aridity_total, aridity_verified))
message(sprintf("→ 冷胁迫TFBS验证集总量：%d | GRN被验证数量：%d", tfbs_cold_total, cold_verified))

# ------------------------------------------------------------------------------
# 2. 构建对比数据集（TFBS总量 vs GRN验证数量）
# ------------------------------------------------------------------------------
comparison_data <- tibble(
  condition = factor(c("Alkalinity", "Aridity", "Cold"), 
                     levels = c("Alkalinity", "Aridity", "Cold")),
  tfbs_total = c(tfbs_alkalinity_total, tfbs_aridity_total, tfbs_cold_total),
  grn_verified = c(alkalinity_verified, aridity_verified, cold_verified)
) %>%
  # 重塑数据为长格式，适配堆叠图
  pivot_longer(-condition, names_to = "data_type", values_to = "count") %>%
  mutate(
    data_type = factor(data_type,
                       levels = c("grn_verified", "tfbs_total"),
                       labels = c("GRN Verified by TFBS", "TFBS"))  # 关键修改：将TFBS Validation Set Total改为TFBS
  )

# 查看对比数据
message("\n→ TFBS总量与GRN验证数量对比表：")
print(comparison_data)

# ------------------------------------------------------------------------------
# 3. 绘制堆叠图（TFBS总量 vs GRN验证数量，红蓝配色 + 边框）
# ------------------------------------------------------------------------------
# 定义红蓝配色（标签与数据集中的labels完全对应，避免不匹配）
color_palette <- c(
  "GRN Verified by TFBS" = "#DC143C",    # GRN验证数量用红色
  "TFBS" = "#2E86AB"                    # TFBS总量用蓝色（与简化后的标签对应）
)

ggplot(comparison_data, aes(x = condition, y = count, fill = data_type)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.4) +
  scale_fill_manual(values = color_palette) +  # 关键：使用自定义红蓝配色
  labs(
    x = "Stress Condition",
    y = "Number of Regulatory Pairs",
    fill = "Data Category",
    title = "TFBS Validation of Gene Regulatory Networks"
  ) +
  theme_minimal() +
  # 核心修改：添加边框 + 优化轴线样式
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    # 坐标轴轴线样式
    axis.line = element_line(color = "black", linewidth = 0.5),  
    axis.ticks = element_line(color = "black", linewidth = 0.5), 
    axis.ticks.length = unit(2, "mm"),                          
    panel.grid = element_blank(),                               
    # 关键：添加图片整体边框（四个边）
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    # 确保边框完整显示
    plot.margin = margin(10, 10, 10, 10, "pt")
  ) +
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black",
    fontface = "bold"
  )

# 保存图片为PDF格式
ggsave("TFBS_Total_vs_GRN_Verified_Stacked_Barplot_RedBlue.pdf", 
       width = 6.5, height = 7, dpi = 300, bg = "white",
       device = "pdf",  # 指定输出为PDF格式
       useDingbats = FALSE)  # 避免PDF中的字体兼容性问题
message("\n→ 红蓝配色堆叠图已保存为PDF格式（含坐标轴轴线+整体边框）：TFBS_Total_vs_GRN_Verified_Stacked_Barplot_RedBlue.pdf")
