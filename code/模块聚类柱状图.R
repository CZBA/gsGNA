# 加载必要的库
library(ggplot2)
library(tidyr)
library(Cairo)  # 确保PDF字体嵌入和兼容性

# 创建数据框（保留原始数值，仅修改标签为英文）
data <- data.frame(
  Algorithm = c("MCODE", "infomap", "prop", "eigen", "louvain", "walktrap", "FN", "MCL"),
  Alkalinity = c(1, 186, 1, 3, 3, 1, 2, 1),
  Aridity = c(1, 150, 1, 3, 3, 3, 3, 1),
  Cold = c(1, 275, 1, 2, 6, 1, 3, 1)
)

# 数据重塑为长格式（英文列名）
data_long <- pivot_longer(
  data, 
  cols = -Algorithm, 
  names_to = "Stress_Type", 
  values_to = "Valid_Module_Count"
)

# 定义科研配色（Nature期刊标准，适配3个胁迫类型）
research_palette <- c(
  "Alkalinity" = "#0072B2",    # 碱性-深海蓝
  "Aridity" = "#D55E00",       # 干旱-橙红
  "Cold" = "#009E73"           # 寒冷-祖母绿
)

# 绘制柱状图（全英文+科研配色）
p <- ggplot(data_long, aes(x = Algorithm, y = Valid_Module_Count, fill = Stress_Type)) +
  # 柱状图主体（科研级宽度/间距）
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  # 数值标签（优化位置和样式）
  geom_text(
    aes(label = Valid_Module_Count),
    position = position_dodge(width = 0.8),
    vjust = -0.3, size = 3.5,
    color = "gray30"  # 标签文字用深灰，更协调
  ) +
  # 英文标题和标签（科研论文标准）
  labs(
    title = "Comparison of Valid Module Counts by Clustering Algorithms",
    x = "Clustering Algorithm",
    y = "Number of Valid Modules",
    fill = "Stress Factor"
  ) +
  # 科研主题优化
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "gray30"),
    axis.text.y = element_text(size = 12, color = "gray30"),
    axis.title = element_text(size = 13, face = "bold", color = "gray30"),
    legend.title = element_text(size = 12, face = "bold", color = "gray30"),
    legend.text = element_text(size = 11, color = "gray30"),
    panel.grid = element_line(color = "gray90"),  # 浅灰网格，不抢主体
    panel.border = element_rect(color = "gray50") # 中灰边框，更柔和
  ) +
  # 应用科研配色
  scale_fill_manual(values = research_palette) +
  # 调整y轴范围，避免标签超出（保留逻辑，适配英文图）
  ylim(0, max(data_long$Valid_Module_Count) * 1.1)  

# 保存为PDF格式（矢量图，嵌入字体确保兼容性）
# 使用CairoPDF避免字体丢失，适配所有系统
ggsave(
  "Valid_Module_Count_by_Algorithm.pdf", 
  plot = p, 
  width = 10, 
  height = 7,
  device = cairo_pdf,  # 嵌入字体，解决跨系统显示问题
  bg = "white"
)

# 可选：保存为高分辨率PNG（用于预览）
# ggsave("Valid_Module_Count_by_Algorithm.png", plot = p, width = 10, height = 7, dpi = 300)