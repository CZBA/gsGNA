library(tidyverse)
library(dplyr)
library(readxl)
library(GenomicRanges)
library(ggplot2)
library(openxlsx)
library(edgeR)
library(org.Osativa.eg.db)
library(VennDiagram)
library(gridExtra)
# 定义函数用于读取 Excel 文件并选择指定列
read_excel_select <- function(file_path, cols) {
  read_excel(file_path) %>% dplyr::select(all_of(cols))
}

# 定义去零值函数，去除所有列值都为 0 的行
remove_zero_rows <- function(data) {
  data %>% filter(rowSums(dplyr::select(., where(is.numeric))) != 0)
}

# 读取并处理 RNA 序列数据
# 碱性
cpm_Alkalinity <- read_excel("GSE104928/GSE104928.xlsx") %>% remove_zero_rows()
# 干旱
cpm_Aridity <- read_excel_select("GSE121303_Processed_data.xlsx", c(2, 21,22,23)) %>% remove_zero_rows()
# 去除干旱数据中的重复基因 ID
cpm_Aridity <- cpm_Aridity %>% distinct(Gene_ID,.keep_all = TRUE)
# 寒冷
cpm_Cold <- read_excel_select("GSE112547.xlsx", 1:6) %>% remove_zero_rows()

# 读取水稻 TF 数据
# 从 TSV 文件读取
rice_TF <- read_tsv("Osj_TF_list.txt", col_names = TRUE) %>% dplyr::select(Gene_ID)
# 从 Excel 文件读取
rice_TF_1 <- read_excel("Osj_TF_list.xlsx") %>% 
  dplyr::select(Gene_ID) %>% 
  filter(Gene_ID != "None") %>% 
  na.omit()

# 定义函数用于筛选 TF 数据
filter_TF <- function(data, gene_col, tf_list) {
  data %>% filter(!!sym(gene_col) %in% tf_list$Gene_ID)
}

# 筛选 TF 数据
TF_Alkalinity <- filter_TF(cpm_Alkalinity, "geneID", rice_TF)
TF_Aridity <- filter_TF(cpm_Aridity, "Gene_ID", rice_TF_1)
TF_Cold <- filter_TF(cpm_Cold, "Gene_ID", rice_TF)

# ========== 提取第一列数据 ==========
# 提取原始RNA数据的第一列
cpm_Alkalinity_col1 <- cpm_Alkalinity %>% dplyr::select(1)  # 碱性RNA数据第一列
cpm_Aridity_col1 <- cpm_Aridity %>% dplyr::select(1)        # 干旱RNA数据第一列
cpm_Cold_col1 <- cpm_Cold %>% dplyr::select(1)              # 寒冷RNA数据第一列

# 提取TF筛选后数据的第一列
TF_Alkalinity_col1 <- TF_Alkalinity %>% dplyr::select(1)    # 碱性TF数据第一列
TF_Aridity_col1 <- TF_Aridity %>% dplyr::select(1)          # 干旱TF数据第一列
TF_Cold_col1 <- TF_Cold %>% dplyr::select(1)                # 寒冷TF数据第一列

# ========== 核心：转换ID并替换原ID ==========
# 定义通用的水稻ID转换函数（新增替换逻辑）
convert_and_replace_id <- function(id_data, from_type = "LOC_Os", to_type = "ENTREZID") {
  # 提取ID向量（处理单列数据框）
  gene_ids <- pull(id_data, 1)
  # 获取原列名（保证替换后列名和原数据一致）
  original_colname <- colnames(id_data)[1]
  
  # 空值检查
  if (length(gene_ids) == 0) {
    stop("输入的ID数据为空，请检查！")
  }
  
  # 使用org.Osativa.eg.db进行ID转换
  converted_ids <- mapIds(
    org.Osativa.eg.db,
    keys = gene_ids,
    keytype = from_type,
    column = to_type,
    multiVals = "first"  # 多个匹配时取第一个，避免返回列表
  )
  
  # 构建替换后的数据框（保留原列名，仅替换值）
  result <- tibble(
    !!original_colname := ifelse(is.na(converted_ids), gene_ids, converted_ids)
    # 转换失败的ID保留原值，也可改为 filter(!is.na(converted_ids)) 过滤掉
  )
  
  return(result)
}

# 1. 转换并替换cpm_Aridity_col1的ID（原始干旱RNA第一列）
cpm_Aridity_col1_replaced <- convert_and_replace_id(
  id_data = cpm_Aridity_col1,
  from_type = "RAP",  # 你的原始ID类型
  to_type = "GID"     # 目标ID类型
)

# 2. 转换并替换TF_Aridity_col1的ID（TF筛选后干旱数据第一列）
TF_Aridity_col1_replaced <- convert_and_replace_id(
  id_data = TF_Aridity_col1,
  from_type = "RAP",
  to_type = "GID"
)

# ========== 可选：查看替换前后对比 ==========
cat("cpm_Aridity_col1 替换前前5行：\n")
print(head(cpm_Aridity_col1, 5))
cat("\ncpm_Aridity_col1 替换后前5行：\n")
print(head(cpm_Aridity_col1_replaced, 5))

cat("\nTF_Aridity_col1 替换前前5行：\n")
print(head(TF_Aridity_col1, 5))
cat("\nTF_Aridity_col1 替换后前5行：\n")
print(head(TF_Aridity_col1_replaced, 5))

# ========== 清理工作环境：保留替换后的数据 ==========
# 更新保留变量列表，保留替换后的数据（替换原变量或新增均可）
keep_vars <- c(
  "cpm_Alkalinity_col1", "cpm_Aridity_col1_replaced", "cpm_Cold_col1",  # 替换后的cpm_Aridity_col1
  "TF_Alkalinity_col1", "TF_Aridity_col1_replaced", "TF_Cold_col1"  # 替换后的TF_Aridity_col1
   # 可选：保留原始ID数据用于对比
)

# 获取当前环境中所有变量名
all_vars <- ls(envir = globalenv())

# 筛选出需要删除的变量
remove_vars <- setdiff(all_vars, keep_vars)

# 排除自定义函数，避免误删
remove_vars <- remove_vars[!remove_vars %in% c("read_excel_select", "remove_zero_rows", "filter_TF", "convert_and_replace_id")]

# 执行删除操作
rm(list = remove_vars)

# ========== 1. 提取各数据集的ID集合（基于转换后的数据） ==========
# RNA数据ID集合
rna_alkalinity_ids <- pull(cpm_Alkalinity_col1, 1) %>% unique()  # 碱性RNA ID
rna_aridity_ids <- pull(cpm_Aridity_col1_replaced, 1) %>% unique()  # 干旱RNA ID（转换后）
rna_cold_ids <- pull(cpm_Cold_col1, 1) %>% unique()  # 寒冷RNA ID

# TF数据ID集合
tf_alkalinity_ids <- pull(TF_Alkalinity_col1, 1) %>% unique()  # 碱性TF ID
tf_aridity_ids <- pull(TF_Aridity_col1_replaced, 1) %>% unique()  # 干旱TF ID（转换后）
tf_cold_ids <- pull(TF_Cold_col1, 1) %>% unique()  # 寒冷TF ID

# ========== 2. 定义维恩图绘制函数（通用版，修改为PDF输出+科研配色） ==========
draw_venn <- function(id_list, title, fill_colors, filename) {
  # id_list: 命名的ID集合列表（如list(Alkalinity=ids1, Aridity=ids2, Cold=ids3)）
  # title: 维恩图标题
  # fill_colors: 填充色向量
  # filename: 输出文件名（无需后缀）
  
  # 绘制维恩图
  venn_plot <- venn.diagram(
    x = id_list,
    filename = NULL,  # 先不输出文件，返回图形对象
    col = "black",    # 边框颜色
    fill = fill_colors,  # 科研配色
    alpha = 0.6,      # 透明度（调整为更适合科研的0.6）
    label.col = "black",  # 数字标签颜色
    cex = 1.2,        # 数字标签大小
    fontfamily = "serif",
    cat.col = "black",  # 分类标签颜色
    cat.cex = 1.2,     # 分类标签大小
    cat.fontfamily = "serif",
    main = title,
    main.cex = 1.5,
    main.fontfamily = "serif"
  )
  
  # 保存为PDF格式（高分辨率，适合科研发表）
  pdf(paste0(filename, ".pdf"), width = 8, height = 8)  # 设置尺寸为8x8英寸，适合期刊要求
  grid.draw(venn_plot)
  dev.off()
  
  # 返回图形对象用于后续组合
  return(venn_plot)
}

# ========== 3. 绘制RNA数据维恩图（科研配色） ==========
# 科研常用配色（Nature/Science风格）：深蓝色、深红色、深绿色
rna_colors <- c("#0072B2", "#D55E00", "#009E73")  
rna_id_list <- list(
  Alkalinity = rna_alkalinity_ids,
  Aridity = rna_aridity_ids,
  Cold = rna_cold_ids
)

rna_venn <- draw_venn(
  id_list = rna_id_list,
  title = "RNA-seq Data (Converted IDs)",
  fill_colors = rna_colors,
  filename = "RNA_Data_Venn"
)

# ========== 4. 绘制TF数据维恩图（科研配色） ==========
# 科研常用配色：深紫色、深橙色、深青色
tf_colors <- c("#CC79A7", "#F0E442", "#56B4E9")  
tf_id_list <- list(
  Alkalinity = tf_alkalinity_ids,
  Aridity = tf_aridity_ids,
  Cold = tf_cold_ids
)

tf_venn <- draw_venn(
  id_list = tf_id_list,
  title = "TF Data (Converted IDs)",
  fill_colors = tf_colors,
  filename = "TF_Data_Venn"
)

# ========== 5. 统计RNA和TF数据的核心指标（用于绘制柱状图） ==========
# 统计各处理组的基因数量
gene_counts <- tibble(
  Group = c("Alkalinity", "Aridity", "Cold"),  # 处理组
  RNA_Count = c(
    length(rna_alkalinity_ids),  # 碱性RNA基因数
    length(rna_aridity_ids),     # 干旱RNA基因数（转换后）
    length(rna_cold_ids)        # 寒冷RNA基因数
  ),
  TF_Count = c(
    length(tf_alkalinity_ids),   # 碱性TF基因数
    length(tf_aridity_ids),      # 干旱TF基因数（转换后）
    length(tf_cold_ids)          # 寒冷TF基因数
  )
)

# 查看统计结果
cat("RNA和TF基因数量统计：\n")
print(gene_counts)

# ========== 6. 绘制RNA vs TF 数量对比柱状图（科研级样式） ==========
# 数据重塑（长格式，适合ggplot2）
gene_counts_long <- gene_counts %>%
  pivot_longer(
    cols = c(RNA_Count, TF_Count),
    names_to = "DataType",
    values_to = "GeneNumber"
  ) %>%
  mutate(
    # 美化标签
    DataType = factor(DataType, levels = c("RNA_Count", "TF_Count"), labels = c("RNA-seq", "TF")),
    Group = factor(Group, levels = c("Alkalinity", "Aridity", "Cold"))
  )

# 绘制分组柱状图
bar_plot <- ggplot(gene_counts_long, aes(x = Group, y = GeneNumber, fill = DataType)) +
  # 柱状图主体（宽度0.7，避免过宽）
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  # 科研配色（与维恩图风格统一）
  scale_fill_manual(values = c("#0072B2", "#CC79A7")) +
  # 添加数值标签（显示在柱子顶部）
  geom_text(
    aes(label = GeneNumber),
    position = position_dodge(width = 0.7),
    vjust = -0.5,  # 标签在柱子上方
    size = 4,      # 字体大小
    fontface = "bold"
  ) +
  # 坐标轴和标题设置（科研图表样式）
  labs(
    x = "Stress Treatment",  # X轴标题
    y = "Number of Genes",   # Y轴标题
    fill = "Data Type",      # 图例标题
    title = "Comparison of Gene Numbers (RNA-seq vs TF)",  # 主标题
    subtitle = "After ID Conversion"  # 副标题
  ) +
  # 主题设置（期刊级样式）
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()  # 去除网格线（更简洁）
  )

# 保存柱状图（PDF格式，高分辨率）
ggsave(
  filename = "RNA_TF_GeneCount_BarPlot.pdf",
  plot = bar_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# 显示图表（RStudio中）
print(bar_plot)

# ========== 7. 可选：绘制TF占RNA比例柱状图（更直观展示TF富集） ==========
# 计算TF占RNA的比例
gene_counts_ratio <- gene_counts %>%
  mutate(
    TF_Ratio = (TF_Count / RNA_Count) * 100,  # 百分比
    Group = factor(Group, levels = c("Alkalinity", "Aridity", "Cold"))
  )

# 绘制比例柱状图
ratio_plot <- ggplot(gene_counts_ratio, aes(x = Group, y = TF_Ratio, fill = Group)) +
  geom_col(width = 0.7, alpha = 0.8) +
  # 按处理组配色（与维恩图一致）
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  # 数值标签（显示百分比）
  geom_text(
    aes(label = sprintf("%.1f%%", TF_Ratio)),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = "Stress Treatment",
    y = "TF Genes Ratio (%)",
    title = "TF Genes Proportion in RNA-seq Data",
    subtitle = "After ID Conversion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none",  # 单组填充无需图例
    panel.grid = element_blank()
  )

# 保存比例图
ggsave(
  filename = "TF_Ratio_BarPlot.pdf",
  plot = ratio_plot,
  width = 8,
  height = 7,
  dpi = 300
)

# 显示比例图
print(ratio_plot)

# ========== 8. 可选：组合两个柱状图（一键输出） ==========
# 合并两个图表（使用gridExtra）
combined_plot <- grid.arrange(bar_plot, ratio_plot, ncol = 2)

# 保存组合图
ggsave(
  filename = "RNA_TF_Combined_BarPlots.pdf",
  plot = combined_plot,
  width = 18,
  height = 7,
  dpi = 300
)
