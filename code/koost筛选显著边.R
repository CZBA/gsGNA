library(vroom)
library(dplyr)
library(tidyr)

# 读取数据
Alkalinity <- vroom("alkalinity_edges.csv", delim = ",", col_names = TRUE)
Aridity <- vroom("aridity_edges.csv", delim = ",", col_names = TRUE)
Cold <- vroom("cold_edges.csv", delim = ",", col_names = TRUE)

# 读取TF文件并提取唯一值
tf_file <- "TFTarget.txt"
TF <- read.delim(tf_file, sep = "\t", header = TRUE)
if (ncol(TF) < 1) stop("TF文件至少需要一列数据")
TF <- select(TF, 1)
tf_values <- unique(TF[[1]])
if (length(tf_values) == 0) stop("未从TF文件中提取到有效值")

# 异常值处理函数（IQR法截断）
handle_outliers <- function(vec) {
  # 计算四分位和IQR
  q1 <- quantile(vec, 0.25, na.rm = TRUE)
  q3 <- quantile(vec, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  # 确定上下边界（1.5*IQR为常用阈值，可根据数据调整）
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  # 截断异常值（用边界值替换）
  vec[vec < lower_bound] <- lower_bound
  vec[vec > upper_bound] <- upper_bound
  return(vec)
}

# 通用Z值计算函数（先处理异常值，再标准化）
transform_zscore <- function(vec) {
  # 第一步：处理异常值
  vec_handled <- handle_outliers(vec)
  
  # 第二步：计算均值和标准差（基于处理后的数据）
  mean_val <- mean(vec_handled, na.rm = TRUE)
  sd_val <- sd(vec_handled, na.rm = TRUE)
  
  # 处理标准差为0的极端情况（所有值相同）
  if (sd_val == 0) {
    warning("数据的标准差为0，可能所有值相同，标准化结果均为0")
    return(rep(0, length(vec)))
  }
  
  # 计算Z值：(处理后的值 - 均值) / 标准差
  z_scores <- (vec_handled - mean_val) / sd_val
  return(z_scores)
}

# P值计算函数（基于Z值的正态分布假设）
calculate_p_values <- function(z_vector) {
  # 双侧检验：2*(1 - 正态分布下Z值的累积概率)
  p_values <- 2 * (1 - pnorm(abs(z_vector)))
  # 处理极端大的Z值（避免P值为0，用最小浮点数替代）
  p_values[p_values == 0] <- .Machine$double.xmin
  return(p_values)
}

# 基于P值的显著性筛选函数
is_significant <- function(p_values, pvalue.cutoff = 0.05) {
  significant <- p_values < pvalue.cutoff
  return(significant)
}

# 数据处理主函数（包含异常值处理步骤）
process_data <- function(data, dataset_name, tf_values, output_file = NULL, 
                         p_cutoff = 0.05) {
  cat(paste0("\n开始处理", dataset_name, "数据\n"))
  
  # 检查数据是否为空
  if(nrow(data) == 0) {
    stop(paste0(dataset_name, "数据为空，请检查输入文件"))
  }
  
  # 检查列数
  if (ncol(data) < 2) {
    stop(paste0(dataset_name, "数据至少需要2列（源、目标），当前列数：", ncol(data)))
  }
  cat("假设第一列为源节点，第二列为目标节点，第三列（或指定列）为权重\n")
  
  # 记录初始行数
  rows_before <- nrow(data)
  
  # 去重处理
  filtered_data <- distinct(data)
  rows_after_dedup <- nrow(filtered_data)
  cat("去重完成：从", rows_before, "行减少到", rows_after_dedup, "行，共删除", 
      rows_before - rows_after_dedup, "行重复数据\n")
  
  # 检查权重列
  if(!"EdgeWeight" %in% colnames(filtered_data)) {
    if(ncol(filtered_data) < 3) {
      stop(paste0(dataset_name, "数据列数不足3列，无法确定权重列位置"))
    }
    colnames(filtered_data)[3] <- "EdgeWeight"
    cat("假设第三列是权重列，并命名为'EdgeWeight'\n")
  }
  
  # 去除权重为0或NA的行
  rows_before_weight <- nrow(filtered_data)
  filtered_data <- filtered_data %>% 
    drop_na(EdgeWeight) %>% 
    filter(EdgeWeight != 0)
  rows_after_weight <- nrow(filtered_data)
  removed_zero <- rows_before_weight - rows_after_weight
  cat("去除权重为0或NA的行：共删除", removed_zero, "行，剩余", rows_after_weight, "行\n")
  
  # 用TF数据筛选第一列
  rows_before_tf <- nrow(filtered_data)
  filtered_data <- filtered_data %>% 
    filter(.[[1]] %in% tf_values)
  rows_after_tf <- nrow(filtered_data)
  cat("TF筛选完成：从", rows_before_tf, "行筛选到", rows_after_tf, "行匹配数据\n")
  
  if (rows_after_tf == 0) {
    warning(paste0(dataset_name, "数据在TF筛选后无剩余记录，返回空数据框"))
    return(filtered_data)
  }
  
  # 保存原始权重（用于后续对比异常值处理效果）
  filtered_data <- filtered_data %>%
    mutate(EdgeWeight_original = EdgeWeight)
  
  # 计算处理异常值后的Z值（核心改进：标准化前已处理异常值）
  filtered_data <- filtered_data %>%
    mutate(EdgeWeight.z = transform_zscore(EdgeWeight))
  
  # 计算P值
  filtered_data <- filtered_data %>%
    mutate(EdgeWeight.p_value = calculate_p_values(EdgeWeight.z))
  
  # 显著性筛选
  significant <- is_significant(filtered_data$EdgeWeight.p_value, pvalue.cutoff = p_cutoff)
  significant_edges <- filtered_data[significant, ]
  
  cat(paste0("显著边筛选完成（P值基于正态分布假设，阈值=", p_cutoff, "）：",
             "保留", nrow(significant_edges), "行显著边数据\n"))
  
  # 保存结果
  if(!is.null(output_file)) {
    dir <- dirname(output_file)
    if (!dir.exists(dir) && dir != "") {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    write.csv(significant_edges, output_file, row.names = FALSE)
    cat(paste0("已成功保存", dataset_name, "处理结果到", output_file, "\n"))
  }
  
  return(significant_edges)
}

# 处理三个数据集（包含异常值处理步骤）
alkalinity_significant <- process_data(
  Alkalinity, 
  "Alkalinity", 
  tf_values, 
  "Alkalinity_significant_edges.csv"
)

aridity_significant <- process_data(
  Aridity, 
  "Aridity", 
  tf_values, 
  "Aridity_significant_edges.csv"
)

cold_significant <- process_data(
  Cold, 
  "Cold", 
  tf_values, 
  "Cold_significant_edges.csv"
)
