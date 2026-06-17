# 加载必要的库
library(readr)
library(dplyr)

# 读取数据
Alkalinity <- read_tsv("network_Alkalinity_expr_output.tsv", col_names = TRUE)
Aridity <- read_tsv("network_Aridity_expr_output.tsv", col_names = TRUE)
Cold <- read_tsv("network_Cold_expr_output.tsv", col_names = TRUE)


# 假设TF文件仍然需要用来筛选
tf_file <- "TFTarget.txt"
TF <- read.delim(tf_file, sep = "\t", header = TRUE) %>% 
  select(1)
tf_values <- unique(TF[[1]])

# Z-score转换函数
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}

# 基于Z-score的显著性筛选函数
biadjacency_matrix <- function(zscore_vector, pvalue.cutoff = 0.05) {
  # 对于Z-score，对应的p值cutoff 0.05约等于±1.96
  z_cutoff <- qnorm(1 - pvalue.cutoff/2)  # 计算双侧检验的Z临界值
  significant <- abs(zscore_vector) > z_cutoff
  return(significant)
}

# 数据处理主函数
process_data <- function(data, dataset_name, tf_values, output_file = NULL) {
  cat(paste0("\n开始处理", dataset_name, "数据\n"))
  
  # 记录初始行数
  rows_before <- nrow(data)
  
  # 去重处理
  filtered_data <- distinct(data)
  rows_after_dedup <- nrow(filtered_data)
  cat("去重完成：从", rows_before, "行减少到", rows_after_dedup, "行，共删除", 
      rows_before - rows_after_dedup, "行重复数据\n")
  
  # 检查是否有名为EdgeWeight的列，如果没有则假设第三列是权重
  if(!"EdgeWeight" %in% colnames(filtered_data)) {
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
  
  # 计算Z-score
  filtered_data <- filtered_data %>%
    mutate(EdgeWeight.zscore = transform_zscore(EdgeWeight))
  
  # 显著性筛选
  significant <- biadjacency_matrix(filtered_data$EdgeWeight.zscore, pvalue.cutoff = 0.05)
  significant_edges <- filtered_data[significant, ]
  cat("显著边筛选完成：保留", nrow(significant_edges), "行显著边数据\n")
  
  # 如果指定了输出文件，则保存结果
  if(!is.null(output_file)) {
    write.csv(significant_edges, output_file, row.names = FALSE)
    cat(paste0("已成功保存", dataset_name, "处理结果到", output_file, "\n"))
  }
  
  return(significant_edges)
}

# 处理三个数据集
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
