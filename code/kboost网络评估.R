# 加载依赖包（仅保留必要包）
library(tidyverse)
library(pROC)      # 计算AUROC
library(PRROC)     # 计算AUPR
library(vroom)     # 快速读取数据
library(ROSE) 

# ------------------------------------------------------------------------------
# 1. 数据加载与预处理
# ------------------------------------------------------------------------------
# 读取三个胁迫条件的GRN边列表（TF-靶基因-权重）并标准化处理
load_grn_data <- function(file_path) {
  vroom(file_path, show_col_types = FALSE)[, 1:3] %>%
    rename(TF = 1, Target = 2, EdgeWeight = 3) %>%  # 统一列名
    drop_na() %>%                                    # 删除含NA的行
    distinct(TF, Target, .keep_all = TRUE)           # 去重：保留唯一TF-靶基因对
}

# 加载数据
Alkalinity_data <- load_grn_data("Alkalinity_significant_edges.csv")
Aridity_data <- load_grn_data("Aridity_significant_edges.csv")
Cold_data <- load_grn_data("Cold_significant_edges.csv")

# 读取ChIP-seq验证数据并生成验证对
chip_seq <- read_tsv("chip_data.txt", col_names = TRUE) %>% 
  dplyr::select(TF = 1, Target = 2) %>%
  distinct()  # 确保验证集中无重复对
chip_pairs <- paste(chip_seq$TF, chip_seq$Target, sep = "_")  # 生成验证对字符串

# ------------------------------------------------------------------------------
# 2. 计算AUROC和AUPR评估指标
# 评估函数（修复K > N的核心问题）
evaluate_regulation <- function(pred_data, true_pairs, oversample = TRUE, p = 0.5, seed = 123) {
  # 添加标签
  pred_data <- pred_data %>%
    mutate(Label = ifelse(paste(TF, Target) %in% paste(true_pairs$TF, true_pairs$Target), 1, 0))
  
  # 计算标签数量
  pos_count <- sum(pred_data$Label == 1, na.rm = TRUE)
  neg_count <- sum(pred_data$Label == 0, na.rm = TRUE)
  
  cat("标签分布 - 正例:", pos_count, " 负例:", neg_count, "\n")
  
  if (pos_count == 0 || neg_count == 0) {
    warning("数据集中缺少正例或负例，无法进行有效评估")
    return(NULL)
  }
  
  set.seed(seed)
  train_data <- pred_data
  test_data <- pred_data
  
  # 过采样处理
  if (oversample) {
    train_data$Label <- as.factor(train_data$Label)
    oversampled_data <- ROSE(
      Label ~ EdgeWeight,
      data = train_data,
      seed = seed,
      p = p
    )$data
    eval_data <- oversampled_data
  } else {
    eval_data <- train_data
  }
  
  # 计算ROC和PR曲线指标
  test_roc_obj <- roc(test_data$Label, test_data$EdgeWeight)
  test_auroc <- pROC::auc(test_roc_obj)
  
  test_pr_obj <- pr.curve(
    scores.class0 = test_data$EdgeWeight[test_data$Label == 1],
    scores.class1 = test_data$EdgeWeight[test_data$Label == 0],
    curve = TRUE
  )
  test_aupr <- test_pr_obj$auc.integral
  
  # 核心修复：重新定义总体范围（解决K > N问题）
  threshold <- quantile(test_data$EdgeWeight, 0.9)
  pred_pairs_test <- paste(test_data$TF[test_data$EdgeWeight > threshold], 
                           test_data$Target[test_data$EdgeWeight > threshold]) %>%
    unique()
  
  # 关键修改：将总体定义为"预测数据和真实数据的并集"
  all_pred_pairs <- unique(paste(pred_data$TF, pred_data$Target))
  all_true_pairs <- unique(paste(true_pairs$TF, true_pairs$Target))
  all_possible_pairs <- unique(c(all_pred_pairs, all_true_pairs))  # 合并去重作为总体
  N <- length(all_possible_pairs)  # 新的总体大小（确保N >= K）
  
  true_pairs_str <- all_true_pairs  # 真实关系
  K <- length(true_pairs_str)       # 真实关系数（现在K <= N）
  
  overlap_pairs <- intersect(pred_pairs_test, true_pairs_str)
  k <- length(overlap_pairs)        # 命中数
  n <- length(pred_pairs_test)      # 预测数
  
  overlap_rate <- ifelse(K > 0, length(overlap_pairs) / K * 100, 0)
  
  # 参数合法性检查
  valid_params <- all(c(
    K >= 0,
    N >= K,          # 现在N是并集，确保N >= K
    n >= 0,
    n <= N,
    k >= 0,
    k <= min(K, n)
  ))
  
  # 计算p值
  if (valid_params && K > 0 && N > 0 && n > 0) {
    p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  } else {
    p_value <- NA
    warning("超几何分布参数无效，参数：K=", K, ", N=", N, ", n=", n, ", k=", k)
  }
  
  return(list(
    overlap_rate = overlap_rate,
    auroc = test_auroc,
    aupr = test_aupr,
    p_value = p_value,
    hyper_params = list(K=K, N=N, n=n, k=k, valid=valid_params)
  ))
}

# 后续评估代码保持不变
stress_list <- list(
  "Alkalinity" = Alkalinity_data,
  "Aridity" = Aridity_data,
  "Cold" = Cold_data
)

results_df <- data.frame(
  Stress = character(),
  Overlap_Rate = numeric(),
  AUROC = numeric(),
  AUPR = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

cat("===== 评估结果 =====\n")
for (stress_name in names(stress_list)) {
  stress_data <- stress_list[[stress_name]]
  
  cat(paste0("正在评估 ", stress_name, " 与 chip_seq 数据 (", Sys.time(), ")\n"))
  
  eval_result <- evaluate_regulation(
    stress_data, chip_seq,
    oversample = TRUE,
    p = 0.5,
    seed = 123
  )
  
  if (!is.null(eval_result)) {
    cat("  超几何分布参数 - K:", eval_result$hyper_params$K, 
        "N:", eval_result$hyper_params$N, 
        "n:", eval_result$hyper_params$n, 
        "k:", eval_result$hyper_params$k, 
        "有效:", eval_result$hyper_params$valid, "\n")
    
    results_df <- rbind(results_df, data.frame(
      Stress = stress_name,
      Overlap_Rate = round(eval_result$overlap_rate, 2),
      AUROC = round(eval_result$auroc, 3),
      AUPR = round(eval_result$aupr, 3),
      P_Value = ifelse(is.na(eval_result$p_value), NA, eval_result$p_value)
    ))
    
    cat(paste0(stress_name, " 评估结果:\n"))
    cat("  重叠率: ", results_df$Overlap_Rate[nrow(results_df)], "%\n", sep = "")
    cat("  AUROC: ", results_df$AUROC[nrow(results_df)], "\n", sep = "")
    cat("  AUPR: ", results_df$AUPR[nrow(results_df)], "\n", sep = "")
    cat("  超几何p值: ", ifelse(is.na(eval_result$p_value), "参数无效", 
                            format(eval_result$p_value, scientific = TRUE)), "\n\n")
  } else {
    cat(paste0(stress_name, " 评估失败，跳过\n\n"))
  }
}

write.csv(results_df, "evaluation_results_fixed_final.csv", row.names = FALSE)
cat("\n===== 汇总结果 =====\n")
print(results_df)