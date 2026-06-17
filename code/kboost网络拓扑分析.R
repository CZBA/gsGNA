# 加载必要的库
library(igraph)
library(ggplot2)
library(tidyr)
library(dplyr)

# 从CSV文件读取边列表的函数
read_edge_list_from_csv <- function(file_path, weight_threshold = NULL) {
  edge_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  edge_data <- edge_data %>% rename(from = source, to =target, weight = EdgeWeight)
  
  if (!is.null(weight_threshold)) {
    edge_data <- edge_data[abs(edge_data$weight) > weight_threshold, ]
  }
  return(edge_data)
}

# 从边列表创建图的函数
create_graph_from_edge_list <- function(edge_list) {
  graph <- graph_from_data_frame(edge_list, directed = TRUE)
  E(graph)$weight <- edge_list$weight
  E(graph)$color <- ifelse(E(graph)$weight < 0, "blue", "red")
  return(graph)
}

# 生成随机网络的函数
generate_random_network <- function(original_graph) {
  num_nodes <- vcount(original_graph)
  num_edges <- ecount(original_graph)
  
  # 创建具有相同节点数和边数的随机网络
  random_graph <- sample_gnm(num_nodes, num_edges, directed = TRUE)
  
  # 保留原始网络的权重分布（正负比例）
  original_weights <- E(original_graph)$weight
  positive_weights <- original_weights[original_weights > 0]
  negative_weights <- original_weights[original_weights < 0]
  
  # 计算正负权重比例
  pos_ratio <- length(positive_weights) / length(original_weights)
  
  # 为随机网络分配权重
  num_pos <- round(num_edges * pos_ratio)
  num_neg <- num_edges - num_pos
  
  # 从原始权重中随机抽样
  random_pos_weights <- sample(positive_weights, num_pos, replace = TRUE)
  random_neg_weights <- sample(negative_weights, num_neg, replace = TRUE)
  
  # 组合正负权重并随机打乱
  random_weights <- c(random_pos_weights, random_neg_weights)
  random_weights <- sample(random_weights)
  
  # 为随机网络添加权重
  E(random_graph)$weight <- random_weights
  E(random_graph)$color <- ifelse(E(random_graph)$weight < 0, "blue", "red")
  
  return(random_graph)
}

# 网络拓扑分析函数
network_topology_analysis <- function(graph) {
  num_nodes <- vcount(graph)
  num_edges <- ecount(graph)
  if (num_nodes == 0 || num_edges == 0) return(NULL)
  
  avg_degree <- mean(degree(graph))
  clustering_coefficient <- transitivity(graph, type = "global")
  avg_path_length <- mean_distance(graph, weights = NA)
  density <- edge_density(graph)
  degree_dist <- degree_distribution(graph)
  
  return(list(
    num_nodes = num_nodes, num_edges = num_edges,
    avg_degree = avg_degree,
    clustering_coefficient = clustering_coefficient,
    avg_path_length = avg_path_length,
    density = density,
    degree_distribution = degree_dist
  ))
}

# 比较原始网络和随机网络的函数
compare_networks <- function(original_graph, name) {
  cat(paste0("Analyzing ", name, " network:\n"))
  
  # 分析原始网络
  original_result <- network_topology_analysis(original_graph)
  
  # 生成并分析随机网络
  random_graph <- generate_random_network(original_graph)
  random_result <- network_topology_analysis(random_graph)
  
  # 打印结果
  if (!is.null(original_result) && !is.null(random_result)) {
    cat("\n=== Original Network ===\n")
    cat("Number of nodes:", original_result$num_nodes, "\n")
    cat("Number of edges:", original_result$num_edges, "\n")
    cat("Average degree:", original_result$avg_degree, "\n")
    cat("Clustering coefficient:", original_result$clustering_coefficient, "\n")
    cat("Average path length:", original_result$avg_path_length, "\n")
    cat("Network density:", original_result$density, "\n\n")
    
    cat("=== Random Network ===\n")
    cat("Number of nodes:", random_result$num_nodes, "\n")
    cat("Number of edges:", random_result$num_edges, "\n")
    cat("Average degree:", random_result$avg_degree, "\n")
    cat("Clustering coefficient:", random_result$clustering_coefficient, "\n")
    cat("Average path length:", random_result$avg_path_length, "\n")
    cat("Network density:", random_result$density, "\n\n")
  }
  
  return(list(original = original_result, random = random_result))
}

# ----------------------
# 拓扑指标分面散点图函数（全英文）
# ----------------------
plot_topology_metrics <- function(original_analysis, random_analysis, plot_name) {
  results <- data.frame(
    Network_Type = c("Original Network", "Random Network"),
    Average_Degree = c(original_analysis$avg_degree, random_analysis$avg_degree),
    Clustering_Coefficient = c(original_analysis$clustering_coefficient, random_analysis$clustering_coefficient),
    Average_Path_Length = c(original_analysis$avg_path_length, random_analysis$avg_path_length)
  )
  
  results_long <- gather(results, key = "Metrics", value = "Values", -Network_Type)
  
  p <- ggplot(results_long, aes(x = Network_Type, y = Values, color = Network_Type)) +
    geom_point(size = 4, alpha = 0.8) +
    facet_wrap(~Metrics, scales = "free_y") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("Original Network" = "#1f77b4", "Random Network" = "#ff7f0e")) +
    labs(
      title = paste(plot_name, "Network Topological Metrics"),
      y = "Value"
    )
  
  print(p)
  ggsave(paste0(plot_name, "_topology_metrics.pdf"), plot = p, width = 10, height = 6, bg = "white")
  return(p)
}

# ----------------------
# ----------------------
# 度分布分析与幂律展示函数（最终修复版）
# ----------------------
plot_degree_distribution <- function(original_analysis, random_analysis, plot_name) {
  if (is.null(original_analysis$degree_distribution) || is.null(random_analysis$degree_distribution)) {
    warning(paste0(plot_name, " degree distribution data is invalid, skipping plotting"))
    return(NULL)
  }
  
  # 处理度分布数据
  original_dist <- original_analysis$degree_distribution
  random_dist <- random_analysis$degree_distribution
  max_degree <- max(length(original_dist), length(random_dist))
  original_dist <- c(original_dist, rep(0, max_degree - length(original_dist)))
  random_dist <- c(random_dist, rep(0, max_degree - length(random_dist)))
  degrees <- 1:max_degree
  original_nonzero <- original_dist > 0
  random_nonzero <- random_dist > 0
  
  # 幂律拟合（仅对原始网络）
  if(sum(original_nonzero) >= 2) {
    log_k <- log(degrees[original_nonzero])
    log_p <- log(original_dist[original_nonzero])
    power_law_fit <- lm(log_p ~ log_k)
    gamma <- -coefficients(power_law_fit)[[2]]
    intercept <- coefficients(power_law_fit)[[1]]
    fit_line <- exp(intercept) * (degrees^(-gamma))
  } else {
    gamma <- NA
    fit_line <- NULL
    warning(paste0(plot_name, " insufficient data for power-law fitting"))
  }
  
  # 输出为PDF设备
  pdf(paste0(plot_name, "_power_law_degree_distribution.pdf"), width = 10, height = 8)
  
  # 绘制双对数坐标图
  plot(degrees[original_nonzero], original_dist[original_nonzero], 
       log = "xy",  
       type = "p", pch = 16, col = "blue", cex = 1.2,
       xlab = "Degree (k)",
       ylab = "Frequency P(k)",
       main = ifelse(is.na(gamma),
                     paste(plot_name, "Degree Distribution (Log-Log Scale)"),
                     paste(plot_name, "Degree Distribution (Log-Log Scale)\nPower-Law Exponent (gamma) =", round(gamma, 2))),
       cex.lab = 1.2, cex.main = 1.2,
       ylim = if(length(c(original_dist[original_nonzero], random_dist[random_nonzero])) > 0) {
         c(min(original_dist[original_nonzero], random_dist[random_nonzero]), 
           max(original_dist[original_nonzero], random_dist[random_nonzero]))
       } else {
         c(1e-6, 1)  # 兜底默认值
       })
  
  # 添加随机网络的度分布
  if(sum(random_nonzero) > 0) {
    points(degrees[random_nonzero], random_dist[random_nonzero], 
           pch = 16, col = "red", cex = 1.2)
  }
  
  # 添加幂律拟合线
  if(!is.null(fit_line) && sum(original_nonzero) >= 2) {
    lines(degrees[original_nonzero], fit_line[original_nonzero], 
          col = "black", lwd = 2, lty = 2)
  }
  
  # 绘制图例（增加完整的鲁棒性检查）
  # 1. 定义图例元素（根据是否有拟合线调整）
  legend_labels <- c("Original Network", "Random Network")
  legend_colors <- c("blue", "red")
  legend_pch <- c(16, 16)
  legend_lty <- c(NA, NA)
  
  if(!is.null(fit_line) && sum(original_nonzero) >= 2) {
    legend_labels <- c(legend_labels, "Power-Law Fit")
    legend_colors <- c(legend_colors, "black")
    legend_pch <- c(legend_pch, NA)
    legend_lty <- c(legend_lty, 2)
  }
  
  # 2. 先绘制图例，再获取其坐标（避免plot=FALSE的空值问题）
  legend_obj <- legend("topright", 
                       legend = legend_labels,
                       col = legend_colors, 
                       pch = legend_pch, 
                       lty = legend_lty, 
                       lwd = 2, 
                       cex = 1.0,
                       bty = "n",  # 先不显示默认边框
                       inset = 0.02,
                       x.intersp = 1.2,
                       y.intersp = 1.5)
  
  # 3. 安全绘制自定义边框（增加空值检查）
  if (!is.null(legend_obj) && !is.null(legend_obj$rect)) {
    # 提取图例坐标并检查长度
    left <- ifelse(length(legend_obj$rect$left) > 0, legend_obj$rect$left, par("usr")[1] + 0.1)
    bottom <- ifelse(length(legend_obj$rect$bottom) > 0, legend_obj$rect$bottom, par("usr")[3] + 0.1)
    right <- ifelse(length(legend_obj$rect$right) > 0, legend_obj$rect$right, par("usr")[2] - 0.1)
    top <- ifelse(length(legend_obj$rect$top) > 0, legend_obj$rect$top, par("usr")[4] - 0.1)
    
    # 绘制边框（确保坐标有效）
    if (left < right && bottom < top) {
      rect(left, bottom, right, top, border = "black", lwd = 2)
    }
  }
  
  dev.off()  # 关闭PDF设备
  
  return(list(gamma = gamma))
}

# ----------------------
# 主程序
# ----------------------
weight_threshold <- NULL

# 读取数据（需确保CSV文件路径正确）
alkalinity_edges <- read_edge_list_from_csv("Alkalinity_significant_edges.csv", weight_threshold)
aridity_edges <- read_edge_list_from_csv("Aridity_significant_edges.csv", weight_threshold)
cold_edges <- read_edge_list_from_csv("Cold_significant_edges.csv", weight_threshold)

# 创建图对象
graph_aridity <- create_graph_from_edge_list(aridity_edges)
graph_cold <- create_graph_from_edge_list(cold_edges)
graph_alkalinity <- create_graph_from_edge_list(alkalinity_edges)

# 比较原始网络与随机网络
aridity_comparison <- compare_networks(graph_aridity, "Aridity")
cold_comparison <- compare_networks(graph_cold, "Cold")
alkalinity_comparison <- compare_networks(graph_alkalinity, "Alkalinity")

# 批量生成PDF格式图表
# 1. 碱度网络
plot_degree_distribution(alkalinity_comparison$original, alkalinity_comparison$random, "Alkalinity")
plot_topology_metrics(alkalinity_comparison$original, alkalinity_comparison$random, "Alkalinity")

# 2. 干旱度网络
plot_degree_distribution(aridity_comparison$original, aridity_comparison$random, "Aridity")
plot_topology_metrics(aridity_comparison$original, aridity_comparison$random, "Aridity")

# 3. 寒冷度网络
plot_degree_distribution(cold_comparison$original, cold_comparison$random, "Cold")
plot_topology_metrics(cold_comparison$original, cold_comparison$random, "Cold")
