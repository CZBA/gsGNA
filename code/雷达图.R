# 加载所需包
library(reshape2)
library(fmsb)
library(Cairo)  # 解决中文渲染

# 安全提取拓扑指标函数
safe_extract <- function(x, default = 0) {
  if (is.null(x) || is.na(x) || is.infinite(x)) return(default)
  return(x)
}

# 提取各网络拓扑指标
topo_metrics <- data.frame(
  Network = c("碱性", "干旱", "寒冷"),
  Nodes = c(
    safe_extract(alkalinity_comparison$original$num_nodes),
    safe_extract(aridity_comparison$original$num_nodes),
    safe_extract(cold_comparison$original$num_nodes)
  ),
  Edges = c(
    safe_extract(alkalinity_comparison$original$num_edges),
    safe_extract(aridity_comparison$original$num_edges),
    safe_extract(cold_comparison$original$num_edges)
  ),
  Avg_Degree = c(
    safe_extract(alkalinity_comparison$original$avg_degree),
    safe_extract(aridity_comparison$original$avg_degree),
    safe_extract(cold_comparison$original$avg_degree)
  ),
  Clustering_Coeff = c(
    safe_extract(alkalinity_comparison$original$clustering_coefficient),
    safe_extract(aridity_comparison$original$clustering_coefficient),
    safe_extract(cold_comparison$original$clustering_coefficient)
  ),
  Avg_Path_Length = c(
    safe_extract(alkalinity_comparison$original$avg_path_length),
    safe_extract(aridity_comparison$original$avg_path_length),
    safe_extract(cold_comparison$original$avg_path_length)
  ),
  check.names = FALSE
)

# ========== 科研配色方案选择 ==========
# 方案1：Nature期刊常用配色（低饱和+高区分）
line_colors <- c("#0072B2", "#D55E00", "#009E73")  # 蓝、橙红、绿
# 方案2：Science期刊常用配色（灰度+彩色搭配）
# line_colors <- c("#332288", "#88CCEE", "#117733")

fill_colors <- adjustcolor(line_colors, alpha.f = 0.2)  # 半透明填充（科研图标准）

# 雷达图绘制函数（科研配色版）
radar_plot <- function(df, metric_order) {
  metric_order <- intersect(metric_order, colnames(df))
  if (length(metric_order) == 0) {
    warning("无有效指标可绘制雷达图")
    return(NULL)
  }
  
  # 安全标准化函数
  safe_normalize <- function(x) {
    x_range <- max(x) - min(x)
    if (x_range == 0 || is.na(x_range)) return(rep(0.5, length(x)))
    return((x - min(x)) / x_range)
  }
  
  # 数据标准化
  df_scaled <- df
  df_scaled[metric_order] <- lapply(df_scaled[metric_order], function(col) {
    col[is.na(col) | is.infinite(col)] <- 0
    safe_normalize(col)
  })
  
  # 准备雷达图数据
  df_radar <- rbind(
    rep(1, length(metric_order)),
    rep(0, length(metric_order)),
    df_scaled[metric_order]
  )
  
  # 绘图参数重置
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # 绘制雷达图（科研配色+中文支持）
  par(mfrow = c(1, 1), mar = c(0, 0, 2, 0), family = "SimHei")
  fmsb::radarchart(
    df_radar,
    axistype = 1,
    pcol = line_colors[1:nrow(df_scaled)],  # 适配网络数量
    pfcol = fill_colors[1:nrow(df_scaled)],
    plwd = 2.5,  # 线条加粗（科研图更清晰）
    plty = 1,
    cglcol = "gray90",  # 浅灰网格（不抢主体色）
    cglty = 1,
    cglwd = 1,
    axislabcol = "gray60",  # 坐标轴标签灰色（低调）
    title = "水稻非生物胁迫网络拓扑特征雷达图",
    vlcex = 0.9,  # 指标标签放大（更易读）
    na.itp = FALSE
  )
  
  # 添加图例（科研图标准位置）
  if (nrow(df_scaled) > 0) {
    legend("bottomright",  # 右下角不遮挡图形
           legend = df_scaled$Network, 
           col = line_colors[1:nrow(df_scaled)], 
           lwd = 2.5, 
           bty = "n",
           cex = 0.9,
           text.col = "gray30")  # 图例文字灰色（更协调）
  }
}

# 定义指标顺序
metric_order <- c("Nodes", "Edges", "Avg_Degree", "Clustering_Coeff", "Avg_Path_Length")

# 输出PDF（Cairo支持中文+科研配色）
CairoPDF(
  "stress_network_topology_radar.pdf",
  width = 10,
  height = 8,
  bg = "white",
  family = "SimHei"
)
radar_plot(topo_metrics, metric_order)
dev.off()

# 清除所有警告
assign("last.warning", NULL, envir = baseenv())