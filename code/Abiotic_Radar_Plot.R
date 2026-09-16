# ====================== Load required packages ======================
library(igraph)
library(fmsb)
library(Cairo)  # Handle Chinese rendering

# ====================== Basic utility functions (read edge table, build graph, topology calculation) ======================
# Read a CSV edge list
read_edge_list_from_csv <- function(file_path, weight_threshold = NULL) {
  edge_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  edge_data <- edge_data %>% rename(from = TF, to = Target, weight = XGB_Score)
  
  if (!is.null(weight_threshold)) {
    edge_data <- edge_data[abs(edge_data$weight) > weight_threshold, ]
  }
  return(edge_data)
}

# Build a directed regulatory network graph
create_graph_from_edge_list <- function(edge_list) {
  graph <- graph_from_data_frame(edge_list, directed = TRUE)
  E(graph)$weight <- edge_list$weight
  E(graph)$color <- ifelse(E(graph)$weight < 0, "blue", "red")
  return(graph)
}

# Single-network topology metrics calculation
network_topology_analysis <- function(graph) {
  num_nodes <- vcount(graph)
  num_edges <- ecount(graph)
  if (num_nodes == 0 || num_edges == 0) return(NULL)
  
  avg_degree <- mean(degree(graph))
  clustering_coefficient <- transitivity(graph, type = "global")
  avg_path_length <- mean_distance(graph, weights = NA)
  density <- edge_density(graph)
  
  return(list(
    num_nodes = num_nodes,
    num_edges = num_edges,
    avg_degree = avg_degree,
    clustering_coefficient = clustering_coefficient,
    avg_path_length = avg_path_length,
    density = density
  ))
}

# ====================== Safe extraction function ======================
safe_extract <- function(x, default = 0) {
  if (is.null(x) || is.na(x) || is.infinite(x)) return(default)
  return(x)
}

# ====================== Radar-chart plotting function (research color scheme) ======================
radar_plot <- function(df, metric_order, title = "Network Topology Radar Chart") {
  metric_order <- intersect(metric_order, colnames(df))
  if (length(metric_order) == 0) {
    warning("No valid metrics to plot the radar chart")
    return(NULL)
  }
  
  # Safe normalization function
  safe_normalize <- function(x) {
    x[is.na(x) | is.infinite(x)] <- 0
    x_range <- max(x) - min(x)
    if (x_range == 0 || is.na(x_range)) return(rep(0.5, length(x)))
    return((x - min(x)) / x_range)
  }
  
  # Data normalization
  df_scaled <- df
  df_scaled[metric_order] <- lapply(df_scaled[metric_order], function(col) {
    safe_normalize(col)
  })
  
  # Prepare the radar-chart data
  df_radar <- rbind(
    rep(1, length(metric_order)),
    rep(0, length(metric_order)),
    df_scaled[metric_order]
  )
  
  # Research color scheme
  n_networks <- nrow(df_scaled)
  line_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442")[1:n_networks]
  fill_colors <- adjustcolor(line_colors, alpha.f = 0.2)
  
  # Reset plotting parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Draw the radar chart
  par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
  
  fmsb::radarchart(
    df_radar,
    axistype = 1,
    pcol = line_colors,
    pfcol = fill_colors,
    plwd = 2.5,
    plty = 1,
    cglcol = "gray90",
    cglty = 1,
    cglwd = 1,
    axislabcol = "gray50",
    title = title,
    vlcex = 1.0,
    na.itp = FALSE,
    caxislabels = c("0.0", "0.25", "0.5", "0.75", "1.0")
  )
  
  # Add a legend
  if (n_networks > 0) {
    legend("bottomright",
           legend = df_scaled$Network,
           col = line_colors,
           lwd = 2.5,
           bty = "n",
           cex = 0.9,
           text.col = "gray30")
  }
}

# ====================== Main program entry ======================
# 1. Parameter settings
weight_threshold <- NULL

# 2. Read the edge lists of the three abiotic-stress XGB fusions
cat("===== Starting to read the data files =====\n")
Alkalinity_edges <- read_edge_list_from_csv("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv", weight_threshold)
Aridity_edges <- read_edge_list_from_csv("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv", weight_threshold)
Cold_edges <- read_edge_list_from_csv("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv", weight_threshold)

# 3. Build the network graph objects
cat("===== Building the network graphs =====\n")
graph_Alkalinity <- create_graph_from_edge_list(Alkalinity_edges)
graph_Aridity <- create_graph_from_edge_list(Aridity_edges)
graph_Cold <- create_graph_from_edge_list(Cold_edges)

# 4. Compute the network topology metrics
cat("===== Computing network topology metrics =====\n")
Alkalinity_metrics <- network_topology_analysis(graph_Alkalinity)
Aridity_metrics <- network_topology_analysis(graph_Aridity)
Cold_metrics <- network_topology_analysis(graph_Cold)

# ====================== Radar-chart data preparation ======================
topo_metrics <- data.frame(
  Network = c("Alkalinity", "Aridity", "Cold"),
  Nodes = c(
    safe_extract(Alkalinity_metrics$num_nodes),
    safe_extract(Aridity_metrics$num_nodes),
    safe_extract(Cold_metrics$num_nodes)
  ),
  Edges = c(
    safe_extract(Alkalinity_metrics$num_edges),
    safe_extract(Aridity_metrics$num_edges),
    safe_extract(Cold_metrics$num_edges)
  ),
  Avg_Degree = c(
    safe_extract(Alkalinity_metrics$avg_degree),
    safe_extract(Aridity_metrics$avg_degree),
    safe_extract(Cold_metrics$avg_degree)
  ),
  Clustering_Coeff = c(
    safe_extract(Alkalinity_metrics$clustering_coefficient),
    safe_extract(Aridity_metrics$clustering_coefficient),
    safe_extract(Cold_metrics$clustering_coefficient)
  ),
  Avg_Path_Length = c(
    safe_extract(Alkalinity_metrics$avg_path_length),
    safe_extract(Aridity_metrics$avg_path_length),
    safe_extract(Cold_metrics$avg_path_length)
  ),
  Density = c(
    safe_extract(Alkalinity_metrics$density),
    safe_extract(Aridity_metrics$density),
    safe_extract(Cold_metrics$density)
  ),
  check.names = FALSE
)

# Print the topology metrics for inspection
cat("\n===== Topology metrics summary =====\n")
print(topo_metrics)

# ====================== Generate the radar-chart PDF (basic version) ======================
cat("\n===== Starting to draw the radar chart =====\n")

# Define the metric order (can be adjusted as needed)
metric_order <- c("Nodes", "Edges", "Avg_Degree", "Clustering_Coeff", "Avg_Path_Length")

# Use CairoPDF output
CairoPDF(
  "abiotic_stress_network_topology_radar.pdf",
  width = 10,
  height = 8,
  bg = "white"
)

# Draw the radar chart
radar_plot(
  topo_metrics, 
  metric_order, 
  title = "Network Topology Radar Chart of Rice under Abiotic Stresses"
)

dev.off()
cat("Radar chart saved as: abiotic_stress_network_topology_radar.pdf\n")

# ====================== Generate the radar chart (version with density) ======================
cat("\n===== Starting to draw the radar chart (with density) =====\n")

metric_order_density <- c("Nodes", "Edges", "Avg_Degree", "Clustering_Coeff", "Avg_Path_Length", "Density")

CairoPDF(
  "abiotic_stress_network_topology_radar_with_density.pdf",
  width = 10,
  height = 8,
  bg = "white"
)

radar_plot(
  topo_metrics, 
  metric_order_density, 
  title = "Network Topology Radar Chart of Rice under Abiotic Stresses (with Density)"
)

dev.off()
cat("Radar chart with density saved as: abiotic_stress_network_topology_radar_with_density.pdf\n")

# ====================== Save the topology metrics table ======================
cat("\n===== Saving the topology metrics table =====\n")
write.csv(topo_metrics, "abiotic_stress_network_topology_metrics.csv", row.names = FALSE)
cat("Topology metrics table saved as: abiotic_stress_network_topology_metrics.csv\n")

# ====================== Completion notice ======================
cat("\n===== All analyses completed! =====\n")
cat("Generated file list:\n")
cat("1. abiotic_stress_network_topology_radar.pdf\n")
cat("2. abiotic_stress_network_topology_radar_with_density.pdf\n")
cat("3. abiotic_stress_network_topology_metrics.csv\n")
