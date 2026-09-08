# ====================== Load required packages ======================
library(igraph)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

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

# Single-network topology metrics calculation (raw degree vector added for plotting)
network_topology_analysis <- function(graph) {
  num_nodes <- vcount(graph)
  num_edges <- ecount(graph)
  if (num_nodes == 0 || num_edges == 0) return(NULL)
  
  avg_degree <- mean(degree(graph))
  clustering_coefficient <- transitivity(graph, type = "global")
  avg_path_length <- mean_distance(graph, weights = NA)
  density <- edge_density(graph)
  degree_dist <- degree_distribution(graph)
  degree_raw <- degree(graph)
  
  return(list(
    num_nodes = num_nodes,
    num_edges = num_edges,
    avg_degree = avg_degree,
    clustering_coefficient = clustering_coefficient,
    avg_path_length = avg_path_length,
    density = density,
    degree_distribution = degree_dist,
    degree_raw = degree_raw
  ))
}

# Generate a random GNM network of the same size (preserving the positive/negative weight distribution)
generate_random_network <- function(original_graph) {
  num_nodes <- vcount(original_graph)
  num_edges <- ecount(original_graph)
  random_graph <- sample_gnm(num_nodes, num_edges, directed = TRUE)
  
  original_weights <- E(original_graph)$weight
  positive_weights <- original_weights[original_weights > 0]
  negative_weights <- original_weights[original_weights < 0]
  pos_ratio <- length(positive_weights) / length(original_weights)
  
  num_pos <- round(num_edges * pos_ratio)
  num_neg <- num_edges - num_pos
  
  random_pos_weights <- sample(positive_weights, num_pos, replace = TRUE)
  random_neg_weights <- sample(negative_weights, num_neg, replace = TRUE)
  random_weights <- c(random_pos_weights, random_neg_weights)
  random_weights <- sample(random_weights)
  
  E(random_graph)$weight <- random_weights
  E(random_graph)$color <- ifelse(E(random_graph)$weight < 0, "blue", "red")
  return(random_graph)
}

# Multiple-permutation random networks, output metric mean + SD (solves the instability of a single random run)
generate_random_network_stats <- function(original_graph, n_perm = 20) {
  res_list <- map(1:n_perm, ~{
    g_rand <- generate_random_network(original_graph)
    network_topology_analysis(g_rand)
  }) %>% discard(is.null)
  
  metrics_df <- map_dfr(res_list, ~tibble(
    avg_degree = .x$avg_degree,
    clustering_coefficient = .x$clustering_coefficient,
    avg_path_length = .x$avg_path_length,
    density = .x$density
  ))
  
  stat_summary <- metrics_df %>% 
    summarise(across(everything(), list(mean = mean, sd = sd)))
  return(stat_summary)
}

# Single-stress network comparison: real network + random network permutation statistics
compare_networks <- function(original_graph, name, n_perm = 20) {
  cat(paste0("===== Analyzing ", name, " network, permutation count: ", n_perm, " =====\n"))
  original_result <- network_topology_analysis(original_graph)
  random_stats <- generate_random_network_stats(original_graph, n_perm)
  
  out <- list(
    name = name,
    original = original_result,
    random_mean = random_stats %>% select(ends_with("_mean")),
    random_sd = random_stats %>% select(ends_with("_sd"))
  )
  return(out)
}

# ====================== Visualization function 1: multi-stress topology metrics error-bar bar chart ======================
plot_all_topology_bar <- function(comparison_list, save_name = "All_Stress_Topology_Bar") {
  df_all <- map_dfr(comparison_list, function(x) {
    obs <- tibble(
      Stress = x$name,
      Network = "Original",
      avg_degree = x$original$avg_degree,
      clustering_coefficient = x$original$clustering_coefficient,
      avg_path_length = x$original$avg_path_length,
      density = x$original$density,
      sd_val = 0
    )
    
    rand <- tibble(
      Stress = x$name,
      Network = "Random",
      avg_degree = x$random_mean$avg_degree_mean,
      clustering_coefficient = x$random_mean$clustering_coefficient_mean,
      avg_path_length = x$random_mean$avg_path_length_mean,
      density = x$random_mean$density_mean,
      sd_val = c(x$random_sd$avg_degree_sd,
                 x$random_sd$clustering_coefficient_sd,
                 x$random_sd$avg_path_length_sd,
                 x$random_sd$density_sd)
    )
    bind_rows(obs, rand)
  })
  
  df_long <- df_all %>%
    pivot_longer(
      cols = c(avg_degree, clustering_coefficient, avg_path_length, density),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = c("avg_degree", "clustering_coefficient", "avg_path_length", "density"),
        labels = c("Average Degree", "Clustering Coefficient", "Average Path Length", "Network Density")
      )
    )
  
  p <- ggplot(df_long, aes(x = Stress, y = Value, fill = Network)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = Value - sd_val, ymax = Value + sd_val),
      position = position_dodge(width = 0.8), width = 0.2
    ) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Original" = "#2E86AB", "Random" = "#A23B72")) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Topological Metrics of GRNs vs Random Networks",
      x = "Stress Condition",
      y = "Metric Value",
      fill = "Network Type"
    )
  
  ggsave(paste0(save_name, ".pdf"), plot = p, width = 12, height = 8, bg = "white")
  return(p)
}

# ====================== Visualization function 2: multi-stress overlaid log-log degree distribution (with power-law fit) ======================
plot_degree_ggplot <- function(comparison_list, save_name = "All_Stress_Degree_Dist") {
  degree_df <- map_dfr(comparison_list, function(x) {
    raw_deg <- x$original$degree_raw
    deg_tab <- table(raw_deg) %>% as.data.frame()
    colnames(deg_tab) <- c("k", "count")
    deg_tab$k <- as.numeric(as.character(deg_tab$k))
    deg_tab$freq <- deg_tab$count / sum(deg_tab$count)
    deg_tab$Stress <- x$name
    return(deg_tab)
  })
  
  fit_df <- map_dfr(comparison_list, function(x) {
    raw_deg <- x$original$degree_raw
    deg_tab <- table(raw_deg) %>% as.data.frame()
    colnames(deg_tab) <- c("k", "count")
    deg_tab$k <- as.numeric(as.character(deg_tab$k))
    deg_tab$freq <- deg_tab$count / sum(deg_tab$count)
    
    log_data <- deg_tab %>% filter(freq > 0) %>% mutate(lk = log10(k), lp = log10(freq))
    if (nrow(log_data) >= 3) {
      fit <- lm(lp ~ lk, data = log_data)
      gamma <- -coef(fit)[2]
      r2 <- round(summary(fit)$r.squared, 2)
      new_x <- tibble(lk = seq(min(log_data$lk), max(log_data$lk), length.out = 50))
      new_x$lp <- predict(fit, new_x)
      new_x$k <- 10^new_x$lk
      new_x$freq <- 10^new_x$lp
      new_x$Stress <- x$name
      new_x$label <- paste0("\u03b3 = ", round(gamma, 2), " R²=", r2)
      return(new_x)
    } else {
      return(tibble())
    }
  })
  
  p <- ggplot() +
    geom_point(data = degree_df, aes(x = k, y = freq, color = Stress), size = 2, alpha = 0.7) +
    geom_line(data = fit_df, aes(x = k, y = freq, color = Stress), linetype = "dashed", linewidth = 1) +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Degree Distribution (Log-Log Scale) of Gene Regulatory Networks",
      x = "Node Degree k (log10)",
      y = "Frequency P(k) (log10)",
      color = "Stress"
    )
  
  ggsave(paste0(save_name, ".pdf"), plot = p, width = 10, height = 7, bg = "white")
  return(p)
}

# ====================== Main program entry ======================
# 1. Parameter settings
weight_threshold <- NULL
perm_times <- 20 # Number of random-network permutation repeats

# 2. Read the edge lists of the three stresses (modify to your local CSV paths)
alkalinity_edges <- read_edge_list_from_csv("Alkalinity_XGBFusion_3DCEMA_GENIE3_GRNBoost_IGEGRNs.csv", weight_threshold)
aridity_edges <- read_edge_list_from_csv("Aridity_XGBFusion_3DCEMA_GRNBoost_IGEGRNs_Kboost.csv", weight_threshold)
cold_edges <- read_edge_list_from_csv("Cold_XGBFusion_3DCEMA_GENIE3_GRNBoost.csv", weight_threshold)

# 3. Build the network graph objects
graph_alkalinity <- create_graph_from_edge_list(alkalinity_edges)
graph_aridity <- create_graph_from_edge_list(aridity_edges)
graph_cold <- create_graph_from_edge_list(cold_edges)

# Consolidate lists for batch plotting
graph_all <- list(Alkalinity = graph_alkalinity, Aridity = graph_aridity, Cold = graph_cold)
edge_all <- list(Alkalinity = alkalinity_edges, Aridity = aridity_edges, Cold = cold_edges)

# 4. Real-network vs random-network topology comparison (with multiple-permutation statistics)
alkalinity_comparison <- compare_networks(graph_alkalinity, "Alkalinity", n_perm = perm_times)
aridity_comparison <- compare_networks(graph_aridity, "Aridity", n_perm = perm_times)
cold_comparison <- compare_networks(graph_cold, "Cold", n_perm = perm_times)
comp_all <- list(alkalinity_comparison, aridity_comparison, cold_comparison)

# 5. Batch generate all 4 categories of PDF charts
cat("\n===== Starting to plot the topology-metrics comparison bar chart =====\n")
p1 <- plot_all_topology_bar(comp_all)
print(p1)

cat("\n===== Starting to plot the log-log degree distribution =====\n")
p2 <- plot_degree_ggplot(comp_all)
print(p2)

cat("\nAll charts plotted; PDF files saved to the working directory!\n")
