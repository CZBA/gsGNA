# ============================================================================
# gsGNA | Cross-stress GRN conservation and regulatory rewiring analysis
# Input: FINAL locked XGBoost selected-combination Top 10% GRNs
# Output: network topology metrics, hub summaries, edge-conservation tests,
#         TF-target rewiring analyses, and publication-ready comparison figures.
# Required R packages: data.table, igraph, ggplot2, VennDiagram
# ============================================================================

required <- c("data.table", "igraph", "ggplot2", "VennDiagram", "fmsb")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install the required R packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(fmsb)
})

# ---- 0. Settings ------------------------------------------------------------
input_root <- Sys.getenv("GSGNA_INPUT_ROOT", unset=".")
ensemble_dir <- file.path(input_root, "gsGNA_FixedTop10")
stresses <- c("Drought", "Alkaline", "Cold")
network_files <- setNames(file.path(
  ensemble_dir, "FINAL_Locked_Test_And_Network", stresses, "Final_Ranked_Network",
  "XGBoost_Ranked_SelectedUnion_AllEdges.csv"
), stresses)
input_candidate_stage <- "top10"
network_edge_fraction <- 0.10
run_selected_union_comparison <- FALSE  # Optional comparison with the selected base-method candidates.
out_dir <- file.path(input_root, "gsGNA_Cross_Stress_Network_Analysis")
random_seed <- 20261027L
data_table_threads <- 8L
hub_fraction <- 0.10    # Top 10% of positive-outdegree TFs; ties kept
n_null <- 1000L         # TF-outdegree-preserving overlap null (pairwise/triple)
n_tf_boot <- 1000L      # TF-level bootstrap for core/periphery differences
min_tf_degree <- 10L    # per-TF stress-specific enrichment minimum outdegree
export_large_edge_lists <- TRUE
weight_threshold <- NULL   # Optional minimum score applied after rank selection.
path_exact_node_limit <- 3000L
path_n_sources <- 200L

stopifnot(is.finite(network_edge_fraction), network_edge_fraction > 0,
          network_edge_fraction <= 1, hub_fraction > 0, hub_fraction <= 1,
          n_null >= 1L, n_tf_boot >= 1L)

missing_final <- network_files[!file.exists(network_files)]
if (length(missing_final)) {
  stop("Missing ranked XGBoost networks:\n",
       paste(names(missing_final), missing_final, sep=" -> ", collapse="\n"))
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
set.seed(random_seed)
setDTthreads(threads=data_table_threads)

stress_code <- c(Drought = 1L, Alkaline = 2L, Cold = 4L)
pair_names <- combn(stresses, 2L, simplify = FALSE)
write_dt <- function(x, filename) {
  sep <- if (grepl("\\.tsv$", filename, ignore.case = TRUE)) "\t" else ","
  fwrite(x, file.path(out_dir, filename), sep = sep)
}
write_ids <- function(x, filename) writeLines(as.character(x), file.path(out_dir, filename))
safe_div <- function(a, b) ifelse(is.finite(b) & b > 0, a / b, NA_real_)
STRESS_COLORS <- c(Drought = "#2C967D", Alkaline = "#D97C6D", Cold = "#6097E7")
GROUP_COLORS <- c("Shared hub TFs" = "#C96A65", "Other TFs" = "#679FBD")


# ---- 1. Load the top-ranked edges --------------------------------------------
load_grn <- function(file) {
  x <- fread(file, select=c("TF", "Target", "EdgeWeight"), showProgress=FALSE)
  if (!all(c("TF", "Target", "EdgeWeight") %in% names(x)) || !nrow(x)) {
    stop("Expected a nonempty ranked network with TF, Target, and EdgeWeight: ", file)
  }
  x[, EdgeWeight := suppressWarnings(as.numeric(EdgeWeight))]
  if (any(!is.finite(x$EdgeWeight))) stop("Non-finite ranked edge weights: ", file)
  export_summary <- file.path(dirname(file), "Final_Network_Export_Summary.csv")
  if (file.exists(export_summary)) {
    audit <- fread(export_summary, showProgress=FALSE)
    required <- c("File", "Exported_Edges", "N_Union_Candidates")
    if (!all(required %in% names(audit)) || nrow(audit) != 1L ||
        !isTRUE(audit$File == basename(file)) || !isTRUE(audit$Exported_Edges == nrow(x)) ||
        !isTRUE(audit$N_Union_Candidates == nrow(x))) {
      stop("Expected the full ranked selected-union network: ", file)
    }
  }
  setorderv(x, c("EdgeWeight", "TF", "Target"), c(-1L, 1L, 1L))
  x <- x[seq_len(ceiling(nrow(x)*network_edge_fraction))]
  x[, `:=`(TF = toupper(trimws(as.character(TF))),
           Target = toupper(trimws(as.character(Target))),
           EdgeWeight = suppressWarnings(as.numeric(EdgeWeight)))]
  x <- x[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) &
           !TF %chin% c("NA", "NAN") & !Target %chin% c("NA", "NAN") &
           TF != Target & is.finite(EdgeWeight)]
  if (!is.null(weight_threshold)) {
    if (!is.numeric(weight_threshold) || length(weight_threshold) != 1L ||
        !is.finite(weight_threshold)) stop("weight_threshold must be a finite number or NULL")
    # XGBoost score is an edge confidence/importance; its sign does NOT show
    # activation/repression. Filter score >= cutoff, not abs(score).
    x <- x[EdgeWeight >= weight_threshold]
  }
  if (!nrow(x)) stop("No valid TF-target links after filtering: ", file)
  setorder(x, -EdgeWeight)
  x <- unique(x, by = c("TF", "Target"))
  setkey(x, TF, Target)
  x
}
read_selected_methods <- function(path) {
  if (!file.exists(path)) stop("Missing locked model file: ", path)
  models <- fread(path, showProgress=FALSE)
  required <- c("Dataset", "Strategy", "Combination", "Input_Stage")
  if (!all(required %in% names(models)) || nrow(models) != length(stresses) ||
      anyDuplicated(models$Dataset) || !setequal(models$Dataset, stresses)) {
    stop("Expected one locked model with strategy, combination and input stage per stress: ", path)
  }
  if (anyNA(models$Strategy) || any(models$Strategy != "XGBoost") ||
      anyNA(models$Input_Stage) || any(models$Input_Stage != input_candidate_stage)) {
    stop("Expected locked XGBoost models at input stage ", input_candidate_stage, ": ", path)
  }
  combinations <- setNames(as.character(models$Combination), models$Dataset)[stresses]
  selected <- lapply(combinations, function(x) trimws(strsplit(x, "+", fixed=TRUE)[[1L]]))
  valid_methods <- c("GENIE3", "Kboost", "GRNBoost2", "3DCEMA", "DeepRIG", "IGEGRNs")
  if (any(!vapply(selected, function(x) {
    length(x) > 0L && all(x %in% valid_methods) && !anyDuplicated(x)
  }, logical(1)))) stop("Invalid base-algorithm combination in: ", path)
  selected
}

selection_file <- file.path(ensemble_dir, "FINAL_PerStress_Models_LOCKED.csv")
selected_methods <- if (file.exists(selection_file) || run_selected_union_comparison) {
  read_selected_methods(selection_file)
}
writeLines(c(
  paste("Ranked input networks:", paste(network_files, collapse=" | ")),
  sprintf("Network edge fraction: %g%%", 100*network_edge_fraction),
  sprintf("Hub fraction: %g%% of positive-outdegree TFs (ties included)", 100*hub_fraction),
  paste("Random seed:", random_seed),
  paste("Overlap null replicates:", n_null),
  paste("TF bootstrap replicates:", n_tf_boot),
  "Shared observed TF and target IDs provide a proxy for candidate eligibility.",
  "The overlap null preserves TF outdegree."
), file.path(out_dir, "00_Run_Configuration.txt"))

networks <- lapply(network_files, load_grn)
input_audit <- data.table(
  Stress = stresses,
  File = unname(network_files),
  Selected_Methods = if (is.null(selected_methods)) NA_character_ else
    vapply(selected_methods, paste, collapse="+", FUN.VALUE=character(1)),
  Input_Stage = input_candidate_stage,
  Output_Stage = "Final_Ranked_Selected_Union",
  Network_Edge_Fraction = network_edge_fraction,
  Unique_TFs = vapply(networks, function(x) uniqueN(x$TF), integer(1)),
  Unique_Targets = vapply(networks, function(x) uniqueN(x$Target), integer(1)),
  Unique_Edges = vapply(networks, nrow, integer(1))
)
fwrite(input_audit, file.path(out_dir, "00_Final_Locked_Input_Audit.csv"))
cat("[Loading] Networks loaded: ",
    paste(sprintf("%s=%s", stresses,
                  format(vapply(networks, nrow, integer(1L)), big.mark = ",")),
          collapse = "; "), "\n", sep="")

# ---- 1b. Original topology and radar, corrected definitions ----------------
# Directed TF -> Target binary regulatory network:
# Nodes: unique TF and Target identifiers; edges: unique directed pairs.
# Average total degree = 2*E/V, for compatibility with the original script.
# Directed density = E/(V*(V-1)), self-loops already excluded.
# Global clustering: computed on simple UNDIRECTED projection; there is no
# signed-activation interpretation of the XGB score.
# Path length: largest WEAKLY connected component, converted to undirected;
# exact for small components, sampling sources for big components.
# These measures are DESCRIPTIVE and do not establish architecture conservation.

network_topology_analysis <- function(edge_df, stress) {
  g <- igraph::graph_from_data_frame(edge_df[, .(TF, Target)], directed = TRUE)
  nv <- igraph::vcount(g); ne <- igraph::ecount(g)
  if (nv < 2L || ne < 1L) stop("Insufficient network: ", stress)
  ug <- igraph::as_undirected(g, mode = "collapse")
  ug <- igraph::simplify(ug, remove.multiple = TRUE, remove.loops = TRUE)
  weak <- igraph::components(ug)
  giant_id <- which.max(weak$csize)
  giant_vertices <- which(weak$membership == giant_id)
  giant <- igraph::induced_subgraph(ug, giant_vertices)
  if (igraph::vcount(giant) < 2L) {
    apl <- NA_real_; path_method <- "Not_estimable"
  } else if (igraph::vcount(giant) <= path_exact_node_limit) {
    apl <- igraph::mean_distance(giant, directed = FALSE,
                                unconnected = FALSE, weights = NA)
    path_method <- "Exact_largest_weak_component"
  } else {
    sample_source <- sample.int(igraph::vcount(giant),
                                min(path_n_sources, igraph::vcount(giant)))
    dd <- igraph::distances(giant, v = sample_source,
                           to = igraph::V(giant), mode = "all", weights = NA)
    dd <- as.matrix(dd)
    apl <- mean(dd[is.finite(dd) & dd > 0])
    path_method <- paste0("Sampled_", length(sample_source), "_source_vertices")
    rm(dd)
  }
  tri <- suppressWarnings(igraph::transitivity(ug, type = "global"))
  if (!is.finite(tri)) tri <- NA_real_
  ans <- data.table(
    Network = stress,
    Nodes = nv,
    Edges = ne,
    TFs_Observed = uniqueN(edge_df$TF),
    Targets_Observed = uniqueN(edge_df$Target),
    Avg_Degree = 2 * ne / nv,
    Avg_TF_OutDegree = ne / uniqueN(edge_df$TF),
    Clustering_Coeff = tri,
    Avg_Path_Length = apl,
    Density = ne / (as.double(nv) * (nv - 1)),
    Largest_Weak_Component_Nodes = length(giant_vertices),
    Largest_Weak_Component_Fraction = length(giant_vertices) / nv,
    Number_Weak_Components = weak$no,
    Path_Length_Method = path_method
  )
  rm(g, ug, giant)
  gc(verbose = FALSE)
  ans
}

topo_metrics <- rbindlist(lapply(stresses, function(s) {
  cat("[Topology] ", s, " ...\n", sep = "")
  network_topology_analysis(networks[[s]], s)
}))
write_dt(topo_metrics, "abiotic_stress_network_topology_metrics.csv")
cat("\n==== Topology metrics ====\n"); print(topo_metrics)

# Radar is a visualization only: with just 3 GRNs, min-max scaling can make
# small differences look dramatic. Each axis is independently scaled 0-1.
# A constant metric becomes 0.5. The values in CSV are the original units.
radar_plot <- function(df, metrics, filename, title) {
  stopifnot(length(metrics) >= 3L, all(metrics %in% names(df)))
  normalize_col <- function(v) {
    v <- as.numeric(v)
    v[!is.finite(v)] <- NA_real_
    if (all(is.na(v))) return(rep(0.5, length(v)))
    v[is.na(v)] <- stats::median(v, na.rm = TRUE)
    r <- range(v)
    if (diff(r) <= .Machine$double.eps * max(1, abs(r)))
      return(rep(0.5, length(v)))
    (v - r[1L]) / diff(r)
  }
  actual <- as.data.frame(lapply(df[, ..metrics], normalize_col))
  radar_input <- as.data.frame(rbind(rep(1, length(metrics)),
                                     rep(0, length(metrics)),
                                     as.matrix(actual)))
  names(radar_input) <- metrics
  cols <- c("#0072B2", "#D55E00", "#009E73")
  n <- nrow(df)
  opened <- FALSE
  if (requireNamespace("Cairo", quietly = TRUE)) {
    try({Cairo::CairoPDF(file.path(out_dir, filename), width = 10, height = 8,
                        bg = "white"); opened <- TRUE}, silent = TRUE)
  }
  if (!opened) {
    grDevices::pdf(file.path(out_dir, filename), width = 10, height = 8,
                   useDingbats = FALSE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); grDevices::dev.off()}, add = TRUE)
  par(mar = c(1, 1, 4, 1))
  fmsb::radarchart(
    radar_input, maxmin = TRUE, axistype = 1,
    pcol = cols[seq_len(n)],
    pfcol = grDevices::adjustcolor(cols[seq_len(n)], alpha.f = .18),
    plwd = 2.5, plty = 1, seg = 4,
    cglcol = "gray85", cglty = 1, cglwd = 1,
    axislabcol = "black", vlcex = .85,
    caxislabels = c("0", "0.25", "0.5", "0.75", "1"),
    title = title
  )
  legend("bottomright", legend = df$Network,
         col = cols[seq_len(n)], lwd = 2.5, bty = "n", cex = .85)
}

radar_plot(
  topo_metrics,
  c("Nodes", "Edges", "Avg_Degree", "Clustering_Coeff", "Avg_Path_Length"),
  "abiotic_stress_network_topology_radar.pdf",
  "Rice abiotic stress GRN topology (min-max scaled)"
)
radar_plot(
  topo_metrics,
  c("Nodes", "Edges", "Avg_Degree", "Clustering_Coeff", "Avg_Path_Length", "Density"),
  "abiotic_stress_network_topology_radar_with_density.pdf",
  "Rice abiotic stress GRN topology with network density"
)
writeLines(c(
  "Radar axes independently min-max scaled over only three networks; descriptive only.",
  "Average total degree=2E/V; average TF outdegree=E/observed TF count.",
  "Global clustering on simple undirected projection (direction ignored).",
  "Shortest path on undirected largest weakly connected component; method in CSV.",
  "XGBoost scores are confidence/ranking values, NOT activation/repression signs."
), file.path(out_dir, "00_Topology_Metric_Definitions.txt"))
cat("[Topology] Topology CSV and two radar PDFs generated.\n")

# ---- 2. Original Top 10% hub calling and outputs -----------------------------
get_degree_hubs <- function(x, fraction) {
  d <- x[, .(degree = .N), by = TF]
  setorder(d, -degree, TF)
  n_top <- max(1L, as.integer(ceiling(nrow(d) * fraction)))
  cutoff <- d$degree[n_top]
  d[, `:=`(Degree_Rank = seq_len(.N),
           Degree_Rank_Tied = frank(-degree, ties.method = "min"),
           Degree_Fold = degree / mean(degree),
           Hub_Top10 = degree >= cutoff)]
  list(degrees = d, hubs = d[Hub_Top10 == TRUE],
       diag = data.table(N_TFs_Observed = nrow(d), N_Top10_Nominal = n_top,
                         N_Hub_Actual = sum(d$Hub_Top10),
                         Actual_Percent = 100 * mean(d$Hub_Top10),
                         Cutoff_Degree = cutoff,
                         Boundary_Tie_Count = sum(d$degree == cutoff),
                         Mean_Degree = mean(d$degree),
                         Median_Degree = median(d$degree),
                         Total_Edges = sum(d$degree)))
}
degree_results <- lapply(networks, get_degree_hubs, fraction = hub_fraction)
hub_sets <- lapply(degree_results, function(z) z$hubs$TF)
common_hubs <- Reduce(intersect, hub_sets)
hub_only <- lapply(stresses, function(s) {
  setdiff(hub_sets[[s]], unlist(hub_sets[setdiff(stresses, s)], use.names = FALSE))
})
names(hub_only) <- stresses
hub_summary <- rbindlist(lapply(stresses, function(s) {
  cbind(data.table(Stress = s), degree_results[[s]]$diag)
}))
write_dt(hub_summary, "Three_Stress_Top10_Threshold_Summary.csv")
for (s in stresses) {
  write_dt(degree_results[[s]]$degrees, paste0(s, "_All_TF_Degrees.csv"))
  write_dt(degree_results[[s]]$hubs, paste0(s, "_HubTF_Top10_Details.csv"))
  write_ids(hub_sets[[s]], paste0(s, "_HubTFs_Top10.txt"))
  write_ids(hub_only[[s]], paste0(s, "_Only_HubTFs_Top10.txt"))
}
write_ids(common_hubs, "Three_Stress_Common_HubTFs_Top10.txt")
for (p in pair_names) {
  write_ids(intersect(hub_sets[[p[1L]]], hub_sets[[p[2L]]]),
            paste0(paste(p, collapse = "_"), "_Shared_HubTFs_Top10.txt"))
}
venn_grob <- venn.diagram(
  x = hub_sets, filename = NULL, category.names = stresses, disable.logging = TRUE,
  fill = c("#82B5E4", "#E8ABBA", "#A7D8AE"),
  alpha = 0.55, col = "transparent", cex = 0.85, cat.cex = 0.95
)
pdf(file.path(out_dir, "Three_Stress_HubTF_Venn_Top10.pdf"), width = 8, height = 8)
grid.newpage(); grid.draw(venn_grob); dev.off()

# Original hub TF --> targets, including shared-hub pooled targets and exact pairs.
hub_edges <- lapply(stresses, function(s) networks[[s]][TF %chin% hub_sets[[s]]])
names(hub_edges) <- stresses
common_hub_edges <- lapply(stresses, function(s) networks[[s]][TF %chin% common_hubs])
names(common_hub_edges) <- stresses
for (s in stresses) {
  write_dt(hub_edges[[s]], paste0(s, "_HubTF_Targets_Top10.tsv"))
  write_dt(common_hub_edges[[s]], paste0(s, "_CommonHubTF_Targets_Top10.tsv"))
  write_dt(networks[[s]][TF %chin% hub_only[[s]]],
           paste0(s, "_StressSpecificHubTF_Targets_Top10.tsv"))
}
hub_target_sets <- lapply(common_hub_edges, function(x) unique(x$Target))
write_ids(Reduce(intersect, hub_target_sets),
          "CommonHubTF_PooledTargets_All3_Top10.txt")
for (s in stresses) {
  others <- setdiff(stresses,s)
  write_ids(setdiff(hub_target_sets[[s]],
                    unique(unlist(hub_target_sets[others],use.names=FALSE))),
            paste0("CommonHubTF_PooledTargets_",s,"_Only_Top10.txt"))
}
common_hub_pair_full <- Reduce(function(a, b) merge(a, b, by = c("TF", "Target")),
  lapply(common_hub_edges, function(x) x[, .(TF, Target)]))
write_dt(common_hub_pair_full, "CommonHubTF_SamePairs_All3_Top10.tsv")

# Ranks of common hubs (originally from complete networks, ties min-ranked).
rank_tables <- lapply(stresses, function(s) {
  x <- copy(degree_results[[s]]$degrees[, .(TF, degree, Degree_Rank_Tied)])
  setnames(x, c("degree", "Degree_Rank_Tied"),
           paste0(s, c("_Degree", "_Rank")))
  x
})
common_hub_ranks <- Reduce(function(a,b) merge(a,b,by="TF"), rank_tables)
common_hub_ranks <- common_hub_ranks[TF %chin% common_hubs]
write_dt(common_hub_ranks, "Three_Stress_Common_Hub_Ranks_Top10.csv")
if (nrow(common_hub_ranks) >= 3L) {
  mat <- as.matrix(common_hub_ranks[, paste0(stresses, "_Degree"), with=FALSE])
  cors <- suppressWarnings(cor(mat, method = "spearman"))
  write_dt(data.table(Stress = rownames(cors), cors, check.names=FALSE),
           "Three_Stress_Common_Hub_Spearman.csv")
}
# Reproduce the original top-10%-hub rank visualization.
plot_degrees <- rbindlist(lapply(stresses,function(s) {
  d <- copy(degree_results[[s]]$degrees)
  d[,Stress:=s]
  d
}))
p_deg <- ggplot(plot_degrees,aes(x=Degree_Rank,y=degree,fill=Hub_Top10)) +
  geom_col(width=1) + facet_wrap(~Stress,nrow=1,scales="free_x") +
  scale_fill_manual(values=c("FALSE"="#80A7C0","TRUE"="#CB665A"),
                    labels=c("FALSE"="Other TFs","TRUE"="Top-10% hubs")) +
  labs(x="TF rank (descending outdegree)",y="Outdegree (unique targets)",
       title="Top-10% hub TFs in three stress-specific GRNs",fill=NULL) +
  theme_bw(base_size=11)+theme(text=element_text(color="black"),
                                plot.title=element_text(hjust=.5,color="black"),
                                axis.text=element_text(color="black"),
                                axis.title=element_text(color="black"),
                                legend.text=element_text(color="black"),
                                legend.position="bottom")
ggsave(file.path(out_dir,"Three_Stress_TF_Degree_Ranking_Top10.pdf"),
       p_deg,width=13,height=4.5)
cat("[Hubs] Hubs selected. All-three common hubs: ", length(common_hubs), "\n", sep="")

# ---- 3. OPTIONAL: locked selected-method TOP10 input union vs final network ---
# This optional descriptive comparison uses ONLY each stress's selected XGBoost
# method subset recorded during model selection.
find_candidate_file <- function(m, stress) {
  dirs <- if (m == "GRNBoost2") c("GRNBoost2", "GRNBoost") else m
  f <- file.path(input_root, dirs, paste0(stress, "_", input_candidate_stage, "_edges.csv"))
  found <- f[file.exists(f)]
  if (length(found)) found[[1L]] else NA_character_
}
get_top_by_degree <- function(d, fraction = hub_fraction) {
  positive <- sort(d[is.finite(d) & d > 0], decreasing=TRUE)
  if (!length(positive)) return(character())
  cutoff <- positive[max(1L, as.integer(ceiling(length(positive)*fraction)))]
  names(d)[is.finite(d) & d >= cutoff & d > 0]
}
if (run_selected_union_comparison) {
  candidates <- setNames(lapply(stresses, function(s)
    vapply(selected_methods[[s]], find_candidate_file, character(1L), stress=s)),
    stresses)
  if (anyNA(unlist(candidates, use.names=FALSE))) {
    warning("Some selected-method Top10 candidate files are missing. Skipping the optional comparison.")
  } else {
    candidate_hubs <- list()
    compare_tables <- list()
    for (s in stresses) {
      raw <- rbindlist(lapply(candidates[[s]], function(f) {
        z <- fread(f, select=1:2, showProgress=FALSE)
        setnames(z, names(z), c("TF", "Target"))
        z[, `:=`(TF=toupper(trimws(as.character(TF))),
                 Target=toupper(trimws(as.character(Target))))]
        z[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) & TF!=Target]
      }))
      raw <- unique(raw, by=c("TF", "Target"))
      setkey(raw, TF, Target)
      missing_edges <- raw[networks[[s]], on=.(TF, Target), which=TRUE]
      if (anyNA(missing_edges))
        warning(s, ": ", sum(is.na(missing_edges)),
                " FINAL edges not in selected Top10 input union; check output/input versions.")
      cand <- raw[,.(Candidate_Degree=.N),by=TF]
      fin <- degree_results[[s]]$degrees[,.(TF,Final_Degree=degree)]
      cmp <- merge(cand,fin,by="TF",all=TRUE)
      setnafill(cmp,fill=0L,cols=c("Candidate_Degree","Final_Degree"))
      cmp[, Stress := s]
      compare_tables[[s]] <- cmp
      candidate_hubs[[s]] <- get_top_by_degree(setNames(cmp$Candidate_Degree,cmp$TF))
      write_dt(cmp, paste0(s,"_Candidate_vs_Final_Degrees.csv"))
      if (any(cmp$Final_Degree > cmp$Candidate_Degree)) {
        warning(s, ": some final edges absent in its locked selected-method Top10 candidate union; verify data versions.")
      }
    }
    cand_common <- Reduce(intersect,candidate_hubs)
    write_ids(cand_common,"Candidate_ThreeStress_Common_Top10_Hubs.txt")
    write_ids(intersect(cand_common,common_hubs),"Candidate_Final_Shared_Common_Hubs.txt")
    write_dt(rbindlist(lapply(stresses,function(s) {
      data.table(Stress=s,Candidate_Hubs=length(candidate_hubs[[s]]),
                 Final_Hubs=length(hub_sets[[s]]),
                 Overlap=length(intersect(candidate_hubs[[s]],hub_sets[[s]])))
    })),"Candidate_vs_Final_Top10_Hub_Summary.csv")
  }
}

# ---- 4. SHARED ELIGIBILITY: identical TF and Target sets for all 3 ----------
shared_tfs <- Reduce(intersect, lapply(networks, function(x) unique(x$TF)))
shared_targets <- Reduce(intersect, lapply(networks, function(x) unique(x$Target)))
if (!length(shared_tfs) || !length(shared_targets))
  stop("No shared observed TF/target universe. Check identifier schemes.")
write_ids(sort(shared_tfs),"Shared_Eligible_TFs.txt")
write_ids(sort(shared_targets),"Shared_Eligible_Targets.txt")
shared_edges <- lapply(networks, function(x) {
  d <- x[TF %chin% shared_tfs & Target %chin% shared_targets, .(TF, Target)]
  setkey(d, TF, Target)
  d
})
if (any(vapply(shared_edges, nrow, integer(1L)) == 0L))
  stop("One stress has zero edges in the common TF/target scope; comparable analysis impossible.")
summary_universe <- rbindlist(lapply(stresses,function(s) {
  data.table(Stress=s, Raw_TFs=uniqueN(networks[[s]]$TF),
             Raw_Targets=uniqueN(networks[[s]]$Target),
             Raw_Edges=nrow(networks[[s]]),
             Shared_Eligible_TFs=length(shared_tfs),
             Shared_Eligible_Targets=length(shared_targets),
             Scoped_Edges=nrow(shared_edges[[s]]),
             Retained_Fraction=nrow(shared_edges[[s]])/nrow(networks[[s]]))
}))
write_dt(summary_universe,"01_Comparison_Universe_Summary.csv")
cat("[Universe] Shared TFs:",length(shared_tfs)," shared targets:",length(shared_targets),"\n")

# ---- 5. Edge conservation & stress-specific occupancy -----------------------
all_scoped <- rbindlist(lapply(stresses,function(s) {
  z <- copy(shared_edges[[s]])
  z[, `:=`(Stress=s, Stress_Code=stress_code[[s]])]
  z
}))
# The 3 networks have no duplicate (TF,Target) by construction.
occ <- all_scoped[, .(N_Stresses=.N, Stress_Mask=sum(Stress_Code)), by=.(TF,Target)]
setkey(occ,TF,Target)
mask_names <- c("Drought", "Alkaline", "Drought+Alkaline",
                "Cold", "Drought+Cold", "Alkaline+Cold", "All_three")
pattern_summary <- occ[, .(N_Edges=.N), by=.(Stress_Mask,N_Stresses)]
pattern_summary[, Stress_Pattern := mask_names[Stress_Mask]]
setcolorder(pattern_summary,c("Stress_Pattern","Stress_Mask","N_Stresses","N_Edges"))
setorder(pattern_summary,Stress_Mask)
write_dt(pattern_summary,"02_Edge_Occupancy_Seven_Patterns.csv")
if (export_large_edge_lists) {
  write_dt(occ[Stress_Mask==7L,.(TF,Target)],"02_Exactly_Conserved_Edges_All3.tsv")
  for (s in stresses) {
    write_dt(occ[Stress_Mask==stress_code[[s]],.(TF,Target)],
             paste0("02_",s,"_Unique_Eligible_Edges.tsv"))
  }
}
# Annotated edges have one row for each stress in which that edge occurs.
annotated <- occ[all_scoped,on=.(TF,Target)]

# Counts, overlap, retention and Jaccard for RAW and SHARED universes.
compare_edges <- function(a, b) {
  setkey(a,TF,Target)
  setkey(b,TF,Target)
  ab <- nrow(a[b, on=.(TF,Target), nomatch=0L])
  u <- nrow(a)+nrow(b)-ab
  c(Edges_A=nrow(a), Edges_B=nrow(b), Conserved_Edges=ab,
    Union_Edges=u, Jaccard=if(u) ab/u else NA_real_,
    Dice=if(nrow(a)+nrow(b)) 2*ab/(nrow(a)+nrow(b)) else NA_real_,
    A_to_B_Retention=if(nrow(a)) ab/nrow(a) else NA_real_,
    B_to_A_Retention=if(nrow(b)) ab/nrow(b) else NA_real_)
}
pairwise_edge <- rbindlist(lapply(pair_names,function(p) {
  rbindlist(lapply(c("Full_Network","Shared_Eligible_Universe"), function(scope) {
    z <- if(scope=="Full_Network")
      compare_edges(networks[[p[1L]]][,.(TF,Target)],networks[[p[2L]]][,.(TF,Target)])
    else compare_edges(shared_edges[[p[1L]]],shared_edges[[p[2L]]])
    cbind(data.table(Stress_A=p[1L],Stress_B=p[2L],Scope=scope),
          as.data.table(as.list(z)))
  }))
}))
write_dt(pairwise_edge,"03_Pairwise_Edge_Conservation.csv")
all3_count <- nrow(occ[Stress_Mask == 7L])
all3_fraction <- data.table(Stress=stresses,
  Scoped_Edges=vapply(shared_edges,nrow,integer(1L)),
  All3_Edges=all3_count)
all3_fraction[, All3_Retention := All3_Edges/Scoped_Edges]
write_dt(all3_fraction,"03_Three_Stress_Edge_Conservation.csv")

# ---- 6. Degree-controlled overlap null: sample targets by TF ---------------
# With a fixed TF, two random target subsets of sizes a and b in N allowed
# targets have overlap X ~ Hypergeometric(m=b, n=N-b, k=a).
# This controls TF outdegree and self-loops; target indegree is NOT preserved.
tf_degree_table <- data.table(TF=shared_tfs)
tf_degree_table[, Allowed_Targets := length(shared_targets)-as.integer(TF %chin% shared_targets)]
for(s in stresses) {
  col <- paste0(s,"_Degree")
  tmp <- shared_edges[[s]][,.(Value=.N),by=TF]
  tf_degree_table[, (col) := 0L]
  tf_degree_table[tmp, (col) := i.Value, on="TF"]
}
write_dt(tf_degree_table,"04_CommonUniverse_PerTF_Degrees.csv")

null_rows <- list()
for (p in pair_names) {
  a <- p[1L]; b <- p[2L]
  dA <- tf_degree_table[[paste0(a,"_Degree")]]
  dB <- tf_degree_table[[paste0(b,"_Degree")]]
  N <- tf_degree_table$Allowed_Targets
  active <- dA>0 & dB>0 & N>0
  exp_overlap <- sum(as.double(dA[active])*dB[active]/N[active])
  sim <- replicate(n_null,
    sum(rhyper(sum(active),m=dB[active],n=N[active]-dB[active],k=dA[active])))
  observed <- nrow(shared_edges[[a]][shared_edges[[b]],on=.(TF,Target),nomatch=0L])
  null_rows[[length(null_rows)+1L]] <- data.table(
    Stress_A=a,Stress_B=b,Overlap_Type="Pairwise",Observed=observed,
    Null_Mean=mean(sim),Null_SD=sd(sim),Analytic_Expected=exp_overlap,
    Enrichment_Fold=safe_div(observed,exp_overlap),
    Null_95_Lower=as.numeric(quantile(sim,.025)),
    Null_95_Upper=as.numeric(quantile(sim,.975)),
    MonteCarlo_P_UpperTail=(1+sum(sim>=observed))/(n_null+1),
    N_Null_Simulations=n_null)
}
N <- tf_degree_table$Allowed_Targets
d1 <- tf_degree_table$Drought_Degree
d2 <- tf_degree_table$Alkaline_Degree
d3 <- tf_degree_table$Cold_Degree
active <- d1>0 & d2>0 & d3>0 & N>0
exp3 <- sum(as.double(d1[active])*d2[active]*d3[active]/(as.double(N[active])^2))
sim3 <- replicate(n_null, {
  ab <- rhyper(sum(active),m=d2[active],n=N[active]-d2[active],k=d1[active])
  sum(rhyper(sum(active),m=ab,n=N[active]-ab,k=d3[active]))
})
null_rows[[length(null_rows)+1L]] <- data.table(
  Stress_A="All_three",Stress_B="All_three",Overlap_Type="Triple",
  Observed=all3_count,Null_Mean=mean(sim3),Null_SD=sd(sim3),
  Analytic_Expected=exp3,Enrichment_Fold=safe_div(all3_count,exp3),
  Null_95_Lower=as.numeric(quantile(sim3,.025)),
  Null_95_Upper=as.numeric(quantile(sim3,.975)),
  MonteCarlo_P_UpperTail=(1+sum(sim3>=all3_count))/(n_null+1),
  N_Null_Simulations=n_null)
null_summary <- rbindlist(null_rows)
write_dt(null_summary,"04_OutdegreePreserving_Null_Overlap.csv")
cat("[Edges] Edge conservation and outdegree-preserving null complete.\n")

# ---- 7. TF-target rewiring (pairwise, direction retained) --------------------
rewiring <- rbindlist(lapply(pair_names,function(p) {
  a <- p[1L]; b <- p[2L]
  z <- tf_degree_table[,.(TF,Allowed_Targets,
    Degree_A=get(paste0(a,"_Degree")),Degree_B=get(paste0(b,"_Degree")))]
  shared <- shared_edges[[a]][shared_edges[[b]],on=.(TF,Target),nomatch=0L]
  count <- shared[,.(Shared_Targets=.N),by=TF]
  z[,Shared_Targets:=0L]
  z[count,Shared_Targets:=i.Shared_Targets,on="TF"]
  z[, `:=`(Stress_A=a,Stress_B=b,
           Union_Targets=Degree_A+Degree_B-Shared_Targets,
           Unique_A=Degree_A-Shared_Targets,
           Unique_B=Degree_B-Shared_Targets,
           Expected_Shared=as.double(Degree_A)*Degree_B/Allowed_Targets,
           Is_Common_Hub=TF %chin% common_hubs)]
  z[, `:=`(Target_Jaccard=safe_div(Shared_Targets,Union_Targets),
           Rewiring_Fraction=safe_div(Unique_A+Unique_B,Union_Targets),
           Shared_Fold=safe_div(Shared_Targets,Expected_Shared),
           P_Overlap_Enrichment=phyper(Shared_Targets-1L,
              Degree_B,Allowed_Targets-Degree_B,Degree_A,lower.tail=FALSE))]
  z[Union_Targets==0L,P_Overlap_Enrichment:=NA_real_]
  z[, FDR_Overlap_BH:=p.adjust(P_Overlap_Enrichment,method="BH")]
  z
}))
setcolorder(rewiring,c("Stress_A","Stress_B","TF","Is_Common_Hub",
                        setdiff(names(rewiring),c("Stress_A","Stress_B","TF","Is_Common_Hub"))))
write_dt(rewiring,"05_PerTF_TF_Target_Rewiring.csv")
rewiring_summary <- rewiring[Union_Targets>0L, .(
  N_TFs=.N,Median_Jaccard=median(Target_Jaccard,na.rm=TRUE),
  Mean_Jaccard=mean(Target_Jaccard,na.rm=TRUE),
  Mean_Rewiring=mean(Rewiring_Fraction,na.rm=TRUE),
  TFs_Overlap_Enriched_FDR05=sum(FDR_Overlap_BH<.05,na.rm=TRUE)),
  by=.(Stress_A,Stress_B,Is_Common_Hub)]
write_dt(rewiring_summary,"05_TF_Rewiring_Core_vs_Other.csv")

# Per-stress uniqueness relative to the same shared universe. Hypergeometric
# tests compare each TF's UNIQUE edge fraction with network-wide uniqueness.
per_tf <- annotated[, .(Degree=.N,
   Unique_Edges=sum(N_Stresses==1L),
   Two_Stress_Edges=sum(N_Stresses==2L),
   All_Three_Edges=sum(N_Stresses==3L)), by=.(Stress,TF)]
per_tf[, `:=`(Unique_Fraction=Unique_Edges/Degree,
             All3_Fraction=All_Three_Edges/Degree,
             Is_Common_Hub=TF %chin% common_hubs)]
for (s in stresses) {
  ix <- which(per_tf$Stress==s)
  Nedges <- nrow(shared_edges[[s]])
  Kunique <- sum(per_tf$Unique_Edges[ix])
  baseline <- if(Nedges>0) Kunique/Nedges else NA_real_
  per_tf[ix,Network_Unique_Baseline:=baseline]
  per_tf[ix,Unique_Enrichment_Fold:=safe_div(Unique_Fraction,baseline)]
  per_tf[ix,P_Unique_Enrichment:=phyper(Unique_Edges-1L,Kunique,
          Nedges-Kunique,Degree,lower.tail=FALSE)]
  per_tf[ix[per_tf$Degree[ix]<min_tf_degree],P_Unique_Enrichment:=NA_real_]
  per_tf[ix,FDR_Unique_BH:=p.adjust(P_Unique_Enrichment,method="BH")]
}
write_dt(per_tf,"06_PerTF_StressSpecific_Edge_Enrichment.csv")
spec_summary <- per_tf[,.(
   N_TFs=.N,TFs_Eligible=sum(Degree>=min_tf_degree),
   TFs_Unique_Enriched_FDR05=sum(FDR_Unique_BH<.05,na.rm=TRUE),
   Total_Edges=sum(Degree),Total_Unique_Edges=sum(Unique_Edges),
   Edge_Weighted_Unique_Fraction=sum(Unique_Edges)/sum(Degree),
   TF_Mean_Unique_Fraction=mean(Unique_Fraction),
   TF_Mean_All3_Fraction=mean(All3_Fraction)),
   by=.(Stress,Is_Common_Hub)]
write_dt(spec_summary,"06_StressSpecific_Core_vs_Periphery.csv")

# Resample TFs separately within hub / nonhub groups for CI of the difference
# in UNWEIGHTED TF mean unique fraction (other minus shared hub).
core_diff <- rbindlist(lapply(stresses,function(s) {
  z <- per_tf[Stress==s]
  c <- z[Is_Common_Hub==TRUE,Unique_Fraction]
  o <- z[Is_Common_Hub==FALSE,Unique_Fraction]
  if (!length(c) || !length(o))
    return(data.table(Stress=s, Core_N=length(c),Other_N=length(o),
      Diff_Other_Minus_Core=NA_real_,CI_Lower=NA_real_,CI_Upper=NA_real_))
  sim <- replicate(n_tf_boot,mean(sample(o,length(o),replace=TRUE))-
                            mean(sample(c,length(c),replace=TRUE)))
  data.table(Stress=s,Core_N=length(c),Other_N=length(o),
    Diff_Other_Minus_Core=mean(o)-mean(c),
    CI_Lower=as.numeric(quantile(sim,.025)),CI_Upper=as.numeric(quantile(sim,.975)))
}))
write_dt(core_diff,"06_TFClusterBootstrap_Periphery_Minus_Core.csv")
# Complementary direct evidence for a "conserved core": fraction of edges
# present in all three stresses, compared as TF-level means (not raw edge P).
core_preservation <- rbindlist(lapply(stresses,function(s) {
  z <- per_tf[Stress==s]
  hc <- z[Is_Common_Hub==TRUE,All3_Fraction]
  ot <- z[Is_Common_Hub==FALSE,All3_Fraction]
  if (!length(hc) || !length(ot))
    return(data.table(Stress=s,Core_N=length(hc),Other_N=length(ot),
      Diff_Core_Minus_Other=NA_real_,CI_Lower=NA_real_,CI_Upper=NA_real_))
  sim <- replicate(n_tf_boot,mean(sample(hc,length(hc),replace=TRUE))-
                            mean(sample(ot,length(ot),replace=TRUE)))
  data.table(Stress=s,Core_N=length(hc),Other_N=length(ot),
    Diff_Core_Minus_Other=mean(hc)-mean(ot),
    CI_Lower=as.numeric(quantile(sim,.025)),CI_Upper=as.numeric(quantile(sim,.975)))
}))
write_dt(core_preservation,"06_TFClusterBootstrap_Core_Minus_Other_Conserved.csv")
cat("[Rewiring] TF rewiring, specific edge enrichment and TF-bootstrap complete.\n")

# ---- 8. Cross-stress comparison figures (A-D) -------------------------------

# Visual defaults: publication palette; ALL TEXT BLACK; FULL CLOSED FRAMES.
PUB_STRESS <- c(Drought = '#267D79', Alkaline = '#D67A68', Cold = '#5A8FD4')
PUB_CLASS  <- c('Shared hub TFs'='#C66A60','Other TFs'='#6B9BB6')
PUB_INK    <- '#000000'
PUB_FRAME  <- '#000000'

publication_theme <- function(base_size = 10.5) {
  ggplot2::theme_bw(base_size=base_size,base_family='sans') +
    ggplot2::theme(
      text=ggplot2::element_text(color=PUB_INK),
      plot.title=ggplot2::element_text(face='bold',color=PUB_INK,size=11.2,
                                       margin=ggplot2::margin(b=3)),
      plot.subtitle=ggplot2::element_text(color=PUB_INK,size=9.0,
                                          margin=ggplot2::margin(b=6)),
      plot.tag=ggplot2::element_text(face='bold',color=PUB_INK,size=13.5),
      plot.tag.position='topleft',
      axis.title=ggplot2::element_text(color=PUB_INK,size=9.4),
      axis.text=ggplot2::element_text(color=PUB_INK,size=9.2),
      axis.ticks=ggplot2::element_line(color=PUB_FRAME,linewidth=.34),
      panel.border=ggplot2::element_rect(color=PUB_FRAME,fill=NA,linewidth=.60),
      panel.grid.major=ggplot2::element_line(color='#E7E9EB',linewidth=.27),
      panel.grid.minor=ggplot2::element_blank(),
      strip.background=ggplot2::element_rect(fill='#F2F4F5',color=PUB_FRAME,
                                             linewidth=.60),
      strip.text=ggplot2::element_text(face='bold',color=PUB_INK,
                                       margin=ggplot2::margin(t=8,b=8)),
      legend.position='bottom',
      legend.title=ggplot2::element_blank(),
      legend.text=ggplot2::element_text(color=PUB_INK,size=8.5),
      legend.key.height=grid::unit(3.5,'mm'),
      plot.margin=ggplot2::margin(t=9,r=11,b=8,l=11)
    )
}

pub_plot_save <- function(plot, stem, width=7.5, height=5.0) {
  ggplot2::ggsave(file.path(out_dir,paste0(stem,'.pdf')),plot=plot,
                  width=width,height=height,units='in',bg='white',
                  device=grDevices::pdf,useDingbats=FALSE)
  ggplot2::ggsave(file.path(out_dir,paste0(stem,'.png')),plot=plot,
                  width=width,height=height,units='in',dpi=400,bg='white')
}

# Panel A: 3x3 pairwise Jaccard heatmap, diagonal neutral (not unity).
# Pairwise edge Jaccard is calculated above on the shared observed TF/Target IDs.
plot_edges <- data.table::copy(pairwise_edge[Scope=='Shared_Eligible_Universe'])
if (nrow(plot_edges)!=length(pair_names))
  stop('Panel A: expected exactly three pairwise comparisons in shared universe')
pair_id <- function(a,b) ifelse(a < b,paste(a,b,sep='|'),paste(b,a,sep='|'))
plot_edges[, PairID := pair_id(Stress_A,Stress_B)]
if (anyDuplicated(plot_edges$PairID)) stop('Panel A: duplicated pairwise comparison')
heat_map <- data.table::CJ(Stress_X=stresses,Stress_Y=stresses,sorted=FALSE)
heat_map[,PairID:=pair_id(Stress_X,Stress_Y)]
heat_map <- merge(heat_map,plot_edges[,.(PairID,Jaccard,Conserved_Edges)],
                  by='PairID',all.x=TRUE,sort=FALSE)
heat_map[,Diagonal:=Stress_X==Stress_Y]
# Explicit logical-column comparison is required for data.table row filtering.
# heat_map[!Diagonal, ...] may be interpreted as a variable in the calling scope
# by current data.table versions; the comparison below is unambiguous.
if (any(!is.finite(heat_map[Diagonal == FALSE, Jaccard])))
  stop('Panel A: missing or nonfinite Jaccard coefficient')
heat_map[,Stress_X:=factor(Stress_X,levels=stresses)]
heat_map[,Stress_Y:=factor(Stress_Y,levels=rev(stresses))]
heat_map[,Cell_Label:=ifelse(Diagonal,'-',sprintf('%.3f',Jaccard))]
upper_fill <- max(.02,heat_map[Diagonal == FALSE, max(Jaccard, na.rm=TRUE)])
p_edge <- ggplot2::ggplot(heat_map,ggplot2::aes(x=Stress_X,y=Stress_Y)) +
  ggplot2::geom_tile(data=heat_map[Diagonal==TRUE],
                     fill='#F0F2F3',color='white',linewidth=1.3,width=.98,height=.98) +
  ggplot2::geom_tile(data=heat_map[Diagonal==FALSE],
                     ggplot2::aes(fill=Jaccard),color='white',
                     linewidth=1.3,width=.98,height=.98) +
  ggplot2::geom_text(ggplot2::aes(label=Cell_Label),color=PUB_INK,
                     fontface='bold',size=4.15) +
  ggplot2::scale_fill_gradient(
    name='Jaccard similarity',
    low='#E9F2F6',high='#4687A9',
    limits=c(0,upper_fill),
    breaks=pretty(c(0,upper_fill),n=4),
    labels=function(v) sprintf('%.2f',v),
    guide=ggplot2::guide_colorbar(
      title.position='top',
      barheight=grid::unit(31,'mm'),
      barwidth=grid::unit(4.5,'mm'),
      ticks=TRUE,frame.colour='#222222'
    )
  ) +
  ggplot2::coord_fixed(clip='off') +
  ggplot2::labs(tag='A',title='Directed-edge conservation',
                subtitle='Pairwise TF-target Jaccard similarity',
                x=NULL,y=NULL) +
  publication_theme() +
  ggplot2::theme(panel.grid=ggplot2::element_blank(),
                 axis.ticks=ggplot2::element_blank(),
                 axis.text.x=ggplot2::element_text(face='bold',color=PUB_INK),
                 axis.text.y=ggplot2::element_text(face='bold',color=PUB_INK),
                 legend.position='right',
                 legend.title=ggplot2::element_text(face='bold',color=PUB_INK,size=8.5),
                 legend.text=ggplot2::element_text(color=PUB_INK,size=8),
                 legend.margin=ggplot2::margin(l=1,r=0))
data.table::fwrite(heat_map[,.(Stress_X,Stress_Y,Jaccard,Conserved_Edges,Diagonal)],
                   file.path(out_dir,'Fig12A_Jaccard_Heatmap_Data.csv'))

# Panel B: per-TF target rewiring. Distribution rather than three aggregate bars.
# Literature-style violin + narrow white box; NO newly inferred P values.
format_pair <- function(a,b) paste(a,b,sep='\nvs\n')
plot_rw <- data.table::copy(rewiring[Union_Targets>0L & is.finite(Rewiring_Fraction),
                                     .(Stress_A,Stress_B,Is_Common_Hub,Rewiring_Fraction)])
plot_rw[,Comparison:=format_pair(Stress_A,Stress_B)]
pair_levels <- vapply(pair_names,function(p) format_pair(p[[1]],p[[2]]),character(1))
plot_rw[,Comparison:=factor(Comparison,levels=pair_levels)]
plot_rw[,Hub_Class:=ifelse(Is_Common_Hub,'Shared hub TFs','Other TFs')]
plot_rw[,Hub_Class:=factor(Hub_Class,levels=names(PUB_CLASS))]
if (!nrow(plot_rw)) stop('Panel B: no evaluable per-TF rewiring fractions')
box_n <- plot_rw[,.(N_TFs=.N),by=.(Comparison,Hub_Class)]
data.table::fwrite(box_n,file.path(out_dir,'Fig12B_Rewiring_PerGroup_N.csv'))
p_rw <- ggplot2::ggplot(plot_rw,
                        ggplot2::aes(x=Comparison,y=Rewiring_Fraction,
                                    fill=Hub_Class)) +
  ggplot2::geom_violin(ggplot2::aes(group=interaction(Comparison,Hub_Class)),
                       position=ggplot2::position_dodge(width=.78),
                       width=.73,scale='width',trim=TRUE,alpha=.68,
                       color='#454545',linewidth=.28,na.rm=TRUE) +
  ggplot2::geom_boxplot(ggplot2::aes(group=interaction(Comparison,Hub_Class)),
                        position=ggplot2::position_dodge(width=.78),
                        width=.14,fill='white',color='#161616',
                        outlier.shape=NA,linewidth=.46,na.rm=TRUE) +
  ggplot2::scale_fill_manual(values=PUB_CLASS,drop=FALSE) +
  ggplot2::scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),
                              labels=function(x) sprintf('%.0f%%',100*x)) +
  ggplot2::labs(tag='B',title='TF-target rewiring',
                subtitle='Distribution across TFs (1 - target Jaccard)',
                x=NULL,y='Rewired target fraction') +
  publication_theme() +
  ggplot2::theme(panel.grid.major.x=ggplot2::element_blank())

# Panels C & D: effect-size forest plots with the ACTUAL TF-bootstrap 95% CIs.
# Both panels use a COMMON x-axis range (in percentage points).
prepare_diff_plot <- function(dt,estimate_col) {
  wanted <- c('Stress',estimate_col,'CI_Lower','CI_Upper')
  if (!all(wanted %in% names(dt)))
    stop('Missing TF-bootstrap fields: ',paste(setdiff(wanted,names(dt)),collapse=', '))
  z <- data.table::copy(dt)
  z[,`:=`(Estimate_pp=100*get(estimate_col),
          Lower_pp=100*CI_Lower,Upper_pp=100*CI_Upper)]
  z[,Y:=as.numeric(length(stresses)-match(Stress,stresses)+1L)]
  if (anyNA(z$Y) || anyDuplicated(z$Stress) || nrow(z)!=length(stresses))
    stop('Forest panels require one result per stress')
  z[,Valid_CI:=is.finite(Lower_pp)&is.finite(Upper_pp)&
                   Lower_pp<=Estimate_pp&Upper_pp>=Estimate_pp]
  z
}
plot_core <- prepare_diff_plot(core_preservation,'Diff_Core_Minus_Other')
plot_divergent <- prepare_diff_plot(core_diff,'Diff_Other_Minus_Core')
data.table::fwrite(plot_core,file.path(out_dir,'Fig12C_Core_Conservation_Forest_Data.csv'))
data.table::fwrite(plot_divergent,file.path(out_dir,'Fig12D_Periphery_Divergence_Forest_Data.csv'))
range_candidates <- c(0,plot_core$Estimate_pp,plot_core$Lower_pp,plot_core$Upper_pp,
                      plot_divergent$Estimate_pp,plot_divergent$Lower_pp,
                      plot_divergent$Upper_pp)
range_candidates <- range_candidates[is.finite(range_candidates)]
shared_max <- max(5,abs(range_candidates),na.rm=TRUE)
forest_xlim <- c(-1.20*shared_max,1.20*shared_max)

forest_panel <- function(dt,tag,heading,subtitle) {
  with_ci <- dt[Valid_CI==TRUE]
  if (any(!dt$Valid_CI)) warning(tag,': missing/invalid bootstrap CI; points only for affected stresses')
  ggplot2::ggplot(dt,ggplot2::aes(x=Estimate_pp,y=Y)) +
    ggplot2::geom_vline(xintercept=0,linetype='dashed',
                        color='#6C7278',linewidth=.65) +
    ggplot2::geom_segment(data=with_ci,
       ggplot2::aes(x=Lower_pp,xend=Upper_pp,y=Y,yend=Y,color=Stress),
       linewidth=.90,lineend='round',show.legend=FALSE) +
    ggplot2::geom_segment(data=with_ci,
       ggplot2::aes(x=Lower_pp,xend=Lower_pp,y=Y-.085,yend=Y+.085,color=Stress),
       linewidth=.72,show.legend=FALSE) +
    ggplot2::geom_segment(data=with_ci,
       ggplot2::aes(x=Upper_pp,xend=Upper_pp,y=Y-.085,yend=Y+.085,color=Stress),
       linewidth=.72,show.legend=FALSE) +
    ggplot2::geom_point(ggplot2::aes(fill=Stress),shape=21,size=3.6,
                        color='black',stroke=.38,show.legend=FALSE) +
    ggplot2::geom_text(ggplot2::aes(x=Estimate_pp,y=Y+.23,
       label=sprintf('%+.1f pp',Estimate_pp)),
       color=PUB_INK,size=3.0,fontface='bold',na.rm=TRUE) +
    ggplot2::scale_color_manual(values=PUB_STRESS,drop=FALSE) +
    ggplot2::scale_fill_manual(values=PUB_STRESS,drop=FALSE) +
    ggplot2::scale_x_continuous(limits=forest_xlim,
        breaks=pretty(forest_xlim,n=5),
        labels=function(x) sprintf('%.0f',x)) +
    ggplot2::scale_y_continuous(breaks=seq_along(stresses),labels=rev(stresses),
        limits=c(.58,length(stresses)+.51),expand=ggplot2::expansion(mult=0)) +
    ggplot2::labs(tag=tag,title=heading,subtitle=subtitle,
                  x='Mean TF-level difference (percentage points)',y=NULL) +
    publication_theme() +
    ggplot2::theme(panel.grid.major.y=ggplot2::element_blank())
}
p_core <- forest_panel(plot_core,'C','Conserved-edge fraction',
                       'Shared hubs minus other TFs, with 95% TF-bootstrap CI')
p_divergent <- forest_panel(plot_divergent,'D','Stress-specific edge fraction',
                            'Other TFs minus shared hubs, with 95% TF-bootstrap CI')

# Export single plots in both vector PDF and high-resolution PNG.
pub_plot_save(p_edge,'Fig12A_Edge_Conservation',7.4,5.15)
pub_plot_save(p_rw,'Fig12B_PerTF_Rewiring',7.4,5.15)
pub_plot_save(p_core,'Fig12C_Conserved_Core_95CI',7.4,5.15)
pub_plot_save(p_divergent,'Fig12D_StressSpecific_Divergence_95CI',7.4,5.15)

# Supplement: descriptive fraction of unique, pairwise and triple-shared edges
# within the union of all predicted edges in the shared observed ID universe.
occ_supp <- pattern_summary[,.(N_Edges=sum(N_Edges)),by=N_Stresses]
occ_supp[,Fraction:=N_Edges/sum(N_Edges)]
occ_supp[,Group:=factor(N_Stresses,levels=c(1,2,3),
                        labels=c('One stress','Two stresses','All three'))]
occ_supp[,Pct_Label:=sprintf('%.1f%%',100*Fraction)]
p_occ <- ggplot2::ggplot(occ_supp,
                         ggplot2::aes(x=Group,y=Fraction,fill=Group)) +
  ggplot2::geom_col(width=.42,color='black',linewidth=.30) +
  ggplot2::geom_text(ggplot2::aes(label=Pct_Label),vjust=-.6,
                     size=3.75,color='black',fontface='bold') +
  ggplot2::scale_fill_manual(values=c('One stress'='#D7A190',
                                     'Two stresses'='#89A7BF',
                                     'All three'='#5E9D92'),guide='none') +
  ggplot2::scale_y_continuous(labels=function(x) sprintf('%.0f%%',100*x),
      expand=ggplot2::expansion(mult=c(0,.13))) +
  ggplot2::labs(title='Edge occupancy across stresses',
                subtitle='Unique directed TF-target pairs in shared observed ID scope',
                x=NULL,y='Fraction of union edges') +
  publication_theme() +
  ggplot2::theme(panel.grid.major.x=ggplot2::element_blank())
pub_plot_save(p_occ,'FigS_Edge_Occupancy_OneTwoThree_Stresses',7.4,4.9)
data.table::fwrite(occ_supp,file.path(out_dir,'FigS_Edge_Occupancy_Plot_Data.csv'))

# 2x2 layout in A4 LANDSCAPE (297 x 210 mm), with single-page vector PDF.
draw_four_panels <- function(plots) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout=grid::grid.layout(2L,2L,
      widths=grid::unit(c(1,1),'null'),
      heights=grid::unit(c(1,1),'null'))))
  on.exit(grid::popViewport(),add=TRUE)
  for (i in seq_along(plots)) {
    rr <- (i-1L)%/%2L+1L
    cc <- (i-1L)%%2L+1L
    print(plots[[i]],vp=grid::viewport(layout.pos.row=rr,layout.pos.col=cc),
          newpage=FALSE)
  }
}
combined_pdf <- file.path(out_dir,'Cross_Stress_Network_Comparison_ABCD.pdf')
if (capabilities('cairo')) {
  grDevices::cairo_pdf(combined_pdf,width=297/25.4,height=210/25.4,
                       bg='white',onefile=TRUE)
} else {
  grDevices::pdf(combined_pdf,width=297/25.4,height=210/25.4,
                 bg='white',useDingbats=FALSE)
}
tryCatch(draw_four_panels(list(p_edge,p_rw,p_core,p_divergent)),
         finally=grDevices::dev.off())
combined_png <- file.path(out_dir,'Cross_Stress_Network_Comparison_ABCD.png')
grDevices::png(combined_png,width=297/25.4,height=210/25.4,units='in',
               res=400,bg='white')
tryCatch(draw_four_panels(list(p_edge,p_rw,p_core,p_divergent)),
         finally=grDevices::dev.off())
