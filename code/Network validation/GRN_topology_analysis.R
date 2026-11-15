# gsGNA network topology and Erdős–Rényi G(n,m) comparison
#
# Run from the repository root: Rscript code/GRN_topology_analysis.R
# Set GSGNA_INPUT_ROOT to the root used for the ensemble comparison.
# Uses the top 10% of the ranked XGBoost selected-union edges per stress.
# Computes directed degree distributions and compares the largest undirected
# connected component with ER graphs matched for node and edge counts.

required <- c('igraph', 'data.table', 'ggplot2')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install the required R packages: ', paste(missing, collapse=', '))
suppressPackageStartupMessages({
  library(igraph)
  library(data.table)
  library(ggplot2)
})

# ------------------------- USER CONFIGURATION --------------------------------
input_root <- Sys.getenv('GSGNA_INPUT_ROOT', unset='.')
network_base <- file.path(input_root, 'gsGNA_FixedTop10', 'FINAL_Locked_Test_And_Network')
# The order determines the random seed used for each stress.
stress_names <- c('Drought', 'Alkaline', 'Cold')
network_files <- setNames(
  file.path(network_base, stress_names, 'Final_Ranked_Network',
            'XGBoost_Ranked_SelectedUnion_AllEdges.csv'),
  stress_names
)
output_dir <- file.path(input_root, 'gsGNA_Topology_ER_Gnm')
network_edge_fraction <- 0.10
n_rand <- 199L                 # Number of ER replicates per stress.
seed <- 123L
exact_path_max_nodes <- 1500L # All-pairs shortest paths only on small graphs
approx_path_sources <- 64L     # Sampled BFS sources on larger connected graphs
path_batch_sources <- 16L
min_connected_null_fraction <- 0.90 # Below this, sigma is not estimated
min_connected_null_count <- 20L
degree_plot_type <- 'All_Total'  # Alternatives: 'TF_Out' or 'Target_In'.
show_descriptive_loglog_fit <- TRUE

stopifnot(is.finite(network_edge_fraction), network_edge_fraction > 0,
          network_edge_fraction <= 1,
          degree_plot_type %in% c('All_Total', 'TF_Out', 'Target_In'),
          n_rand >= 1L, approx_path_sources >= 2L,
          exact_path_max_nodes >= 2L, path_batch_sources >= 1L,
          min_connected_null_fraction > 0, min_connected_null_fraction <= 1,
          min_connected_null_count >= 1L)

missing_inputs <- network_files[!file.exists(network_files)]
if (length(missing_inputs)) {
  stop('Missing ranked XGBoost networks:\n', paste(missing_inputs, collapse='\n'))
}
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# ------------------------- Read and clean ------------------------------------
read_network <- function(path, stress) {
  raw <- fread(path, select=c('TF', 'Target', 'EdgeWeight'), showProgress=FALSE)
  if (!all(c('TF', 'Target', 'EdgeWeight') %in% names(raw)) || !nrow(raw)) {
    stop('Expected a nonempty ranked network with TF, Target, and EdgeWeight: ', path)
  }
  raw[, EdgeWeight := suppressWarnings(as.numeric(EdgeWeight))]
  if (any(!is.finite(raw$EdgeWeight))) stop('Non-finite ranked edge weights: ', path)
  export_summary <- file.path(dirname(path), 'Final_Network_Export_Summary.csv')
  if (file.exists(export_summary)) {
    audit <- fread(export_summary, showProgress=FALSE)
    required <- c('File', 'Exported_Edges', 'N_Union_Candidates')
    if (!all(required %in% names(audit)) || nrow(audit) != 1L ||
        !isTRUE(audit$File == basename(path)) || !isTRUE(audit$Exported_Edges == nrow(raw)) ||
        !isTRUE(audit$N_Union_Candidates == nrow(raw))) {
      stop('Expected the full ranked selected-union network: ', path)
    }
  }

  setorderv(raw, c('EdgeWeight', 'TF', 'Target'), c(-1L, 1L, 1L))
  n_ranked <- nrow(raw)
  n_keep <- ceiling(n_ranked * network_edge_fraction)
  retained <- raw[seq_len(n_keep)]
  dt <- retained[, .(
    TF=toupper(trimws(as.character(TF))),
    Target=toupper(trimws(as.character(Target))),
    Score=EdgeWeight
  )]
  dt <- dt[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target)]
  dropped_invalid <- n_keep - nrow(dt)
  dropped_loops <- sum(dt$TF == dt$Target)
  dt <- dt[TF != Target]
  before_dedup <- nrow(dt)
  setorder(dt, TF, Target, -Score)
  dt <- unique(dt, by=c('TF', 'Target'))
  dropped_duplicates <- before_dedup - nrow(dt)
  if (nrow(dt) < 2L) stop(stress, ': At least two valid edges are required.')
  message(sprintf('[%s] Retained %s of %s ranked edges before cleaning',
                  stress, n_keep, n_ranked))
  if (any(dt$Score < 0 | dt$Score > 1)) {
    warning(stress, ': Some edge scores are outside [0, 1]. Check the XGBoost input.')
  }
  # Only TF and Target are used for topology; confidence is NOT a distance/sign.
  gd <- igraph::graph_from_data_frame(dt[,.(TF,Target)], directed=TRUE)
  stopifnot(igraph::is_simple(gd), igraph::ecount(gd) == nrow(dt))
  gu <- igraph::as_undirected(gd, mode='collapse')
  gu <- igraph::simplify(gu, remove.multiple=TRUE, remove.loops=TRUE)
  cc <- igraph::components(gu)
  ids <- which(cc$membership == which.max(cc$csize))
  gc <- igraph::induced_subgraph(gu, vids=ids)
  if (igraph::vcount(gc) < 3L || igraph::ecount(gc) < 2L) {
    stop(stress, ': The largest connected component is too small for small-world analysis.')
  }
  audit <- data.table(
    Stress=stress,Input=path,Score_Column='EdgeWeight',
    Ranked_Export_Rows=n_ranked,Requested_Edge_Fraction=network_edge_fraction,
    Original_Rows=n_keep,Removed_Invalid=dropped_invalid,
    Removed_Self_Loops=dropped_loops,Removed_Duplicates=dropped_duplicates,
    Clean_Directed_Edges=nrow(dt),Full_Nodes=igraph::vcount(gd),
    Full_Undirected_Edges=igraph::ecount(gu),
    Full_Directed_Density=igraph::edge_density(gd,loops=FALSE),
    Full_Undirected_Density=igraph::edge_density(gu,loops=FALSE),
    Full_Component_Count=cc$no,
    Largest_Component_Nodes=igraph::vcount(gc),
    Largest_Component_Edges=igraph::ecount(gc),
    Largest_Component_Node_Fraction=igraph::vcount(gc)/igraph::vcount(gu),
    Largest_Component_Edge_Fraction=igraph::ecount(gc)/igraph::ecount(gu)
  )
  if (audit$Largest_Component_Node_Fraction < 0.90) {
    warning(stress, ': The largest connected component contains ',
            sprintf('%.1f%%',100*audit$Largest_Component_Node_Fraction),
            ' of the nodes. Sigma describes this component.')
  }
  list(directed=gd, undirected=gu, core=gc, audit=audit, edges=dt)
}

# ------------------------- Degree descriptions -------------------------------
get_degree_tables <- function(g, dt, stress) {
  vn <- igraph::V(g)$name
  tf_ids <- match(unique(dt$TF),vn)
  target_ids <- match(unique(dt$Target),vn)
  k_out <- igraph::degree(g,v=tf_ids,mode='out',loops=FALSE)
  k_in <- igraph::degree(g,v=target_ids,mode='in',loops=FALSE)
  k_total <- igraph::degree(g,mode='all',loops=FALSE)
  res <- rbindlist(list(
    data.table(Stress=stress,DegreeType='TF_Out',k=as.integer(k_out)),
    data.table(Stress=stress,DegreeType='Target_In',k=as.integer(k_in)),
    data.table(Stress=stress,DegreeType='All_Total',k=as.integer(k_total))
  ))
  if (anyNA(res$k) || any(res$k < 1L)) stop(stress, ': Invalid degree values.')
  res
}

# ------------------------- Path computation ---------------------------------
# Input must be a CONNECTED simple undirected graph.
# Uniform sources selected once and reused for real + every ER graph.
path_stats <- function(g, source_ids=NULL) {
  n <- igraph::vcount(g)
  if (is.null(source_ids)) {
    L <- igraph::mean_distance(g,directed=FALSE,unconnected=FALSE,weights=NA)
    return(list(L=as.numeric(L), Source_SE=NA_real_,Sources=n,Method='Exact'))
  }
  row_avg <- numeric(length(source_ids))
  for (lo in seq.int(1L,length(source_ids),by=path_batch_sources)) {
    hi <- min(length(source_ids),lo+path_batch_sources-1L)
    ids <- source_ids[lo:hi]
    dist <- igraph::distances(g,v=ids,to=igraph::V(g),mode='all',weights=NA)
    dist[cbind(seq_along(ids),ids)] <- NA_real_ # Exclude self-distance=0
    row_avg[lo:hi] <- rowMeans(dist,na.rm=TRUE)
  }
  list(L=mean(row_avg),
       Source_SE=if(length(row_avg)>1L) stats::sd(row_avg)/sqrt(length(row_avg)) else NA_real_,
       Sources=length(source_ids),Method='Sampled_BFS')
}

clustering_global <- function(g) {
  x <- suppressWarnings(igraph::transitivity(g,type='global'))
  if(is.finite(x)) as.numeric(x) else NA_real_
}

one_er_test <- function(core, stress, index) {
  n <- igraph::vcount(core)
  m <- igraph::ecount(core)
  if (as.double(m) > as.double(n)*(n-1)/2) stop('The edge count exceeds the ER G(n,m) limit.')
  set.seed(seed + 1000L*index)
  src <- if(n <= exact_path_max_nodes) NULL else {
    sample.int(n,size=min(approx_path_sources,n),replace=FALSE)
  }
  obs_path <- path_stats(core,source_ids=src)
  obs_C <- clustering_global(core)
  obs <- data.table(
    Stress=stress,ER_Model='Undirected_Erdos_Renyi_Gnm',
    Scope='Fixed_Observed_Largest_Connected_Component',
    Nodes=n,Edges=m,Density=igraph::edge_density(core,loops=FALSE),
    Global_Clustering_C=obs_C,Mean_Shortest_Path_L=obs_path$L,
    Path_Source_SE=obs_path$Source_SE,Path_Sources=obs_path$Sources,
    Path_Method=obs_path$Method
  )

  nulls <- vector('list',n_rand)
  message('[',stress,'] ER G(n,m): n=',n,', m=',m,', repeats=',n_rand)
  for (i in seq_len(n_rand)) {
    set.seed(seed + 100000L*index + i)
    # IMPORTANT: ER network matches the ACTUAL UNDIRECTED CORE's n and m.
    # NOT degree-preserving, NOT an ER G(n,p), NOT a directed graph.
    er <- igraph::sample_gnm(n=n,m=m,directed=FALSE,loops=FALSE)
    connected <- igraph::is_connected(er)
    er_C <- clustering_global(er)
    if (connected) {
      er_path <- path_stats(er,source_ids=src)
    } else {
      # Conditional finite-pair mean of disconnected networks is not used.
      er_path <- list(L=NA_real_,Source_SE=NA_real_)
    }
    nulls[[i]] <- data.table(
      Stress=stress,Replicate=i,Nodes=n,Edges=m,Connected=connected,
      Global_Clustering_C=er_C,Mean_Shortest_Path_L=er_path$L,
      Path_Source_SE=er_path$Source_SE,Density=igraph::edge_density(er,loops=FALSE)
    )
    if (i %% 20L == 0L || i == n_rand) message('[',stress,'] random ',i,'/',n_rand)
  }
  nd <- rbindlist(nulls)
  valid <- nd[Connected==TRUE & is.finite(Global_Clustering_C) &
                is.finite(Mean_Shortest_Path_L) & Mean_Shortest_Path_L>0]
  valid_count <- nrow(valid)
  valid_fraction <- valid_count/n_rand
  enough <- valid_count >= min(min_connected_null_count,n_rand) &&
            valid_fraction >= min_connected_null_fraction

  mu_sd_q <- function(x) {
    x <- x[is.finite(x)]
    c(Mean=if(length(x)) mean(x) else NA_real_,
      SD=if(length(x)>1L) stats::sd(x) else NA_real_,
      Q025=if(length(x)) as.numeric(stats::quantile(x,0.025)) else NA_real_,
      Q975=if(length(x)) as.numeric(stats::quantile(x,0.975)) else NA_real_)
  }
  # BOTH C and L use exactly the same connected-null replicates.
  cstats <- mu_sd_q(valid$Global_Clustering_C)
  lstats <- mu_sd_q(valid$Mean_Shortest_Path_L)
  c_all <- mu_sd_q(nd$Global_Clustering_C)
  ratio_ok <- enough && is.finite(obs_C) && is.finite(obs_path$L) &&
              is.finite(cstats['Mean']) && is.finite(lstats['Mean']) &&
              cstats['Mean']>0 && lstats['Mean']>0 && obs_path$L>0
  gamma <- if(ratio_ok) unname(obs_C/cstats['Mean']) else NA_real_
  lambda <- if(ratio_ok) unname(obs_path$L/lstats['Mean']) else NA_real_
  sigma <- if(ratio_ok) gamma/lambda else NA_real_
  p_C_high <- if (valid_count && is.finite(obs_C)) {
    (1 + sum(valid$Global_Clustering_C >= obs_C))/(valid_count+1)
  } else NA_real_
  p_L_2s <- if (valid_count && is.finite(obs_path$L)) {
    hi <- (1 + sum(valid$Mean_Shortest_Path_L >= obs_path$L))/(valid_count+1)
    lo <- (1 + sum(valid$Mean_Shortest_Path_L <= obs_path$L))/(valid_count+1)
    min(1,2*min(hi,lo))
  } else NA_real_

  sw <- data.table(
    Stress=stress,Null_Model='ER_Gnm_undirected_fixed_core_N_M',
    Nodes=n,Edges=m,
    N_ER_Requested=n_rand,N_ER_Connected_And_Valid=valid_count,
    Fraction_ER_Connected_And_Valid=valid_fraction,
    C_Observed=obs_C,C_ER_Connected_Mean=unname(cstats['Mean']),
    C_ER_Connected_SD=unname(cstats['SD']),
    C_ER_Connected_Q025=unname(cstats['Q025']),
    C_ER_Connected_Q975=unname(cstats['Q975']),
    C_ER_All_Mean=unname(c_all['Mean']),
    L_Observed=obs_path$L,L_ER_Connected_Mean=unname(lstats['Mean']),
    L_ER_Connected_SD=unname(lstats['SD']),
    L_ER_Connected_Q025=unname(lstats['Q025']),
    L_ER_Connected_Q975=unname(lstats['Q975']),
    Gamma=gamma,Lambda=lambda,Sigma=sigma,
    SmallWorld_Heuristic=if(is.finite(sigma)) ifelse(sigma>1,'Sigma_gt_1','Sigma_le_1') else 'Not_estimated',
    Null_Clustering_OneSided_P=p_C_high,
    Null_Path_TwoSided_P=p_L_2s,
    Path_Method=obs_path$Method,Path_Sources=obs_path$Sources,
    Status=if(ratio_ok) 'Descriptive_ER_sigma_available'
           else if(!enough) 'Not_reported_insufficient_connected_ER_nulls'
           else 'Not_reported_invalid_C_or_L_denominator'
  )
  message(sprintf('[%s] Connected valid ER = %d/%d (%.1f%%); gamma=%s; lambda=%s; sigma=%s',
          stress,valid_count,n_rand,100*valid_fraction,
          format(signif(gamma,4)),format(signif(lambda,4)),format(signif(sigma,4))))
  list(observed=obs, null=nd, smallworld=sw)
}

# ------------------------- Main loop ----------------------------------------
all_audit <- list(); all_degree <- list(); all_obs <- list()
all_null <- list(); all_sw <- list()
for (ix in seq_along(network_files)) {
  stress <- names(network_files)[ix]
  message('\n==================== ',stress,' ====================')
  path <- network_files[[stress]]
  obj <- read_network(path,stress)
  all_audit[[stress]] <- obj$audit
  all_degree[[stress]] <- get_degree_tables(obj$directed,obj$edges,stress)
  result <- one_er_test(obj$core,stress,ix)
  all_obs[[stress]] <- result$observed
  all_null[[stress]] <- result$null
  all_sw[[stress]] <- result$smallworld
  rm(obj,result); gc(verbose=FALSE)
}
combine_save <- function(xs,fn) {
  d <- rbindlist(xs,use.names=TRUE,fill=TRUE)
  fwrite(d,file.path(output_dir,fn))
  d
}
audit <- combine_save(all_audit,'00_Input_Network_Audit.csv')
degree_all <- combine_save(all_degree,'01_Degree_Vectors.csv')
degree_summary <- degree_all[, .(
  Count  = as.integer(.N),
  Min    = as.numeric(min(k)),
  Mean   = as.numeric(mean(k)),
  Median = as.numeric(stats::median(k)),
  P90    = as.numeric(stats::quantile(k, 0.90)),
  P95    = as.numeric(stats::quantile(k, 0.95)),
  P99    = as.numeric(stats::quantile(k, 0.99)),
  Max    = as.numeric(max(k))
), by = .(Stress, DegreeType)]

fwrite(
  degree_summary,
  file.path(output_dir, "02_Degree_Summary.csv")
)
observed <- combine_save(all_obs,'05_Observed_Core_Topology.csv')
null_all <- combine_save(all_null,'06_ER_Gnm_Null_Replicates.csv')
smallworld <- combine_save(all_sw,'08_SmallWorld_ER_Gnm_Gamma_Lambda_Sigma.csv')

# ------------------------- Publication text style -----------------------------
# Text stays pure black and all original data colours remain unchanged.
# Increase the height of the FACET TITLE BOXES (e.g. "Average Degree") by
# adding vertical padding to horizontal facet strip labels.
# Adjustable: 10 pt per side (top and bottom); use 12 or 14 for taller boxes.
facet_title_padding_pt <- 10
# Applied AFTER each plot's own theme to override ggplot2's default grey text.
black_text_theme <- ggplot2::theme(
  text               = ggplot2::element_text(colour="#000000"),
  plot.title         = ggplot2::element_text(colour="#000000"),
  plot.subtitle      = ggplot2::element_text(colour="#000000"),
  plot.caption       = ggplot2::element_text(colour="#000000"),
  plot.tag           = ggplot2::element_text(colour="#000000"),
  axis.title         = ggplot2::element_text(colour="#000000"),
  axis.title.x       = ggplot2::element_text(colour="#000000"),
  axis.title.y       = ggplot2::element_text(colour="#000000"),
  axis.text          = ggplot2::element_text(colour="#000000"),
  axis.text.x        = ggplot2::element_text(colour="#000000"),
  axis.text.y        = ggplot2::element_text(colour="#000000"),
  strip.text         = ggplot2::element_text(colour="#000000"),
  strip.text.x       = ggplot2::element_text(
    colour="#000000",
    margin=ggplot2::margin(t=facet_title_padding_pt,
                           b=facet_title_padding_pt, unit="pt")
  ),
  strip.text.y       = ggplot2::element_text(colour="#000000"),
  legend.title       = ggplot2::element_text(colour="#000000"),
  legend.text        = ggplot2::element_text(colour="#000000")
)

# ------------------------- Figure A: degree CCDF ------------------------------
degree_freq <- degree_all[,.(Count=.N),by=.(Stress,DegreeType,k)]
setorder(degree_freq,Stress,DegreeType,-k)
degree_freq[,Prob_GE:=cumsum(Count)/sum(Count),by=.(Stress,DegreeType)]
fwrite(degree_freq,file.path(output_dir,'03_Degree_CCDF_Data.csv'))
p_a <- ggplot(degree_freq,aes(x=k,y=Prob_GE,color=Stress)) +
  geom_step(linewidth=0.8) +
  facet_grid(DegreeType~Stress,scales='free') +
  scale_x_log10() + scale_y_log10() +
  labs(title='Empirical degree CCDF of final rice GRNs',
       subtitle='Descriptive, no power-law fitting or scale-free classification',
       x='Degree k (log10)',y='P(K >= k) (log10)') +
  theme_bw(base_size=11) + theme(legend.position='none')
p_a <- p_a + black_text_theme
ggsave(file.path(output_dir,'Figure_A_Degree_CCDF.pdf'),p_a,width=12,height=8,bg='white')
ggsave(file.path(output_dir,'Figure_A_Degree_CCDF.png'),p_a,width=12,height=8,dpi=300,bg='white')

# ------------------------- Figure B: global clustering ------------------------
p_b <- ggplot(null_all[is.finite(Global_Clustering_C)],
              aes(x=Global_Clustering_C,fill=Stress)) +
  geom_histogram(bins=25,alpha=0.75,color='white') +
  geom_vline(data=observed,aes(xintercept=Global_Clustering_C),
             color='#B43333',linewidth=1) +
  facet_wrap(~Stress,scales='free') +
  labs(title='Observed clustering vs ER G(n,m) null networks',
       subtitle='Red line: observed largest connected component; nulls: fixed n and m',
       x='Global clustering coefficient',y='Number of ER graphs') +
  theme_bw(base_size=11) + theme(legend.position='none')
p_b <- p_b + black_text_theme
ggsave(file.path(output_dir,'Figure_B_ER_Clustering_Null.pdf'),p_b,width=11,height=5,bg='white')
ggsave(file.path(output_dir,'Figure_B_ER_Clustering_Null.png'),p_b,width=11,height=5,dpi=300,bg='white')

# ------------------------- Figure C: C and L observed/null --------------------
obs_plot <- rbindlist(list(
  observed[,.(Stress,Metric='Global clustering C',Type='Observed',
               Value=Global_Clustering_C,SD=0)],
  observed[,.(Stress,Metric='Mean shortest path L',Type='Observed',
               Value=Mean_Shortest_Path_L,SD=0)]
))
random_plot <- rbindlist(list(
  smallworld[,.(Stress,Metric='Global clustering C',Type='ER random (connected)',
                 Value=C_ER_Connected_Mean,SD=C_ER_Connected_SD)],
  smallworld[,.(Stress,Metric='Mean shortest path L',Type='ER random (connected)',
                 Value=L_ER_Connected_Mean,SD=L_ER_Connected_SD)]
))
plot_c <- rbindlist(list(obs_plot,random_plot))
plot_c <- plot_c[is.finite(Value)]
p_c <- ggplot(plot_c,aes(x=Stress,y=Value,fill=Type)) +
  geom_col(position=position_dodge(width=0.75),width=0.66) +
  geom_errorbar(data=plot_c[Type=='ER random (connected)' & is.finite(SD)],
       aes(ymin=pmax(0,Value-SD),ymax=Value+SD),
       position=position_dodge(width=0.75),width=0.2) +
  facet_wrap(~Metric,scales='free_y',ncol=2) +
  scale_fill_manual(values=c('Observed'='#24789D','ER random (connected)'='#B05B80')) +
  labs(title='Topology of observed GRNs vs ER random networks',
       subtitle='Random mean +/- SD; C and L use the same connected ER replicates',
       x='Stress',y='Metric value',fill='Network') +
  theme_bw(base_size=11) + theme(legend.position='bottom')
p_c <- p_c + black_text_theme
ggsave(file.path(output_dir,'Figure_C_Observed_vs_ER.pdf'),p_c,width=10,height=5.5,bg='white')
ggsave(file.path(output_dir,'Figure_C_Observed_vs_ER.png'),p_c,width=10,height=5.5,dpi=300,bg='white')

# ------------------------- Figure D: ratios ----------------------------------
ratio_data <- data.table::melt(
  smallworld[,.(Stress,Gamma,Lambda,Sigma)],id.vars='Stress',
  measure.vars=c('Gamma','Lambda','Sigma'),
  variable.name='Ratio',value.name='Value')
ratio_data <- ratio_data[is.finite(Value)]
if (nrow(ratio_data)) {
  p_d <- ggplot(ratio_data,aes(x=Stress,y=Value,fill=Stress)) +
    geom_col(width=0.60) +
    geom_hline(yintercept=1,linetype='dashed',color='#A34141') +
    geom_text(aes(label=sprintf('%.2f',Value)),vjust=-0.4,size=3.5,colour='#000000') +
    facet_wrap(~Ratio,scales='free_y',ncol=3) +
    labs(title='Small-world descriptors against ER G(n,m)',
         subtitle='Gamma=Cobs/C_ER; Lambda=Lobs/L_ER; Sigma=Gamma/Lambda',
         x='Stress',y='Ratio (reference = 1)') +
    theme_bw(base_size=11) + theme(legend.position='none')
  p_d <- p_d + black_text_theme
  ggsave(file.path(output_dir,'Figure_D_ER_Gamma_Lambda_Sigma.pdf'),p_d,
         width=11,height=5,bg='white')
  ggsave(file.path(output_dir,'Figure_D_ER_Gamma_Lambda_Sigma.png'),p_d,
         width=11,height=5,dpi=300,bg='white')
} else {
  warning('Sigma is unavailable for all stresses. Check ER connectivity and clustering.')
}

# ------------------------- Figure C: FOUR-panel topology bars ----------------
# Each row contains the estimate and SD for one stress, network and metric.
# Mean degree and density are fixed when ER G(n,m) uses the same n and m;
# accordingly both bars overlap in magnitude and their SDs are effectively 0.
stress_levels <- c('Alkaline','Drought','Cold')
metric_levels <- c('Average Degree','Clustering Coefficient',
                   'Average Path Length','Network Density')

obs_bars <- data.table::rbindlist(list(
  observed[, .(Stress,Metric='Average Degree', Network='Original',
               Value=as.numeric(2*Edges/Nodes),SD=0.0)],
  observed[, .(Stress,Metric='Clustering Coefficient',Network='Original',
               Value=as.numeric(Global_Clustering_C),SD=0.0)],
  observed[, .(Stress,Metric='Average Path Length',Network='Original',
               Value=as.numeric(Mean_Shortest_Path_L),SD=0.0)],
  observed[, .(Stress,Metric='Network Density',Network='Original',
               Value=as.numeric(Density),SD=0.0)]
))

# C and L: use the same connected ER replicas as sigma.
# Mean degree and density: all generated ER graphs, because fixed n,m by design.
null_nm <- null_all[, .(
  DegreeMean = as.numeric(mean(2*Edges/Nodes)),
  DegreeSD = as.numeric(stats::sd(2*Edges/Nodes)),
  DensityMean = as.numeric(mean(Density)),
  DensitySD = as.numeric(stats::sd(Density))
), by=Stress]

rand_bars <- data.table::rbindlist(list(
  null_nm[, .(Stress,Metric='Average Degree',Network='ER random',
              Value=DegreeMean,SD=DegreeSD)],
  smallworld[, .(Stress,Metric='Clustering Coefficient',Network='ER random',
                 Value=as.numeric(C_ER_Connected_Mean),
                 SD=as.numeric(C_ER_Connected_SD))],
  smallworld[, .(Stress,Metric='Average Path Length',Network='ER random',
                 Value=as.numeric(L_ER_Connected_Mean),
                 SD=as.numeric(L_ER_Connected_SD))],
  null_nm[, .(Stress,Metric='Network Density',Network='ER random',
              Value=DensityMean,SD=DensitySD)]
))

bars <- rbindlist(list(obs_bars,rand_bars),use.names=TRUE)
bars[, Metric := factor(Metric, levels=metric_levels)]
bars[, Network := factor(Network, levels=c('Original','ER random'))]
bars[, Stress := factor(Stress, levels=stress_levels)]
# Fixed x offsets ensure error bars line up with the correct group, even when
# they are drawn only for the random-network group.
bars[, x_pos := as.numeric(Stress) + ifelse(Network=='Original', -0.19, 0.19)]
fwrite(bars[,.(Stress,Network,Metric,Value,SD)],
       file.path(output_dir,'09_Four_Metric_Plot_Data.csv'))

p_c <- ggplot(bars[is.finite(Value)],aes(x=x_pos,y=Value,fill=Network)) +
  geom_col(width=0.36,color='gray35',linewidth=0.18) +
  geom_errorbar(
    data=bars[Network=='ER random' & is.finite(Value) & is.finite(SD)],
    aes(ymin=pmax(0,Value-SD),ymax=Value+SD),
    width=0.13,linewidth=0.5
  ) +
  facet_wrap(~Metric,scales='free_y',ncol=2) +
  scale_x_continuous(breaks=seq_along(stress_levels),labels=stress_levels) +
  scale_fill_manual(values=c('Original'='#2E86AB','ER random'='#A23B72')) +
  labs(title='Topological Metrics of GRNs vs ER Random Networks',
       subtitle='ER mean +/- SD; C/L use connected ER graphs; degree/density match by construction',
       x='Stress Condition',y='Metric Value',fill='Network Type') +
  theme_bw(base_size=12) +
  theme(text=element_text(color='black'),
        plot.title=element_text(hjust=0.5,size=13,face='bold',color='black'),
        strip.text=element_text(face='bold',color='black'),
        axis.text=element_text(color='black'),
        axis.title=element_text(color='black'),
        legend.position='bottom', panel.grid.minor=element_blank())

p_c <- p_c + black_text_theme
ggsave(file.path(output_dir,'Figure_C_Four_Metric_Topology_Bar.pdf'),p_c,
       width=297,height=210,units='mm',bg='white',device=grDevices::pdf,
       useDingbats=FALSE)
ggsave(file.path(output_dir,'Figure_C_Four_Metric_Topology_Bar.png'),p_c,
       width=297,height=210,units='mm',dpi=400,bg='white')

# ------------------------- Figure E: OVERLAID log-log degree plot -----------
# The dashed least-squares line is a VISUAL guide only. Linear regression on
# a log-log degree histogram is NOT a valid power-law/scale-free inference.
degree_pdf <- degree_all[DegreeType == degree_plot_type,
                         .(Count=.N),by=.(Stress,k)]
degree_pdf <- degree_pdf[k>0 & Count>0]
degree_pdf[, Freq := Count/sum(Count),by=Stress]
setorder(degree_pdf,Stress,k)

# Return a plot line and descriptive slope/R2, but never label as alpha/gamma
# of a fitted statistical power-law distribution.
fit_lines <- list()
fit_summary <- list()
for (stress in stress_levels) {
  dd <- degree_pdf[Stress==stress & k>0 & Freq>0]
  if (nrow(dd)<3L || length(unique(dd$k))<3L) next
  dd[,`:=`(logk=log10(k),logfreq=log10(Freq))]
  mod <- stats::lm(logfreq~logk,data=dd)
  grid_logk <- seq(min(dd$logk),max(dd$logk),length.out=150L)
  fitted_logfreq <- as.numeric(stats::predict(mod,
                                              newdata=data.frame(logk=grid_logk)))
  fit_lines[[stress]] <- data.table(Stress=stress,k=10^grid_logk,
                                    Freq=10^fitted_logfreq)
  fit_summary[[stress]] <- data.table(
    Stress=stress,DegreeType=degree_plot_type,
    LogLog_OLS_Slope=as.numeric(stats::coef(mod)[2L]),
    LogLog_OLS_R2=as.numeric(summary(mod)$r.squared),
    Unique_Degree_Bins=nrow(dd),
    Interpretation='Visual OLS guide only; NOT a power-law fit or scale-free test'
  )
}
fit_df <- if(length(fit_lines)) rbindlist(fit_lines) else
  data.table(Stress=character(),k=numeric(),Freq=numeric())
fit_stats <- if(length(fit_summary)) rbindlist(fit_summary) else
  data.table(Stress=character(),DegreeType=character(),
             LogLog_OLS_Slope=numeric(),LogLog_OLS_R2=numeric(),
             Unique_Degree_Bins=integer(),Interpretation=character())
fwrite(degree_pdf,file.path(output_dir,'10_Degree_LogLog_Frequency_Data.csv'))
fwrite(fit_stats,file.path(output_dir,'11_Degree_LogLog_Descriptive_OLS.csv'))

degree_pdf[, Stress := factor(Stress,levels=stress_levels)]
fit_df[, Stress := factor(Stress,levels=stress_levels)]
p_e <- ggplot() +
  geom_point(data=degree_pdf,aes(x=k,y=Freq,color=Stress),
             size=1.65,alpha=0.65) +
  {if(show_descriptive_loglog_fit && nrow(fit_df)>0L)
    geom_line(data=fit_df,aes(x=k,y=Freq,color=Stress),
              linetype='dashed',linewidth=0.9)
    else NULL} +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides='bl') +
  scale_color_manual(values=c('Alkaline'='#FA7065',
                              'Drought'='#08BF43','Cold'='#6198F5'),drop=FALSE) +
  labs(title='Degree Distribution (Log-Log Scale) of Rice GRNs',
       subtitle=if(show_descriptive_loglog_fit)
         'Dashed lines: descriptive log-log OLS only; NOT a power-law or scale-free test'
       else 'Empirical distributions only; no fitted model',
       x='Node Degree k (log10)',y='Frequency P(k) (log10)',color='Stress') +
  theme_bw(base_size=12) +
  theme(text=element_text(color='black'),
        plot.title=element_text(hjust=0.5,face='bold',color='black'),
        axis.text=element_text(color='black'),
        axis.title=element_text(color='black'),
        panel.grid.minor=element_blank(),legend.position='right')
p_e <- p_e + black_text_theme
ggsave(file.path(output_dir,'Figure_E_All_Stress_Degree_LogLog.pdf'),p_e,
       width=270,height=180,units='mm',bg='white',device=grDevices::pdf,
       useDingbats=FALSE)
ggsave(file.path(output_dir,'Figure_E_All_Stress_Degree_LogLog.png'),p_e,
       width=270,height=180,units='mm',dpi=400,bg='white')


message('Results and figures saved to: ',normalizePath(output_dir,winslash='/'))

# Completion flag contains the complete input path list and output summary location.
writeLines(c(
  paste0('Output dir: ', normalizePath(output_dir, winslash='/', mustWork=TRUE)),
  sprintf('Ranked network inputs (top %g%% analyzed):', 100*network_edge_fraction),
  vapply(names(network_files), function(x) paste(x, network_files[[x]], sep=' : '), character(1)),
  paste0('ER repetitions per stress: ', n_rand),
  'Observed topology, ER replicates and figures are complete.'
), file.path(output_dir, 'NEW_TOP10_ER_ANALYSIS_COMPLETE.txt'))
