# gsGNA Infomap stability from saved random-seed runs
#
# Run from the repository root:
#   Rscript code/Infomap_stability_analysis.R
#
# GSGNA_ENSEMBLE_RUN sets the ensemble directory name under the project root
# (default: gsGNA_FixedTop10).
# GSGNA_INFOMAP_CACHE_DIR accepts an absolute path to an existing cache
# (default: results/GRN_Modules/Infomap/MultiSeed within the project).
# The cache contains 00_Inputs_Manifest.csv, 00_Run_Signatures.rds, and
# <Drought|Alkaline|Cold>/Checkpoints/Run_*.rds. Network and module
# checksums identify the inputs used for the cached runs.
# GSGNA_EXPORT_PAIR_DETAIL=1 also exports individual gene-pair statistics.
#
# Calculates ARI, AMI and NMI on jointly assigned nodes for each pair of runs.
# For the original five modules per stress, calculates within-module gene-pair
# coassignment across alternative runs, with the original seed excluded.
# Requires data.table and aricode. ggplot2 enables optional figures.

for (pkg in c("data.table", "aricode")) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    stop("Missing package: ", pkg, ". Run install.packages('", pkg, "') first.")
  }
}
suppressPackageStartupMessages(library(data.table))

# ------------------------- A. Input/output configuration -----------------------
project_root <- Sys.getenv(
  "GSGNA_PROJECT_ROOT", unset = Sys.getenv("GSGNA_INPUT_ROOT", unset = ".")
)
ensemble_run <- Sys.getenv("GSGNA_ENSEMBLE_RUN", unset = "gsGNA_FixedTop10")
infomap_root <- file.path(project_root, "results", "GRN_Modules", "Infomap")
input_root <- Sys.getenv(
  "GSGNA_INFOMAP_CACHE_DIR", unset = file.path(infomap_root, "MultiSeed")
)
output_root <- file.path(infomap_root, "Stability")
stresses <- c("Drought", "Alkaline", "Cold")
network_files <- setNames(file.path(
  project_root, ensemble_run, "FINAL_Locked_Test_And_Network",
  stresses, "Final_Ranked_Network", "XGBoost_SelectedUnion_Top10Percent_Edges.csv"
), stresses)
original_top5 <- list(
  Drought = c("M00064", "M00054", "M00051", "M00058", "M00088"),
  Alkaline = c("M00066", "M00071", "M00060", "M00104", "M00076"),
  Cold = c("M00003", "M00049", "M00009", "M00030", "M00004")
)

# Optional detailed file of every gene pair in the original 15 modules.
# Off by default: it may be large, but is not needed for summary statistics.
export_pair_detail <- identical(Sys.getenv("GSGNA_EXPORT_PAIR_DETAIL", "0"), "1")

if (!dir.exists(input_root)) {
  stop("Cannot find finished Infomap multi-seed output: ", input_root,
       "\nThis code reads existing Run_*.rds checkpoints; it does NOT rerun Infomap.")
}
manifest_file <- file.path(input_root, "00_Inputs_Manifest.csv")
signature_file <- file.path(input_root, "00_Run_Signatures.rds")
if (!file.exists(manifest_file) || !file.exists(signature_file)) {
  stop("Missing run manifest/signatures in: ", input_root)
}
manifest <- fread(manifest_file, showProgress=FALSE)
required_fields <- c("Stress", "Network", "Network_MD5", "Original_RDS")
if (!all(required_fields %in% names(manifest))) {
  stop("Missing manifest columns: ", paste(setdiff(required_fields,names(manifest)), collapse=", "))
}
signatures <- readRDS(signature_file)
if (is.null(names(signatures)) || anyDuplicated(names(signatures))) {
  stop("Expected uniquely named run signatures in: ", signature_file)
}
if (anyDuplicated(manifest$Stress)) stop("Duplicate stresses in 00_Inputs_Manifest.csv")

dir.create(output_root, recursive=TRUE, showWarnings=FALSE)

# ----------------------- B. Exact node order reconstruction --------------------
# MUST exactly match the original 30-run script's read_network_pairs(); otherwise
# cached numeric member IDs cannot be mapped to the correct genes.
get_original_network_nodes <- function(path) {
  if (!file.exists(path)) stop("Original locked network no longer exists: ", path)
  hdr <- names(fread(path, nrows=0L, showProgress=FALSE))
  if (!all(c("TF","Target") %in% hdr)) {
    stop("Network must contain TF and Target: ",path)
  }
  scores <- intersect(c("Weight","EdgeWeight","XGB_Score","Score"),hdr)
  if (!length(scores)) stop("Original network lacks score column: ",path)
  score_col <- scores[[1L]]
  edges <- fread(path, select=c("TF","Target",score_col), showProgress=FALSE)
  setnames(edges, score_col, "Score")
  edges[, `:=`(
    TF=toupper(trimws(as.character(TF))),
    Target=toupper(trimws(as.character(Target))),
    Score=suppressWarnings(as.numeric(Score))
  )]
  edges <- edges[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) &
                   !(TF %chin% c("NA","NAN","NULL")) &
                   !(Target %chin% c("NA","NAN","NULL")) &
                   TF != Target & is.finite(Score)]
  setorder(edges, TF, Target, -Score)
  edges <- unique(edges, by=c("TF","Target"))
  sort(unique(c(edges$TF,edges$Target)))
}

clean_ids <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sort(unique(x[!is.na(x) & nzchar(x) &
    !(x %in% c("NA", "NAN", "NULL", "GENE", "GID", "TF"))]))
}

# Avoid failures/warnings from empty strata and undefined metrics.
safe_metric <- function(fun, a, b) {
  z <- tryCatch(suppressWarnings(as.numeric(fun(a,b))),error=function(e) NA_real_)
  if (length(z)!=1L || !is.finite(z)) NA_real_ else z
}
stat_mean <- function(x) if (length(x) && any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
stat_median <- function(x) if (length(x) && any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
stat_q <- function(x,p) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else as.numeric(quantile(x,p,names=FALSE))
}

# ---------------------- C. Read all saved run memberships ----------------------
load_one_stress <- function(s) {
  path <- network_files[[s]]
  if (!file.exists(path)) stop("Missing final Top 10% GRN for ", s, ": ", path)
  actual_md5 <- unname(tools::md5sum(path))
  # Prefer the canonical label. Use the checksum to identify a row in a
  # historical manifest so the original cache signature stays intact.
  rr <- manifest[Stress == s]
  if (!nrow(rr)) rr <- manifest[tolower(Network_MD5) == tolower(actual_md5)]
  if (nrow(rr) != 1L || !isTRUE(tolower(rr$Network_MD5[[1L]]) == tolower(actual_md5))) {
    stop(s, ": expected one cache-manifest entry matching the input network MD5: ", path)
  }
  cache_key <- rr$Stress[[1L]]
  if (!cache_key %in% names(signatures)) {
    stop("Missing run signature for the cache entry used by ", s)
  }
  nodes <- get_original_network_nodes(path)
  if (!length(nodes)) stop(s,": empty node list")
  cpdir <- file.path(input_root, s, "Checkpoints")
  files <- list.files(cpdir,pattern="^Run_[0-9]{3}\\.rds$",full.names=TRUE)
  if (length(files)<3L) stop(s,": fewer than 3 checkpoint files under ",cpdir)
  runs <- as.integer(sub("^Run_([0-9]{3})\\.rds$","\\1",basename(files)))
  o <- order(runs)
  runs <- runs[o]; files <- files[o]
  if (anyDuplicated(runs) || !identical(runs,seq_len(length(runs)))) {
    stop(s,": Run_001..Run_NNN checkpoints must be consecutive, with no gaps")
  }
  if (length(runs)!=30L) {
    warning(s,": found ",length(runs)," run checkpoints (expected 30); ",
            "the report will accurately show this count")
  }
  mem <- vector("list",length(files))
  seeds <- integer(length(files))
  cover <- vector("list",length(files))
  for (i in seq_along(files)) {
    z <- readRDS(files[[i]])
    if (!is.list(z) || is.null(z$mem) || is.null(z$signature) || is.null(z$seed)) {
      stop("Invalid checkpoint: ",files[[i]])
    }
    if (!identical(z$signature,signatures[[cache_key]])) {
      stop("Checkpoint signature mismatch: ",files[[i]])
    }
    m <- as.integer(z$mem)
    if (length(m)!=length(nodes) || anyNA(m) || any(m<0L)) {
      stop("Node membership length/values mismatch: ",files[[i]],
           "; expected ",length(nodes)," nodes, observed ",length(m))
    }
    mem[[i]] <- m
    seeds[[i]] <- as.integer(z$seed)
    cover[[i]] <- data.table(
      Stress=s,Manuscript_Stress=s,
      Run=runs[[i]],Random_Seed=seeds[[i]],
      Is_Original_Seed=(runs[[i]]==1L),
      Network_Nodes=length(nodes),
      Assigned_Nodes=sum(m>0L),
      Unassigned_Nodes=sum(m==0L),
      Assignment_Fraction=mean(m>0L),
      Detected_Modules=length(unique(m[m>0L]))
    )
  }
  if (anyDuplicated(seeds)) stop(s,": duplicate random seeds; cannot treat runs as independent starts")
  # Check original RDS and exact original seed checkpoint.
  original_rds <- rr$Original_RDS[[1L]]
  if (!file.exists(original_rds)) {
    original_rds <- file.path(infomap_root, s, "Retained_Modules_10_300.rds")
  }
  if (!file.exists(original_rds)) stop("Original Infomap modules RDS missing: ", original_rds)
  if ("Original_RDS_MD5" %in% names(rr)) {
    expected_md5 <- rr$Original_RDS_MD5[[1L]]
    if (!is.na(expected_md5) && nzchar(expected_md5) &&
        !isTRUE(tolower(unname(tools::md5sum(original_rds))) == tolower(expected_md5))) {
      stop(s, ": original Infomap module file does not match the cached checksum: ", original_rds)
    }
  }
  modall <- readRDS(original_rds)
  if (!is.list(modall) || is.null(names(modall))) stop("Invalid original module RDS: ",original_rds)
  ids <- original_top5[[s]]
  if (!all(ids %in% names(modall))) stop(s,": missing original selected modules")
  refs <- lapply(modall[ids],clean_ids)
  for (mid in names(refs)) {
    indices <- match(refs[[mid]],nodes)
    if (anyNA(indices)) stop(s,"/",mid,": original gene not found in original network")
    origlabels <- mem[[1L]][indices]
    if (any(origlabels==0L) || length(unique(origlabels))!=1L) {
      stop(s,"/",mid,": original-seed checkpoint does not reproduce original module")
    }
    # Also exclude extra genes in that checkpoint module.
    if (sum(mem[[1L]] == origlabels[[1L]])!=length(indices)) {
      stop(s,"/",mid,": original checkpoint module includes extra genes")
    }
  }
  list(stress=s,nodes=nodes,memberships=mem,seeds=seeds,refs=refs,
       run_coverage=rbindlist(cover),n_runs=length(mem),
       network=path,original_rds=original_rds)
}

# ----------------------- D. Entire-partition ARI/AMI/NMI -----------------------
whole_partition_agreement <- function(x) {
  M <- x$memberships
  k <- length(M)
  tab <- vector("list",choose(k,2L))
  pos <- 0L
  for (i in seq_len(k-1L)) {
    a <- M[[i]]
    for (j in seq.int(i+1L,k)) {
      b <- M[[j]]
      mask <- a>0L & b>0L
      n <- sum(mask)
      # If original partitions are incomplete, these are conditional on
      # comparable, jointly assigned genes. No phantom 'unassigned cluster'.
      a2 <- a[mask]; b2 <- b[mask]
      ok <- n>=3L && length(unique(a2))>=1L && length(unique(b2))>=1L
      ARI <- if (ok) safe_metric(aricode::ARI,a2,b2) else NA_real_
      AMI <- if (ok) safe_metric(aricode::AMI,a2,b2) else NA_real_
      NMI <- if (ok) safe_metric(aricode::NMI,a2,b2) else NA_real_
      full <- all(mask)
      pos <- pos+1L
      tab[[pos]] <- data.table(
        Stress=x$stress,Manuscript_Stress=x$stress,
        Run_1=i,Run_2=j,Seed_1=x$seeds[[i]],Seed_2=x$seeds[[j]],
        Comparison_Type=if (i==1L) "Original_vs_Alternative" else "Alternative_vs_Alternative",
        Total_Nodes=length(a),Both_Assigned_Nodes=n,
        Both_Assigned_Fraction=n/length(a),
        # These three columns are only valid on the shared assigned-node set.
        ARI_Conditional=ARI,AMI_Conditional=AMI,NMI_Conditional=NMI,
        # Strict full network score stays NA unless every node assigned.
        ARI_Strict_Full=if (full) ARI else NA_real_,
        AMI_Strict_Full=if (full) AMI else NA_real_,
        NMI_Strict_Full=if (full) NMI else NA_real_,
        Full_Network_Partition_Comparison=full
      )
    }
  }
  out <- rbindlist(tab)
  expected <- choose(k,2L)
  stopifnot(nrow(out)==expected)
  out
}

# -------------------- E. Pair coassignment per ORIGINAL module -----------------
# Each gene pair has equal weight in module summary. We do NOT restrict to
# adjacency edges. This is all pairs WITHIN EACH of the 15 original modules.
module_pair_coassignment <- function(x) {
  alts <- x$memberships[-1L] # remove original random seed: avoids self-inflation
  R <- length(alts)
  if (R<2L) stop("Need at least 2 alternative runs for pair coassignment")
  summaries <- vector("list",length(x$refs))
  details <- if (export_pair_detail) vector("list",length(x$refs)) else NULL
  mid_names <- names(x$refs)
  for (idx in seq_along(mid_names)) {
    mid <- mid_names[[idx]]
    genes <- x$refs[[mid]]
    n <- length(genes)
    if (n<2L) stop("Reference module has fewer than two genes: ",mid)
    gi <- match(genes,x$nodes)
    if (anyNA(gi)) stop("Gene mapping mismatch: ",x$stress,"/",mid)
    same <- matrix(0L,n,n)
    observed <- matrix(0L,n,n)
    assigned_by_gene <- integer(n)
    for (r in seq_along(alts)) {
      z <- alts[[r]][gi]
      good <- z>0L
      assigned_by_gene <- assigned_by_gene + as.integer(good)
      both <- outer(good,good,"&")
      observed <- observed + both
      same <- same + (outer(z,z,"==") & both)
    }
    ut <- upper.tri(same)
    num_same <- same[ut]
    num_both <- observed[ut]
    both_valid <- num_both>0L
    conditional <- rep(NA_real_,length(num_same))
    conditional[both_valid] <- num_same[both_valid]/num_both[both_valid]
    strict <- num_same/R  # unassigned runs count as not recovered
    summaries[[idx]] <- data.table(
      Stress=x$stress,Manuscript_Stress=x$stress,
      Original_Module=mid,Module_Size=n,
      Alternative_Random_Seeds=R,
      Total_WithinModule_GenePairs=length(num_same),
      Pairs_With_AtLeast_One_Joint_Assignment=sum(both_valid),
      Pairs_Never_Jointly_Assigned=sum(!both_valid),
      Pair_Observation_Coverage=mean(both_valid),
      Mean_Joint_Assignment_Fraction=mean(num_both/R),
      Mean_Gene_Assignment_Fraction=mean(assigned_by_gene/R),
      Mean_Pair_Coassignment_Conditional=stat_mean(conditional),
      Median_Pair_Coassignment_Conditional=stat_median(conditional),
      Q25_Pair_Coassignment_Conditional=stat_q(conditional,0.25),
      Q75_Pair_Coassignment_Conditional=stat_q(conditional,0.75),
      Fraction_Observed_Pairs_Coassignment_GEQ050=if (any(both_valid))
        mean(conditional[both_valid]>=0.50) else NA_real_,
      Fraction_Observed_Pairs_Coassignment_GEQ080=if (any(both_valid))
        mean(conditional[both_valid]>=0.80) else NA_real_,
      Mean_Pair_Coassignment_AllAlternativeRuns=mean(strict)
    )
    if (export_pair_detail) {
      ij <- which(ut,arr.ind=TRUE)
      details[[idx]] <- data.table(
        Stress=x$stress,Manuscript_Stress=x$stress,
        Original_Module=mid,Gene_1=genes[ij[,1]],Gene_2=genes[ij[,2]],
        Both_Assigned_Runs=num_both,Same_Community_Runs=num_same,
        Coassignment_Conditional=conditional,
        Coassignment_AllAlternativeRuns=strict
      )
    }
    message("  [",x$stress,"] ",mid,": ",length(num_same),
            " gene pairs; mean conditional coassignment=",
            sprintf("%.4f",summaries[[idx]]$Mean_Pair_Coassignment_Conditional))
  }
  list(summary=rbindlist(summaries),
       pair_detail=if (export_pair_detail) rbindlist(details) else NULL)
}

# ---------------------------- F. Execute cached analyses -----------------------
all_pairs <- list()
all_module <- list()
all_coverage <- list()
manifest_out <- list()
for (s in stresses) {
  message("\n=================== ",s," ===================")
  x <- load_one_stress(s)
  message("Loaded ",x$n_runs," existing runs; ",length(x$nodes),
          " nodes; ",length(x$refs)," original modules.")
  pw <- whole_partition_agreement(x)
  mc <- module_pair_coassignment(x)
  all_pairs[[s]] <- pw
  all_module[[s]] <- mc$summary
  all_coverage[[s]] <- x$run_coverage
  manifest_out[[s]] <- data.table(
    Stress=s,Manuscript_Stress=s,
    Network=x$network,Network_MD5=unname(tools::md5sum(x$network)),
    Original_Modules_RDS=x$original_rds,
    Cached_Runs=x$n_runs,
    Pairwise_Partition_Comparisons=nrow(pw),
    Alternative_Seeds_Used_For_Module_Coassignment=x$n_runs-1L
  )
  if (export_pair_detail) {
    fwrite(mc$pair_detail,file.path(output_root,
      paste0(s,"_WithinModule_GenePair_Coassignment_AllPairs.csv.gz")))
  }
  rm(x,pw,mc); invisible(gc(verbose=FALSE))
}
all_pairs <- rbindlist(all_pairs)
all_module <- rbindlist(all_module)
all_coverage <- rbindlist(all_coverage)
manifest_out <- rbindlist(manifest_out)

summary_agreement <- all_pairs[, .(
  Pair_Comparisons=.N,
  Comparisons_With_Valid_ARI=sum(is.finite(ARI_Conditional)),
  Median_Both_Assigned_Fraction=stat_median(Both_Assigned_Fraction),
  Min_Both_Assigned_Fraction=min(Both_Assigned_Fraction),
  Median_ARI=stat_median(ARI_Conditional),
  Q25_ARI=stat_q(ARI_Conditional,.25),Q75_ARI=stat_q(ARI_Conditional,.75),
  Median_AMI=stat_median(AMI_Conditional),
  Q25_AMI=stat_q(AMI_Conditional,.25),Q75_AMI=stat_q(AMI_Conditional,.75),
  Median_NMI=stat_median(NMI_Conditional),
  Q25_NMI=stat_q(NMI_Conditional,.25),Q75_NMI=stat_q(NMI_Conditional,.75)
),by=.(Stress,Manuscript_Stress,Comparison_Type)]

fwrite(manifest_out,file.path(output_root,"00_Input_and_Cache_Manifest.csv"),bom=TRUE)
fwrite(all_coverage,file.path(output_root,"01_Run_Assignment_Coverage.csv"),bom=TRUE)
fwrite(all_pairs,file.path(output_root,"02_AllRunPairs_ARI_AMI_NMI.csv"),bom=TRUE)
fwrite(summary_agreement,file.path(output_root,"03_WholePartition_ARI_AMI_NMI_Summary.csv"),bom=TRUE)
fwrite(all_module,file.path(output_root,"04_Original15_WithinModule_PairCoassignment_Summary.csv"),bom=TRUE)

# Optional figures (no extra packages required for the numerical results).
if (requireNamespace("ggplot2",quietly=TRUE)) {
  dlong <- melt(
    all_pairs[,.(Stress,Comparison_Type,ARI=ARI_Conditional,
                 AMI=AMI_Conditional,NMI=NMI_Conditional)],
    id.vars=c("Stress","Comparison_Type"),
    variable.name="Metric",value.name="Agreement"
  )
  p <- ggplot2::ggplot(dlong,ggplot2::aes(x=Stress,y=Agreement,fill=Comparison_Type)) +
    ggplot2::geom_boxplot(outlier.size=.4,width=.65,na.rm=TRUE) +
    ggplot2::facet_wrap(~Metric,nrow=1,scales="free_y") +
    ggplot2::labs(x=NULL,y="Agreement on jointly assigned nodes",
                  fill="Comparison",title="Infomap partition agreement across random seeds") +
    ggplot2::theme_bw(base_size=11)
  ggplot2::ggsave(file.path(output_root,"Fig1_WholePartition_ARI_AMI_NMI.png"),
                  p,width=12,height=4.8,dpi=300)

  p2 <- ggplot2::ggplot(all_module,
    ggplot2::aes(x=Original_Module,y=Mean_Pair_Coassignment_Conditional)) +
    ggplot2::geom_col(fill="#4E79A7",width=.7) +
    ggplot2::facet_wrap(~Manuscript_Stress,scales="free_x") +
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::labs(x="Original Infomap module",
       y="Mean within-module pair coassignment",
       title="Original 15 modules: gene-pair stability across alternative seeds",
       subtitle="Conditioned on both genes being assigned; original seed excluded") +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
  ggplot2::ggsave(file.path(output_root,"Fig2_Original15_Pair_Coassignment.png"),
                  p2,width=11,height=5,dpi=300)
}

message("\nAll cache-only stability analyses completed: ",output_root)
message("\nWhole partition ARI/AMI/NMI summaries (conditional, coverage reported):")
print(summary_agreement[,.(Stress,Comparison_Type,Pair_Comparisons,
                           Median_Both_Assigned_Fraction,
                           Median_ARI,Median_AMI,Median_NMI)])
message("\nOriginal 15 modules: within-module mean gene-pair coassignment:")
print(all_module[,.(Stress,Original_Module,Alternative_Random_Seeds,
                    Mean_Pair_Coassignment_Conditional,
                    Mean_Pair_Coassignment_AllAlternativeRuns,
                    Mean_Joint_Assignment_Fraction)])
