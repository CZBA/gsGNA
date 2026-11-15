# gsGNA hub TF identification and three-stress overlap
#
# Run from the repository root: Rscript code/hub_TF_identification.R
# Uses the locked XGBoost combinations and ranked selected-union exports.
# Retains the top 10% of the ranked edges for each stress before computing
# TF outdegrees. Hub thresholds are 5%, 10%, and 20% of positive-outdegree
# TFs per stress, with ties included.
#
# The primary null samples edges uniformly from each stress's selected
# base-method union. TF-label permutation keeps the observed outdegree
# distribution fixed. Both assess three-stress shared hub counts.
# Optional directed degree-preserving rewiring assesses TF-target pair overlap.

required_packages <- c('data.table', 'ggplot2')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_packages)) stop('Install the required R packages: ',paste(missing_packages,collapse=', '))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------- USER SETTINGS -------------------------------
input_root <- Sys.getenv('GSGNA_INPUT_ROOT', unset='.')
ensemble_dir <- file.path(input_root, 'gsGNA_FixedTop10')
# Retain the original stress order for reproducible random draws.
stress_names <- c('Drought', 'Alkaline', 'Cold')
network_files <- setNames(file.path(
  ensemble_dir, 'FINAL_Locked_Test_And_Network', stress_names, 'Final_Ranked_Network',
  'XGBoost_Ranked_SelectedUnion_AllEdges.csv'
), stress_names)
selection_file <- file.path(ensemble_dir, 'FINAL_PerStress_Models_LOCKED.csv')
input_stage <- 'top10'           # Input threshold for the base GRN algorithms.
network_edge_fraction <- 0.10   # Fraction of ranked ensemble edges used for hub analysis.
method_names <- c('GENIE3', 'Kboost', 'GRNBoost2', '3DCEMA', 'DeepRIG', 'IGEGRNs')
output_dir <- file.path(input_root, 'gsGNA_Hub_TF_Analysis')

# Complementary cross-stress label permutation: keeps each observed TF degree
# distribution and number of nonzero-degree TFs fixed, relabels eligible TFs.
run_tf_label_permutation <- TRUE
n_label_permutation <- 1000L

# Holm correction is across the three tested hub thresholds WITHIN each null
# family; these are not TF-level adjusted P-values.
threshold_p_adjust <- 'holm'


primary_fraction <- 0.10        # candidate hub TF outdegree top 10%
sensitivity_fractions <- c(0.05, 0.10, 0.20)
n_random <- 1000L
random_seed <- 123L
verify_final_subset <- TRUE     # check every final edge belongs to its OWN locked union
progress_every <- 100L

# Optional directed degree-preserving rewiring. Can be expensive even for the
# common-hub subnetwork. To produce the complementary structural comparison,
# explicitly set run_degree_preserving <- TRUE and assess mixing/convergence.
run_degree_preserving <- FALSE
rewire_scope <- 'common_hubs'  # local test; 'full' tests the whole network
n_rewire <- 199L               # paper: >=999 when computationally feasible
rewire_trials_per_edge <- 2     # ATTEMPTED swaps / edge / replicate; assess mixing at 5-10
# -----------------------------------------------------------------------------

stopifnot(identical(input_stage, 'top10'),
          is.finite(network_edge_fraction), network_edge_fraction > 0,
          network_edge_fraction <= 1,
          threshold_p_adjust %in% c('holm','BH'), n_label_permutation>=99L,
          is.finite(primary_fraction), primary_fraction > 0, primary_fraction <= 1,
          length(sensitivity_fractions) >= 1L,
          all(is.finite(sensitivity_fractions)),
          all(sensitivity_fractions > 0 & sensitivity_fractions <= 1),
          any(abs(sensitivity_fractions - primary_fraction) < 1e-12),
          n_random >= 99L, progress_every >= 1L,
          rewire_scope %in% c('common_hubs','full'))
fractions <- sort(unique(sensitivity_fractions))
frac_key <- function(x) paste0('Top', format(100*x, scientific=FALSE, trim=TRUE), 'pct')
frac_label <- function(x) sprintf('Top %g%%', 100*x)

read_locked_methods <- function(path) {
  if (!file.exists(path)) stop('Missing locked model file: ', path)
  models <- fread(path, showProgress=FALSE)
  required <- c('Dataset', 'Strategy', 'Combination', 'Input_Stage')
  if (!all(required %in% names(models))) {
    stop('Missing columns in the locked model file: ',
         paste(setdiff(required, names(models)), collapse=', '))
  }
  if (nrow(models) != length(stress_names) || anyDuplicated(models$Dataset) ||
      !setequal(models$Dataset, stress_names)) {
    stop('Expected one locked model for each stress: ', paste(stress_names, collapse=', '))
  }
  if (anyNA(models$Strategy) || any(models$Strategy != 'XGBoost') ||
      anyNA(models$Input_Stage) || any(models$Input_Stage != input_stage)) {
    stop('The locked model file must contain XGBoost models at input stage ', input_stage, '.')
  }
  combinations <- setNames(as.character(models$Combination), models$Dataset)
  selected <- lapply(combinations[stress_names], function(combination) {
    if (is.na(combination)) return(character())
    trimws(strsplit(combination, '+', fixed=TRUE)[[1L]])
  })
  invalid <- !vapply(selected, function(methods) {
    length(methods) > 0L && all(methods %in% method_names) && !anyDuplicated(methods)
  }, logical(1))
  if (any(invalid)) {
    stop('Invalid locked base-algorithm combination for: ',
         paste(names(selected)[invalid], collapse=', '))
  }
  selected
}

selected_methods <- read_locked_methods(selection_file)
missing_networks <- network_files[!file.exists(network_files)]
if (length(missing_networks)) {
  stop('Missing ranked XGBoost networks: ', paste(missing_networks, collapse=' | '))
}
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
fwrite(rbindlist(lapply(names(network_files),function(ds) data.table(
  Stress=ds,Stage=input_stage,Selected_Methods=paste(selected_methods[[ds]],collapse='+'),
  Selected_Method_Count=length(selected_methods[[ds]]),
  Network_Edge_Fraction=network_edge_fraction, Ranked_Edge_File=network_files[[ds]]))),
  file.path(output_dir,'00_Input_Manifest.csv'))

# ------------------------ 1. INPUT + PAIR NORMALIZATION ----------------------
clean_pairs <- function(z) {
  if (!all(c('TF','Target') %in% names(z))) stop('The edge table must contain TF and Target columns.')
  z <- z[, .(TF=toupper(trimws(as.character(TF))),
             Target=toupper(trimws(as.character(Target))))]
  z <- z[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) & TF != Target]
  unique(z, by=c('TF','Target'))
}

read_final <- function(file, candidate_count) {
  dt <- fread(file, showProgress=FALSE, select=c('TF', 'Target', 'EdgeWeight'))
  if (!all(c('TF', 'Target', 'EdgeWeight') %in% names(dt))) {
    stop('The ranked network must contain TF, Target, and EdgeWeight: ', file)
  }
  dt[, EdgeWeight := suppressWarnings(as.numeric(EdgeWeight))]
  if (any(!is.finite(dt$EdgeWeight))) stop('Non-finite ranked edge weights: ', file)
  if (!candidate_count || nrow(dt) != candidate_count) {
    stop('Expected a full ranked export of ', candidate_count, ' candidate edges, but found ',
         nrow(dt), ': ', file)
  }
  # Use the same score and gene-ID ordering as the ensemble export.
  setorderv(dt, c('EdgeWeight', 'TF', 'Target'), c(-1L, 1L, 1L))
  n_keep <- ceiling(candidate_count * network_edge_fraction)
  ans <- clean_pairs(dt[seq_len(n_keep)])
  if (nrow(ans) != n_keep) stop('Duplicate or invalid edges in the retained network: ', file)
  ans
}

method_file <- function(ds, method) {
  folders <- if (method=='GRNBoost2') c('GRNBoost2','GRNBoost') else method
  options <- file.path(input_root, folders, paste0(ds, '_', input_stage, '_edges.csv'))
  found <- options[file.exists(options)]
  if (!length(found)) stop('Missing base-algorithm candidate file: ', paste(options,collapse=' | '))
  found[1L]
}

read_candidate_union <- function(ds) {
  methods_here <- selected_methods[[ds]]
  if (is.null(methods_here) || !length(methods_here))
    stop('Missing locked base-algorithm combination for: ',ds)
  pieces <- lapply(methods_here, function(method) {
    p <- method_file(ds, method)
    dt <- fread(p, showProgress=FALSE)
    if (ncol(dt) < 3L) stop('The candidate edge file must have at least three columns: ', p)
    # Exactly like original read_grn, retain only finite weights.
    w <- suppressWarnings(as.numeric(dt[[3L]]))
    tmp <- data.table(TF=dt[[1L]], Target=dt[[2L]])[is.finite(w)]
    rm(dt, w)
    clean_pairs(tmp)
  })
  cand <- unique(rbindlist(pieces, use.names=TRUE), by=c('TF','Target'))
  rm(pieces)
  setkey(cand, TF, Target)
  cand
}

outdegree <- function(ed) ed[, .(Degree=.N), by=TF]

hub_from_counts <- function(tf, deg, pct) {
  ok <- is.finite(deg) & deg > 0
  if (!any(ok)) return(character())
  values <- sort(deg[ok], decreasing=TRUE)
  k <- max(1L, ceiling(length(values)*pct))
  threshold <- values[k]
  tf[ok & deg >= threshold]  # include ALL TFs tied at cutoff
}

summarize_fractions <- function(tf, deg, fractions) {
  lapply(fractions, function(p) hub_from_counts(tf,deg,p))
}

# ---- Read stress data once; verify the final edges are a subset of candidates
cat('\n[Inputs] Reading selected-method candidate unions and ranked XGBoost networks\n')
caps <- vector('list', length(network_files)); names(caps) <- names(network_files)
final_edges <- vector('list', length(network_files)); names(final_edges) <- names(network_files)
candidate_hubs <- vector('list', length(network_files)); names(candidate_hubs) <- names(network_files)
observed_hubs <- vector('list', length(network_files)); names(observed_hubs) <- names(network_files)
input_audit <- vector('list', length(network_files))

for (s in names(network_files)) {
  cand <- read_candidate_union(s)
  final <- read_final(network_files[[s]], nrow(cand))
  if (verify_final_subset) {
    ix <- cand[final, on=.(TF, Target), which=TRUE]
    n_missing <- sum(is.na(ix))
    if (n_missing > 0L) stop(s, ': ',n_missing,
                             ' retained edges are outside the selected-method candidate union.')
  } else n_missing <- NA_integer_

  cd <- cand[, .(Candidate_Degree=.N), by=TF]
  fd <- outdegree(final)
  if (nrow(final) > nrow(cand)) stop(s, ': The retained network is larger than the candidate union.')
  # Retain candidate TF universe and the FULL observed TF degree multiset,
  # including zeros: this is the label-permutation null's reference population.
  observed_deg <- fd$Degree[match(cd$TF,fd$TF)]
  observed_deg[is.na(observed_deg)] <- 0
  if (sum(observed_deg)!=nrow(final)) stop(s, ': The sum of TF outdegrees differs from the retained edge count.')
  caps[[s]] <- list(TF=cd$TF, Capacity=as.numeric(cd$Candidate_Degree),
                    Observed_Degree=as.numeric(observed_deg),
                    N_Candidate=nrow(cand), N_Final=nrow(final))
  observed_hubs[[s]] <- lapply(fractions, function(p)
    hub_from_counts(fd$TF, fd$Degree, p))
  candidate_hubs[[s]] <- lapply(fractions, function(p)
    hub_from_counts(cd$TF, cd$Candidate_Degree, p))
  names(observed_hubs[[s]]) <- names(candidate_hubs[[s]]) <- vapply(fractions, frac_key, character(1))
  final_edges[[s]] <- final
  input_audit[[match(s,names(network_files))]] <- data.table(
    Stress=s, Candidate_Union_Edges=nrow(cand), Final_Edges=nrow(final),
    Locked_Selected_Methods=paste(selected_methods[[s]],collapse='+'),
    Candidate_N_PositiveDegreeTF=uniqueN(cand$TF), Final_N_PositiveDegreeTF=uniqueN(final$TF),
    Final_Edge_Retention=nrow(final)/nrow(cand), Missing_Final_Pairs=n_missing)
  cat(sprintf(' %-10s Candidate=%s, Final=%s, Final_TFs=%d\n', s,
              format(nrow(cand),big.mark=','), format(nrow(final),big.mark=','), uniqueN(final$TF)))
  rm(cand, final, cd, fd, observed_deg); gc(verbose=FALSE)
}
input_audit <- rbindlist(input_audit)
fwrite(input_audit, file.path(output_dir, '01_Input_Candidate_Final_Audit.csv'))

# -------------------------- 2. OBSERVED HUB RESULTS --------------------------
common_hub_sets <- function(hubs_by_stress, key) {
  Reduce(intersect, lapply(hubs_by_stress, `[[`, key))
}
observed_common <- candidate_common <- vector('list', length(fractions))
names(observed_common) <- names(candidate_common) <- vapply(fractions, frac_key, character(1))
obs_rows <- vector('list', length(fractions))
for (j in seq_along(fractions)) {
  key <- frac_key(fractions[j])
  observed_common[[key]] <- common_hub_sets(observed_hubs,key)
  candidate_common[[key]] <- common_hub_sets(candidate_hubs,key)
  obs_rows[[j]] <- data.table(
    Threshold=fractions[j], Threshold_Label=frac_label(fractions[j]),
    Candidate_Common_Hubs=length(candidate_common[[key]]),
    Observed_Common_Hubs=length(observed_common[[key]]),
    Drought_Hubs=length(observed_hubs$Drought[[key]]),
    Alkaline_Hubs=length(observed_hubs$Alkaline[[key]]),
    Cold_Hubs=length(observed_hubs$Cold[[key]]))
  fwrite(data.table(TF=sort(observed_common[[key]])),
         file.path(output_dir,paste0('02_Observed_',key,'_CommonHubTFs.csv')))
}
observed_summary <- rbindlist(obs_rows)
fwrite(observed_summary,file.path(output_dir,'02_Observed_Candidate_Hub_Summary.csv'))

# Export every stress-specific hub TF (not only the three-stress intersection)
# for downstream cross-stress comparisons and functional annotation.
hub_detail <- rbindlist(lapply(names(network_files), function(ds) {
  cc <- caps[[ds]]
  rbindlist(lapply(seq_along(fractions),function(j) {
    key <- frac_key(fractions[j]); selected <- observed_hubs[[ds]][[key]]
    data.table(Stress=ds,Threshold=fractions[j],
               Threshold_Label=frac_label(fractions[j]),TF=selected,
               Observed_Outdegree=as.numeric(cc$Observed_Degree[match(selected,cc$TF)]))
  }))
}))
fwrite(hub_detail,file.path(output_dir,'02_AllStress_Observed_HubTFs_ByThreshold.csv'))
cat('\nObserved and candidate-union hub counts\n'); print(observed_summary)


# ===================== 2b. EXPORT HUB TF + SHARED TARGET LISTS ==================
# This uses ONLY the observed Top10% outdegree hub sets and the FINAL XGBoost
# TF-Target edge lists already loaded above. No null model, P value or training
# is used to select individual TFs or common targets.
export_hub_tf_and_shared_targets <- function(observed_hubs, final_edges,
                                             output_dir, primary_fraction=0.10) {
  if (!all(stress_names %in% names(observed_hubs)) ||
      !all(stress_names %in% names(final_edges))) {
    stop('Missing hub sets or edge tables for one or more stresses.')
  }
  key <- paste0('Top', format(primary_fraction*100,scientific=FALSE,trim=TRUE),'pct')
  export_dir <- file.path(output_dir,'06_HubTF_Exclusive_Shared_Targets',key)
  per_tf_dir <- file.path(export_dir,'Per_Common_TF_Shared_Targets')
  dir.create(per_tf_dir,recursive=TRUE,showWarnings=FALSE)

  write_ids <- function(ids,file) {
    ids <- sort(unique(as.character(ids)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    writeLines(ids,file.path(export_dir,file),useBytes=TRUE)
    invisible(length(ids))
  }
  tf_sets <- lapply(stress_names,function(s) {
    z <- observed_hubs[[s]][[key]]
    if (is.null(z)) stop('Missing ', key, ' hub TF set for ', s)
    sort(unique(as.character(z)))
  })
  names(tf_sets) <- stress_names
  common3 <- Reduce(intersect,tf_sets)
  write_ids(common3,'ThreeStress_Common_HubTF.txt')
  union3 <- Reduce(union,tf_sets)
  write_ids(union3,'ThreeStress_Union_HubTF.txt')

  # Per-stress: all, exclusive to one stress, shared with >=1 other,
  # and shared by all three.  Files contain ONE TF ID per line.
  for (s in stress_names) {
    others <- setdiff(stress_names,s)
    other_union <- union(tf_sets[[others[1L]]],tf_sets[[others[2L]]])
    write_ids(tf_sets[[s]],paste0(s,'_All_HubTF.txt'))
    write_ids(setdiff(tf_sets[[s]],other_union),paste0(s,'_Specific_Only_HubTF.txt'))
    write_ids(intersect(tf_sets[[s]],other_union),paste0(s,'_Shared_AtLeastTwo_HubTF.txt'))
    write_ids(common3,paste0(s,'_Common_AllThree_HubTF.txt'))
    for (other in others) {
      write_ids(intersect(tf_sets[[s]],tf_sets[[other]]),
                paste0(s,'_SharedWith_',other,'_IncludingAllThree_HubTF.txt'))
    }
  }
  # Exactly-two-stress TFs: explicitly exclude those present in the third.
  for (pair in combn(stress_names,2L,simplify=FALSE)) {
    remaining <- setdiff(stress_names,pair)
    ids <- setdiff(intersect(tf_sets[[pair[1L]]],tf_sets[[pair[2L]]]),
                   tf_sets[[remaining]])
    write_ids(ids,paste0(paste(pair,collapse='_'),'_Only_SharedHubTF.txt'))
  }

  # Membership table: one row per unique Hub TF with Venn-region provenance.
  membership <- data.table::data.table(TF=sort(union3))
  for (s in stress_names) membership[,(paste0('In_',s)):=TF %in% tf_sets[[s]]]
  membership[,N_Stresses:=as.integer(In_Drought)+as.integer(In_Alkaline)+as.integer(In_Cold)]
  membership[,Category:=data.table::fifelse(N_Stresses==3L,'Shared_AllThree',
                                            data.table::fifelse(N_Stresses==2L,'Shared_ExactlyTwo','Stress_Specific'))]
  membership[,Stress_Sets:=vapply(TF,function(tf)
    paste(stress_names[vapply(stress_names,function(s)tf %in% tf_sets[[s]],logical(1))],
          collapse='+'),character(1))]
  data.table::fwrite(membership,file.path(export_dir,'HubTF_Stress_Membership.csv'))
  region_counts <- membership[,.(N_TFs=.N),by=.(Stress_Sets,Category)]
  data.table::setorder(region_counts,Category,Stress_Sets)
  data.table::fwrite(region_counts,file.path(export_dir,'HubTF_Venn_Region_Counts.csv'))

  # Strict target conservation: the SAME TF -> SAME TARGET pair must appear in
  # ALL THREE stress GRNs. Intersecting pooled Target IDs is NOT equivalent.
  per_stress_commonTF_edges <- lapply(stress_names,function(s) {
    d <- final_edges[[s]]
    if (!all(c('TF','Target') %in% names(d))) stop(s,' edge table is missing TF or Target.')
    d <- unique(d[TF %in% common3,.(TF,Target)],by=c('TF','Target'))
    data.table::setkeyv(d,c('TF','Target'))
    d
  })
  names(per_stress_commonTF_edges) <- stress_names
  strict_common_pairs <- Reduce(function(a,b)
    merge(a,b,by=c('TF','Target'),all=FALSE,sort=FALSE),
    per_stress_commonTF_edges)
  data.table::setorder(strict_common_pairs,TF,Target)
  # TXT table with TF and Target columns; also CSV for downstream R analysis.
  data.table::fwrite(strict_common_pairs,
                     file.path(export_dir,'ThreeStress_CommonHubTF_Strict_Shared_TF_Target_Pairs.txt'),
                     sep='\t',quote=FALSE)
  data.table::fwrite(strict_common_pairs,
                     file.path(export_dir,'ThreeStress_CommonHubTF_Strict_Shared_TF_Target_Pairs.csv'))
  write_ids(strict_common_pairs$Target,
            'ThreeStress_CommonHubTF_Strict_Shared_TargetIDs.txt')

  # One TXT per common TF, including empty files for TFs without a conserved
  # target. This avoids hiding zero-overlap TFs.
  all_tf_summary <- data.table::data.table(TF=sort(common3))
  for (s in stress_names) {
    x <- per_stress_commonTF_edges[[s]][,.N,by=TF]
    n <- x$N[match(all_tf_summary$TF,x$TF)]; n[is.na(n)] <- 0L
    all_tf_summary[,(paste0('N_',s,'_Targets')):=as.integer(n)]
  }
  per_tf_shared <- strict_common_pairs[,.N,by=TF]
  n_common <- per_tf_shared$N[match(all_tf_summary$TF,per_tf_shared$TF)]
  n_common[is.na(n_common)] <- 0L
  all_tf_summary[,N_Strict_Shared_Targets:=as.integer(n_common)]

  # Pool the predicted edges of shared hub TFs across the three stresses.
  all_common_tf_edges_long <- data.table::rbindlist(lapply(stress_names,function(s) {
    x <- per_stress_commonTF_edges[[s]]
    data.table::data.table(
      Stress=s,
      TF=x$TF, Target=x$Target)
  }), use.names=TRUE)
  all_common_tf_pairs <- unique(all_common_tf_edges_long[,.(TF,Target)],by=c('TF','Target'))
  data.table::setorder(all_common_tf_pairs,TF,Target)
  all_targets_union <- sort(unique(all_common_tf_pairs$Target))

  # Best choice for GO/KEGG: one Target ID per line, no duplicates.
  write_ids(all_targets_union,'ThreeStress_CommonHubTF_All_Targets_Union.txt')
  # Every distinct regulator-target association found in at least one stress.
  data.table::fwrite(all_common_tf_pairs,
                     file.path(export_dir,'ThreeStress_CommonHubTF_All_TF_Target_Pairs_Union.csv'))
  # Full per-stress provenance for every TF-Target edge.
  data.table::fwrite(all_common_tf_edges_long,
                     file.path(export_dir,'ThreeStress_CommonHubTF_All_TF_Target_Pairs_ByStress.csv'))

  pair_frequency <- all_common_tf_edges_long[,.(
    N_Stresses=data.table::uniqueN(Stress),
    Found_In=paste(sort(unique(Stress)),collapse='+')),by=.(TF,Target)]
  data.table::setorder(pair_frequency,TF,Target)
  data.table::fwrite(pair_frequency,
                     file.path(export_dir,'ThreeStress_CommonHubTF_All_Pairs_StressOccurrence.csv'))
  per_tf_all <- all_common_tf_pairs[,.(N_All_Unique_Targets=.N),by=TF]
  all_tf_summary[,N_All_Unique_Targets:=0L]
  all_tf_summary[per_tf_all,N_All_Unique_Targets:=i.N_All_Unique_Targets,on='TF']

  # Optional per-TF files: all targets across any of the three GRNs.
  per_tf_all_dir <- file.path(export_dir,'Per_Common_TF_All_Targets')
  dir.create(per_tf_all_dir,recursive=TRUE,showWarnings=FALSE)
  all_by_tf <- split(all_common_tf_pairs$Target, all_common_tf_pairs$TF)
  for (tf in sort(common3)) {
    safe_tf <- gsub('[^A-Za-z0-9_.-]','_',tf)
    writeLines(sort(unique(all_by_tf[[tf]])),
               file.path(per_tf_all_dir,paste0(safe_tf,'_All_TargetIDs.txt')),useBytes=TRUE)
  }
  message(sprintf('[ALL TARGETS] shared Hub TFs=%d; union Target IDs=%d; unique TF-Target pairs=%d',
                  length(common3),length(all_targets_union),nrow(all_common_tf_pairs)))
  for (tf in all_tf_summary$TF) {
    ids <- strict_common_pairs[TF==tf,Target]
    safe_tf <- gsub('[^A-Za-z0-9_.-]','_',tf)
    # Place inside per-TF folder rather than the export root.
    writeLines(sort(unique(ids)),
               file.path(per_tf_dir,paste0(safe_tf,'_Shared_TargetIDs.txt')),useBytes=TRUE)
  }
  data.table::fwrite(all_tf_summary,
                     file.path(export_dir,'ThreeStress_CommonHubTF_PerTF_Target_Counts.csv'))

  # Supplementary, RELAXED interpretation: a Target hit by ANY common Hub TF in
  # each stress. It may be regulated by DIFFERENT TFs in different stresses.
  pooled_targets <- lapply(per_stress_commonTF_edges,function(d) unique(d$Target))
  relaxed_common_targets <- Reduce(intersect,pooled_targets)
  write_ids(relaxed_common_targets,
            'ThreeStress_CommonHubTF_Pooled_TargetIntersection_DifferentTFsAllowed.txt')

  summary <- data.table::rbindlist(list(
    data.table::data.table(Item='Total_unique_hub_TFs',Count=length(union3)),
    data.table::data.table(Item='Hub_TFs_shared_all_three',Count=length(common3)),
    data.table::data.table(Item='Strict_shared_TF_Target_pairs',Count=nrow(strict_common_pairs)),
    data.table::data.table(Item='Distinct_Targets_in_strict_shared_pairs',
                           Count=data.table::uniqueN(strict_common_pairs$Target)),
    data.table::data.table(Item='Common_TFs_with_at_least_one_strict_target',
                           Count=sum(all_tf_summary$N_Strict_Shared_Targets>0L)),
    data.table::data.table(Item='Relaxed_pooled_Target_intersection',
                           Count=length(relaxed_common_targets)),
    data.table::data.table(Item='All_Targets_Union_Of_Three_Stress_GRNs',
                           Count=length(all_targets_union)),
    data.table::data.table(Item='All_Unique_TF_Target_Pairs_Union',
                           Count=nrow(all_common_tf_pairs))
  ))
  data.table::fwrite(summary,file.path(export_dir,'Export_Summary.csv'))
  writeLines(c(
    sprintf('Input: top %g%% of the ranked stress-specific XGBoost selected unions.',
            100*network_edge_fraction),
    sprintf('Hub definition: top %g%% of positive-outdegree TFs in EACH stress, ties included.',100*primary_fraction),
    'Stress-specific TF = present in the stress hub set, absent from both other hub sets.',
    'Shared at least two = present in the stress and >=1 other stress.',
    'Pair-only = present in exactly the two named stress hub sets; absent from the third.',
    'Three-stress common TF = present in all three stress hub sets.',
    'STRICT common target: identical (TF,Target) prediction present in all three networks; the TF itself must be a common hub.',
    'Per_Common_TF_Shared_Targets has one target-ID TXT for EACH common TF (empty if none).',
    'RELAXED pooled target intersection permits different TFs per stress, and is NOT an edge-conservation statistic.',
    'ALL TARGETS UNION: targets of any of the shared Hub TFs in at least one of the three final GRNs; IDs are deduplicated.',
    'ThreeStress_CommonHubTF_All_TF_Target_Pairs_Union.csv retains one row per distinct predicted TF-Target pair.',
    'ThreeStress_CommonHubTF_All_TF_Target_Pairs_ByStress.csv additionally preserves which stress contains each pair.',
    'Per_Common_TF_All_Targets contains one target-ID file for each shared Hub TF, pooling targets across the three GRNs.',
    'Lists are predicted regulatory associations; not experimentally confirmed edges.'),
    file.path(export_dir,'README_Definitions.txt'),useBytes=TRUE)
  message('[Hub TF and target export] ',normalizePath(export_dir,winslash='/',mustWork=TRUE))
  print(summary)
  invisible(list(directory=export_dir,summary=summary))
}
# Export immediately after hub identification, before time-consuming null sims.
export_result <- export_hub_tf_and_shared_targets(
  observed_hubs=observed_hubs, final_edges=final_edges,
  output_dir=output_dir, primary_fraction=primary_fraction)
# ========================= END EXTRA EXPORT ==================================

# ------ 3. MULTIVARIATE HYPERGEOMETRIC RANDOM EDGE-RETENTION NULL ------------
# For each TF i, Candidate_Degree_i specifies the number of eligible edges.
# We sample exactly N_Final of N_Candidate edges uniformly WITHOUT replacement.
# A sequential hypergeometric draw gives exactly the same TF degree distribution
# as explicit uniform edge sampling, while avoiding millions of sampled edge IDs.
rmvhyper_counts <- function(capacity, n_keep) {
  if (any(!is.finite(capacity)) || any(capacity < 0) ||
      any(abs(capacity-round(capacity)) > 1e-8)) stop('Invalid TF candidate outdegree.')
  total <- sum(capacity)
  if (n_keep < 0 || n_keep > total) stop('The retained edge count falls outside the candidate universe.')
  n <- length(capacity)
  if (!n) return(numeric())
  if (n_keep == total) return(as.numeric(capacity))
  if (n_keep == 0) return(rep(0,n))
  res <- numeric(n)
  remains_total <- total
  remains_keep <- n_keep
  if (n > 1L) for (i in seq_len(n-1L)) {
    m <- capacity[i]
    if (remains_keep == 0) {
      draw <- 0
    } else if (remains_keep == remains_total) {
      draw <- m
    } else if (m == 0) {
      draw <- 0
    } else {
      draw <- rhyper(1L, m=m, n=remains_total-m, k=remains_keep)
    }
    res[i] <- draw
    remains_total <- remains_total-m
    remains_keep <- remains_keep-draw
  }
  res[n] <- remains_keep
  if (any(res < 0) || any(res > capacity) || sum(res) != n_keep)
    stop('Multivariate hypergeometric draw violates its count constraints.')
  res
}

cat('\n[Edge-retention null] Replicates: ',n_random,'\n',sep='')
set.seed(random_seed)
null_mat <- matrix(NA_integer_,nrow=n_random,ncol=length(fractions),
                   dimnames=list(NULL,vapply(fractions, frac_key, character(1))))
for (b in seq_len(n_random)) {
  simulated <- lapply(caps, function(one) {
    deg <- rmvhyper_counts(one$Capacity,one$N_Final)
    h <- summarize_fractions(one$TF,deg,fractions)
    names(h) <- vapply(fractions,frac_key,character(1))
    h
  })
  for (j in seq_along(fractions)) {
    null_mat[b,j] <- length(common_hub_sets(simulated,frac_key(fractions[j])))
  }
  if (b %% progress_every == 0L || b == n_random)
    cat(sprintf(' Simulation %d / %d\n',b,n_random))
}

null_long <- rbindlist(lapply(seq_along(fractions), function(j)
  data.table(Iteration=seq_len(n_random),Threshold=fractions[j],
             Threshold_Label=frac_label(fractions[j]),
             Null_Common_Hubs=null_mat[,j])))
fwrite(null_long,file.path(output_dir,'03_Null_Iterations_AllThresholds.csv'))

null_stats <- rbindlist(lapply(seq_along(fractions), function(j) {
  k <- frac_key(fractions[j]); null <- null_mat[,j]
  obs <- length(observed_common[[k]])
  mu <- mean(null)
  q <- unname(quantile(null,probs=c(.025,.5,.975),type=1))
  data.table(Threshold=fractions[j],Threshold_Label=frac_label(fractions[j]),
             Observed_Common_Hubs=obs,Candidate_Common_Hubs=length(candidate_common[[k]]),
             Null_Mean=mu,Null_Median=q[2L],Null_L95=q[1L],Null_U95=q[3L],
             Null_Min=min(null),Null_Max=max(null),
             Enrichment_Fold=if(mu>0) obs/mu else NA_real_,
             Difference_From_Null_Mean=obs-mu,
             Exceed_Count=sum(null>=obs),
             Empirical_P_Upper=(1+sum(null>=obs))/(n_random+1),
             N_Random=n_random,
             Null_Model='Selected-union uniform edge retention; TF degrees random')
}))
null_stats[,Empirical_P_Holm_Across_Thresholds:=stats::p.adjust(Empirical_P_Upper,
                                                                method=threshold_p_adjust)]
null_stats[,Analysis_Role:=ifelse(abs(Threshold-primary_fraction)<1e-12,
                                  'Primary_Top10', 'Sensitivity_Exploratory')]
fwrite(null_stats,file.path(output_dir,'03_Null_Enrichment_Summary.csv'))
cat('\n[Edge-retention results]\n'); print(null_stats)

# -------- COMPLEMENT: TF-LABEL PERMUTATION, DEGREE DISTRIBUTION FIXED ----------
# This is NOT an exact-degree-preserving graph rewiring. It preserves the
# observed outdegree multiset within each stress (including TFs with zero
# final outdegree) but randomizes which candidate TF gets which degree.
# Consequently, the number of hubs per stress and ties are held fixed.
label_stats <- NULL
label_long <- NULL
if (run_tf_label_permutation) {
  cat('\n[TF-label permutation] Permuting TF identities within each stress\n')
  set.seed(random_seed+7219L)
  label_null <- matrix(NA_integer_,nrow=n_label_permutation,ncol=length(fractions))
  for (b in seq_len(n_label_permutation)) {
    sim_sets <- lapply(caps,function(cc) {
      permuted_degrees <- cc$Observed_Degree[sample.int(length(cc$Observed_Degree))]
      ans <- summarize_fractions(cc$TF,permuted_degrees,fractions)
      names(ans) <- vapply(fractions,frac_key,character(1))
      ans
    })
    for (j in seq_along(fractions))
      label_null[b,j] <- length(common_hub_sets(sim_sets,frac_key(fractions[j])))
    if (b %% progress_every==0L || b==n_label_permutation)
      cat(sprintf(' label perm %d / %d\n',b,n_label_permutation))
  }
  label_long <- rbindlist(lapply(seq_along(fractions),function(j)
    data.table(Iteration=seq_len(n_label_permutation),Threshold=fractions[j],
               Threshold_Label=frac_label(fractions[j]),Null_Common_Hubs=label_null[,j])))
  label_stats <- rbindlist(lapply(seq_along(fractions),function(j) {
    z <- label_null[,j]; key <- frac_key(fractions[j]); obs <- length(observed_common[[key]])
    data.table(Threshold=fractions[j],Threshold_Label=frac_label(fractions[j]),
               Observed_Common_Hubs=obs,Null_Mean=mean(z),
               Null_L95=as.numeric(quantile(z,.025,type=1)),
               Null_U95=as.numeric(quantile(z,.975,type=1)),
               Enrichment_Fold=if(mean(z)>0) obs/mean(z) else NA_real_,
               Difference_From_Null_Mean=obs-mean(z),
               Exceed_Count=sum(z>=obs),
               Empirical_P_Upper=(1+sum(z>=obs))/(n_label_permutation+1),
               N_Random=n_label_permutation,
               Null_Model='Observed outdegree multiset preserved; candidate TF labels permuted')
  }))
  label_stats[,Empirical_P_Holm_Across_Thresholds:=stats::p.adjust(
    Empirical_P_Upper,method=threshold_p_adjust)]
  fwrite(label_long,file.path(output_dir,'03b_TFLabelPermutation_Null_Iterations.csv'))
  fwrite(label_stats,file.path(output_dir,'03b_TFLabelPermutation_Null_Summary.csv'))
  print(label_stats)
}

# ------------------------------- 4. FIGURES ---------------------------------
# All data are taken from the existing analysis objects / output CSVs.
# Observed values are RED, null expectations and intervals BLUE.
# The P values are Holm-adjusted across the hub cutoffs *within each null family*.
# They are NOT P values for individual TFs.
cat('\n[Figures] Redraw with point-adjacent numeric labels and Holm P...\n')

OBS_RED <- '#BF3E46'
NULL_BLUE <- '#387DA9'
NULL_LIGHT <- '#A9CCDD'
HUB_GREEN <- '#D7E9DE'

plot_theme <- function() {
  ggplot2::theme_bw(base_size=11.5) +
    ggplot2::theme(
      text=ggplot2::element_text(color='black'),
      axis.title=ggplot2::element_text(color='black',size=11),
      axis.text=ggplot2::element_text(color='black',size=10),
      plot.title=ggplot2::element_text(face='bold',color='black',hjust=0.5,size=13),
      plot.subtitle=ggplot2::element_text(color='black',size=9.5),
      plot.caption=ggplot2::element_text(color='black',size=8,hjust=0),
      legend.text=ggplot2::element_text(color='black'),
      legend.position='bottom',
      strip.text=ggplot2::element_text(face='bold',color='black'),
      strip.background=ggplot2::element_rect(fill='#F3F5F6',color='#CED3D6'),
      panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(color='#E7E9EC',linewidth=0.36)
    )
}
plot_save <- function(plot,name,w=8.5,h=5.5) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(file.path(output_dir,paste0(name,'.pdf')),plot=plot,
                  width=w,height=h,bg='white',device=grDevices::pdf,useDingbats=FALSE)
  ggplot2::ggsave(file.path(output_dir,paste0(name,'.png')),plot=plot,
                  width=w,height=h,bg='white',dpi=400)
}
format_p <- function(p) {
  vapply(p,function(z) {
    if (!is.finite(z)) return('NA')
    if (z < 0.001) return('< 0.001')
    sprintf('= %.4f',z)
  },character(1))
}
threshold_levels <- vapply(fractions,frac_label,character(1))
primary_key <- frac_key(primary_fraction)
primary_label <- frac_label(primary_fraction)
primary_obs <- null_stats[abs(Threshold-primary_fraction)<1e-12]
main_sim <- null_long[abs(Threshold-primary_fraction)<1e-12]
if (nrow(primary_obs)!=1L || !nrow(main_sim))
  stop('Expected one result for the primary hub threshold.')

# ---- FIGURE 1: count histogram; each label is located by its OWN reference line.
hist1 <- main_sim[,.(Frequency=.N),by=Null_Common_Hubs]
setorder(hist1,Null_Common_Hubs)
max_f1 <- max(hist1$Frequency)
mean_x <- primary_obs$Null_Mean[[1L]]
obs_x <- primary_obs$Observed_Common_Hubs[[1L]]
fold_x <- primary_obs$Enrichment_Fold[[1L]]
p1_val <- primary_obs$Empirical_P_Holm_Across_Thresholds[[1L]]
x1_lim <- c(min(c(hist1$Null_Common_Hubs,mean_x,obs_x))-1.0,
            max(c(hist1$Null_Common_Hubs,mean_x,obs_x))+1.7)
# Blue mean above its dotted line; red observed and Holm P beside red dashed line.
p1 <- ggplot(hist1,aes(x=Null_Common_Hubs,y=Frequency))+
  geom_col(width=0.80,fill=NULL_LIGHT,color='white',linewidth=0.25)+
  geom_vline(xintercept=mean_x,linetype='dotted',color=NULL_BLUE,linewidth=0.85)+
  geom_vline(xintercept=obs_x,linetype='dashed',color=OBS_RED,linewidth=0.98)+
  annotate('text',x=mean_x,y=max_f1*1.37,
           label=sprintf('Null mean\n%.2f',mean_x),
           color=NULL_BLUE,fontface='bold',size=3.85,hjust=0.5,vjust=1,lineheight=0.95)+
  annotate('text',x=obs_x,y=max_f1*1.37,
           label=sprintf('Observed\n%d',obs_x),
           color=OBS_RED,fontface='bold',size=3.85,hjust=0.5,vjust=1,lineheight=0.95)+
  annotate('text',x=obs_x,y=max_f1*1.14,
           label=sprintf('Holm-adjusted P %s\nFold = %.2fx',format_p(p1_val),fold_x),
           color=OBS_RED,fontface='plain',size=3.3,hjust=0.5,vjust=1,lineheight=1.0)+
  coord_cartesian(xlim=x1_lim,ylim=c(0,max_f1*1.5),clip='off')+
  labs(title=paste0(primary_label,' common hub TFs: observed vs. candidate-edge null'),
       x='Three-stress shared hub TFs under the candidate-edge null',
       y='Simulation frequency',
       caption='Dashed red: observed; dotted blue: null mean. P adjusted across tested hub cutoffs (Holm).')+
  plot_theme()
plot_save(p1,'Fig1_Top10_Null_Distribution',w=9.7,h=6.0)

# ---- FIGURE 2: all labels anchored to the relevant triangle, circle, or CI.
fig2dat <- data.table::copy(null_stats)
fig2dat[,Threshold_Label:=factor(Threshold_Label,levels=threshold_levels)]
obs_scale <- max(fig2dat$Observed_Common_Hubs,fig2dat$Null_U95,na.rm=TRUE)
ygap <- max(2.5,obs_scale*0.033)
fig2dat[,`:=`(
  ObsY = Observed_Common_Hubs + ygap*0.65,
  PY = Observed_Common_Hubs + ygap*2.1,
  NullY = Null_L95 - ygap*0.75,
  Observed_Label=sprintf('%d',Observed_Common_Hubs),
  P_Label=paste0('Holm P ',format_p(Empirical_P_Holm_Across_Thresholds)),
  Null_Label=sprintf('%.2f\n[%.0f, %.0f]',Null_Mean,Null_L95,Null_U95)
)]
p2 <- ggplot(fig2dat,aes(x=Threshold_Label))+
  geom_rect(data=fig2dat[abs(Threshold-primary_fraction)<1e-12],
            aes(xmin=as.numeric(Threshold_Label)-0.49,
                xmax=as.numeric(Threshold_Label)+0.49,
                ymin=-Inf,ymax=Inf),
            inherit.aes=FALSE,fill='#EFF2F5',alpha=0.65)+
  geom_errorbar(aes(ymin=Null_L95,ymax=Null_U95),
                color=NULL_BLUE,width=0.10,linewidth=0.85)+
  geom_point(aes(y=Null_Mean),color=NULL_BLUE,size=3.7,shape=16)+
  geom_point(aes(y=Observed_Common_Hubs),color=OBS_RED,size=4.6,shape=17)+
  geom_text(aes(y=ObsY,label=Observed_Label),color=OBS_RED,
            size=4.05,fontface='bold',vjust=0)+
  geom_text(aes(y=PY,label=P_Label),color='black',fontface='italic',
            size=3.25,vjust=0)+
  geom_text(aes(y=NullY,label=Null_Label),color=NULL_BLUE,
            fontface='bold',size=3.3,vjust=1.02,lineheight=0.96)+
  scale_y_continuous(limits=c(0,max(fig2dat$PY,na.rm=TRUE)+ygap*2.3),
                     expand=expansion(mult=c(0.02,0.01)))+
  labs(title='Cross-stress common hub TFs across cutoffs',
       subtitle='Red: observed counts and Holm-adjusted P; blue: null means and 2.5–97.5% intervals',
       x='TF outdegree hub threshold',y='Number of shared hub TFs',
       caption='Grey shading denotes the primary Top 10% cutoff; intervals describe null simulations, not replicate CIs.')+
  plot_theme()+
  theme(panel.grid.major.x=element_blank())
plot_save(p2,'Fig2_Observed_vs_Null_ByThreshold',w=9.5,h=6.6)

# ---- FIGURE 3: directly label bar (fold), and inside bar (obs & null).
plot_fold <- data.table::copy(null_stats)
plot_fold[,Threshold_Label:=factor(Threshold_Label,levels=threshold_levels)]
plot_fold[,`:=`(
  Fold_Label=sprintf('%.2fx',Enrichment_Fold),
  Obs_Label=sprintf('Observed: %d',Observed_Common_Hubs),
  Null_Label=sprintf('Null: %.2f',Null_Mean),
  P_Label=paste0('Holm P ',format_p(Empirical_P_Holm_Across_Thresholds))
)]
plot_fold[,`:=`(FoldY=Enrichment_Fold+0.065,PY=Enrichment_Fold+0.20)]
p3 <- ggplot(plot_fold,aes(x=Threshold_Label,y=Enrichment_Fold))+
  geom_hline(yintercept=1,linetype='dashed',color='#6B737A',linewidth=0.65)+
  geom_col(width=0.62,fill=HUB_GREEN,color='#B1C6BA',linewidth=0.45)+
  geom_text(aes(y=FoldY,label=Fold_Label),color=OBS_RED,
            fontface='bold',size=4.2,vjust=0)+
  geom_text(aes(y=PY,label=P_Label),color='black',
            fontface='italic',size=3.45,vjust=0)+
  geom_text(aes(y=0.65,label=Obs_Label),color=OBS_RED,
            fontface='bold',size=3.6)+
  geom_text(aes(y=0.47,label=Null_Label),color=NULL_BLUE,
            fontface='bold',size=3.6)+
  scale_y_continuous(limits=c(0,max(plot_fold$PY,na.rm=TRUE)+0.16),
                     expand=expansion(mult=c(0,0)))+
  labs(title='Shared hub TF enrichment relative to candidate-edge null',
       subtitle='Effect size = observed shared hubs / mean(null)',
       x='TF outdegree hub threshold',y='Observed / null mean',
       caption='P values are Holm-adjusted across the three tested cutoffs; no per-TF P values are implied.')+
  plot_theme()+theme(panel.grid.major.x=element_blank())
plot_save(p3,'Fig3_Shared_Hub_Enrichment_Fold',w=9.0,h=6.0)

hub_counts_table <- rbindlist(lapply(stress_names, function(stress) {
  data.table(Stress=stress,
             Group=c('Final XGBoost GRN', 'Selected-method candidate union'),
             Hubs=c(length(observed_hubs[[stress]][[primary_key]]),
                    length(candidate_hubs[[stress]][[primary_key]])))
}))
point_table <- rbindlist(lapply(stress_names, function(stress) {
  one <- caps[[stress]]
  data.table(Stress=stress, TF=one$TF, Candidate_Degree=one$Capacity,
             Final_Degree=one$Observed_Degree)
}))
fwrite(hub_counts_table, file.path(output_dir, '02_HubCounts_Candidate_vs_Final.csv'))
fwrite(point_table, file.path(output_dir, '04_Candidate_Final_TF_Outdegree_Data.csv'))

# ---- FIGURE 4: direct red values on final bars, blue values on candidate bars.
hub_comp <- data.table::copy(hub_counts_table)
hub_comp[,Stress:=factor(Stress,levels=names(network_files))]
group_levels <- c('Final XGBoost GRN','Selected-method candidate union')
hub_comp[,Group:=factor(Group,levels=group_levels)]
p4 <- ggplot(hub_comp,aes(x=Stress,y=Hubs,fill=Group))+
  geom_col(position=position_dodge(width=0.75),width=0.66)+
  geom_text(aes(label=Hubs,color=Group),position=position_dodge(width=0.75),
            vjust=-0.45,fontface='bold',size=4.2)+
  scale_fill_manual(values=c('Final XGBoost GRN'='#CF876B',
                             'Selected-method candidate union'='#8CAAC5'),drop=FALSE)+
  scale_color_manual(values=c('Final XGBoost GRN'=OBS_RED,
                              'Selected-method candidate union'=NULL_BLUE),guide='none')+
  scale_y_continuous(expand=expansion(mult=c(0,0.16)))+
  labs(title=paste0('Number of ',primary_label,' hub TFs by stress'),
       x=NULL,y='Hub TF count',fill=NULL)+
  plot_theme()
plot_save(p4,'Fig4_Top10_Candidate_vs_Final_HubCounts',w=9.0,h=5.6)

# ---- FIGURE 5: summarize scatter in-panel, without labeling individual TFs.
point_df <- data.table::copy(point_table)
cor_df <- point_df[,.(
  Spearman_Rho=if(uniqueN(Candidate_Degree)>1L && uniqueN(Final_Degree)>1L)
    suppressWarnings(stats::cor(Candidate_Degree,Final_Degree,method='spearman'))
  else NA_real_,
  N_TFs=.N,
  X_Annot=min(Candidate_Degree+1)*1.12,
  Y_Annot=max(Final_Degree+1)/1.30),by=Stress]
point_df[,Stress:=factor(Stress,levels=names(network_files))]
cor_df[,Stress:=factor(Stress,levels=names(network_files))]
cor_df[,Annot:=sprintf('Spearman rho = %.3f\nN(TFs) = %d',Spearman_Rho,N_TFs)]
p5 <- ggplot(point_df,aes(x=Candidate_Degree+1,y=Final_Degree+1))+
  geom_point(alpha=0.4,size=0.75,color=NULL_BLUE)+
  geom_label(data=cor_df,aes(x=X_Annot,y=Y_Annot,label=Annot),
             inherit.aes=FALSE,hjust=0,vjust=0,fill='white',
             color='black',size=3.05,label.size=0.18)+
  scale_x_log10()+scale_y_log10()+facet_wrap(~Stress,nrow=1)+
  labs(title='Candidate versus final TF outdegree',
       x='Candidate-union TF outdegree + 1 (log10)',
       y='Final XGBoost TF outdegree + 1 (log10)')+
  plot_theme()+theme(legend.position='none')
plot_save(p5,'Fig5_Candidate_vs_Final_TF_Outdegree',w=11.7,h=4.8)

# ---- FIGURE 5b: visibly distinguish near-zero random mean from far-right obs.
if (run_tf_label_permutation && !is.null(label_stats) && nrow(label_stats)) {
  lp <- label_stats[abs(Threshold-primary_fraction)<1e-12]
  lps <- label_long[abs(Threshold-primary_fraction)<1e-12]
  if (nrow(lp)!=1L || !nrow(lps)) stop('Missing TF-label permutation results for the primary hub threshold.')
  h2 <- lps[,.(Frequency=.N),by=Null_Common_Hubs]
  n2max <- max(h2$Frequency)
  obs2 <- lp$Observed_Common_Hubs[[1L]]
  null2 <- lp$Null_Mean[[1L]]
  p2_adj <- lp$Empirical_P_Holm_Across_Thresholds[[1L]]
  fe2 <- lp$Enrichment_Fold[[1L]]
  p5b <- ggplot(h2,aes(x=Null_Common_Hubs,y=Frequency))+
    geom_col(width=0.78,fill=NULL_LIGHT,color='white')+
    geom_vline(xintercept=null2,color=NULL_BLUE,linetype='dotted',linewidth=0.85)+
    geom_vline(xintercept=obs2,color=OBS_RED,linetype='dashed',linewidth=1.0)+
    annotate('text',x=max(2,null2+1),y=n2max*1.34,
             label=sprintf('Null mean = %.2f',null2),hjust=0,
             color=NULL_BLUE,fontface='bold',size=3.9)+
    annotate('text',x=obs2-.6,y=n2max*1.32,
             label=sprintf('Observed = %d',obs2),hjust=1,
             color=OBS_RED,fontface='bold',size=4.0)+
    annotate('text',x=obs2-.6,y=n2max*1.14,
             label=sprintf('Holm-adjusted P %s\nFold = %.2fx',format_p(p2_adj),fe2),
             hjust=1,color=OBS_RED,size=3.4,lineheight=1.0)+
    coord_cartesian(xlim=c(-0.8,max(obs2,max(h2$Null_Common_Hubs))+1.5),
                    ylim=c(0,n2max*1.48),clip='off')+
    labs(title='Three-stress shared hub TFs under TF-label permutation',
         x='Number of shared hub TFs under permuted TF labels',
         y='Simulation frequency',
         caption='TF-label permutation preserves the outdegree multiset, not each TF degree or the TF–target edge structure.')+
    plot_theme()
  plot_save(p5b,'Fig5b_TFLabelPermutation_Top10_SharedHub_Null',w=9.7,h=5.8)
}

# ---- FIGURE 6: publication palette; all seven Venn-region counts retained.
# Vector PDF is generated DIRECTLY by grid (not converted from PNG).
if (requireNamespace('VennDiagram',quietly=TRUE)) {
  sets <- lapply(observed_hubs,function(x) x[[primary_key]])
  if (!identical(sort(names(sets)),sort(stress_names)))
    warning('The Venn diagram expects one hub set for each stress.')
  venn_obj <- VennDiagram::venn.diagram(
    x=sets,filename=NULL,category.names=names(sets),
    fill=c('#387DA9','#DC795D','#43A18D'),alpha=0.57,
    col=c('#2A658D','#BD634D','#317E6B'),lwd=1.0,
    cex=1.28,cat.cex=1.20,
    fontface='plain',cat.fontface='bold',
    fontfamily='sans',cat.fontfamily='sans',
    cat.col=c('black','black','black'),margin=0.07)
  grDevices::pdf(file.path(output_dir,'Fig6_ThreeStress_Top10_Hub_Venn.pdf'),
                 width=7.33,height=7.33,paper='special',onefile=FALSE,
                 useDingbats=FALSE,bg='white')
  grid::grid.newpage();grid::grid.draw(venn_obj);grDevices::dev.off()
  grDevices::png(file.path(output_dir,'Fig6_ThreeStress_Top10_Hub_Venn.png'),
                 width=2600,height=2600,res=400,bg='white')
  grid::grid.newpage();grid::grid.draw(venn_obj);grDevices::dev.off()
} else {
  warning('Install VennDiagram to generate the hub overlap diagram (Fig. 6).')
}
message('[Figures] Saved to: ',normalizePath(output_dir,winslash='/'))

# ------------------------- 5. OPTIONAL DEGREE-PRESERVING NULL ---------------
if (run_degree_preserving) {
  if (!requireNamespace('igraph',quietly=TRUE))
    stop('Install igraph to run degree-preserving rewiring.')
  if (n_rewire < 19L) stop('n_rewire must be at least 19.')
  common_tf <- observed_common[[primary_key]]
  if (!length(common_tf)) stop('The rewiring test requires at least one hub shared by all three stresses.')
  cat('\n[Degree-preserving rewiring] Scope: ',rewire_scope,', B=',n_rewire,'\n',sep='')
  # Restrict only if user opts into the faster, explicitly LOCAL scope.
  # When scope='full', the entire final graph's in/out degree is preserved.
  edges_rewire <- lapply(final_edges,function(e) {
    if (rewire_scope == 'common_hubs') e[TF %in% common_tf] else e
  })
  gr <- lapply(edges_rewire,function(e)
    igraph::graph_from_data_frame(e,directed=TRUE))
  del <- lapply(gr, function(g) list(out=igraph::degree(g,mode='out'),
                                     indegree=igraph::degree(g,mode='in')))
  extract_pairs <- function(g, relevant_tfs) {
    e <- igraph::as_data_frame(g,what='edges')
    if (!nrow(e)) return(character())
    e <- e[e$from %in% relevant_tfs,,drop=FALSE]
    unique(paste(e$from,e$to,sep='\t'))
  }
  obs_edge_overlap <- length(Reduce(intersect,lapply(gr,extract_pairs,
                                                     relevant_tfs=common_tf)))
  null_edges <- numeric(n_rewire)
  set.seed(random_seed+111L)
  for (b in seq_len(n_rewire)) {
    randomized <- lapply(seq_along(gr),function(i) {
      g <- gr[[i]]
      n_trials <- max(100, ceiling(igraph::ecount(g)*rewire_trials_per_edge))
      igraph::rewire(g,with=igraph::keeping_degseq(loops=FALSE,niter=n_trials))
    })
    for (i in seq_along(randomized)) {
      if (!identical(as.numeric(igraph::degree(randomized[[i]],mode='out')),
                     as.numeric(del[[i]]$out)) ||
          !identical(as.numeric(igraph::degree(randomized[[i]],mode='in')),
                     as.numeric(del[[i]]$indegree)))
        stop('Rewiring changed the node indegree or outdegree.')
    }
    null_edges[b] <- length(Reduce(intersect,lapply(randomized,extract_pairs,
                                                    relevant_tfs=common_tf)))
    if (b %% progress_every == 0L || b == n_rewire)
      cat(sprintf(' Rewiring replicate %d / %d\n',b,n_rewire))
  }
  rewire_summary <- data.table(
    Scope=rewire_scope,Number_Of_Common_Hubs=length(common_tf),
    Observed_ThreeStress_SamePairs=obs_edge_overlap,
    Null_Mean=mean(null_edges),Null_L95=unname(quantile(null_edges,.025,type=1)),
    Null_U95=unname(quantile(null_edges,.975,type=1)),
    Fold=if(mean(null_edges)>0)obs_edge_overlap/mean(null_edges) else NA_real_,
    P_Upper=(1+sum(null_edges>=obs_edge_overlap))/(n_rewire+1),
    N_Rewire=n_rewire,Trials_Per_Edge=rewire_trials_per_edge,
    Degree_Preservation=if (rewire_scope=='full')
      'Whole GRN in/out degree preserved' else
        'Common-hub edge subgraph in/out degree preserved ONLY')
  fwrite(rewire_summary,file.path(output_dir,'05_DegreePreserving_Overlap_Summary.csv'))
  fwrite(data.table(Iteration=seq_len(n_rewire),Null_ThreeStress_SamePairs=null_edges),
         file.path(output_dir,'05_DegreePreserving_Overlap_Null.csv'))
  plot_dp <- ggplot(data.table(Null_Pairs=null_edges),aes(x=Null_Pairs)) +
    geom_histogram(bins=30,fill='#6B96BA',color='white') +
    geom_vline(xintercept=obs_edge_overlap,color='#B44142',linetype='dashed',linewidth=1) +
    annotate('label',x=mean(null_edges),y=Inf,
             label=sprintf('Observed = %d\nNull mean = %.2f\nFold = %.2fx\nP = %.4f',
                           obs_edge_overlap,mean(null_edges),rewire_summary$Fold,rewire_summary$P_Upper),
             hjust=0,vjust=1.2,size=3.2,color='black',fill='white') +
    labs(title='Shared TF-target pairs under degree-preserving rewiring',
         subtitle=sprintf('Scope=%s | Observed=%d | Null mean=%.1f | P=%.4f',
                          rewire_scope,obs_edge_overlap,mean(null_edges),rewire_summary$P_Upper),
         x='Three-stress identical TF-target pairs',y='Simulation frequency') + plot_theme()
  plot_save(plot_dp,'Fig7_DegreePreserving_SamePairs_Null',w=9,h=5.8)
  print(rewire_summary)
}
