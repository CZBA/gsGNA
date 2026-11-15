# gsGNA ensemble comparison for rice abiotic stress GRNs
#
# Run from the repository root: Rscript code/gsGNA_Ensemble_Comparison.R
# Install the R packages listed in 'required' before running.
#
# Input: <method>/<stress>_top10_edges.csv (TF, Target, Weight in columns 1-3)
# and chip_seq.csv (TF and Target in columns 1-2). Set GSGNA_INPUT_ROOT to
# use an input directory other than the current working directory.
# If chip_seq.csv is absent from the input root, use data/chip_seq.csv
#
# Six fusion strategies each evaluate all 57 subsets of two to six algorithms.
# TF-wise development CV selects one fusion strategy across the three stresses
# and a base-algorithm subset per stress. Nested CV evaluates this selection
# procedure on the development TFs. TF-disjoint testing follows model lock.
#
# Models are fitted on their selected candidate unions. Comparisons use the
# fixed six-method Top 10% union, with score zero for edges outside a subset's
# union. ChIP-supported edges are positives and the rest are unverified.
# Final networks include all edges in the selected union. The scores rank
# candidate interactions and do not represent P values or FDR estimates.

required <- c("data.table", "pROC", "PRROC", "xgboost", "ada", "rpart", "glmnet")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) {
  stop("Install the required R packages: ", paste(missing_pkgs, collapse=", "))
}
suppressPackageStartupMessages(library(data.table))

# =============================== CONFIG =====================================
seed <- 123L
input_root <- Sys.getenv("GSGNA_INPUT_ROOT", unset=".")
# This stress order determines which random seed is assigned to each analysis.
stress_names <- c("Drought", "Alkaline", "Cold")
method_names <- c("GENIE3", "Kboost", "GRNBoost2", "3DCEMA", "DeepRIG", "IGEGRNs")
# The input threshold is fixed before development model selection.
threshold_stages <- "top10"
reference_stage <- "top10"    # fixed six-method Top10 union
input_stage <- threshold_stages[1L] # changed only by orchestrator
threshold_sweep_development <- TRUE
# Nested CV evaluates the FULL inner selection process on development TFs ONLY.
nested_cv_enabled <- TRUE
nested_outer_folds <- 5L  # existing development fivefold assignment
nested_inner_folds <- 4L
nested_bootstrap_reps <- 1000L
threshold_study_root <- file.path(input_root, "gsGNA_FixedTop10")
# Selection MUST use a shared evaluation universe for meaningful AUROC/AUPR comparisons.
# The alternative is supplied for descriptive/exploratory within-combination ranking only.
selection_universe <- "Common_All6"  # Fixed: model selection always uses a common evaluation universe
# Each combination is ALWAYS FIT on its own input-method union, independent of the above.
# CV/test report BOTH common-universe MAIN metrics and combination-union
# descriptive metrics. Only COMMON AUROC/AUPR are used for model selection.
# Each subset's feature training, candidate scoring and final full selected-union export use
# that subset's own union, not the six-method union.
chip_path <- file.path(input_root, "chip_seq.csv")
repo_chip_path <- file.path("data", "chip_seq.csv")
if (!file.exists(chip_path) && file.exists(repo_chip_path)) {
  chip_path <- repo_chip_path
}
chip_tf_col <- 1L      # TF column in chip_seq.csv
chip_target_col <- 2L  # Target column in chip_seq.csv

# Keep TFs disjoint and stratify by the number of ChIP-positive edges.
split_modes <- "stratified_random"
root_base <- file.path(threshold_study_root,paste0("Stage_",input_stage))
# Resume completed combinations only when the input and settings match.
resume_existing <- TRUE
# Use a fresh output directory when changing inputs or settings.
allow_overwrite <- FALSE

run_mode <- "FULL"            # Evaluate the complete subset search.
include_singletons <- FALSE   # Evaluate single algorithms as separate baselines.
only_all6 <- FALSE            # TRUE limits each strategy to the six-method subset.

# Evaluate held-out TFs after locking the development-selected models.
evaluate_locked_test <- FALSE # forced: sweep development phase; final test is a separate locked function
outer_test_fraction <- 0.20
outer_min_test_tfs <- 3L
cv_folds <- 5L
train_unlabeled_ratio <- 1.0
cv_unlabeled_ratio <- Inf    # Inf retains all unverified candidates. A finite ratio samples them.
outer_unlabeled_ratio <- Inf
# Select models using AUROC pooled across all five out-of-fold predictions.
selection_metric <- "Pooled_AUROC" # fixed consistently for both subset and overall strategy selection
# Shared-universe evaluation (CoReGRN-style zero padding):
# Unpredicted edges receive score 0, including on the common all-six background.
# Genuine zero-scored nominated edges tie with unpredicted edges. No random ranking.
# This constant is NEVER used in model fitting, own-union ranking or network export.
outside_union_eval_score <- 0.0
normalize_weights <- TRUE
weight_modes <- setNames(rep("nonnegative",length(method_names)), method_names)
weight_modes["3DCEMA"] <- "absolute" # Rank 3DCEMA edges by absolute weight.
score_directions <- setNames(rep("higher_better",length(method_names)),method_names)
score_directions["Kboost"] <- "higher_better"

# Bagging resamples algorithm columns with replacement within each subset,
# then averages edge scores across bootstrap replicates.
bagging_run_mode <- "publication"  # debug=25, mid=200, final=500, publication=2000
bagging_B <- switch(bagging_run_mode, debug=25L, mid=200L, final=500L, publication=2000L,
                    stop("bagging_run_mode must be debug, mid, final, or publication"))

xgb_rounds <- 150L
xgb_depth <- 3L
xgb_eta <- 0.10
xgb_nthread <- 4L
ada_iter <- 100L
ada_nu <- 0.5
ridge_lambda <- 0.01
bootstrap_tf_reps <- 2000L
# Repeat TF-wise CV on the fixed development set to assess stability.
reviewer_cv_repeats <- 10L
reviewer_ci_reps <- 2000L
reviewer_repeat_enabled <- TRUE
# Interpret held-out results as exploratory if they informed model selection.
heldout_evidence_status <- "Exploratory_if_previously_used_for_model_selection"
# Compare XGBoost against each other strategy using paired TF-level metrics.
superiority_reference <- "XGBoost"
superiority_alpha <- 0.05
superiority_permutation_reps <- 19999L  # Monte Carlo replicates for the one-sided sign-flip test.
superiority_bootstrap_reps <- 2000L
superiority_adjustment <- "holm"    # Five comparisons per stress and metric.
# Optional held-out benchmark of all six development-locked strategies.
run_locked_six_comparison <- FALSE # final evaluation is selected strategy only
# Optional: independent TF-regulation / other experimentally supported edges,
# not used by any model fit, tuning or final test. NULL means no such validation.
external_validation_path <- NULL
external_tf_col <- 1L
external_target_col <- 2L
bootstrap_min_tfs <- 3L
network_top_pct <- 1.0 # export all edges in selected union
predict_batch_size <- 200000L
save_test_predictions <- FALSE  # TRUE writes edge-level test predictions to disk.
save_winner_networks <- TRUE    # only the FINAL cross-stage locked model is exported
# The exported TopRank edges are not statistically significant unless an
# independent permutation/null-model FDR analysis is performed separately.
# ============================================================================

valid_modes <- "stratified_random"
if (!length(split_modes) || anyNA(split_modes) || anyDuplicated(split_modes) ||
    any(!split_modes %in% valid_modes)) {
  stop("The split mode must be stratified_random: ", paste(valid_modes,collapse=", "))
}
stopifnot(identical(threshold_stages,"top10"),
          reference_stage == "top10",
          input_stage %in% threshold_stages,
          identical(selection_universe,"Common_All6"),
          identical(run_mode,"FULL"),
          identical(selection_metric,"Pooled_AUROC"), # overall selection below also uses AUROC
          nested_outer_folds==cv_folds,nested_inner_folds>=3L,
          nested_inner_folds < length(method_names)+cv_folds,
          nested_bootstrap_reps>=100L,
          is.finite(outside_union_eval_score), outside_union_eval_score == 0,
          outer_test_fraction > 0, outer_test_fraction < 1,
          outer_min_test_tfs >= 3L, cv_folds >= 3L,
          train_unlabeled_ratio > 0, cv_unlabeled_ratio > 0, outer_unlabeled_ratio > 0,
          identical(names(weight_modes),method_names),
          all(weight_modes %in% c("nonnegative","signed","absolute")),
          identical(names(score_directions),method_names),
          all(score_directions %in% c("higher_better","lower_better")),
          ridge_lambda > 0, network_top_pct > 0, network_top_pct <= 1,
          reviewer_cv_repeats >= 2L, reviewer_ci_reps >= 1001L,
          superiority_permutation_reps >= 999L, superiority_bootstrap_reps >= 1001L,
          superiority_adjustment %in% c("holm","BH"),
          bagging_run_mode %in% c("debug","mid","final","publication"),
          bagging_B >= 1L, bagging_B == as.integer(bagging_B))

# Each split mode writes to its own output directory.
root_out <- NA_character_
split_mode <- NA_character_

# --------------------------- DATA INTEGRATION -------------------------------
canon <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Robust ChIP column resolution:
# (1) honor explicitly specified names or 1-based indices;
# (2) match role-specific names, not the ambiguous name "Gene" as Target;
# (3) if there are exactly two columns, infer the other column after finding one;
# (4) never silently select the same column for TF and Target.
resolve_chip_columns <- function(z, tf_spec=chip_tf_col, target_spec=chip_target_col) {
  hdr <- names(z)
  keys <- canon(hdr)
  n <- ncol(z)

  specified_index <- function(x, role) {
    if (is.null(x)) return(NA_integer_)
    if (length(x) != 1L || is.na(x)) stop("ChIP ", role, " column must be specified by one valid name or index")
    if (is.numeric(x)) {
      if (!is.finite(x) || x != as.integer(x) || x < 1 || x > n)
        stop("ChIP ", role, " column index is outside 1..", n, ": ", x)
      return(as.integer(x))
    }
    if (!is.character(x) || !nzchar(x)) stop("ChIP ", role, " column name is invalid")
    hits <- which(hdr == x)
    if (!length(hits)) hits <- which(keys == canon(x))
    if (length(hits) != 1L)
      stop("ChIP ", role, " column name does not match exactly one column: ", x,
           ". Available columns: ", paste(hdr, collapse=" | "))
    hits
  }
  detect <- function(exact, prefixes, role) {
    matches <- which(keys %in% exact)
    if (!length(matches)) {
      pattern <- paste0("^(?:", paste(prefixes, collapse="|"), ")")
      matches <- grep(pattern, keys, perl=TRUE)
    }
    if (length(matches) > 1L)
      stop("ChIP ", role, " matches multiple columns: ", paste(hdr[matches], collapse=", "),
           ". Set chip_", if (role=="TF") "tf" else "target", "_col")
    if (length(matches)) as.integer(matches) else NA_integer_
  }

  i <- specified_index(tf_spec, "TF")
  j <- specified_index(target_spec, "Target")
  if (is.na(i)) {
    i <- detect(c("tf", "tfid", "tfname", "tfgene", "tfgeneid",
                  "regulator", "regulatorid", "regulatorgene",
                  "source", "sourcegene", "sourceid", "transcriptionfactor",
                  "transcriptionfactorid"),
                c("tf", "regulator", "sourcegene", "transcriptionfactor"), "TF")
  }
  if (is.na(j)) {
    j <- detect(c("target", "targetid", "targetgene", "targetgeneid",
                  "targetname", "targetlocus", "targetlocusid",
                  "regulatedgene", "downstreamgene"),
                c("target", "regulatedgene", "downstreamgene"), "Target")
  }

  # Bare 'Gene'/'GeneID' is not intrinsically a target: use it only after TF
  # has been unambiguously assigned and there is just one remaining gene column.
  if (is.na(j) && !is.na(i) && n > 2L) {
    weak <- which(keys %in% c("gene", "geneid", "genename", "genesymbol", "locus", "locusid"))
    weak <- setdiff(weak, i)
    if (length(weak) == 1L) j <- as.integer(weak)
  }

  if (n == 2L) {
    if (is.na(i) && is.na(j)) {
      stop("ChIP headers are unrecognized: ", paste(hdr, collapse=" | "),
           ". Set chip_tf_col and chip_target_col to verified names or indices.")
    } else if (is.na(i)) {
      i <- setdiff(1:2, j)
      message("Using the other ChIP column for TF: ", hdr[i])
    } else if (is.na(j)) {
      j <- setdiff(1:2, i)
      message("Using the other ChIP column for Target: ", hdr[j])
    }
  }

  if (is.na(i) || is.na(j)) {
    stop("Could not identify the ChIP columns. Set chip_tf_col and chip_target_col in the configuration.",
         "\nAvailable columns: ", paste(sprintf("[%d] %s", seq_along(hdr), hdr), collapse=" | "))
  }
  if (i == j) stop("ChIP TF and Target resolve to the same column [", i, "] ", hdr[i],
                   ". Check chip_tf_col and chip_target_col.")
  message(sprintf("ChIP column mapping: TF=[%d] %s | Target=[%d] %s",
                  i, hdr[i], j, hdr[j]))
  c(TF=as.integer(i), Target=as.integer(j))
}

read_chip <- function(path) {
  if (!file.exists(path)) stop("Missing ChIP file: ", path)
  z <- fread(path, showProgress=FALSE, check.names=FALSE)
  if (ncol(z) < 2L) stop("The ChIP file must contain at least two columns: TF and Target.")
  loc <- resolve_chip_columns(z)
  i <- unname(loc["TF"]); j <- unname(loc["Target"])
  out <- unique(data.table(TF=toupper(trimws(as.character(z[[i]]))),
                           Target=toupper(trimws(as.character(z[[j]])))))
  out <- out[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) & TF!=Target]
  if (!nrow(out)) stop("No valid ChIP TF-target pairs after filtering. Check column order and gene IDs.")
  setkey(out, TF, Target)
  cat("ChIP unique positive pairs: ", nrow(out), "\n", sep="")
  cat("ChIP preview (first 3 pairs):\n")
  print(utils::head(out, 3L))
  out
}
get_grn_file <- function(ds, method, stage=input_stage) {
  folder_candidates <- switch(method,
                              GENIE3="GENIE3", Kboost="Kboost", GRNBoost2=c("GRNBoost2", "GRNBoost"),
                              `3DCEMA`="3DCEMA", DeepRIG="DeepRIG", IGEGRNs="IGEGRNs")
  # Prefer the current labels across all method directories, then legacy filenames.
  legacy_label <- switch(ds, Drought="Aridity", Alkaline="Alkalinity", NULL)
  if (stage %in% c("top5", "top10", "top20")) {
    suffix <- paste0("_", stage, "_edges.csv")
  } else if (stage == "significant") {
    suffix <- if (method == "GENIE3") "_significant_edges.txt" else "_significant_edges.csv"
  } else {
    stop("input_stage must be top5, top10, top20, or significant.")
  }
  filenames <- paste0(c(ds, legacy_label), suffix)
  opts <- unlist(lapply(filenames, function(filename) {
    file.path(input_root, folder_candidates, filename)
  }), use.names=FALSE)
  matched <- opts[file.exists(opts)]
  if (!length(matched)) {
    stop("Missing ", method, " input for ", ds, ". Checked: ", paste(opts, collapse="; "))
  }
  matched[1L]
}
read_grn <- function(ds, method) {
  p <- get_grn_file(ds, method)
  z <- fread(p, showProgress=FALSE, check.names=FALSE)
  if (ncol(z) < 3L) stop("The GRN file must contain TF, Target, and Weight columns: ", p)
  message("GRN column mapping: ", p, " | TF=", names(z)[1L],
          " | Target=", names(z)[2L], " | Weight=", names(z)[3L])
  dt <- data.table(TF=toupper(trimws(as.character(z[[1L]]))),
                   Target=toupper(trimws(as.character(z[[2L]]))),
                   Weight=suppressWarnings(as.numeric(z[[3L]])))
  dt <- dt[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) &
             TF!=Target & is.finite(Weight)]
  if (!nrow(dt)) stop("No valid GRN edges after filtering: ", p)
  mode <- unname(weight_modes[[method]])
  if (is.null(mode) || !mode %in% c("nonnegative", "signed", "absolute"))
    stop("Invalid or missing weight mode for: ", method)
  direction <- unname(score_directions[[method]])
  if (is.null(direction) || !direction %in% c("higher_better", "lower_better"))
    stop("Invalid or missing score direction for: ", method)
  n_neg <- sum(dt$Weight < 0)
  raw_min <- min(dt$Weight)
  raw_max <- max(dt$Weight)
  raw_med <- stats::median(dt$Weight)
  diag <- data.table(
    Dataset=ds, Method=method, Filename=p, Weight_Mode=mode,
    Score_Direction=direction,
    Input_TF_Column=names(z)[1L], Input_Target_Column=names(z)[2L],
    Input_Weight_Column=names(z)[3L], N_Input_Rows=nrow(z),
    N_Discarded_Rows=nrow(z)-nrow(dt), N_Valid_Rows=nrow(dt),
    N_Duplicate_Rows=nrow(dt)-uniqueN(dt, by=c("TF", "Target")), N_Negative=n_neg,
    N_Zero=sum(dt$Weight == 0), N_Positive=sum(dt$Weight > 0),
    Fraction_Negative=n_neg / nrow(dt),
    Raw_Min=raw_min, Raw_Median=raw_med, Raw_Max=raw_max
  )
  cat(sprintf(" %-10s weights: min=%.10e, median=%.10e, max=%.10e; negative=%s/%s (%.2f%%); mode=%s; direction=%s\n",
              method, raw_min, raw_med, raw_max,
              format(n_neg, big.mark=","), format(nrow(dt), big.mark=","),
              100*n_neg/nrow(dt), mode, direction))
  # Tiny Kboost weights are normal: scientific notation prevents misleading 0.0000 displays.
  if (method == "Kboost") {
    cat(sprintf(" Kboost diagnostic: strictly-positive=%d, exact-zero=%d, unique-scores=%d\n",
                sum(dt$Weight > 0), sum(dt$Weight == 0),
                data.table::uniqueN(dt$Weight)))
    if (data.table::uniqueN(dt$Weight) < 2L)
      warning("Kboost scores are constant; check for lost precision. Only membership information remains.")
  }
  if (mode == "nonnegative" && n_neg > 0L) {
    stop("Negative weights in ", p, " with weight mode ", mode,
         ". Check the original score definition and set weight_modes['", method,
         "'] to 'signed' to retain signs or 'absolute' to rank by magnitude.")
  }
  if (mode == "absolute") dt[, Weight := abs(Weight)]
  # All six methods are configured with higher_better, including Kboost.
  # No value reversal is applied to Kboost despite its small magnitude.
  # This branch is retained for future methods with genuinely inverse scores.
  flip_max <- if (direction == "lower_better") max(dt$Weight) else NA_real_
  if (direction == "lower_better") {
    dt[, Weight := flip_max - Weight]
    message(sprintf(" %s score direction reversed, max=%.10e", method, flip_max))
  }
  diag[, `:=`(Flip_Max=flip_max,
              Processed_Min=min(dt$Weight),
              Processed_Median=stats::median(dt$Weight),
              Processed_Max=max(dt$Weight))]
  # Duplicate edges retain the highest confidence score.
  dt <- dt[, .(Weight=max(Weight)), by=.(TF, Target)]
  dt[, Method := method]
  setcolorder(dt, c("TF", "Target", "Method", "Weight"))
  attr(dt, "weight_diagnostics") <- diag
  cat(sprintf(" %-10s %s distinct edges (%s)\n", method,
              format(nrow(dt), big.mark=","), basename(p)))
  dt
}
# Reference is built from Top10 for all six input methods (label-free).
# Strict nestedness is checked when stage edges are read. If nonnested files
# are detected we STOP rather than silently dropping stage-specific edges.
# This avoids a background that changes between Top5/Top10/Top20.
read_reference_universe <- function(ds) {
  parts <- lapply(method_names, function(m) {
    path <- get_grn_file(ds,m,stage=reference_stage)
    z <- fread(path,select=1:2,showProgress=FALSE,check.names=FALSE)
    if (ncol(z)!=2L) stop("Invalid reference TF/Target columns: ",path)
    out <- data.table(TF=toupper(trimws(as.character(z[[1L]]))),
                      Target=toupper(trimws(as.character(z[[2L]]))))
    out <- out[!is.na(TF)&!is.na(Target)&nzchar(TF)&nzchar(Target)&TF!=Target]
    unique(out,by=c("TF","Target"))
  })
  universe <- unique(rbindlist(parts,use.names=TRUE),by=c("TF","Target"))
  setkey(universe,TF,Target)
  if (!nrow(universe)) stop("No reference edges for: ",ds)
  universe
}

build_data <- function(ds, chip) {
  parts <- lapply(method_names, function(m) read_grn(ds, m))
  weight_report <- rbindlist(lapply(parts, function(z) attr(z, "weight_diagnostics")),
                             use.names=TRUE, fill=TRUE)
  ds_diag_dir <- file.path(root_out, ds)
  dir.create(ds_diag_dir, recursive=TRUE, showWarnings=FALSE)
  fwrite(weight_report, file.path(ds_diag_dir, "Input_Weight_Diagnostics.csv"))
  long <- rbindlist(parts, use.names=TRUE)
  w <- dcast(long, TF + Target ~ Method, value.var="Weight", fill=NA_real_)
  rm(long, parts)
  for (m in method_names) if (!m %in% names(w)) w[, (m) := NA_real_]
  setkey(w, TF, Target)
  reference <- read_reference_universe(ds)
  # Input Top10 edges and fixed reference should match exactly.
  extra <- w[!reference,on=.(TF,Target),.(TF,Target)]
  if (nrow(extra)) stop(ds, ": ", nrow(extra),
    " Top10 input edges lie outside the fixed Top10 reference; check files and gene IDs.")
  # The i-join retains every reference edge and includes missing stage scores as NA.
  w <- w[reference,on=.(TF,Target)]
  rm(reference)
  # Fix row order before extracting matrices; never reorder the edge table later.
  setkey(w, TF, Target)
  w[, Label := 0L]
  w[chip, on=.(TF, Target), Label := 1L]
  edges <- w[, .(TF, Target, Label)]
  presence <- as.matrix(w[, ..method_names])
  presence <- 1.0 * !is.na(presence)
  # Membership bits are defined exclusively by original input edges.
  # Zero-valued weights still count as PRESENT; missing values do not.
  # Integer values 1,2,4,8,16,32 encode the six input methods.
  method_mask <- as.integer(presence %*% as.integer(2^(seq_along(method_names)-1L)))
  raw <- as.matrix(w[, ..method_names]); storage.mode(raw) <- "double"
  raw[!is.finite(raw)] <- 0
  colnames(raw) <- method_names; colnames(presence) <- method_names
  rm(w)
  # Bagging scales each zero-filled column using its global min and max.
  # Store the scaling coefficients to avoid allocating another score matrix.
  if (normalize_weights) {
    bag_min <- apply(raw, 2L, min)
    bag_max <- apply(raw, 2L, max)
    bag_rng <- bag_max - bag_min
    bag_scale <- ifelse(is.finite(bag_rng) & bag_rng > 0, 1 / bag_rng, 0)
    bag_offset <- ifelse(bag_scale != 0, -bag_min * bag_scale, 0)
  } else {
    bag_scale <- setNames(rep(1, ncol(raw)), colnames(raw))
    bag_offset <- setNames(rep(0, ncol(raw)), colnames(raw))
  }
  # Transductive, GRN scores only. Supervised scaling is separately train-only.
  normalized <- raw
  if (normalize_weights) {
    for (j in seq_len(ncol(raw))) {
      idx <- presence[,j] == 1
      if (!any(idx)) next
      r <- range(raw[idx,j])
      normalized[idx,j] <- if (diff(r) <= 0) 1 else
        (raw[idx,j] - r[1L]) / (r[2L] - r[1L])
      normalized[!idx,j] <- 0
    }
  }
  ranks <- matrix(0, nrow=nrow(edges), ncol=length(method_names),
                  dimnames=list(NULL, method_names))
  for (j in seq_along(method_names)) {
    idx <- which(presence[,j] == 1)
    if (!length(idx)) next
    rr <- rank(-raw[idx,j], ties.method="average")
    ranks[idx,j] <- if (length(idx) == 1L) 1 else
      (length(idx) - rr) / (length(idx) - 1L)
  }
  matched <- sum(edges$Label == 1L)
  positive_tfs <- unique(edges$TF[edges$Label == 1L])
  excluded_tfs <- setdiff(unique(edges$TF), positive_tfs)
  coverage <- data.table(
    Dataset=ds, Input_Stage=input_stage, Candidate_Universe="Fixed_Top10_All6_Union",
    N_Union_Candidates=nrow(edges), N_Stage_Candidates=sum(method_mask!=0L),
     N_Stage_Covered_Positive=sum(edges$Label==1L & method_mask!=0L),
     Stage_ChIP_Edge_Coverage=sum(edges$Label==1L & method_mask!=0L)/nrow(chip),
     Reference_Stage=reference_stage,Evaluation_Universe="Top10_AllSix_Union",
     N_ChIP_Positives=nrow(chip),
    N_Covered_ChIP_Positives=matched, ChIP_Edge_Coverage=matched/nrow(chip),
    N_ChIP_TFs=uniqueN(chip$TF), N_Covered_Positive_TFs=length(positive_tfs),
    ChIP_TF_Coverage=length(positive_tfs)/uniqueN(chip$TF),
    N_Candidate_TFs=uniqueN(edges$TF), N_Excluded_NoPositive_TFs=length(excluded_tfs),
    N_Excluded_NoPositive_Candidates=sum(edges$TF %in% excluded_tfs))
  # Coverage is not exported before model selection: it includes test TF labels.
  if (matched < 20L) stop("Fewer than 20 covered ChIP positives; check gene IDs and inputs.")
  stopifnot(nrow(edges)==nrow(raw), nrow(raw)==nrow(presence), nrow(raw)==nrow(ranks))
  list(edges=edges, raw=raw, norm=normalized, presence=presence,
       method_mask=method_mask, ranks=ranks,
       bagging_scale=as.numeric(bag_scale), bagging_offset=as.numeric(bag_offset),
       chip_total=nrow(chip), coverage=coverage)
}

# -------------------- PER-COMBINATION CANDIDATE UNIVERSE -----------------------
# The fixed Top10 all-six union defines TF folds and common evaluation background.
# A subset S has its OWN universe U_S = union(E_m, m in S) on the six-union rows.
# No ChIP labels are used to build membership (only the presence matrix).
method_union_flag <- function(data, cols) {
  stopifnot(length(cols)>=1L, all(cols %in% seq_along(method_names)))
  bits <- as.integer(sum(2^(cols-1L)))
  bitwAnd(data$method_mask, bits) != 0L
}
filter_union_rows <- function(rows, union_flag) {
  rows[union_flag[rows]]
}

# Training negative sampling is performed AFTER selecting the subset union.
# This preserves train_unlabeled_ratio for each of the 57 model combinations.
get_own_training <- function(data, training_tfs, union_flag, local_seed) {
  if (!length(training_tfs)) stop("Empty training TF set")
  eligible <- which(union_flag & data$edges$TF %in% training_tfs)
  p <- eligible[data$edges$Label[eligible] == 1L]
  u <- eligible[data$edges$Label[eligible] == 0L]
  if (!length(p) || !length(u))
    stop("Selected combination union lacks positive or unverified training edges")
  # sample_group is defined later and called only after script setup finishes.
  sample_group(p, u, train_unlabeled_ratio, local_seed)
}

# COMMON-EVALUATION ONLY: edges omitted by a combination receive a zero score.
# Real zero-scored candidates remain zero and therefore tie with omitted edges.
# Candidate scores are nonnegative for the six implemented strategies.
# Zero padding is not used for model fitting or network export.
apply_common_eval_floor <- function(scores, common_rows, union_flag) {
  if (length(scores) != length(common_rows) || anyNA(common_rows) ||
      any(common_rows < 1L) || any(common_rows > length(union_flag)))
    stop("Common-evaluation score/row alignment failure")
  inside <- union_flag[common_rows]
  if (anyNA(inside) || any(!is.finite(scores[inside])) ||
      any(scores[inside] < 0))
    stop("Invalid inside-union score: expected finite nonnegative values")
  scores[!inside] <- outside_union_eval_score
  scores
}

# Supervised models predict only their own candidate union. All other
# common-background edges are assigned ZERO (ties allowed with real zero scores).
predict_own_on_common <- function(model, data, common_rows, union_flag) {
  out <- rep(outside_union_eval_score, length(common_rows))
  eligible <- which(union_flag[common_rows])
  if (length(eligible))
    out[eligible] <- predict_in_batches(model, data, rows=common_rows[eligible])
  apply_common_eval_floor(out, common_rows, union_flag)
}

# Fixed Top10-union AUROC/AUPR is the principal comparison across stages/subsets;
# own-union AUROC/AUPR answers the distinct question of within-subset ranking.
# These are NOT interchangeable. For AUPR report the baseline and ChIP coverage.
# ----------------------- TF-DISJOINT SPLITS ---------------------------------
# For all modes, TFs with no observed ChIP positive OR no unverified candidate
# cannot yield defined per-TF AUROC/AUPR. Such TFs are excluded consistently.
get_eligible_tf_stats <- function(edges) {
  s <- edges[, .(N_Pos=sum(Label==1L), N_Unverified=sum(Label==0L)), by=TF]
  out <- s[N_Pos>0L & N_Unverified>0L]
  if (nrow(out) < cv_folds + outer_min_test_tfs)
    stop("Too few TFs with positive and unverified candidate edges: ",nrow(out),
         ". Required: ",cv_folds + outer_min_test_tfs)
  out
}
choose_outer_tfs <- function(edges, ds, seed_local) {
  stats <- get_eligible_tf_stats(edges)
  n_tf <- nrow(stats)
  n_test <- max(outer_min_test_tfs, as.integer(round(n_tf * outer_test_fraction)))
  if (n_tf-n_test < cv_folds) stop("Too few development TFs for cross-validation after assigning the test TFs.")

  if (split_mode == "stratified_random") {
    # Sample test TFs independently within strata defined by ChIP-positive edge counts.
    set.seed(seed_local)
    stats <- stats[sample.int(.N)]             # Break ties reproducibly using the seeded shuffle.
    setorder(stats,-N_Pos)
    stats[,Stratum := cut(seq_len(.N), breaks=3L, labels=FALSE,
                          include.lowest=TRUE)]
    stratum_ids <- sort(unique(stats$Stratum))
    sizes <- vapply(stratum_ids, function(g) sum(stats$Stratum==g), integer(1))
    desired <- n_test * sizes / sum(sizes)
    quota <- pmin(sizes, floor(desired))
    # Allocate at least one test TF to each nonempty stratum where possible.
    if (n_test >= length(stratum_ids)) quota <- pmax(quota,1L)
    while (sum(quota) < n_test) {
      cap <- which(quota < sizes)
      if (!length(cap)) stop("The stratum quotas cannot accommodate the requested number of test TFs.")
      j <- cap[which.max(desired[cap]-quota[cap])]
      quota[j] <- quota[j] + 1L
    }
    if (sum(quota) != n_test) stop("The sum of the stratum quotas does not equal the requested test size.")
    test <- unlist(lapply(seq_along(stratum_ids), function(i) {
      pool <- stats[Stratum==stratum_ids[i],TF]
      sample(pool, quota[i], replace=FALSE)
    }), use.names=FALSE)
    return(list(train=setdiff(stats$TF,test),test=test))
  }
  stop("The split mode must be stratified_random: ",split_mode)
}

# Balance the number of TFs first, then the positive edge counts across folds.
make_tf_folds <- function(edges, tf_train, k, seed_local) {
  pos <- edges[Label==1L & TF %in% tf_train, .N, by=TF]
  if (!setequal(pos$TF,tf_train) || nrow(pos)<k) stop("Too few development TFs or missing ChIP-positive edges.")
  set.seed(seed_local)
  pos <- pos[sample.int(nrow(pos))]
  setorder(pos,-N)
  cap <- rep(length(tf_train) %/% k,k)
  extra <- length(tf_train) %% k
  if (extra>0L) cap[seq_len(extra)] <- cap[seq_len(extra)] + 1L
  count <- integer(k); load <- numeric(k); assignment <- integer(nrow(pos))
  for (i in seq_len(nrow(pos))) {
    f_avail <- which(count<cap)
    f <- f_avail[which.min(load[f_avail])]
    assignment[i] <- f
    count[f] <- count[f]+1L
    load[f] <- load[f]+pos$N[i]
  }
  names(assignment) <- pos$TF
  stopifnot(setequal(unique(assignment),seq_len(k)),max(count)-min(count)<=1L)
  assignment
}
sample_group <- function(pos_indices, unlabeled_indices, ratio, seed_local) {
  if (!length(pos_indices) || !length(unlabeled_indices))
    stop("A partition lacks positive or unverified candidate edges.")
  if (is.infinite(ratio)) return(c(pos_indices,unlabeled_indices))
  set.seed(seed_local)
  n <- max(1L,as.integer(min(length(unlabeled_indices),ceiling(length(pos_indices)*ratio))))
  c(pos_indices,unlabeled_indices[sample.int(length(unlabeled_indices),n,replace=FALSE)])
}
get_indices <- function(edges,tfs,ratio,seed_local) {
  eligible <- which(edges$TF %in% tfs)
  p <- eligible[edges$Label[eligible]==1L]
  u <- eligible[edges$Label[eligible]==0L]
  sample_group(p,u,ratio,seed_local)
}
make_split_plan <- function(edges,ds,ds_seed) {
  outer <- choose_outer_tfs(edges,ds,ds_seed)
  folds <- make_tf_folds(edges,outer$train,cv_folds,ds_seed+1L)
  fold_list <- lapply(seq_len(cv_folds),function(k) {
    valid_tfs <- names(folds)[folds==k]
    train_tfs <- setdiff(outer$train,valid_tfs)
    tr <- get_indices(edges,train_tfs,train_unlabeled_ratio,ds_seed+10L+k)
    va <- get_indices(edges,valid_tfs,cv_unlabeled_ratio,ds_seed+100L+k)
    stopifnot(length(intersect(edges$TF[tr],edges$TF[va]))==0L)
    list(fold=k,train=tr,valid=va,train_tfs=train_tfs,valid_tfs=valid_tfs)
  })
  outer_train <- get_indices(edges,outer$train,train_unlabeled_ratio,ds_seed+500L)
  outer_test <- get_indices(edges,outer$test,outer_unlabeled_ratio,ds_seed+600L)
  stopifnot(length(intersect(edges$TF[outer_train],edges$TF[outer_test]))==0L)
  list(outer=outer,fold_map=folds,folds=fold_list,train=outer_train,test=outer_test)
}

# ---------------------------- EVALUATION ------------------------------------
metric <- function(label, score) {
  # Never silently drop an invalid prediction and change the shared background.
  if (length(label) != length(score) || !length(label))
    stop("Labels and scores must have the same positive length.")
  if (anyNA(label) || any(!label %in% c(0L,1L))) stop("Labels must be 0 or 1, without NA.")
  if (any(!is.finite(score))) stop("Non-finite predictions: shared rows cannot be dropped.")
  y <- as.integer(label); p <- as.numeric(score)
  np <- sum(y==1L); nu <- sum(y==0L); b <- np/length(y)
  out <- c(AUROC=NA_real_, AUPR=NA_real_, Baseline=b, Lift=NA_real_,
           N_Pos=np, N_Unverified=nu)
  if (!np || !nu) return(out)
  out["AUROC"] <- as.numeric(pROC::auc(pROC::roc(y, p, levels=c(0,1), direction="<", quiet=TRUE)))
  out["AUPR"] <- PRROC::pr.curve(scores.class0=p[y==1L],
                                 scores.class1=p[y==0L], curve=FALSE)$auc.integral
  out["Lift"] <- out["AUPR"] / b
  out
}
# A missing class (or an empty combination union) is not a numeric AUROC/AUPR.
# Preserve counts and baseline, and allow the common-background comparison to continue.
metric_safe <- function(label, score) {
  if (length(label) != length(score)) stop("score/label length mismatch")
  if (!length(label))
    return(c(AUROC=NA_real_,AUPR=NA_real_,Baseline=NA_real_,
             Lift=NA_real_,N_Pos=0,N_Unverified=0))
  metric(label,score)
}

per_tf_metric <- function(tf, label, score) {
  if (length(tf)!=length(label) || length(tf)!=length(score))
    stop("TF, labels and scores must align")
  z <- data.table(TF=as.character(tf), Label=as.integer(label), Score=as.numeric(score))
  z[, {
    met <- metric(Label, Score)
    .(AUROC=unname(met["AUROC"]), AUPR=unname(met["AUPR"]),
      Baseline=unname(met["Baseline"]), Lift=unname(met["Lift"]),
      N_Pos=as.integer(met["N_Pos"]), N_Unverified=as.integer(met["N_Unverified"]))
  }, by=TF]
}
cluster_ci <- function(tf, y, score, reps, local_seed) {
  out <- c(AUROC_L95=NA_real_, AUROC_U95=NA_real_,
           AUPR_L95=NA_real_, AUPR_U95=NA_real_, Valid_Boot=0)
  if (reps <= 0) return(out)
  groups <- split(seq_along(tf), tf)
  if (length(groups) < bootstrap_min_tfs) return(out)
  set.seed(local_seed)
  vals <- matrix(NA_real_, nrow=reps, ncol=2L)
  for (b in seq_len(reps)) {
    rows <- unlist(groups[sample.int(length(groups), length(groups), replace=TRUE)], use.names=FALSE)
    m <- metric(y[rows], score[rows]); vals[b,] <- m[c("AUROC", "AUPR")]
  }
  valid <- complete.cases(vals); out["Valid_Boot"] <- sum(valid)
  if (sum(valid) < max(15L, ceiling(reps*0.5))) return(out)
  out[1:2] <- unname(quantile(vals[valid,1], c(.025,.975)))
  out[3:4] <- unname(quantile(vals[valid,2], c(.025,.975)))
  out
}

# ------------------------- UNSUPERVISED SCORES ------------------------------
# Bagging: algorithm-column Bootstrap with replacement, NOT CART/bootstrap of edges.
# For B rounds and M selected algorithms: draw B*M integer column indices,
# aggregate their frequencies, then use ONE matrix multiply; this is
# mathematically equivalent to averaging rowMeans over all B bootstrap replicates.
bagging_column_coefficients <- function(M, B, local_seed) {
  if (M < 1L || B < 1L) stop("Bagging needs >=1 algorithm and >=1 bootstrap round")
  set.seed(local_seed)
  sampled <- sample.int(M, size=as.integer(M * B), replace=TRUE)
  tabulate(sampled, nbins=M) / (M * B)
}
unsup_scores <- function(kind, data, cols, local_seed) {
  if (kind == "Bagging") {
    coeff <- bagging_column_coefficients(length(cols), bagging_B, local_seed)
    # norm_j(edge) = raw_j(edge)*scale_j + offset_j.
    # Coefficients apply to each method's globally Min-Max-normalized scores.
    effective_weights <- coeff * data$bagging_scale[cols]
    intercept <- sum(coeff * data$bagging_offset[cols])
    ans <- as.numeric(data$raw[,cols,drop=FALSE] %*% effective_weights) + intercept
    if (normalize_weights) ans <- pmin(1, pmax(0, ans))
    if (length(ans) != nrow(data$edges) || any(!is.finite(ans)))
      stop("Bagging produced invalid scores")
    return(ans)
  }
  if (kind == "HardVoting") {
    votes <- rowSums(data$presence[,cols,drop=FALSE])
    # A strict majority; all nonmajority edges receive score 0 IN COMMON universe.
    score <- rowMeans(data$norm[,cols,drop=FALSE])
    score[votes <= length(cols)/2] <- 0
    return(score)
  }
  if (kind == "QRF") return(rowMeans(data$ranks[,cols,drop=FALSE]))
  stop("Unknown unsupervised method: ", kind)
}

# --------------------------- SUPERVISED SCORES ------------------------------
# Fit/transform data-driven min-max using ONLY each training fold. Presence
# indicators are kept unchanged; they distinguish absent from zero weight.
features_for <- function(data, rows, cols) {
  cbind(data$raw[rows,cols,drop=FALSE], data$presence[rows,cols,drop=FALSE])
}
# glmnet needs >=2 columns. For one varying predictor, fit the same ridge-logistic
# objective directly, with an unpenalized intercept and a penalized slope.
fit_ridge_1d <- function(x, y, lambda) {
  loss <- function(beta) {
    eta <- beta[1L] + beta[2L]*x
    mean(pmax(eta, 0) + log1p(exp(-abs(eta))) - y*eta) + lambda*beta[2L]^2/2
  }
  gradient <- function(beta) {
    residual <- plogis(beta[1L] + beta[2L]*x) - y
    c(mean(residual), mean(residual*x) + lambda*beta[2L])
  }
  fit <- optim(c(qlogis(mean(y)), 0), fn=loss, gr=gradient,
               method="BFGS", control=list(maxit=2000L, reltol=1e-10))
  if (fit$convergence != 0L || any(!is.finite(fit$par)))
    stop("Single-predictor ridge fit did not converge.")
  list(coefficients=fit$par, objective=fit$value, convergence=fit$convergence)
}
fit_supervised <- function(kind, data, cols, train_rows, local_seed) {
  y <- data$edges$Label[train_rows]
  if (!all(c(0L,1L) %in% y)) stop("Supervised training requires both positive and unverified edges.")
  X <- features_for(data, train_rows, cols)
  names_x <- c(paste0("W_",method_names[cols]), paste0("P_",method_names[cols]))
  colnames(X) <- names_x
  if (normalize_weights) {
    center <- apply(X[,seq_along(cols),drop=FALSE],2,min)
    scales <- apply(X[,seq_along(cols),drop=FALSE],2,max) - center
    scales[!is.finite(scales) | scales <= 0] <- 1
    X[,seq_along(cols)] <- sweep(sweep(X[,seq_along(cols),drop=FALSE],2,center,"-"),2,scales,"/")
  } else { center <- rep(0,length(cols)); scales <- rep(1,length(cols)) }
  keep <- apply(X, 2L, function(z) length(unique(z)) > 1L)
  if (!any(keep)) return(list(kind="Constant", constant=mean(y)))
  X <- X[,keep,drop=FALSE]
  set.seed(local_seed)
  if (kind == "XGBoost") {
    fit <- xgboost::xgb.train(
      params=list(objective="binary:logistic", eval_metric="auc", max_depth=xgb_depth,
                  eta=xgb_eta, subsample=.8, colsample_bytree=.8, nthread=xgb_nthread,
                  seed=local_seed),
      data=xgboost::xgb.DMatrix(X,label=y), nrounds=xgb_rounds, verbose=0)
  } else if (kind == "AdaBoost") {
    train_df <- as.data.frame(X, check.names=FALSE)
    train_df$label <- factor(y, levels=c(0,1))
    fit <- ada::ada(label ~ ., data=train_df, iter=ada_iter, nu=ada_nu,
                    type="real", control=rpart::rpart.control(maxdepth=2L,
                                                              minsplit=20L, cp=0.001))
  } else if (kind == "Stacking") {
    if (ncol(X) == 1L) {
      fit <- fit_ridge_1d(as.numeric(X[,1L]), y, ridge_lambda)
      kind <- "Ridge1D"
    } else {
      fit <- glmnet::glmnet(X, y, family="binomial", alpha=0,
                            lambda=ridge_lambda, standardize=FALSE)
    }
  } else stop("Unknown supervised model: ", kind)
  list(kind=kind, fit=fit, cols=cols, center=center, scales=scales,
       keep=keep, names_x=names_x)
}
predict_supervised <- function(model, data, rows) {
  if (model$kind == "Constant") return(rep(model$constant, length(rows)))
  cols <- model$cols
  X <- features_for(data, rows, cols)
  colnames(X) <- model$names_x
  if (normalize_weights) {
    X[,seq_along(cols)] <- sweep(sweep(X[,seq_along(cols),drop=FALSE],
                                       2,model$center,"-"),2,model$scales,"/")
    X[,seq_along(cols)] <- pmax(0, pmin(1, X[,seq_along(cols),drop=FALSE]))
  }
  X <- X[,model$keep,drop=FALSE]
  if (model$kind == "XGBoost")
    return(as.numeric(predict(model$fit, xgboost::xgb.DMatrix(X))))
  if (model$kind == "Stacking")
    return(as.numeric(predict(model$fit, newx=X, type="response", s=ridge_lambda)))
  if (model$kind == "Ridge1D")
    return(as.numeric(plogis(model$fit$coefficients[1L] +
                               model$fit$coefficients[2L]*X[,1L])))
  if (model$kind == "AdaBoost") {
    probs <- predict(model$fit, newdata=as.data.frame(X,check.names=FALSE), type="prob")
    if (is.null(dim(probs))) stop("AdaBoost prediction did not return a probability matrix.")
    if (!is.null(colnames(probs)) && "1" %in% colnames(probs))
      return(as.numeric(probs[,"1"]))
    if (ncol(probs) != 2L) stop("AdaBoost prediction did not return two probability columns.")
    return(as.numeric(probs[,2L]))
  }
  stop("Unknown model type")
}
predict_in_batches <- function(model, data, rows=seq_len(nrow(data$edges)),
                               size=predict_batch_size) {
  if (!length(rows)) return(numeric())
  ans <- numeric(length(rows))
  for (start in seq.int(1L, length(rows), by=size)) {
    end <- min(length(rows), start+size-1L)
    positions <- seq.int(start, end)
    ans[positions] <- predict_supervised(model, data, rows[positions])
  }
  if (any(!is.finite(ans))) stop("Non-finite scores in batched prediction.")
  ans
}

# ------------------------ CROSS-VALIDATED SELECTION -------------------------
strategies <- c("Bagging","HardVoting","QRF","XGBoost","AdaBoost","Stacking")
supervised <- c("XGBoost","AdaBoost","Stacking")
combinations <- function() {
  if (only_all6) return(list(seq_along(method_names)))
  first <- if (include_singletons) 1L else 2L
  unlist(lapply(seq.int(first, length(method_names)), function(k)
    combn(seq_along(method_names), k, simplify=FALSE)), recursive=FALSE)
}
rank_cv <- function(tab) {
  out <- copy(tab)
  prefix <- if (selection_universe == "Common_All6") "CV_" else "CV_OwnUnion_"
  primary <- paste0(prefix, selection_metric)
  secondary <- if (selection_metric == "Pooled_AUROC")
    paste0(prefix, "Pooled_AUPR") else paste0(prefix, "Pooled_AUROC")
  setorderv(out, c(primary,secondary,"NumMethods","Strategy","Combination"),
            c(-1L,-1L,1L,1L,1L), na.last=TRUE)
  out
}
# Compare each strategy's CV-selected subset with the six-algorithm union.
# FULL uses all 57/63 input subsets for EACH of the six strategies.
candidate_indices <- function(ds, strategy, combination_names) {
  # STRICT: the search space is fixed without looking at heldout results.
  seq_along(combination_names)
}

cv_one <- function(strategy, data, cols, plan, comb_idx, ds_seed) {
  own_flag <- method_union_flag(data, cols)
  fold_rows <- vector("list",length(plan$folds))
  tf_rows <- vector("list",length(plan$folds))
  own_tf_rows <- vector("list",length(plan$folds))
  oof_labels <- vector("list",length(plan$folds))
  oof_scores <- vector("list",length(plan$folds))
  own_labels <- vector("list",length(plan$folds))
  own_scores <- vector("list",length(plan$folds))
  if (!strategy %in% supervised)
    full_score <- unsup_scores(strategy,data,cols,ds_seed+comb_idx*37L)
  for (fold in plan$folds) {
    valid <- fold$valid  # fixed common union for every candidate combination
    own_valid <- filter_union_rows(valid,own_flag)
    own_train <- if (strategy %in% supervised)
      get_own_training(data, setdiff(plan$outer$train, fold$valid_tfs), own_flag,
                       ds_seed + 10L + fold$fold) else
      which(own_flag & data$edges$TF %in% setdiff(plan$outer$train, fold$valid_tfs))
    attempt <- tryCatch({
      if (strategy %in% supervised) {
        if (!length(own_train)) stop("No combination-union training edges")
        model <- fit_supervised(strategy,data,cols,own_train,
                                ds_seed+comb_idx*100L+fold$fold)
        pred <- predict_own_on_common(model,data,valid,own_flag)
        rm(model)
      } else {
        pred <- apply_common_eval_floor(full_score[valid], valid, own_flag)
      }
      # Outside own union is always 0 for COMMON evaluation; genuine zeros tie.
      y <- data$edges$Label[valid]
      met <- metric(y,pred)
      tfm <- per_tf_metric(data$edges$TF[valid],y,pred)
      pos_in_valid <- match(own_valid, valid)
      oy <- y[pos_in_valid]; op <- pred[pos_in_valid]
      own_met <- metric_safe(oy,op)  # returns NA if no positive or no unverified edges
      own_tfm <- if (length(own_valid))
        per_tf_metric(data$edges$TF[own_valid],oy,op) else data.table()
      if (any(!is.finite(met[c("AUROC","AUPR")])) ||
          any(!is.finite(tfm$AUROC)) || any(!is.finite(tfm$Lift)))
        stop("Undefined pooled/per-TF COMMON-union metric in validation fold")
      list(fold=data.table(
             Fold=fold$fold,AUROC=met["AUROC"],AUPR=met["AUPR"],
             Baseline=met["Baseline"],Lift=met["Lift"],N_Pos=met["N_Pos"],
             N_Unverified=met["N_Unverified"],N_TF=nrow(tfm),
             OwnUnion_AUROC=own_met["AUROC"],OwnUnion_AUPR=own_met["AUPR"],
             OwnUnion_Baseline=own_met["Baseline"],OwnUnion_Lift=own_met["Lift"],
             OwnUnion_N_Pos=own_met["N_Pos"],
             OwnUnion_N_Unverified=own_met["N_Unverified"],
             OwnUnion_N_TF=if (nrow(own_tfm)) uniqueN(own_tfm$TF) else 0L,
             Train_N_Common=length(fold$train),Train_N_OwnUnion=length(own_train),
             Status="OK",Error=""),
           tf=tfm[,Fold:=fold$fold],
           own_tf=if (nrow(own_tfm)) own_tfm[,Fold:=fold$fold] else own_tfm,
           label=as.integer(y),score=as.numeric(pred),
           own_label=as.integer(oy),own_score=as.numeric(op))
    },error=function(e) {
      y <- data$edges$Label[valid]
      list(fold=data.table(Fold=fold$fold,AUROC=NA_real_,AUPR=NA_real_,
                           Baseline=mean(y),Lift=NA_real_,N_Pos=sum(y==1L),
                           N_Unverified=sum(y==0L),N_TF=length(fold$valid_tfs),
                           OwnUnion_AUROC=NA_real_,OwnUnion_AUPR=NA_real_,
                           OwnUnion_Baseline=NA_real_,OwnUnion_Lift=NA_real_,
                           OwnUnion_N_Pos=NA_integer_,OwnUnion_N_Unverified=NA_integer_,
                           OwnUnion_N_TF=NA_integer_,Train_N_Common=length(fold$train),
                           Train_N_OwnUnion=length(own_train),
                           Status="FAILED",Error=conditionMessage(e)),
           tf=data.table(),own_tf=data.table(),label=integer(),score=numeric(),
           own_label=integer(),own_score=numeric())
    })
    fold_rows[[fold$fold]] <- attempt$fold
    tf_rows[[fold$fold]] <- attempt$tf
    own_tf_rows[[fold$fold]] <- attempt$own_tf
    oof_labels[[fold$fold]] <- attempt$label
    oof_scores[[fold$fold]] <- attempt$score
    own_labels[[fold$fold]] <- attempt$own_label
    own_scores[[fold$fold]] <- attempt$own_score
  }
  fold_out <- rbindlist(fold_rows,fill=TRUE)
  empty <- c(AUROC=NA_real_,AUPR=NA_real_,Baseline=NA_real_,
             Lift=NA_real_,N_Pos=NA_real_,N_Unverified=NA_real_)
  pooled <- own_pooled <- empty
  if (nrow(fold_out)==cv_folds && all(fold_out$Status=="OK")) {
    pooled <- metric(unlist(oof_labels,use.names=FALSE),
                     unlist(oof_scores,use.names=FALSE))
    own_pooled <- metric_safe(unlist(own_labels,use.names=FALSE),
                         unlist(own_scores,use.names=FALSE))
  }
  list(folds=fold_out,per_tf=rbindlist(tf_rows,fill=TRUE),
       own_per_tf=rbindlist(own_tf_rows,fill=TRUE),
       pooled=pooled, own_pooled=own_pooled,
       n_own_all=sum(own_flag),
       n_own_dev=sum(own_flag & data$edges$TF %in% plan$outer$train),
       n_own_dev_pos=sum(own_flag & data$edges$Label==1L &
                         data$edges$TF %in% plan$outer$train))
}


locked_result <- function(ds, strategy, combination, scores, data, plan, chip,
                          ds_seed, ci_index, result_type,
                          cv_auroc=NA_real_, cv_aupr=NA_real_, global_winner=FALSE,
                           cols=NULL) {
  y <- data$edges$Label[plan$test]
  test_tf <- data$edges$TF[plan$test]
  met <- metric(y, scores)
  if (any(!is.finite(met[c("AUROC", "AUPR")]))) stop("Undefined locked-test metric.")
  ci <- cluster_ci(test_tf, y, scores, bootstrap_tf_reps, ds_seed+12000L+ci_index)
  n_tf <- uniqueN(test_tf)
  tfm <- per_tf_metric(test_tf,y,scores)
  own_idx <- if (is.null(cols)) seq_along(plan$test) else
    which(method_union_flag(data,cols)[plan$test])
  own <- metric_safe(y[own_idx],scores[own_idx])
  own_tfm <- if (length(own_idx))
    per_tf_metric(test_tf[own_idx],y[own_idx],scores[own_idx]) else data.table()
  own_ci <- if (length(own_idx) && is.finite(own["AUROC"]) &&
                uniqueN(test_tf[own_idx])>=bootstrap_min_tfs)
    cluster_ci(test_tf[own_idx],y[own_idx],scores[own_idx],
               bootstrap_tf_reps,ds_seed+800000L+ci_index) else
    c(AUROC_L95=NA_real_,AUROC_U95=NA_real_,
      AUPR_L95=NA_real_,AUPR_U95=NA_real_,Valid_Boot=0)
  n_chip_test <- nrow(chip[TF %in% plan$outer$test])
  macro_auc <- mean(tfm$AUROC,na.rm=TRUE)
  macro_lift <- mean(tfm$Lift,na.rm=TRUE)
  ci_status <- if (bootstrap_tf_reps <= 0L) "Skipped" else
    if (n_tf < bootstrap_min_tfs) "TooFewTFs" else
      if (any(!is.finite(ci[1:4]))) "TooFewValidReplicates" else "ConditionalOnFittedModel"
  data.table(Dataset=ds, Result_Type=result_type, Strategy=strategy,
             Combination=combination, Selection_Metric=selection_metric,
             Selection_Universe=selection_universe,
             Train_Universe=if (is.null(cols)) "Not_Applicable" else "Combination_Union",
             CV_Pooled_AUROC=cv_auroc, CV_Pooled_AUPR=cv_aupr,
             Is_Global_CV_Winner=global_winner,
             Evaluation="ObservedChIP_Positive_vs_Unverified_PU_Proxy",
             Candidate_Universe="Union_of_supplied_GRN_files", Input_Stage=input_stage,
             Test_AUROC=met["AUROC"], Test_AUPR=met["AUPR"],
             Test_Macro_AUROC=macro_auc, Test_Macro_AUPR_Lift=macro_lift,
             Test_Evaluable_TFs=sum(is.finite(tfm$AUROC)),
             Test_AUPR_Baseline=met["Baseline"], Test_AUPR_Lift=met["Lift"],
                           Test_Pos=met["N_Pos"], Test_Unverified=met["N_Unverified"], N_Test_TFs=n_tf,
              OwnUnion_Test_AUROC=unname(own["AUROC"]),
              OwnUnion_Test_AUPR=unname(own["AUPR"]),
              OwnUnion_Test_AUPR_Baseline=unname(own["Baseline"]),
              OwnUnion_Test_AUPR_Lift=unname(own["Lift"]),
              OwnUnion_Test_N_Pos=as.integer(own["N_Pos"]),
              OwnUnion_Test_N_Unverified=as.integer(own["N_Unverified"]),
              OwnUnion_Test_N_TFs=if (nrow(own_tfm)) uniqueN(own_tfm$TF) else 0L,
              OwnUnion_Test_ChIP_Coverage=if (n_chip_test) unname(own["N_Pos"])/n_chip_test else NA_real_,
              OwnUnion_AUROC_L95=unname(own_ci["AUROC_L95"]),
              OwnUnion_AUROC_U95=unname(own_ci["AUROC_U95"]),
              OwnUnion_AUPR_L95=unname(own_ci["AUPR_L95"]),
              OwnUnion_AUPR_U95=unname(own_ci["AUPR_U95"]),
              Common_Test_N_Candidates=length(plan$test),
              Matched_ChIP_Pos_Coverage=sum(data$edges$Label==1L)/data$chip_total,
             Matched_ChIP_TF_Coverage=uniqueN(data$edges$TF[data$edges$Label==1L])/uniqueN(chip$TF),
             Locked_TF_Pos_Coverage=met["N_Pos"]/nrow(chip[TF %in% plan$outer$test]),
             CV_Unlabeled_Ratio=cv_unlabeled_ratio, Test_Unlabeled_Ratio=outer_unlabeled_ratio,
             AUROC_L95=ci["AUROC_L95"], AUROC_U95=ci["AUROC_U95"],
             AUPR_L95=ci["AUPR_L95"], AUPR_U95=ci["AUPR_U95"],
             CI_Status=ci_status, Valid_Boot=ci["Valid_Boot"])
}

export_network <- function(ds, strategy, full_score, data, ds_out, cols=NULL) {
  if (length(full_score) != nrow(data$edges) || any(!is.finite(full_score)))
    stop("Invalid full-network scores.")
  union_flag <- if (is.null(cols)) rep(TRUE,length(full_score)) else
    method_union_flag(data,cols)
  # Network selection uses ONLY the selected methods' reported candidate union.
  # This is DIFFERENT from the common-universe zero padding used in evaluation.
  candidate_rows <- which(union_flag)
  candidate_scores <- full_score[candidate_rows]
  n_positive <- sum(candidate_scores > 0)
  status <- "OK"; threshold <- NA_real_; requested <- 0L; at_boundary <- 0L
  n_export_zero <- 0L
  top <- data.table(TF=character(), Target=character(), EdgeWeight=numeric(),
                    N_Methods_Support=integer(), Supporting_Methods=character(),
                    Edge_Selection=character())
  if (!length(candidate_rows)) {
    status <- "NoSelectedUnionCandidates"
  } else {
    requested <- max(1L, ceiling(length(candidate_rows)*network_top_pct))
    # Export 100% of the selected combination candidate union. Deterministically resolve
    # boundary ties using gene identifiers; zeros may be selected if needed.
    ord <- order(-candidate_scores,
                 data$edges$TF[candidate_rows],data$edges$Target[candidate_rows],
                 method="radix")
    selected <- candidate_rows[ord[seq_len(requested)]]
    threshold <- min(full_score[selected])
    at_boundary <- sum(candidate_scores == threshold)
    n_export_zero <- sum(full_score[selected] == 0)
    top <- data$edges[selected, .(TF, Target)]
    top[, EdgeWeight := full_score[selected]]
    # Support refers ONLY to the base algorithms in the selected combination.
    support_cols <- if (is.null(cols)) seq_along(method_names) else cols
    support <- data$presence[selected,support_cols,drop=FALSE]
    top[, N_Methods_Support := as.integer(rowSums(support))]
    support_mask <- as.integer(support %*% (2^(seq_len(ncol(support))-1L)))
    support_lookup <- vapply(0:(2^ncol(support)-1L), function(k) {
      present <- as.integer(intToBits(k)[seq_len(ncol(support))]) == 1L
      paste(method_names[support_cols][which(present)], collapse="+")
    },character(1))
    top[, Supporting_Methods := support_lookup[support_mask+1L]]
    top[, Edge_Selection := "SelectedUnion_TopRank_Not_Statistical_Pvalue"]
    setorderv(top, c("EdgeWeight", "TF", "Target"), c(-1L, 1L, 1L))
  }
  pct_name <- format(100*network_top_pct, trim=TRUE, scientific=FALSE)
  filename <- paste0(strategy, "_Ranked_SelectedUnion_AllEdges.csv")
  fwrite(top, file.path(ds_out, filename))
  data.table(Dataset=ds, Strategy=strategy, File=filename, Status=status,
             N_Union_Candidates=sum(union_flag),
             N_Common_All6_Candidates=length(full_score),
             Candidate_Universe=if (is.null(cols)) "Common_All6" else "Combination_Union",
             N_Positive_Scores=n_positive,
             Requested_Fraction_Of_Selected_Union=network_top_pct,
             Requested_Edges=requested, Exported_Edges=nrow(top),
             Exported_ZeroScore_Edges=n_export_zero,
             Score_Threshold=threshold, Edges_Tied_At_Threshold=at_boundary,
             Boundary_Tie_Policy="FixedCount_ScoreDesc_ThenTF_ThenTarget",
             Selection_Interpretation="All edges ranked in selected union; no empirical p-value or FDR")
}

# ============================================================================
# Per-strategy reporting and resumable per-combination checkpoints
# ============================================================================

# Save each checkpoint to a fresh temporary file in the destination directory,
# then rename. Checkpoints never overwrite existing files.
atomic_rds_new <- function(object, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  if (file.exists(path)) stop("Refusing to overwrite an existing checkpoint: ",path)
  tmp <- tempfile(pattern=".tmp_checkpoint_", tmpdir=dirname(path), fileext=".rds")
  on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
  saveRDS(object, tmp, compress="gzip")
  if (!file.rename(tmp,path)) stop("Atomic checkpoint rename failed: ",path)
  invisible(path)
}

atomic_csv <- function(z,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(pattern=".tmp_csv_",tmpdir=dirname(path),fileext=".csv")
  on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
  fwrite(z,tmp)
  if (!file.rename(tmp,path)) {
    # Windows cannot rename over an existing file: safe fallback to a copy.
    if (!file.copy(tmp,path,overwrite=TRUE)) stop("Failed to write CSV: ",path)
    unlink(tmp)
  }
  invisible(path)
}

atomic_lines <- function(txt,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(pattern=".tmp_txt_",tmpdir=dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
  writeLines(as.character(txt), tmp, useBytes=TRUE)
  if (!file.rename(tmp,path)) {
    if (!file.copy(tmp,path,overwrite=TRUE)) stop("Failed to write: ",path)
    unlink(tmp)
  }
  invisible(path)
}

atomic_new_marker <- function(txt,path) {
  if (!file.exists(path)) atomic_lines(txt,path)
}

# Include both parameters and file metadata in a run fingerprint so that
# checkpoints from another input stage, split, seed or model settings are not
# accidentally reused. Store a separate fingerprint for each stress/mode.
make_run_fingerprint <- function(ds,active_mode) {
  in_paths <- vapply(method_names,function(m) get_grn_file(ds,m),character(1))
  ref_paths <- vapply(method_names,function(m) get_grn_file(ds,m,stage=reference_stage),character(1))
  files <- unique(c(chip_path,unname(in_paths),unname(ref_paths)))
  if (any(!file.exists(files))) stop("Missing fingerprint input file: ",
                                     paste(files[!file.exists(files)],collapse=", "))
  fi <- file.info(files)
  input_info <- data.frame(Path=normalizePath(files,winslash="/",mustWork=TRUE),
                           Bytes=as.numeric(fi$size),Modified=as.numeric(fi$mtime),stringsAsFactors=FALSE)
  cfg <- list(
    Version="gsGNA_ThresholdSweep_FixedTop10_v1", Dataset=ds, SplitMode=active_mode,
    ThresholdStages=threshold_stages,ReferenceStage=reference_stage,
    ThresholdSweepDevelopment=threshold_sweep_development,
    SelectionUniverse=selection_universe,TrainingUniverse="Combination_Union",
    OutsideUnionEvalScore=outside_union_eval_score,
    Seed=seed,InputStage=input_stage,InputFiles=input_info,
    MethodNames=method_names, WeightModes=weight_modes,
    ScoreDirections=score_directions, Normalize=normalize_weights,
    ModelStrategies=strategies,RunMode=run_mode,IncludeSingletons=include_singletons,
    OnlyAll6=only_all6,
    OuterFraction=outer_test_fraction,OuterMinTF=outer_min_test_tfs,
    CVFolds=cv_folds,TrainRatio=train_unlabeled_ratio,
    CVRatio=cv_unlabeled_ratio,TestRatio=outer_unlabeled_ratio,
    Metric=selection_metric, XGB=c(xgb_rounds,xgb_depth,xgb_eta,xgb_nthread),
    ADA=c(ada_iter,ada_nu),Ridge=ridge_lambda,
    BaggingMode=bagging_run_mode,BaggingB=bagging_B,
    ReviewerRepeats=reviewer_cv_repeats,ReviewerCI=reviewer_ci_reps,
    ReviewerEnabled=reviewer_repeat_enabled, HeldoutStatus=heldout_evidence_status,
    Superiority=c(superiority_reference,superiority_adjustment),
    SuperiorityBootstrap=superiority_bootstrap_reps,
    SuperiorityPermutation=superiority_permutation_reps,
    RunLockedSix=run_locked_six_comparison,
    ExternalValidationPath=external_validation_path,
    EvaluateHeldout=evaluate_locked_test,
    BootstrapTF=bootstrap_tf_reps,BootstrapMinTF=bootstrap_min_tfs,
    ExportTopPct=network_top_pct,SavePredictions=save_test_predictions,
    SaveNetworks=save_winner_networks, Batch=predict_batch_size)
  tmp <- tempfile(fileext=".rds")
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(cfg,tmp,version=2)
  list(Fingerprint=unname(tools::md5sum(tmp)), Config=cfg)
}

verify_or_create_signature <- function(ds_dir, fp) {
  dir.create(ds_dir,recursive=TRUE,showWarnings=FALSE)
  sig_path <- file.path(ds_dir,"Run_Fingerprint.txt")
  if (file.exists(sig_path)) {
    old <- readLines(sig_path,warn=FALSE)
    if (length(old)!=1L || !identical(old,fp$Fingerprint))
      stop("Inputs or settings changed since the previous run: ", ds_dir,
           "\nRestore the prior configuration or choose a fresh threshold_study_root.")
    if (!resume_existing) stop("The run directory already exists and resume_existing=FALSE: ",ds_dir)
  } else {
    existing <- setdiff(list.files(ds_dir,all.files=TRUE,no..=TRUE),
                        c(".DS_Store"))
    if (length(existing)) stop("The run directory has existing files but no Run_Fingerprint.txt: ",ds_dir)
    atomic_rds_new(fp$Config,file.path(ds_dir,"Run_Config_Exact.rds"))
    atomic_lines(fp$Fingerprint,sig_path)
  }
  invisible(TRUE)
}

# Persist a single completed (even if FAILED) combination. All OOF labels and
# predictions are only in memory while the combination is running; the saved
# records include all resulting metrics and failure diagnostics.
one_combination <- function(ds,s,j,comb_list,comb_names,data,plan,ds_seed,fp) {
  cols <- comb_list[[j]]; name <- comb_names[j]
  res <- tryCatch(cv_one(s,data,cols,plan,j,ds_seed),error=function(e) {
    list(folds=rbindlist(lapply(plan$folds,function(f) {
      y <- data$edges$Label[f$valid]
      data.table(Fold=f$fold,AUROC=NA_real_,AUPR=NA_real_,Baseline=mean(y),
                 Lift=NA_real_,N_Pos=sum(y==1L),N_Unverified=sum(y==0L),
                 N_TF=length(f$valid_tfs),Status="FAILED",Error=conditionMessage(e))
    })),per_tf=data.table(),own_per_tf=data.table(),
    pooled=c(AUROC=NA_real_,AUPR=NA_real_,Baseline=NA_real_,
             Lift=NA_real_,N_Pos=NA_real_,N_Unverified=NA_real_),
    own_pooled=c(AUROC=NA_real_,AUPR=NA_real_,Baseline=NA_real_,
                 Lift=NA_real_,N_Pos=NA_real_,N_Unverified=NA_real_),
    n_own_all=NA_integer_,n_own_dev=NA_integer_,n_own_dev_pos=NA_integer_)
  })
  folds <- res$folds; tfm <- res$per_tf; own_tfm <- res$own_per_tf
  pooled <- res$pooled; own_pooled <- res$own_pooled
  folds[,c("Dataset","Strategy","Combination","NumMethods") :=
          list(ds,s,name,length(cols))]
  valid_folds <- sum(folds$Status=="OK" & is.finite(folds$AUROC) & is.finite(folds$AUPR))
  complete_common <- nrow(folds)==cv_folds && valid_folds==cv_folds &&
    nrow(tfm)==length(plan$outer$train) &&
    uniqueN(tfm$TF)==length(plan$outer$train) &&
    all(is.finite(tfm$AUROC)) && all(is.finite(tfm$Lift)) &&
    all(is.finite(pooled[c("AUROC","AUPR","Lift")]))
  complete_own <- complete_common &&
    all(is.finite(own_pooled[c("AUROC","AUPR","Lift")])) &&
    all(is.finite(folds$OwnUnion_AUROC)) &&
    all(is.finite(folds$OwnUnion_AUPR))
  if (nrow(tfm)) {
    tfm[,c("Dataset","Strategy","Combination","NumMethods") :=
          list(ds,s,name,length(cols))]
  }
  if (nrow(own_tfm)) {
    own_tfm[,c("Dataset","Strategy","Combination","NumMethods") :=
                list(ds,s,name,length(cols))]
  }
  own_dev_pos_total <- data$edges[Label==1L & TF %in% plan$outer$train,.N]
  summary <- data.table(
    Dataset=ds,Strategy=s,Combination=name,NumMethods=length(cols),
    Selection_Metric=selection_metric,Selection_Universe=selection_universe,
    Train_Universe="Combination_Union",Evaluation_Common="Common_All6",
    Outside_Union_Evaluation_Score=outside_union_eval_score,
    Evaluation_Own="Combination_Union",
    N_OwnUnion_AllTF_Candidates=res$n_own_all,
    N_OwnUnion_Development_Candidates=res$n_own_dev,
    N_OwnUnion_Development_ChIP_Pos=res$n_own_dev_pos,
    OwnUnion_Development_Pos_Coverage_Within_All6=if (own_dev_pos_total>0)
      res$n_own_dev_pos/own_dev_pos_total else NA_real_,
    CV_Pooled_AUROC=if (complete_common) unname(pooled["AUROC"]) else NA_real_,
    CV_Pooled_AUPR=if (complete_common) unname(pooled["AUPR"]) else NA_real_,
    CV_Pooled_AUPR_Baseline=if (complete_common) unname(pooled["Baseline"]) else NA_real_,
    CV_Pooled_AUPR_Lift=if (complete_common) unname(pooled["Lift"]) else NA_real_,
    CV_Pooled_N_Pos=if (complete_common) as.integer(pooled["N_Pos"]) else NA_integer_,
    CV_Pooled_N_Unverified=if (complete_common) as.integer(pooled["N_Unverified"]) else NA_integer_,
    CV_OwnUnion_Pooled_AUROC=if (complete_own) unname(own_pooled["AUROC"]) else NA_real_,
    CV_OwnUnion_Pooled_AUPR=if (complete_own) unname(own_pooled["AUPR"]) else NA_real_,
    CV_OwnUnion_Pooled_AUPR_Baseline=if (complete_own) unname(own_pooled["Baseline"]) else NA_real_,
    CV_OwnUnion_Pooled_AUPR_Lift=if (complete_own) unname(own_pooled["Lift"]) else NA_real_,
    CV_OwnUnion_Pooled_N_Pos=if (complete_own) as.integer(own_pooled["N_Pos"]) else NA_integer_,
    CV_OwnUnion_Pooled_N_Unverified=if (complete_own) as.integer(own_pooled["N_Unverified"]) else NA_integer_,
    CV_OwnUnion_N_TF_Evaluable=if (nrow(own_tfm))
      sum(is.finite(own_tfm$AUROC)) else 0L,
    CV_Macro_AUROC=if (complete_common) mean(tfm$AUROC) else NA_real_,
    CV_Macro_AUPR_Lift=if (complete_common) mean(tfm$Lift) else NA_real_,
    CV_SD_TF_AUROC=if (complete_common) sd(tfm$AUROC) else NA_real_,
    CV_N_TF=if (complete_common) nrow(tfm) else uniqueN(tfm$TF),
    CV_TF_Above_05=if (complete_common) sum(tfm$AUROC>0.5) else NA_integer_,
    CV_AUROC_Mean=if (complete_common) mean(folds$AUROC) else NA_real_,
    CV_AUROC_SD=if (complete_common) sd(folds$AUROC) else NA_real_,
    CV_AUPR_Mean=if (complete_common) mean(folds$AUPR) else NA_real_,
    CV_AUPR_SD=if (complete_common) sd(folds$AUPR) else NA_real_,
    CV_Lift_Mean=if (complete_common) mean(folds$Lift) else NA_real_,
    CV_Baseline_Mean=mean(folds$Baseline),Valid_Folds=valid_folds,
    Status=if (complete_common &&
               (selection_universe=="Common_All6" || complete_own)) "OK" else "FAILED",
    OwnUnion_Status=if (complete_own) "OK" else "INCOMPLETE",
    Error=paste(unique(folds$Error[folds$Status!="OK"]),collapse=" | "),
    Split_Mode=split_mode)
  list(Fingerprint=fp$Fingerprint, Dataset=ds,Strategy=s,
       Combination=name,CombinationIndex=j,Summary=summary,Folds=folds,
       PerTF=tfm,OwnPerTF=own_tfm)
}


load_strategy_checkpoints <- function(ds,s,indices,comb_list,comb_names,data,plan,ds_seed,fp,ds_out) {
  cp_dir <- file.path(ds_out,"Checkpoints",s)
  dir.create(cp_dir,recursive=TRUE,showWarnings=FALSE)
  blocks <- vector("list",length(indices))
  for (i in seq_along(indices)) {
    j <- indices[i]
    path <- file.path(cp_dir,sprintf("Combo_%03d.rds",j))
    if (file.exists(path)) {
      ck <- tryCatch(readRDS(path),error=function(e)
        stop("Checkpoint damaged: ",path," — ",conditionMessage(e)))
      if (!identical(ck$Fingerprint,fp$Fingerprint) ||
          !identical(ck$Dataset,ds) || !identical(ck$Strategy,s) ||
          !identical(ck$Combination,comb_names[j]) ||
          !identical(as.integer(ck$CombinationIndex),as.integer(j)))
        stop("Checkpoint metadata mismatch: ",path)
      cat(sprintf(" [%s / %s] resume %d/%d: %s\n",ds,s,i,length(indices),comb_names[j]))
    } else {
      cat(sprintf(" [%s / %s] run %d/%d: %s\n",ds,s,i,length(indices),comb_names[j]))
      ck <- one_combination(ds,s,j,comb_list,comb_names,data,plan,ds_seed,fp)
      atomic_rds_new(ck,path)
    }
    blocks[[i]] <- ck
    # Optional human-readable, incremental progress within each strategy.
    sdir <- file.path(ds_out,"Strategies",s)
    progress <- rbindlist(lapply(blocks[seq_len(i)],`[[`,"Summary"),fill=TRUE)
    atomic_csv(progress,file.path(sdir,"CV_Progress_Partial.csv"))
  }
  list(
    Summary=rank_cv(rbindlist(lapply(blocks,`[[`,"Summary"),fill=TRUE)),
    Folds=rbindlist(lapply(blocks,`[[`,"Folds"),fill=TRUE),
    PerTF=rbindlist(lapply(blocks,`[[`,"PerTF"),fill=TRUE),
    OwnPerTF=rbindlist(lapply(blocks,`[[`,"OwnPerTF"),fill=TRUE),
    CompletedCombos=length(indices))
}

# As soon as one strategy finishes all combinations, persist its full results
# and selected best subset; then update the per-stress all-strategies reports.
write_strategy_cv <- function(ds,s,sc,ds_out,comb_names,ds_seed) {
  sdir <- file.path(ds_out,"Strategies",s)
  dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
  atomic_csv(sc$Summary,file.path(sdir,"CV_Subsets.csv"))
  atomic_csv(sc$Folds,file.path(sdir,"CV_PerFold.csv"))
  atomic_csv(sc$PerTF,file.path(sdir,"CV_PerTF.csv"))
  atomic_csv(sc$OwnPerTF,file.path(sdir,"CV_OwnUnion_PerTF.csv"))
  failures <- sc$Folds[Status!="OK"]
  atomic_csv(failures,file.path(sdir,"CV_Failures.csv"))
  eligible <- rank_cv(sc$Summary[Status=="OK" & Valid_Folds==cv_folds])
  if (!nrow(eligible)) {
    atomic_lines(paste0("NO VALID CV COMBINATION: ",ds," / ",s),
                 file.path(sdir,"NO_VALID_WINNER.txt"))
    warning(ds," / ",s,": no completely valid CV combination; next strategy continues")
    return(NULL)
  }
  winner <- copy(eligible[1L])
  winner[, Is_Global_CV_Winner:=FALSE]
  atomic_csv(winner,file.path(sdir,"Best_Combination_CV.csv"))
  atomic_csv(winner,file.path(ds_out,paste0(s,"_Best_Combination_CV.csv")))
  if (s=="Bagging") {
    cols <- match(strsplit(winner$Combination,"+",fixed=TRUE)[[1L]],method_names)
    j <- match(winner$Combination,comb_names)
    coeff <- bagging_column_coefficients(length(cols),bagging_B,ds_seed+j*37L)
    atomic_csv(data.table(Dataset=ds,Split_Mode=split_mode,Strategy=s,
                          Combination=winner$Combination,Method=method_names[cols],
                          Bootstrap_Selections=as.integer(round(coeff*length(cols)*bagging_B)),
                          Relative_Weight=coeff,B=bagging_B),
               file.path(sdir,"Bagging_Winner_Algorithm_Column_Frequencies.csv"))
  }
  atomic_new_marker("CV results and winner saved.",file.path(sdir,"CV_DONE.txt"))
  selected_auc <- if (selection_universe=="Common_All6")
    winner$CV_Pooled_AUROC else winner$CV_OwnUnion_Pooled_AUROC
  selected_pr <- if (selection_universe=="Common_All6")
    winner$CV_Pooled_AUPR else winner$CV_OwnUnion_Pooled_AUPR
  cat(sprintf(" [SAVED] %s / %s best=%s  selection=%s  AUROC=%.5f  AUPR=%.5f\n",
              ds,s,winner$Combination,selection_universe,selected_auc,selected_pr))
  winner
}

write_strategy_test <- function(ds,s,winner,data,plan,chip,ds_seed,ds_out,comb_names) {
  sdir <- file.path(ds_out,"Strategies",s)
  done <- file.path(sdir,"TEST_DONE.txt")
  result_path <- file.path(sdir,"Locked_Test.csv")
  if (file.exists(done) && file.exists(result_path)) {
    cat(" [RESUME] test already saved for ",ds," / ",s,"\n",sep="")
    return(fread(result_path))
  }
  if (file.exists(done) && !file.exists(result_path))
    stop("TEST_DONE marker exists but result is missing: ",sdir)
  cols <- match(strsplit(winner$Combination,"+",fixed=TRUE)[[1L]],method_names)
  if (anyNA(cols)) stop("Cannot parse winner methods: ",winner$Combination)
  idx <- match(s,strategies)
  out <- tryCatch({
    if (s %in% supervised) {
      own_flag <- method_union_flag(data,cols)
      own_train <- get_own_training(data,plan$outer$train,own_flag,ds_seed+500L)
      model <- fit_supervised(s,data,cols,own_train,ds_seed+9000L+idx)
      test_scores <- predict_own_on_common(model,data,plan$test,own_flag)
      rm(model)
    } else {
      j <- match(winner$Combination,comb_names)
      full_score <- unsup_scores(s,data,cols,ds_seed+j*37L)
      test_scores <- full_score[plan$test]
    }
    test_scores <- apply_common_eval_floor(
      test_scores, plan$test, method_union_flag(data,cols))
    test_res <- locked_result(ds,s,winner$Combination,test_scores,data,plan,chip,
                              ds_seed,ci_index=idx,result_type="CV_Selected_Fusion",
                              cv_auroc=winner$CV_Pooled_AUROC,cv_aupr=winner$CV_Pooled_AUPR,
                              global_winner=TRUE,cols=cols)
    test_res[,`:=`(Split_Mode=split_mode,
                   Evidence_Status=heldout_evidence_status)]
    atomic_csv(test_res,result_path)
    atomic_csv(test_res,file.path(ds_out,paste0(s,"_CVWinner_Locked_Test.csv")))
    if (save_test_predictions) {
      pred <- data$edges[plan$test,.(TF,Target,Label)]
      pred[,Score:=test_scores]
      atomic_csv(pred,file.path(sdir,"Locked_Test_Predictions.csv"))
      atomic_csv(per_tf_metric(pred$TF,pred$Label,pred$Score),
                 file.path(sdir,"Locked_Test_PerTF.csv"))
    }
    # Final network export is executed separately AFTER global model lock,
    # never as a side effect of strategy-wise holdout evaluation.
    atomic_new_marker("TF-disjoint holdout result; model locked by CV before evaluation.",done)
    test_res
  },error=function(e) {
    msg <- paste0(ds," / ",s,": ",conditionMessage(e))
    atomic_lines(msg,file.path(sdir,"TEST_ERROR.txt"))
    warning("[TEST FAILED] ",msg,"; CV results are intact; try again after fix")
    NULL
  })
  out
}

# Network export when holdout is intentionally skipped.
write_network_no_test <- function(ds,s,winner,data,plan,ds_seed,ds_out,comb_names) {
  if (!save_winner_networks) return(invisible(NULL))
  sdir <- file.path(ds_out,"Strategies",s)
  done <- file.path(sdir,"NETWORK_DONE.txt")
  if (file.exists(done)) return(invisible(NULL))
  cols <- match(strsplit(winner$Combination,"+",fixed=TRUE)[[1L]],method_names)
  if (s %in% supervised) {
    own_flag <- method_union_flag(data,cols)
    own_train <- get_own_training(data,plan$outer$train,own_flag,ds_seed+500L)
    mdl <- fit_supervised(s,data,cols,own_train,ds_seed+9000L+match(s,strategies))
    scores <- numeric(nrow(data$edges))
    own_rows <- which(own_flag)
    scores[own_rows] <- predict_in_batches(mdl,data,rows=own_rows)
    rm(mdl)
  } else {
    j <- match(winner$Combination,comb_names)
    scores <- unsup_scores(s,data,cols,ds_seed+j*37L)
  }
  # Export-only: exclude all unselected edges; membership is determined by union_flag,
  # not by whether a candidate's score equals zero.
  scores[!method_union_flag(data,cols)] <- 0
  exp <- export_network(ds,s,scores,data,sdir,cols=cols)
  atomic_csv(exp,file.path(sdir,"Network_Export_Summary.csv"))
  atomic_new_marker("Network exported.",done)
  invisible(exp)
}

# Incrementally assemble all finished strategy reports, so stopping before
# strategy 6 does not lose strategy 1..5 statistics.
rollup_stress <- function(ds_out) {
  read_parts <- function(relpath) {
    paths <- file.path(ds_out,"Strategies",strategies,relpath)
    hits <- paths[file.exists(paths)]
    if (!length(hits)) return(data.table())
    rbindlist(lapply(hits,fread),use.names=TRUE,fill=TRUE)
  }
  cvs <- read_parts("CV_Subsets.csv")
  folds <- read_parts("CV_PerFold.csv")
  tf <- read_parts("CV_PerTF.csv")
  tf_own <- read_parts("CV_OwnUnion_PerTF.csv")
  fails <- read_parts("CV_Failures.csv")
  winners <- read_parts("Best_Combination_CV.csv")
  tests <- read_parts("Locked_Test.csv")
  if (nrow(cvs)) atomic_csv(rank_cv(cvs),file.path(ds_out,"All_6Strategies_CV_Subsets.csv"))
  if (nrow(folds)) {
    atomic_csv(folds,file.path(ds_out,"CV_PerFold.csv"))
    atomic_csv(fails,file.path(ds_out,"CV_Failures.csv"))
    all6_name <- paste(method_names,collapse="+")
    base <- folds[Combination==all6_name,
                  .(Dataset,Strategy,Fold,AUROC_All6=AUROC,AUPR_All6=AUPR)]
    if (nrow(base)) {
      paired <- merge(folds,base,by=c("Dataset","Strategy","Fold"),all.x=TRUE)
      # The following are comparable because BOTH methods are scored on the SAME edges.
      paired[,c("AUROC_Diff_Vs_All6","AUPR_Diff_Vs_All6") :=
               list(AUROC-AUROC_All6,AUPR-AUPR_All6)]
      atomic_csv(paired,file.path(ds_out,"Paired_CV_Fold_Differences_Vs_All6.csv"))
    }
  }
  if (nrow(tf)) atomic_csv(tf,file.path(ds_out,"CV_PerTF.csv"))
  if (nrow(tf_own)) atomic_csv(tf_own,file.path(ds_out,"CV_OwnUnion_PerTF.csv"))
  if (nrow(winners)) {
    ranked <- rank_cv(winners)
    ranked[,Is_Global_CV_Winner := .I==1L]
    atomic_csv(ranked,file.path(ds_out,"CV_Selected_Winners.csv"))
    atomic_csv(ranked,file.path(ds_out,"SixStrategies_Best_Combinations_CV.csv"))
  }
  if (nrow(tests)) {
    atomic_csv(tests,file.path(ds_out,"Locked_Test_SixStrategies.csv"))
    if (nrow(winners)) {
      common <- c("Dataset","Strategy","Combination")
      fields <- c(common,"Test_AUROC","Test_AUPR","Test_AUPR_Lift",
                  "Test_Macro_AUROC","Test_Macro_AUPR_Lift")
      joined <- merge(winners,tests[,..fields],by=common,all.x=TRUE)
      atomic_csv(joined,file.path(ds_out,"SixStrategies_Best_CV_With_Locked_Test.csv"))
    }
  }
  progress <- data.table(Strategy=strategies,
                         CV_Done=file.exists(file.path(ds_out,"Strategies",strategies,"CV_DONE.txt")),
                         Test_Done=if (evaluate_locked_test)
                           file.exists(file.path(ds_out,"Strategies",strategies,"TEST_DONE.txt")) else NA,
                         Best_Exists=file.exists(file.path(ds_out,"Strategies",strategies,"Best_Combination_CV.csv")))
  atomic_csv(progress,file.path(ds_out,"Strategy_Progress_Status.csv"))
  invisible(progress)
}

# Update the summary across stresses and split modes after each strategy.
rollup_global <- function() {
  get_parts <- function(pathfun) {
    paths <- unlist(lapply(split_modes,function(mode)
      vapply(stress_names,function(ds) pathfun(mode,ds),character(1))),use.names=FALSE)
    hits <- paths[file.exists(paths)]
    if (!length(hits)) return(data.table())
    rbindlist(lapply(hits,fread),use.names=TRUE,fill=TRUE)
  }
  winners <- get_parts(function(mode,ds) file.path(root_base,mode,ds,"CV_Selected_Winners.csv"))
  tests <- get_parts(function(mode,ds) file.path(root_base,mode,ds,"Locked_Test_SixStrategies.csv"))
  if (nrow(winners)) {
    atomic_csv(winners,file.path(root_base,"Progress_AllModes_Best_Combinations.csv"))
    for (s in strategies) {
      atomic_csv(winners[Strategy==s],file.path(root_base,paste0("Progress_",s,"_Best_Combinations.csv")))
    }
  }
  if (nrow(tests)) atomic_csv(tests,file.path(root_base,"Progress_AllModes_Locked_Tests.csv"))
  invisible(TRUE)
}

write_baselines <- function(ds,ds_out,data,plan,chip,ds_seed) {
  done <- file.path(ds_out,"BASELINES_DONE.txt")
  if (file.exists(done)) return(invisible(TRUE))
  baseline_cv <- vector("list",length(method_names)+1L)
  baseline_cv_tf <- vector("list",length(method_names)+1L)
  test_rows <- list()
  oof <- unlist(lapply(plan$folds,function(f) f$valid),use.names=FALSE)
  for (b in seq_len(length(method_names)+1L)) {
    bn <- if (b<=length(method_names)) paste0("Single_",method_names[b]) else "Mean_All6"
    bs <- if (b<=length(method_names)) data$norm[,b] else rowMeans(data$norm)
    # Apply the SAME zero padding rule used by fusion models on the common universe.
    # For each single algorithm, the own-union scores are unchanged.
    if (b <= length(method_names)) {
      bs <- apply_common_eval_floor(
        bs, seq_len(nrow(data$edges)), method_union_flag(data, b))
    }
    bf <- rbindlist(lapply(plan$folds,function(f) {
      va <- f$valid; met <- metric(data$edges$Label[va],bs[va])
      data.table(Fold=f$fold,Pooled_AUROC=unname(met["AUROC"]),
                 Pooled_AUPR=unname(met["AUPR"]),Pooled_Lift=unname(met["Lift"]))
    }))
    bt <- rbindlist(lapply(plan$folds,function(f) {
      va <- f$valid
      z <- per_tf_metric(data$edges$TF[va],data$edges$Label[va],bs[va])
      z[,Fold:=f$fold]
      z
    }))
    bt[,`:=`(Dataset=ds,Baseline_Method=bn)]
    own_flag <- if (b<=length(method_names)) method_union_flag(data,b) else
      rep(TRUE,nrow(data$edges))
    own_oof <- filter_union_rows(oof,own_flag)
    own_pm <- metric_safe(data$edges$Label[own_oof],bs[own_oof])
    own_tf <- if (length(own_oof))
      per_tf_metric(data$edges$TF[own_oof],data$edges$Label[own_oof],bs[own_oof])
      else data.table()
    baseline_cv_tf[[b]] <- bt
    pm <- metric(data$edges$Label[oof],bs[oof])
    baseline_cv[[b]] <- data.table(Dataset=ds,Baseline=bn,
                                   Pooled_AUROC=unname(pm["AUROC"]),Pooled_AUPR=unname(pm["AUPR"]),
                                   Pooled_AUPR_Baseline=unname(pm["Baseline"]),
                                   Pooled_AUPR_Lift=unname(pm["Lift"]),
                                   Macro_AUROC=mean(bt$AUROC,na.rm=TRUE),
                                   Macro_AUPR_Lift=mean(bt$Lift,na.rm=TRUE),
                                   Mean_Pooled_AUROC=mean(bf$Pooled_AUROC,na.rm=TRUE),
                                   Mean_Pooled_AUPR=mean(bf$Pooled_AUPR,na.rm=TRUE),
                                   Mean_Pooled_Lift=mean(bf$Pooled_Lift,na.rm=TRUE),
                                   N_TF=sum(is.finite(bt$AUROC)),
                                   OwnUnion_Pooled_AUROC=unname(own_pm["AUROC"]),
                                   OwnUnion_Pooled_AUPR=unname(own_pm["AUPR"]),
                                   OwnUnion_Pooled_AUPR_Baseline=unname(own_pm["Baseline"]),
                                   OwnUnion_Pooled_AUPR_Lift=unname(own_pm["Lift"]),
                                   OwnUnion_N_Pos=as.integer(own_pm["N_Pos"]),
                                   OwnUnion_N_Unverified=as.integer(own_pm["N_Unverified"]),
                                   OwnUnion_N_Evaluable_TF=if (nrow(own_tf))
                                     sum(is.finite(own_tf$AUROC)) else 0L,
                                   OwnUnion_Dev_ChIP_Coverage_Within_All6=if (sum(data$edges$Label[oof])>0)
                                     sum(data$edges$Label[own_oof])/sum(data$edges$Label[oof]) else NA_real_,
                                   Training_Universe="Not_Trained_Original_Single_Method")
    if (evaluate_locked_test) {
      scores <- bs[plan$test]
      tr <- locked_result(ds,bn,
                          if (b<=length(method_names)) method_names[b] else paste(method_names,collapse="+"),
                          scores,data,plan,chip,ds_seed,ci_index=100L+b,
                          result_type="Fixed_Baseline",
                          cols=if (b<=length(method_names)) b else seq_along(method_names))
      tr[,`:=`(Split_Mode=split_mode,
               Evidence_Status=heldout_evidence_status)]
      test_rows[[b]] <- tr
    }
  }
  baseline_dt <- rbindlist(baseline_cv)
  atomic_csv(baseline_dt,file.path(ds_out,"CV_Fixed_Baselines.csv"))
  atomic_csv(baseline_dt[,.(Dataset,Baseline,OwnUnion_Pooled_AUROC,
    OwnUnion_Pooled_AUPR,OwnUnion_Pooled_AUPR_Baseline,OwnUnion_Pooled_AUPR_Lift,
    OwnUnion_N_Pos,OwnUnion_N_Unverified,OwnUnion_N_Evaluable_TF,
    OwnUnion_Dev_ChIP_Coverage_Within_All6)],file.path(ds_out,"CV_Fixed_Baselines_OwnUnion.csv"))
  atomic_csv(rbindlist(baseline_cv_tf),file.path(ds_out,"CV_Fixed_Baselines_PerTF.csv"))
  if (evaluate_locked_test) {
    test_dt <- rbindlist(test_rows,fill=TRUE)
    atomic_csv(test_dt,file.path(ds_out,"Locked_Test_Fixed_Baselines.csv"))
    atomic_csv(test_dt[,.(Dataset,Strategy,Combination,
      OwnUnion_Test_AUROC,OwnUnion_Test_AUPR,OwnUnion_Test_AUPR_Baseline,
      OwnUnion_Test_AUPR_Lift,OwnUnion_Test_N_Pos,OwnUnion_Test_N_Unverified,
      OwnUnion_Test_ChIP_Coverage)],
      file.path(ds_out,"Locked_Test_Fixed_Baselines_OwnUnion.csv"))
    fus <- file.path(ds_out,"Locked_Test_SixStrategies.csv")
    if (file.exists(fus)) atomic_csv(rbindlist(list(fread(fus),test_dt),
                                               use.names=TRUE,fill=TRUE),file.path(ds_out,"Locked_Test_Fusions_And_Baselines.csv"))
  }
  atomic_new_marker("Fixed baselines completed.",done)
}

# ============================================================================
# Development CV diagnostics and model lock
# ============================================================================
# Randomized TF-wise folds for repeated, post-selection stability analysis.
# Positives stratify TF groups; an entire TF remains in exactly one fold.
make_tf_folds_replicate <- function(edges, tf_train, k, local_seed) {
  stat <- edges[TF %in% tf_train, .(N_Pos=sum(Label==1L)), by=TF]
  stopifnot(setequal(stat$TF, tf_train), nrow(stat) >= k)
  set.seed(local_seed)
  stat <- stat[sample.int(nrow(stat))]
  setorder(stat, -N_Pos)
  stat[, Stratum := cut(seq_len(.N), breaks=min(3L, .N),
                        labels=FALSE, include.lowest=TRUE)]
  ans <- integer(nrow(stat))
  # randomized round-robin, balanced TF counts and stratum representation
  for (g in unique(stat$Stratum)) {
    ids <- which(stat$Stratum == g)
    ids <- ids[sample.int(length(ids))]
    fold_cycle <- sample(seq_len(k), k)
    ans[ids] <- rep(fold_cycle, length.out=length(ids))
  }
  if (!setequal(unique(ans),seq_len(k))) {
    # Extremely small TF sets: fallback to random balanced assignment.
    ans <- sample(rep(seq_len(k),length.out=nrow(stat)))
  }
  names(ans) <- stat$TF
  stopifnot(!anyDuplicated(names(ans)), setequal(unique(ans),seq_len(k)))
  ans
}

reviewer_oof_repeat <- function(ds, s, cols, data, plan, ds_seed, repeat_id,
                                comb_idx, keep_edge=FALSE) {
  if (repeat_id == 1L) {
    fold_list <- plan$folds
    fold_map <- plan$fold_map
  } else {
    fold_map <- make_tf_folds_replicate(data$edges,plan$outer$train,
                                        cv_folds,ds_seed + 30000L + repeat_id*101L)
    fold_list <- lapply(seq_len(cv_folds),function(k) {
      valid_tfs <- names(fold_map)[fold_map == k]
      tr <- get_indices(data$edges,setdiff(plan$outer$train,valid_tfs),
                        train_unlabeled_ratio,ds_seed+31000L+repeat_id*100L+k)
      va <- get_indices(data$edges,valid_tfs,cv_unlabeled_ratio,
                        ds_seed+32000L+repeat_id*100L+k)
      list(fold=k,train=tr,valid=va,valid_tfs=valid_tfs)
    })
  }
  own_flag <- method_union_flag(data,cols)
  unsup <- if (!s %in% supervised)
    unsup_scores(s,data,cols,ds_seed+comb_idx*37L) else NULL
  preds <- vector("list",cv_folds)
  folds <- vector("list",cv_folds)
  for (fold in fold_list) {
    va <- fold$valid
    if (s %in% supervised) {
      own_train <- get_own_training(
        data,setdiff(plan$outer$train,fold$valid_tfs),own_flag,
        if (repeat_id==1L) ds_seed+10L+fold$fold
        else ds_seed+31000L+repeat_id*100L+fold$fold)
      model <- fit_supervised(s,data,cols,own_train,
                              ds_seed+comb_idx*100L+fold$fold+if (repeat_id==1L) 0L else repeat_id*100000L)
      p <- predict_own_on_common(model,data,va,own_flag)
      rm(model)
    } else {
      p <- apply_common_eval_floor(unsup[va], va, own_flag)
    }
    y <- data$edges$Label[va]
    met <- metric(y,p)
    oi <- which(own_flag[va])
    own_met <- metric_safe(y[oi],p[oi])
    folds[[fold$fold]] <- data.table(Dataset=ds,Strategy=s,Repeat=repeat_id,
                                     Fold=fold$fold,AUROC=unname(met["AUROC"]),AUPR=unname(met["AUPR"]),
                                     Baseline=unname(met["Baseline"]),Lift=unname(met["Lift"]),
                                     N_TF=length(fold$valid_tfs),N_Pos=as.integer(met["N_Pos"]),
                                     N_Unverified=as.integer(met["N_Unverified"]),
                                     OwnUnion_AUROC=unname(own_met["AUROC"]),
                                     OwnUnion_AUPR=unname(own_met["AUPR"]),
                                     OwnUnion_Baseline=unname(own_met["Baseline"]),
                                     OwnUnion_N_Pos=as.integer(own_met["N_Pos"]),
                                     OwnUnion_N_Unverified=as.integer(own_met["N_Unverified"]))
    z <- data$edges[va,.(TF,Target,Label)]
    z[,`:=`(Score=as.numeric(p),Fold=fold$fold,OwnUnion=own_flag[va])]
    preds[[fold$fold]] <- z
  }
  joined <- rbindlist(preds)
  if (anyDuplicated(joined,by=c("TF","Target")) ||
      !setequal(unique(joined$TF),plan$outer$train))
    stop("OOF validation must predict every development TF edge exactly once")
  pm <- metric(joined$Label,joined$Score)
  own_joined <- joined[OwnUnion==TRUE]
  own_pm <- metric_safe(own_joined$Label,own_joined$Score)
  out <- data.table(Dataset=ds,Strategy=s,Repeat=repeat_id,
                    Pooled_AUROC=unname(pm["AUROC"]),Pooled_AUPR=unname(pm["AUPR"]),
                    AUPR_Baseline=unname(pm["Baseline"]),AUPR_Lift=unname(pm["Lift"]),
                    N_Pos=as.integer(pm["N_Pos"]),N_Unverified=as.integer(pm["N_Unverified"]),
                    OwnUnion_Pooled_AUROC=unname(own_pm["AUROC"]),
                    OwnUnion_Pooled_AUPR=unname(own_pm["AUPR"]),
                    OwnUnion_AUPR_Baseline=unname(own_pm["Baseline"]),
                    OwnUnion_AUPR_Lift=unname(own_pm["Lift"]),
                    OwnUnion_N_Pos=as.integer(own_pm["N_Pos"]),
                    OwnUnion_N_Unverified=as.integer(own_pm["N_Unverified"]))
  tfm <- per_tf_metric(joined$TF,joined$Label,joined$Score)
  tfm[,`:=`(Repeat=repeat_id,Strategy=s,Dataset=ds)]
  list(summary=out,folds=rbindlist(folds),per_tf=tfm,
       per_edge=if (keep_edge) joined else NULL,
       foldmap=data.table(TF=names(fold_map),Fold=as.integer(fold_map),
                          Repeat=repeat_id,Strategy=s,Dataset=ds))
}

# Paired TF-cluster bootstrap on identical TF-target edges; uncertainty is
# conditional on these fitted OOF predictions and on the 17–22 development TFs.
paired_tf_boot <- function(a,b,reps,local_seed) {
  stopifnot(identical(a$TF,b$TF),identical(a$Target,b$Target),
            identical(a$Label,b$Label))
  ma <- metric(a$Label,a$Score); mb <- metric(b$Label,b$Score)
  delta <- unname(ma[c("AUROC","AUPR")] - mb[c("AUROC","AUPR")])
  names(delta) <- c("AUROC","AUPR")
  if (reps < 2L) return(data.table(Metric=names(delta),Delta=delta,
                                   Lower95=NA_real_,Upper95=NA_real_,N_Boot=0L))
  groups <- split(seq_len(nrow(a)), a$TF)
  set.seed(local_seed)
  z <- matrix(NA_real_,nrow=reps,ncol=2L)
  for (i in seq_len(reps)) {
    idx <- unlist(groups[sample.int(length(groups),length(groups),replace=TRUE)],
                  use.names=FALSE)
    m1 <- metric(a$Label[idx],a$Score[idx]); m2 <- metric(b$Label[idx],b$Score[idx])
    z[i,] <- unname(m1[c("AUROC","AUPR")] - m2[c("AUROC","AUPR")])
  }
  data.table(Metric=c("AUROC","AUPR"),Delta=delta,
             Lower95=apply(z,2,quantile,probs=.025,na.rm=TRUE),
             Upper95=apply(z,2,quantile,probs=.975,na.rm=TRUE),
             N_Boot=colSums(is.finite(z)))
}

# Fast paired TF-level inference. Repeated CV partitions share the SAME TFs:
# average the repeated OOF metrics within TF BEFORE statistical testing.
# Never treat 10*5 folds or millions of edges as 50 independent observations.
# H1 (one-sided): mean_TF(metric_XGBoost - metric_other) > 0.
# The sign-flip test assumes TF-level paired differences are sign-exchangeable
# under H0; conclusions remain POST-SELECTION EXPLORATORY for development CV.
paired_tf_inference <- function(per_tf_average, ds, seed_local, context="Development_CV_Postselection") {
  reference <- per_tf_average[Strategy==superiority_reference]
  stopifnot(nrow(reference)>0L, !anyDuplicated(reference$TF))
  others <- setdiff(strategies,superiority_reference)
  res <- list()
  for (i in seq_along(others)) {
    other <- others[i]
    x <- merge(reference,per_tf_average[Strategy==other],by="TF",suffixes=c("_XGB","_Other"))
    if (nrow(x)!=nrow(reference)) stop("Missing TF in paired comparison: ",ds," / ",other)
    for (metric_name in c("AUROC","AUPR")) {
      differences <- x[[paste0(metric_name,"_XGB")]]-x[[paste0(metric_name,"_Other")]]
      differences <- as.numeric(differences)
      differences <- differences[is.finite(differences)]
      n <- length(differences)
      if (n<3L) stop("Too few paired TFs for inference: ",ds," / ",other)
      observed <- mean(differences)
      set.seed(seed_local+1000L*i+if (metric_name=="AUROC") 1L else 2L)
      boot_mean <- replicate(superiority_bootstrap_reps,
                             mean(sample(differences,n,replace=TRUE)))
      ci <- as.numeric(quantile(boot_mean,c(.025,.975),names=FALSE))
      if (n <= 15L) {
        # Enumerate ALL 2^n TF-level sign flips for small locked test sets.
        # With n=5, min one-sided exact P=1/32; Holm over five comparisons
        # cannot produce P_adjusted<0.05. Do not claim significance.
        exact_signs <- as.matrix(expand.grid(rep(list(c(-1,1)),n)))
        null_stats <- as.numeric(exact_signs %*% differences)/n
        raw_p <- mean(null_stats >= observed-1e-12)
        permutations_used <- length(null_stats)
        permutation_type <- "Exact"
      } else {
        signs <- matrix(sample(c(-1,1), n*superiority_permutation_reps, replace=TRUE),
                        nrow=superiority_permutation_reps,ncol=n)
        null_stats <- as.numeric(signs %*% differences)/n
        raw_p <- (1+sum(null_stats >= observed-1e-12))/(superiority_permutation_reps+1)
        permutations_used <- superiority_permutation_reps
        permutation_type <- "MonteCarlo"
      }
      res[[length(res)+1L]] <- data.table(
        Dataset=ds, Evidence_Context=context, Reference=superiority_reference,
        Compared_To=other,Metric=metric_name,N_TFs=n,
        XGBoost_MeanTF=mean(x[[paste0(metric_name,"_XGB")]],na.rm=TRUE),
        Other_MeanTF=mean(x[[paste0(metric_name,"_Other")]],na.rm=TRUE),
        Delta_MeanTF=observed,Lower95=ci[1L],Upper95=ci[2L],
        P_OneSided_Raw=raw_p,Bootstrap_Reps=superiority_bootstrap_reps,
        Permutation_Reps=permutations_used,Permutation_Type=permutation_type)
    }
  }
  z <- rbindlist(res)
  # Within each stress, correct five XGBoost-vs-other comparisons PER metric.
  z[,P_Holm:=p.adjust(P_OneSided_Raw,method="holm"),by=.(Dataset,Metric)]
  z[,P_BH:=p.adjust(P_OneSided_Raw,method="BH"),by=.(Dataset,Metric)]
  z[,Supported_After_Holm:=Delta_MeanTF>0 & P_Holm<superiority_alpha]
  z[,CI_Note:="Nominal paired TF-bootstrap CI; not a simultaneous adjusted CI"]
  z[,Inference_Limitation:=if (grepl("Development",context))
    "Post-selection exploratory; model subsets selected on development CV"
    else "Locked TF-holdout; number of independent TF clusters limits significance"]
  setorder(z,Metric,P_Holm,Compared_To)
  z
}

# Fast 95% percentile bootstrap on TF-level mean metrics of ten repeated CVs.
tf_macro_bootstrap_ci <- function(per_tf_average, ds, seed_local) {
  z <- list()
  for (i in seq_along(strategies)) {
    s <- strategies[i]
    d <- per_tf_average[Strategy==s]
    for (metric_name in c("AUROC","AUPR")) {
      vals <- d[[metric_name]]
      vals <- vals[is.finite(vals)]
      if (length(vals)<3L) stop("Too few TFs for TF-macro CI: ", ds, " / ", s)
      set.seed(seed_local+10000L*i+if (metric_name=="AUROC") 1L else 2L)
      bt <- replicate(reviewer_ci_reps,mean(sample(vals,length(vals),replace=TRUE)))
      ci <- quantile(bt,c(.025,.975),names=FALSE)
      z[[length(z)+1L]] <- data.table(
        Dataset=ds,Strategy=s,Metric=metric_name,N_TFs=length(vals),
        TF_Macro_Mean=mean(vals),Lower95=ci[1],Upper95=ci[2],
        Bootstrap_Reps=reviewer_ci_reps,
        Interpretation="TF-macro repeated-CV; each TF appears once after averaging 10 OOF repeats; post-selection")
    }
  }
  rbindlist(z)
}

# ============================================================================
# TF-cluster 95% confidence intervals for pooled fivefold OOF AUROC and AUPR.
#
# This is the SAME statistic and the SAME common evaluation universe used in
# development model selection. The older TF-macro CI is supplementary only.
# Resample WHOLE development TFs with replacement; keep every TF-target edge
# belonging to a sampled TF (including multiplicity). On each resample,
# re-calculate pooled AUROC and PRROC AUPR from the concatenated OOF scores.
# Each strategy uses repeat 1, the exact fivefold split used to select its
# winning subset; the 10 repeats are reported separately for stability.
# CIs are post-selection / conditional on frozen OOF predictions; they are NOT
# nested-CV or model-search adjusted CIs. NO holdout predictions are read.
# ============================================================================
reviewer_pooled_oof_ci <- function(ds, ds_out, ds_seed, winners=NULL) {
  outfile <- file.path(ds_out,"Reviewer_CV_Pooled_OOF_TFCluster_95CI.csv")
  outlong <- file.path(ds_out,"Reviewer_CV_Pooled_OOF_TFCluster_95CI_Long.csv")
  outmerged <- file.path(ds_out,"CV_Selected_Winners_with_Pooled95CI.csv")
  if (is.null(winners)) winners <- fread(file.path(ds_out,"CV_Selected_Winners.csv"))
  if (nrow(winners)!=length(strategies) || !setequal(winners$Strategy,strategies))
    stop("Incomplete six fusion winners for pooled OOF CI: ",ds)

  # A full existing output is reusable. Incomplete/stale outputs fail rather
  # than getting silently treated as valid CIs.
  if (file.exists(outfile) && file.exists(outlong) && file.exists(outmerged)) {
    old <- fread(outfile)
    required_ci_cols <- c("Strategy","Combination","Bootstrap_Reps",
                          "Bootstrap_Valid_Reps","AUROC_L95","AUROC_U95",
                          "AUPR_L95","AUPR_U95")
    idok <- all(required_ci_cols %in% names(old)) &&
      nrow(old)==length(strategies) &&
      setequal(old$Strategy,strategies) &&
      all(old$Bootstrap_Reps==reviewer_ci_reps) &&
      all(is.finite(old$AUROC_L95) & is.finite(old$AUROC_U95) &
          is.finite(old$AUPR_L95) & is.finite(old$AUPR_U95)) &&
      all(old$Bootstrap_Valid_Reps >= ceiling(.5*reviewer_ci_reps))
    matched <- merge(old[,.(Strategy,Combination)],
                     winners[,.(Strategy,Combination)],by="Strategy",
                     suffixes=c("_CI","_Winner"))
    if (isTRUE(idok) && nrow(matched)==length(strategies) &&
        all(matched$Combination_CI==matched$Combination_Winner)) {
      message(" [RESUME] pooled OOF TF-cluster 95% CI already complete: ",ds)
      return(invisible(old))
    }
    stop("Existing pooled CI output is incompatible with current winners/config: ",outfile,
         "; do not mix different CI protocols in a single results directory")
  }
  result <- vector("list",nrow(winners))
  for (i in seq_len(nrow(winners))) {
    w <- winners[i]
    s <- w$Strategy[1L]
    p <- file.path(ds_out,"Strategies",s,"CV_Winner_OOF_Predictions.csv")
    if (!file.exists(p))
      stop("Missing saved repeat-1 winner OOF predictions: ",p,
           "; complete the development CV diagnostics before calculating CIs")
    oof <- fread(p,select=c("TF","Target","Label","Score"),showProgress=FALSE)
    if (!nrow(oof) || anyNA(oof$TF) || anyNA(oof$Target) ||
        anyNA(oof$Label) || any(!oof$Label %in% c(0L,1L)) ||
        any(!is.finite(oof$Score)) ||
        anyDuplicated(oof,by=c("TF","Target")))
      stop("Invalid repeat-1 OOF predictions for ",ds," / ",s)
    n_tf <- uniqueN(oof$TF)
    if (n_tf < bootstrap_min_tfs)
      stop("Too few development TF clusters for pooled CI: ",ds," / ",s)
    met <- metric(oof$Label,oof$Score)
    if (!is.finite(met["AUROC"]) || !is.finite(met["AUPR"]) ||
        abs(unname(met["AUROC"])-w$CV_Pooled_AUROC[1L])>1e-8 ||
        abs(unname(met["AUPR"])-w$CV_Pooled_AUPR[1L])>1e-8)
      stop("POOLED OOF metric mismatch with DEVELOPMENT model-selection metric: ",
           ds," / ",s,". Do not draw mismatched 95% CI.")
    message(sprintf(" [POOLED CI] %s / %s: %s TFs, %s OOF edges, %s bootstrap resamples",
                    ds,s,n_tf,format(nrow(oof),big.mark=","),reviewer_ci_reps))
    ci <- cluster_ci(oof$TF,oof$Label,oof$Score,
                     reviewer_ci_reps, ds_seed+200000L+1009L*i)
    if (unname(ci["Valid_Boot"]) < ceiling(.5*reviewer_ci_reps) ||
        any(!is.finite(ci[c("AUROC_L95","AUROC_U95","AUPR_L95","AUPR_U95")])))
      stop("Insufficient valid TF-cluster bootstrap replicates: ",ds," / ",s)
    result[[i]] <- data.table(
      Dataset=ds,Strategy=s,Combination=w$Combination[1L],
      Evaluation_Universe="Common_All6_Top10",
      OOF_Repeat=1L,N_TF=n_tf,N_Pos=as.integer(met["N_Pos"]),
      N_Unverified=as.integer(met["N_Unverified"]),
      CV_Pooled_AUROC=unname(met["AUROC"]),
      AUROC_L95=unname(ci["AUROC_L95"]),AUROC_U95=unname(ci["AUROC_U95"]),
      CV_Pooled_AUPR=unname(met["AUPR"]),
      AUPR_L95=unname(ci["AUPR_L95"]),AUPR_U95=unname(ci["AUPR_U95"]),
      CV_Pooled_AUPR_Baseline=unname(met["Baseline"]),
      CV_Pooled_AUPR_Lift=unname(met["Lift"]),
      Bootstrap_Reps=as.integer(reviewer_ci_reps),
      Bootstrap_Valid_Reps=as.integer(ci["Valid_Boot"]),
      CI_Type="TF-cluster percentile bootstrap on fivefold pooled OOF predictions",
      CI_Interpretation="Descriptive after subset selection; fixed OOF predictions; not nested or selection-adjusted"
    )
    rm(oof); gc(verbose=FALSE)
  }
  ans <- rbindlist(result)
  setorder(ans,Strategy)
  lng <- rbindlist(list(
    ans[,.(Dataset,Strategy,Combination,Metric="AUROC",
            Estimate=CV_Pooled_AUROC,Lower95=AUROC_L95,Upper95=AUROC_U95,
            N_TF,OOF_Repeat,Evaluation_Universe,Bootstrap_Reps,Bootstrap_Valid_Reps)],
    ans[,.(Dataset,Strategy,Combination,Metric="AUPR",
            Estimate=CV_Pooled_AUPR,Lower95=AUPR_L95,Upper95=AUPR_U95,
            N_TF,OOF_Repeat,Evaluation_Universe,Bootstrap_Reps,Bootstrap_Valid_Reps)]
  ),use.names=TRUE)
  setorder(lng,Metric,Strategy)
  keep_ci <- c("Strategy","Combination","AUROC_L95","AUROC_U95",
               "AUPR_L95","AUPR_U95","Bootstrap_Reps","Bootstrap_Valid_Reps")
  merged <- merge(winners,ans[,..keep_ci],by=c("Strategy","Combination"),all.x=TRUE)
  if (nrow(merged)!=length(strategies) || anyNA(merged$AUROC_L95))
    stop("Pooled CI cannot be joined to CV winners for ",ds)
  atomic_csv(ans,outfile)
  atomic_csv(lng,outlong)
  atomic_csv(merged,outmerged)
  invisible(ans)
}

reviewer_development_stats <- function(ds,data,plan,ds_seed,ds_out,comb_names) {
  done <- file.path(ds_out,"REVIEWER_CV_DIAGNOSTICS_DONE.txt")
  if (file.exists(done)) {
    # An older completed run may have only TF-macro CIs. Reuse its saved OOF
    # predictions and add the correctly matched pooled CIs without refitting.
    reviewer_pooled_oof_ci(ds,ds_out,ds_seed)
    return(invisible(TRUE))
  }
  winners <- fread(file.path(ds_out,"CV_Selected_Winners.csv"))
  stopifnot(nrow(winners)==length(strategies),
            all(winners$Status=="OK"),all(!is.na(winners$CV_Pooled_AUROC)))
  summaries <- list(); per_folds <- list(); all_maps <- list(); per_tf_all <- list()
  for (i in seq_len(nrow(winners))) {
    s <- winners$Strategy[i]
    cols <- match(strsplit(winners$Combination[i],"+",fixed=TRUE)[[1L]],method_names)
    ci <- match(winners$Combination[i],comb_names)
    stopifnot(!anyNA(cols),!is.na(ci))
    sdir <- file.path(ds_out,"Strategies",s)
    for (r in seq_len(if (reviewer_repeat_enabled) reviewer_cv_repeats else 1L)) {
      cat(sprintf(" [DEVELOPMENT CV] %s / %s / repeat %d\n",ds,s,r))
      o <- reviewer_oof_repeat(ds,s,cols,data,plan,ds_seed,r,ci,keep_edge=r==1L)
      summaries[[length(summaries)+1L]] <- o$summary
      per_folds[[length(per_folds)+1L]] <- o$folds
      all_maps[[length(all_maps)+1L]] <- o$foldmap
      per_tf_all[[length(per_tf_all)+1L]] <- o$per_tf
      if (r==1L) {
        setorderv(o$per_edge,c("TF","Target"))
        # Keep 1st-fold OOF data on disk for audit; avoid 2000 expensive
        # TF-cluster bootstrap rescoring passes through millions of edges.
        atomic_csv(o$per_edge,file.path(sdir,"CV_Winner_OOF_Predictions.csv"))
        atomic_csv(o$per_tf,file.path(sdir,"CV_Selected_Winner_PerTF.csv"))
        if (abs(o$summary$Pooled_AUROC-winners$CV_Pooled_AUROC[i])>1e-8 ||
            abs(o$summary$Pooled_AUPR-winners$CV_Pooled_AUPR[i])>1e-8)
          stop("CV winner OOF metric mismatch: ",s)
      }
    }
  }
  pooled <- rbindlist(summaries)
  tfreps <- rbindlist(per_tf_all,fill=TRUE)
  atomic_csv(pooled,file.path(ds_out,"Reviewer_CV_Repeats_Pooled.csv"))
  atomic_csv(rbindlist(per_folds),file.path(ds_out,"Reviewer_CV_Repeats_PerFold.csv"))
  atomic_csv(rbindlist(all_maps),file.path(ds_out,"Reviewer_CV_Repeated_TF_Assignments.csv"))
  atomic_csv(tfreps,file.path(ds_out,"Reviewer_CV_Repeats_PerTF.csv"))
  summarize_metric <- function(v) c(Mean=mean(v),SD=sd(v),Min=min(v),Max=max(v))
  pooled_summary <- pooled[,as.list(c(summarize_metric(Pooled_AUROC),
                                      setNames(summarize_metric(Pooled_AUPR),
                                               paste0("AUPR_",names(summarize_metric(Pooled_AUPR)))))),
                           by=.(Dataset,Strategy)]
  atomic_csv(pooled_summary,file.path(ds_out,"Reviewer_CV_Repeats_Descriptive_Summary.csv"))
  # Each TF is a statistical unit. Repeats of one TF are aggregated first.
  tfavg <- tfreps[,.(AUROC=mean(AUROC),AUPR=mean(AUPR),
                     Lift=mean(Lift)),by=.(Dataset,Strategy,TF)]
  atomic_csv(tfavg,file.path(ds_out,"Reviewer_CV_Repeated_TF_Mean.csv"))
  ci <- tf_macro_bootstrap_ci(tfavg,ds,ds_seed+200000L)
  atomic_csv(ci,file.path(ds_out,"Reviewer_CV_Repeated_TFMacro_Bootstrap_CI.csv"))
  compare <- paired_tf_inference(tfavg,ds,ds_seed+500000L)
  atomic_csv(compare,file.path(ds_out,"Reviewer_XGBoost_vs_Five_Holm_Development.csv"))
  # Useful explicit summary; do not infer victory from absence of a P-value alone.
  proof <- compare[,.(N_Comparisons=.N, N_Positive=sum(Delta_MeanTF>0),
                      N_Holm_Supported=sum(Supported_After_Holm)),by=.(Dataset,Metric)]
  atomic_csv(proof,file.path(ds_out,"Reviewer_XGBoost_Superiority_Check.csv"))
  # IMPORTANT: the MAIN fivefold CV AUROC/AUPR CIs must match pooled OOF
  # selection metrics. TF-macro repeated-CV CIs above remain supplementary.
  reviewer_pooled_oof_ci(ds,ds_out,ds_seed,winners=winners)
  atomic_new_marker("Main pooled OOF TF-cluster CI; supplementary 10x5 TF-macro stability; no holdout labels.",done)
  invisible(TRUE)
}

# Freeze the full CV-selected strategy + feature subset BEFORE held-out test.
freeze_global_winner <- function(ds,ds_out) {
  file <- file.path(ds_out,"CV_Selected_Global_Model_LOCKED.csv")
  z <- fread(file.path(ds_out,"CV_Selected_Winners.csv"))
  z <- rank_cv(z[Status=="OK" & Valid_Folds==cv_folds])
  if (nrow(z)!=length(strategies)) stop("Not all strategy winners are available: ",ds)
  selected <- copy(z[1L]); selected[,`:=`(
    Selection_Source=paste0("Development_TF_Disjoint_FiveFold_",selection_universe,"_",selection_metric),
    Holdout_Used_For_Selection=FALSE,
    Fixed_Input_Stage=input_stage,
    Holdout_Evidence=heldout_evidence_status)]
  if (file.exists(file)) {
    old <- fread(file)
    # Treat matching NA values as equal; an NA comparison must not cause if(NA).
    same_metric <- function(a,b) isTRUE(all.equal(a,b,tolerance=1e-10,
                                                  check.attributes=FALSE))
    if (!identical(old$Strategy,selected$Strategy) ||
        !identical(old$Combination,selected$Combination) ||
        !same_metric(old$CV_Pooled_AUROC,selected$CV_Pooled_AUROC) ||
        !same_metric(old$CV_OwnUnion_Pooled_AUROC,selected$CV_OwnUnion_Pooled_AUROC))
      stop("Global model lock differs from current CV; refusing test: ",file)
  } else atomic_csv(selected,file)
  selected
}


# External evidence may provide an ADDITIONAL, non-independent-of-transcriptome
# verification. It is never used to select strategy/params or train weights.
# Shared ChIP edges are excluded to avoid claiming a second independent source.
external_evidence_validation <- function(ds,score,data,chip,ds_out) {
  if (is.null(external_validation_path) ||
      !is.character(external_validation_path) ||
      !nzchar(external_validation_path)) return(invisible(NULL))
  if (!file.exists(external_validation_path)) stop("Missing external validation file: ",external_validation_path)
  z <- fread(external_validation_path,showProgress=FALSE)
  if (ncol(z)<2L) stop("External validation needs at least 2 columns TF/Target")
  vals <- unique(data.table(TF=toupper(trimws(as.character(z[[external_tf_col]]))),
                            Target=toupper(trimws(as.character(z[[external_target_col]])))))
  vals <- vals[!is.na(TF)&!is.na(Target)&TF!=""&Target!=""&TF!=Target]
  n0 <- nrow(vals)
  vals <- vals[!chip,on=.(TF,Target)]
  excluded <- n0 - nrow(vals)
  if (!nrow(vals)) stop("No external validation edges remain after ChIP-overlap exclusion")
  scored <- data$edges[,.(TF,Target)]
  scored[,Score:=as.numeric(score)]
  scored[,Ext_Label:=0L]
  scored[vals,on=.(TF,Target),Ext_Label:=1L]
  m <- metric(scored$Ext_Label,scored$Score)
  stats <- data.table(Dataset=ds,External_File=external_validation_path,
                      External_N_Input=n0,External_ChIP_Overlap_Excluded=excluded,
                      External_N_Remaining=nrow(vals),N_Candidate_Matched=sum(scored$Ext_Label==1L),
                      N_Candidate_Unverified=sum(scored$Ext_Label==0L),
                      External_AUROC=unname(m["AUROC"]),External_AUPR=unname(m["AUPR"]),
                      External_AUPR_Baseline=unname(m["Baseline"]),
                      External_AUPR_Lift=unname(m["Lift"]),
                      Caveat="PU proxy; verify truly independent experiment provenance; never used for model selection")
  atomic_csv(stats,file.path(ds_out,"External_Validation_Independent_Resource.csv"))
  stats
}

# One global winner, not the best of the six test results.
finalize_global_model <- function(ds,data,plan,chip,ds_seed,ds_out,comb_names) {
  selected <- freeze_global_winner(ds,ds_out)
  strategy <- selected$Strategy[1L]
  cols <- match(strsplit(selected$Combination[1L],"+",fixed=TRUE)[[1L]],method_names)
  sdir <- file.path(ds_out,"Strategies",strategy)
  if (evaluate_locked_test) {
    final_done <- file.path(ds_out,"FINAL_EVAL_DONE.txt")
    if (!file.exists(final_done)) {
      write_strategy_test(ds,strategy,selected,data,plan,chip,ds_seed,ds_out,comb_names)
      path <- file.path(sdir,"Locked_Test.csv")
      if (!file.exists(path)) stop("Locked global winner test failed: ",ds)
      rs <- fread(path)
      rs[,`:=`(Is_Global_CV_Winner=TRUE,
               Model_Locked_Before_Test=TRUE,
               Evidence_Status=heldout_evidence_status)]
      atomic_csv(rs,file.path(ds_out,"CV_LOCKED_GlobalWinner_Heldout.csv"))
      atomic_new_marker("Only global CV winner evaluated; no test-based selection.",final_done)
    }
  }
  if (save_winner_networks) {
    final_done <- file.path(ds_out,"FINAL_NETWORK_DONE.txt")
    if (!file.exists(final_done)) {
      # Supervised fit uses DEVELOPMENT TRAINING TF LABELS ONLY.
      # Any ChIP edges from held-out TFs are never fitted here.
      if (strategy %in% supervised) {
        own_flag <- method_union_flag(data,cols)
        own_train <- get_own_training(data,plan$outer$train,own_flag,ds_seed+500L)
        model <- fit_supervised(strategy,data,cols,own_train,
                                ds_seed+9000L+match(strategy,strategies))
        score <- numeric(nrow(data$edges))
        own_rows <- which(own_flag)
        score[own_rows] <- predict_in_batches(model,data,rows=own_rows)
        rm(model)
      } else {
        j <- match(selected$Combination,comb_names)
        score <- unsup_scores(strategy,data,cols,ds_seed+j*37L)
      }
      export_dir <- file.path(ds_out,"Final_Ranked_Network")
      dir.create(export_dir,recursive=TRUE,showWarnings=FALSE)
      # Export-only: the selected union is enforced again by export_network().
      score[!method_union_flag(data,cols)] <- 0
      exp <- export_network(ds,strategy,score,data,export_dir,cols=cols)
      external_evidence_validation(ds,score,data,chip,ds_out)
      atomic_csv(exp,file.path(export_dir,"Final_Network_Export_Summary.csv"))
      csv <- file.path(export_dir,exp$File[1L])
      if (!file.exists(csv)) stop("Final high-confidence network export missing")
      atomic_csv(data.table(Dataset=ds,Strategy=strategy,
                            Combination=selected$Combination[1L],Input_Stage=input_stage,
                            Final_Ranking_Percent=network_top_pct,
                            Candidate_Universe="Union_of_Selected_Base_Methods",
                            Selected_Input_Methods=selected$Combination[1L],
                            Final_Edge_File=normalizePath(csv,winslash="/",mustWork=TRUE),
                            Significance_Claim="NO_PVALUE_NO_FDR_These_are_empirical_high_ranked_edges",
                            Trained_Only_On="Development_TF_ChIP_edges",
                            Holdout_Edge_Labels_Used=FALSE),
                 file.path(export_dir,"Final_Network_Methodology.csv"))
      rm(score);gc(verbose=FALSE)
      atomic_new_marker("Top score edges exported; high-confidence NOT p<0.05.",final_done)
    }
  }
  invisible(selected)
}

# --------------------- LOCKED SIX-STRATEGY BENCHMARK ---------------------------
# Executed ONLY AFTER all six development CV winners and the global winner
# have been frozen. This is a single prespecified comparison on held-out TFs;
# the result MUST NOT be used to change any model, hyperparameter or subset.
locked_six_strategy_benchmark <- function(ds,data,plan,chip,ds_seed,ds_out,comb_names) {
  if (!evaluate_locked_test || !run_locked_six_comparison) return(invisible(NULL))
  done <- file.path(ds_out,"LOCKED_SIX_BENCHMARK_DONE.txt")
  pooled_file <- file.path(ds_out,"Locked_SixStrategies_Pooled.csv")
  needed_ci <- c("AUROC_L95","AUROC_U95","AUPR_L95","AUPR_U95")
  if (file.exists(done) && file.exists(pooled_file)) {
    old_result <- fread(pooled_file)
    if (all(needed_ci %in% names(old_result)) && nrow(old_result)==length(strategies) &&
        setequal(old_result$Strategy,strategies) &&
        all(vapply(old_result[,..needed_ci],function(z) all(is.finite(z)),logical(1)))) {
      return(invisible(TRUE))
    }
    message("[CI UPGRADE] " ,ds,": six-strategy scores exist but six sets of CIs are missing.")
    # Refit only the six development-locked models and recompute their heldout CIs.
    # Never redo development CV selection or use test labels to tune models.
  }
  lockfile <- file.path(ds_out,"CV_Selected_Global_Model_LOCKED.csv")
  if (!file.exists(lockfile)) stop("Global CV lock is required before test comparison")
  win <- fread(file.path(ds_out,"CV_Selected_Winners.csv"))
  stopifnot(setequal(win$Strategy,strategies),all(win$Status=="OK"))
  if (length(intersect(plan$outer$train,plan$outer$test))) stop("TF leakage in locked test")
  rows <- list(); per_tfs <- list()
  for (i in seq_along(strategies)) {
    s <- strategies[i]
    w <- win[Strategy==s]
    cols <- match(strsplit(w$Combination[1L],"+",fixed=TRUE)[[1L]],method_names)
    if (anyNA(cols)) stop("Invalid combination for ",s)
    if (s %in% supervised) {
      own_flag <- method_union_flag(data,cols)
      own_train <- get_own_training(data,plan$outer$train,own_flag,ds_seed+500L)
      model <- fit_supervised(s,data,cols,own_train,ds_seed+9000L+i)
      score <- predict_own_on_common(model,data,plan$test,own_flag)
      rm(model)
    } else {
      j <- match(w$Combination[1L],comb_names)
      score <- unsup_scores(s,data,cols,ds_seed+j*37L)[plan$test]
    }
    score <- apply_common_eval_floor(
      score, plan$test, method_union_flag(data,cols))
    m <- metric(data$edges$Label[plan$test],score)
    own_idx <- which(method_union_flag(data,cols)[plan$test])
    mo <- metric_safe(data$edges$Label[plan$test[own_idx]],score[own_idx])
    ci <- cluster_ci(data$edges$TF[plan$test],data$edges$Label[plan$test],
                     score,bootstrap_tf_reps,ds_seed+900000L+i)
    rows[[i]] <- data.table(Dataset=ds,Strategy=s,Combination=w$Combination[1L],
                            Test_AUROC=unname(m["AUROC"]),
                            Test_AUPR=unname(m["AUPR"]),
                            Test_Baseline=unname(m["Baseline"]),
                            OwnUnion_Test_AUROC=unname(mo["AUROC"]),
                            OwnUnion_Test_AUPR=unname(mo["AUPR"]),
                            OwnUnion_Test_Baseline=unname(mo["Baseline"]),
                            OwnUnion_Test_Pos=as.integer(mo["N_Pos"]),
                            OwnUnion_Test_Unverified=as.integer(mo["N_Unverified"]),
                            OwnUnion_Test_ChIP_Coverage=if (sum(data$edges$Label[plan$test])>0)
                              unname(mo["N_Pos"])/sum(data$edges$Label[plan$test]) else NA_real_,
                            N_Test_TFs=length(plan$outer$test),
                            AUROC_L95=unname(ci["AUROC_L95"]),
                            AUROC_U95=unname(ci["AUROC_U95"]),
                            AUPR_L95=unname(ci["AUPR_L95"]),
                            AUPR_U95=unname(ci["AUPR_U95"]),
                            N_Valid_Boot=as.integer(ci["Valid_Boot"]),
                            CI_Method="TF-cluster percentile bootstrap, conditional on fitted model",
                            Post_Lock=TRUE,Used_For_Model_Selection=FALSE)
    tfm <- per_tf_metric(data$edges$TF[plan$test],data$edges$Label[plan$test],score)
    tfm[,`:=`(Dataset=ds,Strategy=s)]
    per_tfs[[i]] <- tfm
  }
  atomic_csv(rbindlist(rows),file.path(ds_out,"Locked_SixStrategies_Pooled.csv"))
  test_tfs <- rbindlist(per_tfs)
  atomic_csv(test_tfs,file.path(ds_out,"Locked_SixStrategies_PerTF.csv"))
  inference <- paired_tf_inference(test_tfs,ds,ds_seed+700000L,
                                   context="Locked_TF_Disjoint_Test_After_Model_Freeze")
  # With ~5 test TFs, Holm-adjusted randomization tests have very limited power.
  atomic_csv(inference,file.path(ds_out,"Locked_XGBoost_vs_Five_Holm.csv"))
  atomic_new_marker("All six development-selected combinations tested after global lock; no retuning.",done)
  invisible(TRUE)
}

# ---------------------------- MAIN ------------------------------------------
run_reviewer456_stratified <- function(active_mode) {
  split_mode <<- active_mode
  root_out <<- file.path(root_base,active_mode)

  dir.create(root_out,recursive=TRUE,showWarnings=FALSE)
  chip <- read_chip(chip_path)
  comb_list <- combinations()
  comb_names <- vapply(comb_list,function(ix) paste(method_names[ix],collapse="+"),character(1))
  cat("\n===== SPLIT MODE: ",split_mode," | ",length(comb_list),
      " subsets per strategy | ",length(strategies)," strategies =====\n",sep="")
  for (ds_idx in seq_along(stress_names)) {
    ds <- stress_names[ds_idx]; ds_seed <- seed+10000L*ds_idx
    ds_out <- file.path(root_out,ds)
    fp <- make_run_fingerprint(ds,active_mode)
    verify_or_create_signature(ds_out,fp)
    all_done <- vapply(strategies,function(s) {
      sdir <- file.path(ds_out,"Strategies",s)
      cv_done <- file.exists(file.path(sdir,"CV_DONE.txt")) &&
        file.exists(file.path(sdir,"Best_Combination_CV.csv")) &&
        file.exists(file.path(sdir,"CV_Subsets.csv"))
      cv_done
    },logical(1))
    if (all(all_done) && file.exists(file.path(ds_out,"BASELINES_DONE.txt")) &&
        (threshold_sweep_development ||
           file.exists(file.path(ds_out,"REVIEWER_CV_DIAGNOSTICS_DONE.txt"))) &&
        (!evaluate_locked_test || !run_locked_six_comparison ||
         (file.exists(file.path(ds_out,"LOCKED_SIX_BENCHMARK_DONE.txt")) &&
          file.exists(file.path(ds_out,"Locked_SixStrategies_Pooled.csv")) &&
          {check_ci <- fread(file.path(ds_out,"Locked_SixStrategies_Pooled.csv"))
          ci_names <- c("AUROC_L95","AUROC_U95","AUPR_L95","AUPR_U95")
          all(ci_names %in% names(check_ci)) &&
            nrow(check_ci)==length(strategies) &&
            all(vapply(check_ci[,..ci_names],function(z) all(is.finite(z)),logical(1))) } )) &&
        (!evaluate_locked_test || file.exists(file.path(ds_out,"FINAL_EVAL_DONE.txt"))) &&
        (threshold_sweep_development || !save_winner_networks ||
          file.exists(file.path(ds_out,"FINAL_NETWORK_DONE.txt")))) {
      cat(" [RESUME] entire dataset already finished: ",ds," / ",split_mode,"\n",sep="")
      rollup_stress(ds_out)
      rollup_global()
      next
    }
    cat("\n===== DATASET: ",ds," / ",split_mode," =====\n",sep="")
    data <- build_data(ds,chip)
    plan <- make_split_plan(data$edges,ds,ds_seed)
    stopifnot(!length(intersect(plan$outer$train,plan$outer$test)),
              setequal(names(plan$fold_map),plan$outer$train))
    if (length(plan$outer$test)<10L) message(ds,
                                             ": only ",length(plan$outer$test)," held-out TFs; TF-bootstrap intervals may be unstable.")
    split_stat <- data$edges[TF %in% c(plan$outer$train,plan$outer$test),
                             .(Pos=sum(Label==1L),Unverified=sum(Label==0L)),by=TF]
    split_stat[,`:=`(Split=ifelse(TF %in% plan$outer$train,"Development","LockedTest"),
                     Split_Mode=split_mode)]
    # Do not expose held-out label counts in a pre-selection audit file.
    # Save them only after CV model lock and held-out evaluation.
    fixed_split <- data.table(TF=c(plan$outer$train,plan$outer$test),
                              Split=c(rep("Development",length(plan$outer$train)),
                                      rep("LockedTest",length(plan$outer$test))),
                              Fold=c(as.integer(plan$fold_map[plan$outer$train]),
                                     rep(NA_integer_,length(plan$outer$test))))
    setorder(fixed_split,TF)
    check_path <- file.path(threshold_study_root,"Fixed_TF_Splits",paste0(ds,"_TF_Splits.csv"))
    dir.create(dirname(check_path),recursive=TRUE,showWarnings=FALSE)
    if (file.exists(check_path)) {
      expected <- fread(check_path)
      if (!isTRUE(all.equal(expected,fixed_split,check.attributes=FALSE)))
        stop("Cross-threshold TF splits/folds differ for: ",ds)
    } else atomic_csv(fixed_split,check_path)
    atomic_csv(fixed_split,file.path(ds_out,"TF_Splits.csv"))
    atomic_csv(data$edges[plan$test,.(TF,Target)],
               file.path(ds_out,"Locked_Test_Candidate_IDs_Unlabeled.csv"))
    fold_balance <- rbindlist(lapply(plan$folds,function(f) {
      y <- data$edges$Label[f$valid]
      data.table(Dataset=ds,Fold=f$fold,N_TF=length(f$valid_tfs),
                 N_Pos=sum(y==1L),N_Unverified=sum(y==0L),Baseline=mean(y))
    }))
    atomic_csv(fold_balance,file.path(ds_out,"Balanced_Fold_Check.csv"))
    # Audit TF and target overlap across the data partitions.
    audit <- rbindlist(lapply(plan$folds,function(f) {
      tr_tf <- unique(data$edges$TF[f$train])
      va_tf <- unique(data$edges$TF[f$valid])
      tr_target <- unique(data$edges$Target[f$train])
      va_target <- unique(data$edges$Target[f$valid])
      data.table(Dataset=ds,Fold=f$fold,
                 N_Train_TF=length(tr_tf),N_Validation_TF=length(va_tf),
                 TF_Overlap=length(intersect(tr_tf,va_tf)),
                 N_Train_Target=length(tr_target),N_Validation_Target=length(va_target),
                 Target_Overlap=length(intersect(tr_target,va_target)),
                 Target_Disjoint=length(intersect(tr_target,va_target))==0L,
                 Protocol="TF-disjoint, NOT target-disjoint")
    }))
    audit[,N_Locked_Test_TF:=length(plan$outer$test)]
    audit[,Locked_vs_Development_TF_Overlap:=length(intersect(plan$outer$train,
                                                              plan$outer$test))]
    stopifnot(all(audit$TF_Overlap==0L),
              all(audit$Locked_vs_Development_TF_Overlap==0L))
    atomic_csv(audit,file.path(ds_out,"Reviewer_TF_Target_Leakage_Audit.csv"))
    data_usage <- data.table(
      Stage=c("Training","Model_selection_validation","Locked_holdout_testing",
              "Optional_external_validation","Final_network_prioritization"),
      Data_TF_Role=c("Development TF excluding each CV fold",
                     "Development TF in held-out fivefold CV fold",
                     "Locked test TF only; zero shared TF",
                     "Independent external edge resource if supplied",
                     "Selected-method union for exported edges; fitted supervised labels from development only"),
      Uses_ChIP_Label_For_Fitting=c(TRUE,FALSE,FALSE,FALSE,FALSE),
      Uses_ChIP_Label_For_Metric=c(FALSE,TRUE,TRUE,FALSE,FALSE),
      Uses_For_Model_Selection=c(FALSE,TRUE,FALSE,FALSE,FALSE),
      Comment=c("Training PU proxy labels; sampling limited to training TF",
                paste0("Development CV selects by ",selection_metric," on ",selection_universe),
                "Held-out TFs are evaluated after model lock.",
                "Optional experimental evidence; ChIP overlap excluded; NOT used for selection",
                "Score-based Top-ranked edges, NOT p-value-significant")
    )
    atomic_csv(data_usage,file.path(ds_out,"Reviewer_Data_Usage_Separation.csv"))
    if (evaluate_locked_test) message(ds," / ",active_mode,
                                      ": TF-disjoint holdout is evaluated only after locking development-selected models.")
    for (s in strategies) {
      sdir <- file.path(ds_out,"Strategies",s)
      cvdone <- file.path(sdir,"CV_DONE.txt")
      bestpath <- file.path(sdir,"Best_Combination_CV.csv")
      if (!file.exists(cvdone) || !file.exists(bestpath)) {
        indexes <- candidate_indices(ds,s,comb_names)
        sc <- load_strategy_checkpoints(ds,s,indexes,comb_list,comb_names,
                                        data,plan,ds_seed,fp,ds_out)
        winner <- write_strategy_cv(ds,s,sc,ds_out,comb_names,ds_seed)
        rm(sc)
      } else {
        cat(" [RESUME] CV strategy completed: ",ds," / ",s,"\n",sep="")
        winner <- fread(bestpath)
      }
      # NO holdout evaluation here: finish all six CV strategies first.
      # Once every strategy is known, freeze one global winner and test it once.
      # Update consolidated files immediately after EACH strategy.
      rollup_stress(ds_out)
      rollup_global()
      gc(verbose=FALSE)
      cat(" [CHECKPOINT] ",ds," / ",s,
          " — CV report, best combination and available test summary saved.\n",sep="")
    }
    # Lock the development-selected models before evaluating the held-out TFs.
    rollup_stress(ds_out)
    # Stress-universal strategy is selected only after all three stresses finish.
    if (!threshold_sweep_development)
      stop("Threshold study must keep threshold_sweep_development=TRUE here")
    # Prespecified single-method references are evaluated without selecting on their test results.
    write_baselines(ds,ds_out,data,plan,chip,ds_seed)
    rollup_stress(ds_out)
    rollup_global()
    if (all(file.exists(file.path(ds_out,"Strategies",strategies,"CV_DONE.txt"))))
      atomic_new_marker("All six strategy CV reports saved.",
                        file.path(ds_out,"DATASET_CV_COMPLETE.txt"))
    rm(data,plan); gc(verbose=FALSE)
  }
  # Stage-level summary from completed saved files, without relying on RAM.
  gather <- function(filename) {
    paths <- file.path(root_out,stress_names,filename)
    paths <- paths[file.exists(paths)]
    if (!length(paths)) return(data.table())
    rbindlist(lapply(paths,fread),use.names=TRUE,fill=TRUE)
  }
  cv <- gather("All_6Strategies_CV_Subsets.csv")
  winners <- gather("CV_Selected_Winners.csv")
  tests <- gather("Locked_Test_SixStrategies.csv")
  basecv <- gather("CV_Fixed_Baselines.csv")
  basetest <- gather("Locked_Test_Fixed_Baselines.csv")
  if (nrow(cv)) atomic_csv(cv,file.path(root_out,"ThreeStress_All_CV_Subsets.csv"))
  if (nrow(winners)) {
    atomic_csv(winners,file.path(root_out,"ThreeStress_CV_Selected_Winners.csv"))
    atomic_csv(winners,file.path(root_out,"ThreeStress_SixStrategies_Best_Combinations_CV.csv"))
    for (s in strategies) atomic_csv(winners[Strategy==s],
                                     file.path(root_out,paste0(s,"_ThreeStress_Best_Combinations_CV.csv")))
  }
  if (nrow(tests)) atomic_csv(tests,
                              file.path(root_out,"ThreeStress_LockedHoldout_SixStrategies.csv"))
  if (nrow(basecv)) atomic_csv(basecv,file.path(root_out,"ThreeStress_CV_Fixed_Baselines.csv"))
  if (nrow(basetest)) atomic_csv(basetest,
                                 file.path(root_out,"ThreeStress_LockedHoldout_Fixed_Baselines.csv"))
  if (all(file.exists(file.path(root_out,stress_names,"DATASET_CV_COMPLETE.txt"))))
    atomic_new_marker("CV completed for all requested stresses.",file.path(root_out,"CV_COMPLETE.txt"))
  rollup_global()
  invisible(TRUE)
}

# ============================================================================
# Nested TF-wise cross-validation on development TFs
#
# Outer: five fixed development-TF folds (the locked test TFs are excluded).
# Inner: four TF-disjoint folds on the remaining development TFs. In EVERY
# outer fold, the whole 57-subset x 6-strategy search is rerun separately
# within each stress. One universal fusion STRATEGY is selected by the MEAN
# inner pooled AUROC of three stresses, then worst-stress pooled AUROC, then
# mean pooled AUPR, then the smallest total subset size. Within a strategy,
# each stress can select a different best base-method subset.
# The selected strategy + stress-specific subset are fitted on that outer
# training TF set; prediction is made ONLY on untouched outer validation TFs.
# Outer-fold scores are concatenated by stress to get unbiased-with-respect-to
# inner MODEL-SELECTION pooled OOF AUROC/AUPR, plus TF-bootstrap intervals.
# Conditional on the Top10 input choice; this does NOT validate the threshold
# selection itself if Top10 was selected using data beyond development TFs.
# Checkpoints: every inner configuration, every outer prediction, and locks.
# ============================================================================

nested_root <- file.path(threshold_study_root,
                         sprintf("Nested_Development_%dx%d",nested_outer_folds,nested_inner_folds))

nested_signature <- function() {
  # Bind nested checkpoints to the same underlying input files and model config
  # as standard development CV, plus this new nested-CV protocol.
  fps <- lapply(seq_along(stress_names),function(i)
    make_run_fingerprint(stress_names[i],"stratified_random")$Fingerprint)
  z <- list(Version="Nested_Universal_Selection_v1", Dataset_Fingerprints=fps,
            Outer_Folds=nested_outer_folds,Inner_Folds=nested_inner_folds,
            Number_Of_Subsets=length(combinations()), Strategies=strategies,
            Selection_Metric=selection_metric, Nested_Bootstrap=nested_bootstrap_reps,
            Seed=seed, EvalUniverse=selection_universe,
            LabelPolicy="Development_TFs_Only")
  tmp <- tempfile(fileext=".rds")
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(z,tmp,version=2)
  unname(tools::md5sum(tmp))
}

nested_init <- function() {
  dir.create(nested_root,recursive=TRUE,showWarnings=FALSE)
  key <- nested_signature()
  path <- file.path(nested_root,"NESTED_FINGERPRINT.txt")
  if (file.exists(path)) {
    if (!identical(readLines(path,warn=FALSE),key))
      stop("Nested-CV run fingerprint mismatch. Keep prior results intact and use a new nested_root.")
  } else {
    old <- list.files(nested_root,all.files=TRUE,no..=TRUE)
    if (length(old)) stop("Nested directory is not empty but has no run fingerprint: ",nested_root)
    atomic_lines(key,path)
  }
  key
}

nested_fold_context <- function(data,plan,fold,ds_seed) {
  # Outer validation TFs cannot appear in any inner fit or inner validation.
  outer_validation_tfs <- names(plan$fold_map)[plan$fold_map==fold]
  inner_pool_tfs <- setdiff(plan$outer$train,outer_validation_tfs)
  if (length(inner_pool_tfs)<nested_inner_folds)
    stop("Too few development TFs for nested inner folds")
  fmap <- make_tf_folds(data$edges,inner_pool_tfs,nested_inner_folds,
                        ds_seed+400000L+fold*1009L)
  inner_folds <- lapply(seq_len(nested_inner_folds), function(j) {
    va_tfs <- names(fmap)[fmap==j]
    tr_tfs <- setdiff(inner_pool_tfs,va_tfs)
    va <- get_indices(data$edges,va_tfs,cv_unlabeled_ratio,
                      ds_seed+410000L+fold*1009L+j)
    if (length(intersect(c(tr_tfs,va_tfs),outer_validation_tfs)) ||
        length(intersect(tr_tfs,va_tfs))) stop("Nested inner/outer TF leakage")
    list(Fold=j, Train_TFs=tr_tfs, Validation_TFs=va_tfs, Validation_Rows=va)
  })
  list(OuterFold=fold, OuterValidationTFs=outer_validation_tfs,
       InnerPoolTFs=inner_pool_tfs, InnerFolds=inner_folds)
}

nested_inner_one <- function(ds,strategy,cols,comb_id,data,context,ds_seed) {
  own_flag <- method_union_flag(data,cols)
  unsup <- if (!strategy %in% supervised)
    unsup_scores(strategy,data,cols,ds_seed+comb_id*37L) else NULL
  labs <- vector("list",nested_inner_folds)
  preds <- vector("list",nested_inner_folds)
  folds <- vector("list",nested_inner_folds)
  for (f in context$InnerFolds) {
    va <- f$Validation_Rows
    if (strategy %in% supervised) {
      tr <- get_own_training(data,f$Train_TFs,own_flag,
                             ds_seed+420000L+context$OuterFold*10000L+f$Fold)
      model <- fit_supervised(strategy,data,cols,tr,
               ds_seed+430000L+context$OuterFold*100000L+comb_id*100L+f$Fold)
      p <- predict_own_on_common(model,data,va,own_flag)
      rm(model)
    } else {
      p <- apply_common_eval_floor(unsup[va],va,own_flag)
    }
    y <- data$edges$Label[va]
    m <- metric_safe(y,p)
    if (!all(is.finite(m[c("AUROC","AUPR")])))
      stop("Missing outcome class in nested inner validation fold")
    labs[[f$Fold]] <- y
    preds[[f$Fold]] <- p
    folds[[f$Fold]] <- m
  }
  m <- metric(unlist(labs,use.names=FALSE),unlist(preds,use.names=FALSE))
  if (!all(is.finite(m[c("AUROC","AUPR","Lift")])))
    stop("Undefined nested inner pooled metric")
  data.table(Dataset=ds,OuterFold=context$OuterFold,Strategy=strategy,
             Combination=paste(method_names[cols],collapse="+"),
             CombinationIndex=comb_id,NumMethods=length(cols),
             Inner_Pooled_AUROC=unname(m["AUROC"]),
             Inner_Pooled_AUPR=unname(m["AUPR"]),
             Inner_Pooled_AUPR_Baseline=unname(m["Baseline"]),
             Inner_Pooled_AUPR_Lift=unname(m["Lift"]),
             N_Inner_Training_TFs=length(context$InnerPoolTFs),Status="OK",Error="")
}

nested_search_stress <- function(ds,fold,ds_index,key) {
  dst <- file.path(nested_root,sprintf("OuterFold_%02d",fold),ds)
  selected <- file.path(dst,"Nested_Inner_SixStrategy_Winners.csv")
  if (file.exists(selected)) {
    prev <- fread(selected)
    if (nrow(prev)!=length(strategies) ||
        !setequal(prev$Strategy,strategies) || any(prev$Status!="OK"))
      stop("Nested inner winner file incomplete: ",selected)
    message(" [NESTED RESUME] ",ds," outer=",fold," inner selection finished")
    return(prev)
  }
  # Load only one stress data set at a time; outer/test labels never participate
  # in inner model fitting or ranking. This avoids holding all three in RAM.
  root_out <<- file.path(root_base,"stratified_random")
  split_mode <<- "stratified_random"
  chip <- read_chip(chip_path)
  dat <- build_data(ds,chip)
  ds_seed <- seed+10000L*ds_index
  plan <- make_split_plan(dat$edges,ds,ds_seed)
  context <- nested_fold_context(dat,plan,fold,ds_seed)
  combo <- combinations()
  cnames <- vapply(combo,function(z) paste(method_names[z],collapse="+"),character(1))
  results <- vector("list",length(strategies)*length(combo))
  cursor <- 0L
  for (strategy in strategies) {
    for (j in seq_along(combo)) {
      cursor <- cursor+1L
      ckpt <- file.path(dst,"Checkpoints",strategy,sprintf("Combo_%03d.rds",j))
      if (file.exists(ckpt)) {
        v <- readRDS(ckpt)
        if (!identical(v$Fingerprint,key) || !identical(v$Dataset,ds) ||
            !identical(v$Strategy,strategy) || !identical(v$Combination,cnames[j]) ||
            !identical(v$OuterFold,fold))
          stop("Nested inner checkpoint metadata mismatch: ",ckpt)
        r <- v$Result
      } else {
        cat(sprintf(" [NESTED INNER] %s / outer %d / %s / combo %d/%d\n",
                    ds,fold,strategy,j,length(combo)))
        r <- tryCatch(nested_inner_one(ds,strategy,combo[[j]],j,dat,context,ds_seed),
                      error=function(e) data.table(Dataset=ds,OuterFold=fold,
                        Strategy=strategy,Combination=cnames[j],CombinationIndex=j,
                        NumMethods=length(combo[[j]]),Inner_Pooled_AUROC=NA_real_,
                        Inner_Pooled_AUPR=NA_real_,Inner_Pooled_AUPR_Baseline=NA_real_,
                        Inner_Pooled_AUPR_Lift=NA_real_,
                        N_Inner_Training_TFs=length(context$InnerPoolTFs),
                        Status="FAILED",Error=conditionMessage(e)))
        atomic_rds_new(list(Fingerprint=key,Dataset=ds,Strategy=strategy,
                            OuterFold=fold,Combination=cnames[j],Result=r),ckpt)
      }
      results[[cursor]] <- r
    }
    # Persist within each completed strategy for progress inspection.
    partial <- rbindlist(results[seq_len(cursor)],fill=TRUE)
    atomic_csv(partial,file.path(dst,"Nested_Inner_Partial.csv"))
  }
  z <- rbindlist(results,fill=TRUE)
  winners <- z[Status=="OK" & is.finite(Inner_Pooled_AUROC) &
               is.finite(Inner_Pooled_AUPR)]
  setorder(winners,Strategy,-Inner_Pooled_AUROC,-Inner_Pooled_AUPR,
           NumMethods,Combination)
  winners <- winners[,head(.SD,1L),by=Strategy]
  atomic_csv(z,file.path(dst,"Nested_Inner_AllConfigs.csv"))
  if (nrow(winners)!=length(strategies))
    stop("No complete nested inner winner for at least one strategy: ",ds,
         ". Inspect Nested_Inner_AllConfigs.csv and its failures.")
  atomic_csv(winners,selected)
  rm(dat,plan,context,chip); gc(verbose=FALSE)
  winners
}

nested_choose_universal <- function(fold,all_stress_winners,key) {
  dest <- file.path(nested_root,sprintf("OuterFold_%02d",fold))
  rows <- rbindlist(all_stress_winners,use.names=TRUE,fill=TRUE)
  rank <- rows[,.(Mean_Inner_AUROC=mean(Inner_Pooled_AUROC),
                  Worst_Inner_AUROC=min(Inner_Pooled_AUROC),
                  Mean_Inner_AUPR=mean(Inner_Pooled_AUPR),
                  Sum_NumMethods=sum(NumMethods),N_Stresses=.N),by=Strategy]
  if (nrow(rank)!=length(strategies) || any(rank$N_Stresses!=length(stress_names)))
    stop("Incomplete three-stress inner rankings: outer fold ",fold)
  setorder(rank,-Mean_Inner_AUROC,-Worst_Inner_AUROC,-Mean_Inner_AUPR,
           Sum_NumMethods,Strategy)
  rank[,Rank:=seq_len(.N)]
  choice <- rows[Strategy==rank$Strategy[1L]]
  setorder(choice,Dataset)
  if (nrow(choice)!=length(stress_names)) stop("Incomplete universal inner selection")
  choice[,`:=`(SelectedUniversalStrategy=rank$Strategy[1L],
                InnerMeanThreeStressAUROC=rank$Mean_Inner_AUROC[1L],
                Selection_Used_Outer_Validation_Labels=FALSE,
                Selection_Used_Locked_Test_Labels=FALSE)]
  freeze <- file.path(dest,"Nested_Inner_Universal_Selection_LOCKED.csv")
  if (file.exists(freeze)) {
    old <- fread(freeze)
    if (!identical(old$Strategy,choice$Strategy) ||
        !identical(old$Combination,choice$Combination) ||
        !identical(old$Dataset,choice$Dataset))
      stop("Nested outer fold selection changed: ",fold)
  } else atomic_csv(choice,freeze)
  atomic_csv(rank,file.path(dest,"Nested_Universal_Inner_Strategy_Rankings.csv"))
  choice
}

nested_outer_predict_stress <- function(ds,ds_index,fold,choice,key) {
  dest <- file.path(nested_root,sprintf("OuterFold_%02d",fold),ds)
  ckpt <- file.path(dest,"Nested_Outer_Universal_Predictions.rds")
  if (file.exists(ckpt)) {
    z <- readRDS(ckpt)
    if (!identical(z$Fingerprint,key) || !identical(z$Dataset,ds) ||
        !identical(z$OuterFold,fold) ||
        !identical(z$Strategy,choice$Strategy[1L]) ||
        !identical(z$Combination,choice$Combination[1L]))
      stop("Nested outer prediction checkpoint mismatch: ",ckpt)
    message(" [NESTED RESUME] ",ds," / outer fold ",fold," predictions")
    return(z$Summary)
  }
  root_out <<- file.path(root_base,"stratified_random")
  split_mode <<- "stratified_random"
  chip <- read_chip(chip_path)
  dat <- build_data(ds,chip)
  ds_seed <- seed+10000L*ds_index
  plan <- make_split_plan(dat$edges,ds,ds_seed)
  context <- nested_fold_context(dat,plan,fold,ds_seed)
  outer_train_tfs <- context$InnerPoolTFs
  outer_valid_tfs <- context$OuterValidationTFs
  va <- get_indices(dat$edges,outer_valid_tfs,cv_unlabeled_ratio,
                    ds_seed+450000L+fold)
  s <- choice$Strategy[1L]
  chosen_names <- strsplit(choice$Combination[1L],"+",fixed=TRUE)[[1L]]
  cols <- match(chosen_names,method_names)
  if (anyNA(cols)) stop("Nested selected subset cannot be parsed")
  own_flag <- method_union_flag(dat,cols)
  if (s %in% supervised) {
    tr <- get_own_training(dat,outer_train_tfs,own_flag,
                           ds_seed+460000L+fold)
    model <- fit_supervised(s,dat,cols,tr,
                            ds_seed+470000L+fold*100L+match(s,strategies))
    p <- predict_own_on_common(model,dat,va,own_flag)
    rm(model)
  } else {
    comb_names <- vapply(combinations(),function(z)
      paste(method_names[z],collapse="+"),character(1))
    j <- match(choice$Combination[1L],comb_names)
    if (is.na(j)) stop("Nested combination index not found")
    unsup <- unsup_scores(s,dat,cols,ds_seed+j*37L)
    p <- apply_common_eval_floor(unsup[va],va,own_flag)
    rm(unsup)
  }
  tf <- dat$edges$TF[va]
  y <- dat$edges$Label[va]
  if (!setequal(unique(tf),outer_valid_tfs) ||
      any(tf %in% plan$outer$test) ||
      length(intersect(outer_train_tfs,tf)))
    stop("Nested outer predictions contain forbidden TFs")
  m <- metric(y,p)
  if (!all(is.finite(m[c("AUROC","AUPR")])))
    stop("Undefined nested outer pooled metric")
  summary <- data.table(Dataset=ds,OuterFold=fold,Strategy=s,
                        Combination=choice$Combination[1L],
                        AUROC=unname(m["AUROC"]),AUPR=unname(m["AUPR"]),
                        AUPR_Baseline=unname(m["Baseline"]),
                        AUPR_Lift=unname(m["Lift"]),
                        N_Pos=as.integer(m["N_Pos"]),
                        N_Unverified=as.integer(m["N_Unverified"]),
                        N_Outer_Train_TFs=length(outer_train_tfs),
                        N_Outer_Validation_TFs=length(outer_valid_tfs),
                        TF_Overlap=length(intersect(outer_train_tfs,outer_valid_tfs)))
  atomic_rds_new(list(Fingerprint=key,Dataset=ds,OuterFold=fold,
                      Strategy=s,Combination=choice$Combination[1L],
                      Summary=summary,TF=as.character(tf),Label=as.integer(y),
                      Score=as.numeric(p)),ckpt)
  atomic_csv(summary,file.path(dest,"Nested_Outer_Universal_Performance.csv"))
  rm(dat,plan,context,chip); gc(verbose=FALSE)
  summary
}

run_nested_cv_three_stresses <- function() {
  if (!nested_cv_enabled) return(invisible(NULL))
  if (nested_outer_folds!=cv_folds)
    stop("Nested outer folds must equal the existing development fold count")
  if (!identical(selection_metric,"Pooled_AUROC") ||
      !identical(selection_universe,"Common_All6"))
    stop("Nested search/outer validation must use the same main selection criterion")
  key <- nested_init()
  atomic_csv(data.table(
    InputThreshold="top10",Outer_Development_TF_Folds=nested_outer_folds,
    Inner_TF_Folds=nested_inner_folds,
    Inner_Subsets=length(combinations()),Inner_Strategies=length(strategies),
    Universal_Selection="Mean inner pooled AUROC across three stresses; then minimum AUROC; then mean AUPR; then simpler subsets",
    Inner_Training="Each subset own union; only inner-train TF ChIP labels",
    Inner_Evaluation="Fixed all-six Top10 union; held-out inner TFs",
    Outer_Evaluation="Fixed all-six Top10 union; held-out outer development TFs",
    Independent_TF_Test_Used="No",
    Confidence_Interval="TF-cluster bootstrap on concatenated outer-fold OOF prediction scores",
    Scope="Selection-aware for strategy and subset, conditional on fixed Top10; not a re-evaluation of threshold choice"),
    file.path(nested_root,"Nested_CV_Protocol.csv"))
  done_path <- file.path(nested_root,"NESTED_CV_COMPLETE.txt")
  summary_path <- file.path(nested_root,"Nested_Pooled_OuterOOF_95CI_PerStress.csv")
  if (file.exists(done_path) && file.exists(summary_path)) {
    old <- fread(summary_path)
    if (nrow(old)==length(stress_names) && setequal(old$Dataset,stress_names) &&
        all(is.finite(old$Nested_Pooled_AUROC)) &&
        all(is.finite(old$Nested_Pooled_AUPR))) {
      message(" [NESTED RESUME] all outer-fold results and pooled CIs completed")
      return(invisible(old))
    }
    stop("Nested completion marker is incompatible with saved summary: ",summary_path)
  }
  inner_selection <- list()
  outer_performance <- list()
  rank_tables <- list()
  for (fold in seq_len(nested_outer_folds)) {
    all_win <- lapply(seq_along(stress_names),function(i)
      nested_search_stress(stress_names[i],fold,i,key))
    chosen <- nested_choose_universal(fold,all_win,key)
    inner_selection[[fold]] <- chosen
    rank_tables[[fold]] <- fread(file.path(nested_root,
        sprintf("OuterFold_%02d",fold),"Nested_Universal_Inner_Strategy_Rankings.csv"))[
        ,OuterFold:=fold]
    for (i in seq_along(stress_names)) {
      ds <- stress_names[i]
      w <- chosen[Dataset==ds]
      if (nrow(w)!=1L) stop("Missing stress in nested universal choice: ",ds)
      outer_performance[[length(outer_performance)+1L]] <-
        nested_outer_predict_stress(ds,i,fold,w,key)
    }
  }
  all_choices <- rbindlist(inner_selection,fill=TRUE)
  all_fold_results <- rbindlist(outer_performance,fill=TRUE)
  atomic_csv(all_choices,file.path(nested_root,"Nested_All_OuterFold_Selected_Models.csv"))
  atomic_csv(rbindlist(rank_tables),
             file.path(nested_root,"Nested_Universal_Inner_Rankings_ByOuterFold.csv"))
  atomic_csv(all_fold_results,file.path(nested_root,"Nested_OuterFold_Performance.csv"))
  freq <- unique(all_choices[,.(OuterFold,Strategy)])[,.(N_Selected_OuterFolds=.N),by=Strategy]
  freq[,Frequency:=N_Selected_OuterFolds/nested_outer_folds]
  atomic_csv(freq,file.path(nested_root,"Nested_Strategy_Selection_Frequency.csv"))
  comb_freq <- all_choices[,.(N_Selected_OuterFolds=.N),by=.(Dataset,Strategy,Combination)]
  atomic_csv(comb_freq,file.path(nested_root,"Nested_Subset_Selection_Frequency.csv"))
  ans <- vector("list",length(stress_names))
  for (i in seq_along(stress_names)) {
    ds <- stress_names[i]
    pred <- lapply(seq_len(nested_outer_folds),function(k) {
      z <- readRDS(file.path(nested_root,sprintf("OuterFold_%02d",k),ds,
                                 "Nested_Outer_Universal_Predictions.rds"))
      data.table(TF=z$TF,Label=z$Label,Score=z$Score,OuterFold=k)
    })
    pooled <- rbindlist(pred)
    if (anyNA(pooled$Score) || any(!is.finite(pooled$Score)))
      stop("Nested pooled OOF scores incomplete for ",ds)
    if (uniqueN(pooled$TF) < nested_outer_folds)
      stop("Nested outer OOF has fewer TFs than outer folds for ",ds)
    cross <- pooled[,.(N_Outer_Folds=uniqueN(OuterFold)),by=TF]
    if (any(cross$N_Outer_Folds!=1L))
      stop("Nested OOF duplicates the same TF in multiple outer folds")
    pm <- metric(pooled$Label,pooled$Score)
    ci <- cluster_ci(pooled$TF,pooled$Label,pooled$Score,
                     nested_bootstrap_reps, seed+610000L+i*101L)
    if (any(!is.finite(ci[c("AUROC_L95","AUROC_U95","AUPR_L95","AUPR_U95")])))
      stop("Nested pooled TF-cluster CI could not be calculated: ",ds)
    ans[[i]] <- data.table(Dataset=ds,
      Method="Fivefold outer TF-wise CV with fourfold inner full strategy/subset search",
      Input_Stage=input_stage,Evaluation_Universe="Top10_AllSix_Union",
      N_Outer_Folds=nested_outer_folds,N_Inner_Folds=nested_inner_folds,
      N_Outer_TFs=uniqueN(pooled$TF),
      N_Pos=as.integer(pm["N_Pos"]),N_Unverified=as.integer(pm["N_Unverified"]),
      Nested_Pooled_AUROC=unname(pm["AUROC"]),
      Nested_AUROC_L95=unname(ci["AUROC_L95"]),
      Nested_AUROC_U95=unname(ci["AUROC_U95"]),
      Nested_Pooled_AUPR=unname(pm["AUPR"]),
      Nested_AUPR_L95=unname(ci["AUPR_L95"]),
      Nested_AUPR_U95=unname(ci["AUPR_U95"]),
      Nested_AUPR_Baseline=unname(pm["Baseline"]),
      Nested_AUPR_Lift=unname(pm["Lift"]),
      Bootstrap_TF_Reps=nested_bootstrap_reps,
      Bootstrap_Valid_Reps=unname(ci["Valid_Boot"]),
      Interpretation="Selection-aware for six strategies x 57 subsets, conditional on preset Top10; development TFs only")
    rm(pooled,pred); gc(verbose=FALSE)
  }
  summary <- rbindlist(ans)
  atomic_csv(summary,file.path(nested_root,"Nested_Pooled_OuterOOF_95CI_PerStress.csv"))
  atomic_csv(summary[,.(Mean_Nested_AUROC=mean(Nested_Pooled_AUROC),
                        Worst_Nested_AUROC=min(Nested_Pooled_AUROC),
                        Mean_Nested_AUPR=mean(Nested_Pooled_AUPR),
                        Mean_Nested_AUPR_Lift=mean(Nested_AUPR_Lift),
                        N_Stresses=.N)],
             file.path(nested_root,"Nested_ThreeStress_Overall_Summary.csv"))
  atomic_new_marker("Nested development TF-wise outer/inner evaluation finished; no independent test read.",
                    file.path(nested_root,"NESTED_CV_COMPLETE.txt"))
  message(" [NESTED COMPLETE] ",nested_root)
  invisible(summary)
}

# ======================== FIXED TOP 10% DEVELOPMENT ========================
# Use development CV to select one shared fusion strategy and a base-algorithm
# subset per stress at the fixed Top 10% candidate threshold.
atomic_csv(data.table(Stage=threshold_stages,Reference=reference_stage,
                      Evaluation_Universe="Top10_AllSix_Union",
                      Method="Per-combination own-union fitting, common-universe evaluation",
                      MissingScore=outside_union_eval_score),
           file.path(threshold_study_root,"Threshold_Protocol.csv"))

for (stage in threshold_stages) {
  input_stage <- stage
  root_base <- file.path(threshold_study_root,paste0("Stage_",stage))
  threshold_sweep_development <- TRUE
  evaluate_locked_test <- FALSE
  cat("\n======== DEVELOPMENT THRESHOLD: ",stage," ========\n",sep="")
  for (mode in split_modes) run_reviewer456_stratified(mode)
  rollup_global()
  gc(verbose=FALSE)
}

# ---------------------- DEVELOPMENT-ONLY SELECTION ----------------------------
# Read the CV results once all three stresses have completed the search.
threshold_candidates <- list()
for (stage in threshold_stages) for (ds in stress_names) {
  f <- file.path(threshold_study_root,paste0("Stage_",stage),
                 "stratified_random",ds,"CV_Selected_Winners.csv")
  if (!file.exists(f)) stop("Incomplete threshold development CV: ",f)
  z <- fread(f)
  if (nrow(z)!=length(strategies) || !setequal(z$Strategy,strategies) ||
      any(z$Status!="OK") || any(z$Valid_Folds!=cv_folds) ||
      any(!is.finite(z$CV_Pooled_AUROC)) || any(!is.finite(z$CV_Pooled_AUPR)))
    stop("Invalid/incomplete six-strategy CV winners: ",f)
  z[, `:=`(Input_Stage=stage, Reference_Stage=reference_stage,
           Eval_Reference="Top10_AllSix_Union")]
  threshold_candidates[[length(threshold_candidates)+1L]] <- z
}
per_stress <- rbindlist(threshold_candidates,fill=TRUE)
atomic_csv(per_stress,file.path(threshold_study_root,
                                "Threshold_AllStrategies_PerStress_CV.csv"))

# Strict predeclared universal criterion: mean pooled AUROC; tie break by
# worst-stress pooled AUROC, then mean pooled AUPR, then fewer input methods.
# Ties use the stated secondary sort order.
universal <- per_stress[,.(Mean_CV_AUROC=mean(CV_Pooled_AUROC),
                          Worst_CV_AUROC=min(CV_Pooled_AUROC),
                          Mean_CV_AUPR=mean(CV_Pooled_AUPR),
                          Mean_CV_AUPR_Lift=mean(CV_Pooled_AUPR_Lift),
                          Sum_NumMethods=sum(NumMethods),
                          N_Stresses=.N),by=.(Input_Stage,Strategy)]
if (any(universal$N_Stresses != length(stress_names)) ||
    nrow(universal) != length(threshold_stages)*length(strategies))
  stop("Universal comparison is incomplete")
setorder(universal,-Mean_CV_AUROC,-Worst_CV_AUROC,-Mean_CV_AUPR,
         Sum_NumMethods,Input_Stage,Strategy)
universal[,Rank:=seq_len(.N)]
atomic_csv(universal,file.path(threshold_study_root,
                               "Threshold_Strategy_ThreeStress_Comparison.csv"))
selected_group <- universal[1L]
selected_per_stress <- per_stress[Input_Stage==selected_group$Input_Stage &
                                 Strategy==selected_group$Strategy]
setorder(selected_per_stress,Dataset)
if (nrow(selected_per_stress)!=length(stress_names)) stop("Missing selected stress")
selected_per_stress[,`:=`(Selection_Level="Universal_Threshold_And_Strategy",
                         Selection_Basis="ThreeStress_Mean_Development_Pooled_AUROC",
                         Holdout_Used_For_Selection=FALSE,
                         Model_Locked_Before_Holdout=TRUE)]

# Nested CV is an independent DEVELOPMENT-ONLY evaluation of the full search
# pipeline. Its outer validation metrics are NEVER used to select the final
# universal strategy; the latter is selected on all development TFs below.
# This preserves separation between model search and its unbiased assessment.
if (nested_cv_enabled) run_nested_cv_three_stresses()

# Immutability protection: never overwrite a different frozen selection with
# a later CV run, especially after someone might have seen the test outcome.
lock_path <- file.path(threshold_study_root,"FINAL_Universal_Selection_LOCKED.csv")
stress_lock_path <- file.path(threshold_study_root,"FINAL_PerStress_Models_LOCKED.csv")
if (file.exists(lock_path) || file.exists(stress_lock_path)) {
  if (!file.exists(lock_path) || !file.exists(stress_lock_path))
    stop("Incomplete frozen selection: both lock files must exist")
  old_group <- fread(lock_path)
  old_stress <- fread(stress_lock_path)
  if (!identical(old_group$Input_Stage,selected_group$Input_Stage) ||
      !identical(old_group$Strategy,selected_group$Strategy) ||
      !isTRUE(all.equal(old_group$Mean_CV_AUROC,selected_group$Mean_CV_AUROC,
                        tolerance=1e-9,check.attributes=FALSE)) ||
      !identical(old_stress$Combination,selected_per_stress$Combination) ||
      !identical(old_stress$Dataset,selected_per_stress$Dataset))
    stop("Frozen threshold/model selection changed; refusing to evaluate holdout")
} else {
  atomic_csv(selected_group,lock_path)
  atomic_csv(selected_per_stress,stress_lock_path)
}

# Report each stress's top strategy separately as a descriptive result.
per_stress_best <- per_stress[order(-CV_Pooled_AUROC,-CV_Pooled_AUPR,NumMethods),.SD[1],by=Dataset]
atomic_csv(per_stress_best,file.path(threshold_study_root,
                                     "Threshold_PerStress_Exploratory_Best.csv"))

# ----------------- TEST ONCE AFTER UNIVERSAL FREEZE ---------------------------
# Lock files above exist before ANY test labels are evaluated or reported.
input_stage <- selected_group$Input_Stage[1L]
root_base <- file.path(threshold_study_root,paste0("Stage_",input_stage))
# build_data() saves diagnostic files using root_out; make sure it never writes
# selected-stage diagnostics into the previously processed stage's directory.
root_out <- file.path(root_base,"stratified_random")
threshold_sweep_development <- FALSE
for (ds_idx in seq_along(stress_names)) {
  ds <- stress_names[ds_idx]
  w <- selected_per_stress[Dataset==ds]
  s <- w$Strategy[1L]
  ds_seed <- seed+10000L*ds_idx
  ds_dev <- file.path(root_base,"stratified_random",ds)
  ds_final <- file.path(threshold_study_root,"FINAL_Locked_Test_And_Network",ds)
  dir.create(ds_final,recursive=TRUE,showWarnings=FALSE)
  cat("\n======== LOCKED FINAL: ",ds," | ",input_stage," | ",s," ========\n",sep="")
  split_mode <- "stratified_random"
  chip <- read_chip(chip_path)
  data <- build_data(ds,chip)
  plan <- make_split_plan(data$edges,ds,ds_seed)
  fixed <- fread(file.path(threshold_study_root,"Fixed_TF_Splits",
                           paste0(ds,"_TF_Splits.csv")))
  current <- data.table(TF=c(plan$outer$train,plan$outer$test),
                        Split=c(rep("Development",length(plan$outer$train)),
                                rep("LockedTest",length(plan$outer$test))),
                        Fold=c(as.integer(plan$fold_map[plan$outer$train]),
                               rep(NA_integer_,length(plan$outer$test))))
  setorder(current,TF)
  if (!isTRUE(all.equal(fixed,current,check.attributes=FALSE)))
    stop("Frozen test/CV TF partition differs across phases: ",ds)
  comb_list <- combinations()
  comb_names <- vapply(comb_list,function(ix) paste(method_names[ix],collapse="+"),character(1))
  if (!w$Combination[1L] %in% comb_names) stop("Invalid selected combination: ",ds)
  # The selected threshold + strategy + subset are immutable for this stress.
  frozen_file <- file.path(ds_final,"LOCKED_Chosen_Model.csv")
  if (file.exists(frozen_file)) {
    old <- fread(frozen_file)
    if (!identical(old$Input_Stage,w$Input_Stage) ||
        !identical(old$Strategy,w$Strategy) ||
        !identical(old$Combination,w$Combination))
      stop("Per-stress frozen selection changed: ",ds)
  } else atomic_csv(w,frozen_file)
  # Compute repeated-CV diagnostics at the selected threshold.
  reviewer_development_stats(ds,data,plan,ds_seed,ds_dev,comb_names)
  # Final heldout: exactly the one universally selected fusion per stress.
  evaluate_locked_test <- TRUE
  if (!file.exists(file.path(ds_final,"Strategies",s,"TEST_DONE.txt"))) {
    write_strategy_test(ds,s,w,data,plan,chip,ds_seed,ds_final,comb_names)
  }
  test_file <- file.path(ds_final,"Strategies",s,"Locked_Test.csv")
  if (!file.exists(test_file)) stop("Final locked-test result missing: ",ds)
  test_result <- fread(test_file)
  test_result[,`:=`(Input_Stage=input_stage,Reference_Stage=reference_stage,
                   CrossStage_Selection_Locked=TRUE)]
  atomic_csv(test_result,file.path(ds_final,"FINAL_Selected_Model_TF_Disjoint_Test.csv"))
  # Export scores for the SELECTED combination union only, not all-six union.
  if (save_winner_networks && !file.exists(file.path(ds_final,"FINAL_NETWORK_DONE.txt"))) {
    cols <- match(strsplit(w$Combination[1L],"+",fixed=TRUE)[[1L]],method_names)
    own_flag <- method_union_flag(data,cols)
    score <- numeric(nrow(data$edges))
    if (s %in% supervised) {
      own_train <- get_own_training(data,plan$outer$train,own_flag,ds_seed+500L)
      mdl <- fit_supervised(s,data,cols,own_train,
                            ds_seed+9000L+match(s,strategies))
      own_rows <- which(own_flag)
      score[own_rows] <- predict_in_batches(mdl,data,rows=own_rows)
      rm(mdl)
    } else {
      j <- match(w$Combination[1L],comb_names)
      score <- unsup_scores(s,data,cols,ds_seed+j*37L)
      score[!own_flag] <- 0
    }
    export_dir <- file.path(ds_final,"Final_Ranked_Network")
    dir.create(export_dir,recursive=TRUE,showWarnings=FALSE)
    exp <- export_network(ds,s,score,data,export_dir,cols=cols)
    atomic_csv(exp,file.path(export_dir,"Final_Network_Export_Summary.csv"))
    atomic_csv(data.table(Dataset=ds,Input_Stage=input_stage,Strategy=s,
                          Combination=w$Combination[1L],
                          Evaluation_Background="Top10_AllSix_Union",
                          Network_Candidate_Background="Selected_Combination_Union",
                          Network_Top_Fraction=network_top_pct,
                          Selection_Used_Holdout_Labels="Not used by this run; historical use not verified",
                          Network_Statistical_FDR_Control=FALSE),
               file.path(export_dir,"Final_Network_Methodology.csv"))
    external_evidence_validation(ds,score,data,chip,ds_final)
    atomic_new_marker("Final network selected after universal CV lock.",
                      file.path(ds_final,"FINAL_NETWORK_DONE.txt"))
  }
  evaluate_locked_test <- FALSE
  rm(data,plan,chip);gc(verbose=FALSE)
}

heldout_files <- vapply(stress_names,function(ds)
  file.path(threshold_study_root,"FINAL_Locked_Test_And_Network",ds,
            "FINAL_Selected_Model_TF_Disjoint_Test.csv"),character(1))
if (any(!file.exists(heldout_files))) stop("Missing final test results")
final_test <- rbindlist(lapply(heldout_files,fread),fill=TRUE)
atomic_csv(final_test,file.path(threshold_study_root,
                                "FINAL_ThreeStress_Locked_TF_Disjoint_Tests.csv"))
atomic_new_marker("Fixed Top10: development selection and exploratory heldout evaluation complete.",
                  file.path(threshold_study_root,"STUDY_COMPLETE.txt"))
cat("\nCOMPLETE: fixed Top10 model selection, exploratory test and ranked export: ",
    threshold_study_root,"\n",sep="")
