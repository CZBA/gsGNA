# =============================================================================
# gsGNA | Biotic-stress application and evaluation
#
# Scope:
#   Apply the gsGNA ensemble workflow to three additional biotic-stress datasets:
#   Magnaporthe oryzae (MO), Nephotettix cincticeps feeding (NC), and rice stripe
#   virus (RSV).
#
# Workflow:
#   - Compare all 57 combinations of two to six base GRN inference algorithms.
#   - Select each stress-specific combination using mean repeated-CV pooled AUROC.
#   - Use mean pooled AUPR and then fewer algorithms as tie-breakers.
#   - Compare the selected XGBoost ensemble with six individual algorithms on
#     the same All6 common candidate-edge universe.
#   - Fit and score each ensemble configuration only on its own base-method
#     candidate union. Edges outside that union receive a score of zero during
#     common-universe evaluation and are not used for model training.
#
# Interpretation:
#   - No independent TF-disjoint holdout is performed in this analysis.
#   - Post-selection performance on the same data is descriptive and should not
#     be interpreted as independent or unbiased confirmation of model optimality.
#   - ChIP-supported positives are contrasted with unverified edges. Unverified
#     edges are not confirmed negatives; this is a positive-unlabeled (PU) proxy.
# =============================================================================

required <- c("data.table", "xgboost", "pROC", "PRROC", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Please install the required packages with: install.packages(c(", paste(sprintf("'%s'", missing), collapse=","), "))")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
data.table::setDTthreads(4L)

# --------------------------------- CONFIG -----------------------------------
input_root <- "."
output_root <- file.path(input_root, "gsGNA_Biotic_Application_Evaluation")
stresses <- c("MO", "NC", "RSV")
methods <- c("GENIE3", "Kboost", "GRNBoost2", "3DCEMA", "DeepRIG", "IGEGRNs")
stage <- "top10"                 # Input filenames use the pattern <stress>_<stage>_edges.csv.
chip_file <- file.path(input_root, "chip_seq.csv")
chip_tf_col <- 1L                # Column index for TF identifiers in chip_seq.csv.
chip_target_col <- 2L            # Column index for target identifiers in chip_seq.csv.
seed <- 20261105L
cv_folds <- 5L
cv_repeats <- 10L
train_unverified_ratio <- 1.0    # Training-only sampling ratio; all unverified candidates are retained for evaluation.
xgb_rounds <- 150L
xgb_depth <- 3L
xgb_eta <- 0.10
xgb_threads <- 4L
bootstrap_B <- 2000L             # Number of TF-cluster bootstrap replicates.
prediction_batch <- 200000L
network_top_fraction <- 0.1      # Fraction of the selected-combination union to export; 0.1 exports the top 10%.
export_final_network <- TRUE
resume_subsets <- TRUE           # Resume completed subset runs when inputs and settings are unchanged.
reuse_finished_bootstrap <- TRUE # Reuse compatible completed bootstrap outputs; set FALSE after changing bootstrap settings.
# ----------------------------------------------------------------------------
stopifnot(length(methods) == 6L, cv_folds >= 3L, cv_repeats >= 2L,
          train_unverified_ratio > 0, bootstrap_B >= 100L,
          network_top_fraction > 0, network_top_fraction <= 1)
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)

get_path <- function(stress, method) {
  folders <- if (method == "GRNBoost2") c("GRNBoost2", "GRNBoost") else method
  possible <- file.path(input_root, folders, paste0(stress, "_", stage, "_edges.csv"))
  hits <- possible[file.exists(possible)]
  if (!length(hits)) stop("Missing input file(s): ", paste(possible, collapse=" / "))
  hits[1L]
}

read_chip <- function() {
  if (!file.exists(chip_file)) stop("Required file not found: ", chip_file)
  z <- fread(chip_file, showProgress=FALSE)
  if (max(chip_tf_col, chip_target_col) > ncol(z) || chip_tf_col == chip_target_col)
    stop("Invalid chip_tf_col or chip_target_col setting.")
  out <- data.table(TF=toupper(trimws(as.character(z[[chip_tf_col]]))),
                    Target=toupper(trimws(as.character(z[[chip_target_col]]))))
  out <- unique(out[!is.na(TF) & !is.na(Target) & nzchar(TF) &
                    nzchar(Target) & TF != Target])
  if (!nrow(out)) stop("No valid ChIP interaction labels remain after filtering.")
  setkey(out, TF, Target)
  out
}

read_one_grn <- function(stress, method) {
  path <- get_path(stress, method)
  z <- fread(path, select=1:3, showProgress=FALSE, check.names=FALSE)
  if (ncol(z) != 3L) stop("Expected the first three columns to be TF, Target and Weight: ", path)
  out <- data.table(TF=toupper(trimws(as.character(z[[1L]]))),
                    Target=toupper(trimws(as.character(z[[2L]]))),
                    Weight=suppressWarnings(as.numeric(z[[3L]])))
  out <- out[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) &
               TF != Target & is.finite(Weight)]
  if (!nrow(out)) stop("No valid edges remain in: ", path)
  if (method == "3DCEMA") out[, Weight := abs(Weight)]
  if (any(out$Weight < 0))
    stop(method, " contains negative edge weights. Check the score direction and preprocessing: ", path)
  out <- out[, .(Weight=max(Weight)), by=.(TF, Target)]
  out[, Method := method]
  cat(sprintf("  %-10s %s edges | %s\n", method,
              format(nrow(out), big.mark=","), path))
  out
}

# Build matrices once. Min-max normalization is label-free, within present edges;
# missing-method edges are score 0 and presence indicator 0.
load_stress <- function(stress, chip) {
  long <- rbindlist(lapply(methods, function(m) read_one_grn(stress, m)))
  wide <- dcast(long, TF + Target ~ Method, value.var="Weight", fill=NA_real_)
  rm(long)
  if (!all(methods %in% names(wide))) stop("Inputs are incomplete for one or more of the six base algorithms: ", stress)
  setkey(wide, TF, Target)
  wide[, Label := 0L]
  wide[chip, on=.(TF,Target), Label := 1L]
  N <- nrow(wide)
  X <- matrix(0, nrow=N, ncol=12L,
              dimnames=list(NULL, c(paste0("W_", methods), paste0("P_", methods))))
  mask <- integer(N)
  for (j in seq_along(methods)) {
    m <- methods[j]
    v <- wide[[m]]
    present <- !is.na(v)
    X[, j+6L] <- as.numeric(present)
    mask <- mask + as.integer(present) * as.integer(2^(j-1L))
    if (any(present)) {
      rng <- range(v[present])
      X[present,j] <- if (diff(rng) > 0) (v[present]-rng[1L])/diff(rng) else 1
    }
    wide[, (m) := NULL]
  }
  wide[, Mask := mask]
  eligible_stats <- wide[, .(N_Pos=sum(Label == 1L),
                              N_Unverified=sum(Label == 0L)), by=TF]
  eligible_stats <- eligible_stats[N_Pos > 0L & N_Unverified > 0L]
  if (nrow(eligible_stats) < cv_folds)
    stop(stress, ": fewer than five TFs have both ChIP-positive and unverified candidate edges.")
  eligible_idx <- which(wide$TF %in% eligible_stats$TF)
  eval_data <- wide[eligible_idx]
  eval_X <- X[eligible_idx, , drop=FALSE]
  coverage <- data.table(Stress=stress, Candidate_N=nrow(wide),
                         Candidate_TFs=uniqueN(wide$TF),
                         Eligible_TFs=nrow(eligible_stats),
                         Eligible_N=nrow(eval_data),
                         ChIP_All_N=nrow(chip),
                         ChIP_Covered_N=sum(wide$Label == 1L),
                         ChIP_Eligible_N=sum(eval_data$Label == 1L),
                         Positive_Fraction=mean(eval_data$Label == 1L))
  cat(sprintf("  Eligible TFs: %d | common candidate edges: %s | eligible edges: %s\n",
              nrow(eligible_stats), format(nrow(wide), big.mark=","),
              format(nrow(eval_data), big.mark=",")))
  list(stress=stress, all=wide, Xall=X, eval=eval_data, Xeval=eval_X,
       stats=eligible_stats, coverage=coverage)
}

# Balance TF counts first, then approximate per-fold positive-label load.
make_tf_assignment <- function(stats, k, local_seed) {
  if (nrow(stats) < k) stop("The number of eligible TFs is insufficient for ", k, "-fold cross-validation.")
  z <- copy(stats)
  set.seed(local_seed)
  z <- z[sample.int(nrow(z))]
  setorder(z, -N_Pos)
  quota <- rep(nrow(z) %/% k, k)
  if (nrow(z) %% k) quota[seq_len(nrow(z) %% k)] <- quota[seq_len(nrow(z) %% k)] + 1L
  counts <- integer(k); loads <- numeric(k); fold <- integer(nrow(z))
  for (i in seq_len(nrow(z))) {
    available <- which(counts < quota)
    chosen <- available[which.min(loads[available])]
    fold[i] <- chosen
    counts[chosen] <- counts[chosen] + 1L
    loads[chosen] <- loads[chosen] + z$N_Pos[i]
  }
  setNames(fold, z$TF)
}

make_plans <- function(d) {
  lapply(seq_len(cv_repeats), function(r) {
    tf_to_fold <- make_tf_assignment(d$stats, cv_folds, seed + 10000L*r)
    row_fold <- unname(tf_to_fold[d$eval$TF])
    if (anyNA(row_fold)) stop("One or more TFs could not be assigned to a cross-validation fold.")
    row_fold
  })
}

all_combos <- unlist(lapply(2:6, function(k)
  combn(seq_along(methods), k, simplify=FALSE)), recursive=FALSE)
stopifnot(length(all_combos) == 57L)
comb_name <- function(cols) paste(methods[cols], collapse="+")
comb_bit <- function(cols) as.integer(sum(2^(cols-1L)))

metric <- function(y, score) {
  # Do not silently discard nonfinite predictions or mismatched rows.
  if (length(y) != length(score) || !length(y) ||
      anyNA(y) || any(!y %in% c(0L, 1L)) || any(!is.finite(score))) {
    stop(sprintf(paste0(
      "Invalid evaluation data: number of labels=%s, number of predictions=%s, ",
      "missing labels=%s, non-binary labels=%s, non-finite predictions=%s"),
      length(y), length(score), sum(is.na(y)),
      sum(!is.na(y) & !y %in% c(0L, 1L)), sum(!is.finite(score))))
  }
  npos <- sum(y == 1L); nunv <- sum(y == 0L)
  baseline <- npos / length(y)
  if (!npos || !nunv)
    return(c(AUROC=NA_real_, AUPR=NA_real_, Baseline=baseline,
             Lift=NA_real_, N_Pos=npos, N_Unverified=nunv))
  roc <- as.numeric(pROC::auc(pROC::roc(y, score, levels=c(0,1),
                                          direction="<", quiet=TRUE)))
  pr <- PRROC::pr.curve(scores.class0=score[y == 1L],
                        scores.class1=score[y == 0L], curve=FALSE)$auc.integral
  c(AUROC=roc, AUPR=pr, Baseline=baseline, Lift=pr/baseline,
    N_Pos=npos, N_Unverified=nunv)
}

fit_xgb <- function(X, y, local_seed) {
  if (!all(c(0L,1L) %in% y)) stop("The training data must contain both positive and unverified edges.")
  set.seed(local_seed)
  xgboost::xgb.train(
    params=list(objective="binary:logistic", eval_metric="auc",
                max_depth=xgb_depth, eta=xgb_eta,
                subsample=0.8, colsample_bytree=0.8,
                nthread=xgb_threads, seed=local_seed),
    data=xgboost::xgb.DMatrix(X, label=y),
    nrounds=xgb_rounds, verbose=0)
}

predict_batches <- function(fit, X, ids) {
  if (!length(ids)) return(numeric())
  result <- numeric(length(ids))
  starts <- seq.int(1L, length(ids), by=prediction_batch)
  for (b in starts) {
    end <- min(length(ids), b + prediction_batch - 1L)
    ix <- seq.int(b, end)
    result[ix] <- as.numeric(predict(fit, xgboost::xgb.DMatrix(X[ids[ix], , drop=FALSE])))
  }
  result
}

sample_train <- function(rows, labels, local_seed) {
  pos <- rows[labels[rows] == 1L]
  unv <- rows[labels[rows] == 0L]
  if (!length(pos) || !length(unv)) stop("The training candidate union is missing one of the two label classes for this TF split.")
  amount <- min(length(unv), ceiling(length(pos)*train_unverified_ratio))
  set.seed(local_seed)
  c(pos, if (amount == length(unv)) unv else sample(unv, amount))
}

# Every validation edge has a TF that is completely absent from training.
# Predict only own-union candidates; everything else is 0 for common evaluation.
oof_one_repeat <- function(d, cols, foldmap, local_seed) {
  flag <- bitwAnd(d$eval$Mask, comb_bit(cols)) != 0L
  feature_cols <- c(cols, cols+6L)
  scores <- numeric(nrow(d$eval))
  for (fold in seq_len(cv_folds)) {
    train_pool <- which(flag & foldmap != fold)
    test_rows <- which(flag & foldmap == fold)
    if (!length(test_rows)) next
    tr <- sample_train(train_pool, d$eval$Label,
                       local_seed + fold*101L)
    fit <- fit_xgb(d$Xeval[tr, feature_cols, drop=FALSE],
                   d$eval$Label[tr], local_seed + fold*101L + 7L)
    # rows of the reduced feature matrix index the original test_rows position
    scores[test_rows] <- predict_batches(fit,
      d$Xeval[, feature_cols, drop=FALSE], test_rows)
    rm(fit)
  }
  scores
}

# Durable per-combination CSVs allow an interrupted 57-combination sweep to resume.
run_subset_selection <- function(d, plans, stress_dir, manifest) {
  cache_dir <- file.path(stress_dir, "Subset_CV_Cache")
  dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
  manifest_path <- file.path(cache_dir, "Input_And_Settings_Manifest.rds")
  if (file.exists(manifest_path) && resume_subsets) {
    old <- readRDS(manifest_path)
    if (!identical(old, manifest)) stop("The cached inputs or settings have changed. Use a new output_root or remove the old cache.")
  } else {
    if (length(list.files(cache_dir, pattern="^Combo_[0-9]+\\.csv$")))
      stop("An existing cache was found without a matching manifest. Remove Subset_CV_Cache or use a new output_root before rerunning.")
    saveRDS(manifest, manifest_path)
  }
  all_rows <- vector("list", length(all_combos))
  for (i in seq_along(all_combos)) {
    cols <- all_combos[[i]]
    path <- file.path(cache_dir, sprintf("Combo_%02d.csv", i))
    if (file.exists(path) && resume_subsets) {
      tab <- fread(path)
      if (nrow(tab) != cv_repeats || !identical(sort(tab$Repeat), seq_len(cv_repeats)))
        stop("The number of repeats in the cached result does not match the current settings: ", path)
      all_rows[[i]] <- tab
      next
    }
    cat(sprintf("[%s] XGBoost subset %02d/57: %s\n",
                d$stress, i, comb_name(cols)))
    repetitions <- vector("list", cv_repeats)
    for (r in seq_len(cv_repeats)) {
      score <- oof_one_repeat(d, cols, plans[[r]],
                              seed + i*100000L + r*1000L)
      met <- metric(d$eval$Label, score)
      repetitions[[r]] <- data.table(Stress=d$stress, Combination=comb_name(cols),
        N_Methods=length(cols), Combo_ID=i, Repeat=r,
        AUROC=unname(met["AUROC"]), AUPR=unname(met["AUPR"]),
        AUPR_Baseline=unname(met["Baseline"]), AUPR_Lift=unname(met["Lift"]),
        N_Pos=as.integer(met["N_Pos"]),
        N_Unverified=as.integer(met["N_Unverified"]))
    }
    tab <- rbindlist(repetitions)
    fwrite(tab, path)
    all_rows[[i]] <- tab
  }
  repeat_dt <- rbindlist(all_rows)
  summary_dt <- repeat_dt[, .(
     Mean_AUROC=mean(AUROC), SD_AUROC=sd(AUROC),
     Mean_AUPR=mean(AUPR), SD_AUPR=sd(AUPR),
     Mean_AUPR_Lift=mean(AUPR_Lift), AUPR_Baseline=mean(AUPR_Baseline)),
     by=.(Stress, Combo_ID, Combination, N_Methods)]
  setorder(summary_dt, -Mean_AUROC, -Mean_AUPR, N_Methods, Combo_ID)
  summary_dt[, Rank := seq_len(.N)]
  fwrite(repeat_dt, file.path(stress_dir, "XGBoost_57Subsets_10x5_Repeats.csv"))
  fwrite(summary_dt, file.path(stress_dir, "XGBoost_57Subsets_Ranking.csv"))
  list(winner=summary_dt[1L], summary=summary_dt, repeats=repeat_dt)
}

selected_repeated_oof <- function(d, cols, plans, combo_id) {
  avg <- numeric(nrow(d$eval))
  for (r in seq_len(cv_repeats)) {
    avg <- avg + oof_one_repeat(d, cols, plans[[r]],
                                seed + combo_id*100000L + r*1000L)/cv_repeats
  }
  avg
}

# TF-cluster bootstrap, same resampled TF clusters for every model and
# every paired difference. CI is conditional on the already-selected model.
bootstrap_comparison <- function(d, score_list, stress_dir) {
  labels <- d$eval$Label
  names_models <- names(score_list)
  groups <- split(seq_len(nrow(d$eval)), d$eval$TF)
  G <- length(groups)
  if (G < 3L) stop("TF-cluster bootstrap requires at least three TFs.")

  summary_path <- file.path(stress_dir, "CommonUniverse_Pooled_Performance_TFbootstrap.csv")
  paired_path <- file.path(stress_dir, "CommonUniverse_Paired_Deltas_TFbootstrap.csv")
  if (reuse_finished_bootstrap && file.exists(summary_path) && file.exists(paired_path)) {
    saved_summary <- fread(summary_path)
    saved_paired <- fread(paired_path)
    compatible <- all(c("Stress", "Model", "N_Pos", "N_Unverified", "N_TFs") %in% names(saved_summary)) &&
      all(c("Stress", "Comparator") %in% names(saved_paired)) &&
      nrow(saved_summary) == length(names_models) &&
      setequal(saved_summary$Model, names_models) &&
      nrow(saved_paired) == length(names_models) - 1L &&
      setequal(saved_paired$Comparator, setdiff(names_models, "XGBoost ensemble")) &&
      all(saved_summary$Stress == d$stress) && all(saved_paired$Stress == d$stress) &&
      all(saved_summary$N_Pos == sum(labels == 1L)) &&
      all(saved_summary$N_Unverified == sum(labels == 0L)) &&
      all(saved_summary$N_TFs == G)
    if (isTRUE(compatible)) {
      cat("  Reused completed TF-bootstrap outputs for ", d$stress,
          " (set reuse_finished_bootstrap=FALSE after changing bootstrap settings).\n", sep="")
      return(list(summary=saved_summary, paired=saved_paired))
    }
    cat("  Existing bootstrap CSVs failed compatibility checks; recomputing.\n")
  }

  full <- lapply(score_list, function(sc) metric(labels, sc))
  A <- matrix(NA_real_, nrow=bootstrap_B, ncol=length(score_list),
              dimnames=list(NULL, names_models))
  P <- A
  set.seed(seed + 89011L + match(d$stress, stresses)*1000L)
  for (b in seq_len(bootstrap_B)) {
    cluster_ids <- sample.int(G, G, replace=TRUE)
    sampled_rows <- unlist(groups[cluster_ids], use.names=FALSE)
    y <- labels[sampled_rows]
    if (!any(y==1L) || !any(y==0L)) next
    for (j in seq_along(score_list)) {
      v <- metric(y, score_list[[j]][sampled_rows])
      A[b,j] <- unname(v["AUROC"])
      P[b,j] <- unname(v["AUPR"])
    }
    if (b %% 100L == 0L) cat("  TF bootstrap: ", b, "/", bootstrap_B, "\n", sep="")
  }
  ci <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < bootstrap_B*0.8) return(c(NA_real_, NA_real_))
    as.numeric(quantile(x, c(0.025, 0.975), names=FALSE))
  }
  summary <- rbindlist(lapply(seq_along(score_list), function(j) {
    cm <- full[[j]]
    roc_ci <- ci(A[,j]); pr_ci <- ci(P[,j])
    data.table(Stress=d$stress, Model=names_models[j],
      AUROC=unname(cm["AUROC"]), AUROC_L95=roc_ci[1L], AUROC_U95=roc_ci[2L],
      AUPR=unname(cm["AUPR"]), AUPR_L95=pr_ci[1L], AUPR_U95=pr_ci[2L],
      AUPR_Baseline=unname(cm["Baseline"]), AUPR_Lift=unname(cm["Lift"]),
      N_Pos=as.integer(cm["N_Pos"]), N_Unverified=as.integer(cm["N_Unverified"]),
      N_TFs=G)
  }))
  idx <- match("XGBoost ensemble", names_models)
  paired <- rbindlist(lapply(which(names_models != "XGBoost ensemble"), function(j) {
    ar <- A[,idx] - A[,j]
    ap <- P[,idx] - P[,j]
    dar <- ci(ar); dap <- ci(ap)
    data.table(Stress=d$stress, Comparator=names_models[j],
      Delta_AUROC=unname(full[[idx]]["AUROC"] - full[[j]]["AUROC"]),
      Delta_AUROC_L95=dar[1L], Delta_AUROC_U95=dar[2L],
      Delta_AUPR=unname(full[[idx]]["AUPR"] - full[[j]]["AUPR"]),
      Delta_AUPR_L95=dap[1L], Delta_AUPR_U95=dap[2L],
      Method="Paired_TF_Cluster_Bootstrap_Descriptive_PostSelection")
  }))
  fwrite(summary, file.path(stress_dir, "CommonUniverse_Pooled_Performance_TFbootstrap.csv"))
  fwrite(paired, file.path(stress_dir, "CommonUniverse_Paired_Deltas_TFbootstrap.csv"))
  list(summary=summary, paired=paired)
}

# Supplementary: OWN-UNION scores are not directly compared across models.
# Some single methods may have zero edges among TFs eligible on the All6 union.
# Zero candidates, or candidates containing only one label class, are not
# estimable; report NA with an explicit status rather than crashing or silently
# changing the evaluation population.
write_own_universe <- function(d, scores, win_cols, stress_dir) {
  n <- nrow(d$eval)
  if (nrow(d$Xeval) != n || length(d$eval$Mask) != n ||
      length(d$eval$Label) != n)
    stop("OWN-UNION: Xeval, Mask and Label have inconsistent row counts.")

  rows <- rbindlist(lapply(names(scores), function(name) {
    sc_all <- scores[[name]]
    if (!is.numeric(sc_all) || length(sc_all) != n)
      stop(sprintf("OWN-UNION [%s/%s]: score length=%d, evaluation rows=%d",
                   d$stress, name, length(sc_all), n))

    if (identical(name, "XGBoost ensemble")) {
      flag <- bitwAnd(d$eval$Mask, comb_bit(win_cols)) != 0L
    } else {
      j <- match(name, methods)
      if (is.na(j)) stop("OWN-UNION: unknown model name: ", name)
      flag <- d$Xeval[, j + 6L] == 1
    }
    if (length(flag) != n || anyNA(flag))
      stop(sprintf("OWN-UNION [%s/%s]: flag length=%d, missing flags=%d, evaluation rows=%d",
                   d$stress, name, length(flag), sum(is.na(flag)), n))

    idx <- which(flag)  # Integer row indices avoid NA logical subsetting.
    y <- d$eval$Label[idx]
    sc <- sc_all[idx]
    if (anyNA(y) || any(!is.finite(sc)))
      stop(sprintf("OWN-UNION [%s/%s]: missing labels after filtering=%d, non-finite predictions=%d",
                   d$stress, name, sum(is.na(y)), sum(!is.finite(sc))))

    npos <- sum(y == 1L)
    nunv <- sum(y == 0L)
    status <- if (length(idx) == 0L) "No_candidate_edges" else
              if (npos == 0L) "No_ChIP_positive" else
              if (nunv == 0L) "No_unverified_edges" else "OK"

    if (length(idx) == 0L) {
      # No scored candidates: AUROC, AUPR and baseline are undefined.
      v <- c(AUROC=NA_real_, AUPR=NA_real_, Baseline=NA_real_,
             Lift=NA_real_, N_Pos=0, N_Unverified=0)
    } else {
      v <- metric(y, sc)
    }
    cat(sprintf("  OWN-UNION [%s] %-17s candidates=%s, pos=%s, unverified=%s | %s\n",
                d$stress, name, format(length(idx), big.mark=","),
                format(npos, big.mark=","), format(nunv, big.mark=","), status))

    data.table(Stress=d$stress, Model=name,
      Candidate_N=length(idx), Pos_N=npos, Unverified_N=nunv,
      Evaluation_Status=status,
      Positive_Coverage_Of_Common=npos / sum(d$eval$Label == 1L),
      OwnUnion_AUROC=unname(v["AUROC"]), OwnUnion_AUPR=unname(v["AUPR"]),
      OwnUnion_AUPR_Baseline=unname(v["Baseline"]),
      OwnUnion_AUPR_Lift=unname(v["Lift"]))
  }))
  fwrite(rows, file.path(stress_dir, "OwnUnion_Descriptive_Performance.csv"))
  invisible(rows)
}

export_network <- function(d, cols, stress_dir) {
  if (!export_final_network) return(invisible(NULL))
  flag <- bitwAnd(d$all$Mask, comb_bit(cols)) != 0L
  rows <- which(flag)
  if (!length(rows)) stop("The selected algorithm combination has no candidate edges.")
  train_pool <- which(bitwAnd(d$eval$Mask, comb_bit(cols)) != 0L)
  train_rows <- sample_train(train_pool, d$eval$Label, seed + 75001L)
  feats <- c(cols, cols+6L)
  fit <- fit_xgb(d$Xeval[train_rows, feats, drop=FALSE],
                 d$eval$Label[train_rows], seed + 75002L)
  # Predict ALL TFs in the selected subset's own union; never use All6 union as output.
  score <- predict_batches(fit, d$Xall[, feats, drop=FALSE], rows)
  output <- data.table(TF=d$all$TF[rows], Target=d$all$Target[rows],
                       XGB_Score=score)
  if (network_top_fraction < 1.0) {
    setorder(output, -XGB_Score, TF, Target)
    amount <- max(1L, ceiling(nrow(output)*network_top_fraction))
    output <- output[seq_len(amount)]
  }
  path <- file.path(stress_dir, if (network_top_fraction < 1.0)
    "XGBoost_SelectedUnion_TopRank_Edges.csv" else
    "XGBoost_SelectedUnion_ALL_Edges.csv")
  fwrite(output, path)
  cat("  Exported network: ", path, " (", nrow(output), " edges)\n", sep="")
  invisible(path)
}

make_figure <- function(all_stats) {
  plot_dt <- rbindlist(list(
    all_stats[, .(Stress, Model, Metric="AUROC", Score=AUROC,
                  L95=AUROC_L95, U95=AUROC_U95)],
    all_stats[, .(Stress, Model, Metric="AUPR", Score=AUPR,
                  L95=AUPR_L95, U95=AUPR_U95)]))
  plot_dt[, Model := factor(Model, levels=rev(c("XGBoost ensemble", methods)))]
  plot_dt[, Stress := factor(Stress, levels=stresses)]
  plot_dt[, Type := ifelse(Model == "XGBoost ensemble", "XGBoost ensemble", "Single GRN algorithm")]
  base_lines <- unique(all_stats[, .(Stress, AUPR_Baseline)])
  base_lines[, `:=`(Metric="AUPR", Stress=factor(Stress, levels=stresses))]
  p <- ggplot(plot_dt, aes(x=Score, y=Model, color=Type)) +
    geom_vline(data=base_lines, aes(xintercept=AUPR_Baseline),
               inherit.aes=FALSE, linetype="dashed", linewidth=0.4,
               color="grey50") +
    geom_segment(aes(x=L95, xend=U95, yend=Model),
                 linewidth=0.7, na.rm=TRUE) +
    geom_point(size=2.4) +
    facet_grid(Stress ~ Metric, scales="free_x") +
    scale_color_manual(values=c("XGBoost ensemble"="#C0392B",
                                "Single GRN algorithm"="#286B91")) +
    labs(x="Pooled OOF performance", y=NULL, color=NULL,
         title="gsGNA application in additional biotic stress scenarios",
         subtitle="TF-wise 5-fold CV × 10 repeats; 95% TF-cluster bootstrap CI") +
    theme_bw(base_size=10) +
    theme(text=element_text(color="black"),
          axis.text=element_text(color="black"),
          strip.text=element_text(color="black", face="bold"),
          plot.title=element_text(face="bold", color="black", size=12),
          plot.subtitle=element_text(color="black", size=9),
          legend.position="bottom", panel.grid.minor=element_blank())
  ggsave(file.path(output_root, "MO_NC_RSV_AUROC_AUPR_A4.pdf"),
         p, width=11.69, height=8.27, units="in", device="pdf")
  ggsave(file.path(output_root, "MO_NC_RSV_AUROC_AUPR_A4.png"),
         p, width=11.69, height=8.27, units="in", dpi=300)
}

# ---------------------------------- RUN -------------------------------------
chip <- read_chip()
cat("ChIP positive labels: ", format(nrow(chip), big.mark=","), "\n", sep="")
all_results <- list(); all_pairs <- list(); all_winners <- list(); all_coverage <- list()
for (s in stresses) {
  cat("\n================ STRESS: ", s, " ================\n", sep="")
  stress_dir <- file.path(output_root, s)
  dir.create(stress_dir, recursive=TRUE, showWarnings=FALSE)
  # A light manifest guards against accidental reuse after changing input files.
  files <- c(chip_file, vapply(methods, function(m) get_path(s,m), character(1)))
  info <- file.info(files)
  manifest <- list(files=normalizePath(files, winslash="/"),
                   file_sizes=info$size, file_mtimes=as.numeric(info$mtime),
                   methods=methods, stage=stage, seed=seed,
                   cv_folds=cv_folds, cv_repeats=cv_repeats,
                   train_unverified_ratio=train_unverified_ratio,
                   xgb_rounds=xgb_rounds, xgb_depth=xgb_depth,
                   xgb_eta=xgb_eta, xgb_threads=xgb_threads)
  d <- load_stress(s, chip)
  fwrite(d$coverage, file.path(stress_dir, "Candidate_And_ChIP_Coverage.csv"))
  plans <- make_plans(d)
  selected <- run_subset_selection(d, plans, stress_dir, manifest)
  best <- selected$winner
  cols <- all_combos[[best$Combo_ID]]
  cat("  Winner: ", best$Combination, " | mean AUROC=",
      round(best$Mean_AUROC,4), " | mean AUPR=", round(best$Mean_AUPR,4), "\n", sep="")
  score_xgb <- selected_repeated_oof(d, cols, plans, best$Combo_ID)
  model_scores <- setNames(lapply(seq_along(methods), function(j) d$Xeval[,j]), methods)
  model_scores <- c(list("XGBoost ensemble"=score_xgb), model_scores)
  # Mean of 10 out-of-fold predictions for each edge: no TF's labels were used
  # to train any of the component predictions for that same TF.
  # Validate each model's own-union coverage before the expensive bootstrap.
  write_own_universe(d, model_scores, cols, stress_dir)
  bench <- bootstrap_comparison(d, model_scores, stress_dir)
  fwrite(data.table(Stress=s, Strategy="XGBoost", Winning_Combination=best$Combination,
                    N_Methods=best$N_Methods, Selection_Criterion="Mean_10_Repeated_Pooled_AUROC",
                    Selection_Mean_AUROC=best$Mean_AUROC,
                    Selection_Mean_AUPR=best$Mean_AUPR,
                    Assessment="Same_Data_PostSelection_Repeated_TF_CV"),
         file.path(stress_dir, "Selected_XGBoost_Combination.csv"))
  export_network(d, cols, stress_dir)
  all_results[[s]] <- bench$summary
  all_pairs[[s]] <- bench$paired
  all_winners[[s]] <- data.table(Stress=s, Combo_ID=best$Combo_ID,
                                  Winning_Combination=best$Combination,
                                  Selection_Mean_AUROC=best$Mean_AUROC,
                                  Selection_Mean_AUPR=best$Mean_AUPR)
  all_coverage[[s]] <- d$coverage
  rm(d, plans, selected, score_xgb, model_scores, bench)
  invisible(gc())
}
combined <- rbindlist(all_results)
fwrite(combined, file.path(output_root, "MO_NC_RSV_Main_Performance.csv"))
fwrite(rbindlist(all_pairs), file.path(output_root, "MO_NC_RSV_Paired_Differences.csv"))
fwrite(rbindlist(all_winners), file.path(output_root, "MO_NC_RSV_Best_Subsets.csv"))
fwrite(rbindlist(all_coverage), file.path(output_root, "MO_NC_RSV_Coverage.csv"))
make_figure(combined)
cat("\nDONE. Results: ", normalizePath(output_root, winslash="/"), "\n", sep="")
