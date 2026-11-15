# =============================================================================
# gsGNA | TF perturbation DEG validation
# Changes:
#   1. Match GRN and DEG by EXACT TF + Target pairs (not DEG targets pooled across TFs).
#   2. Each DEG file corresponds to exactly ONE TF; validate this and perform ONE test per stress.
#   3. Fix overlap == 0 handling; Expected is still n*K/N.
#   4. Intersect BOTH sets with the same evaluable target background BEFORE testing.
#   5. Report numerator/denominators, two overlap fractions, fold enrichment, OR,
#      Fisher exact 95% CI for OR, and exploratory gene-bootstrap 95% CI for FE.
#   6. Export tables, Venn diagrams, and effect-size overview plots.
#
# INPUTS:
#   - Final XGBoost selected-union GRNs from gsGNA_Ensemble_Comparison.R.
#   - DroughtDEG.xlsx, AlkalineDEG.xlsx, and ColdDEG.xlsx in input_root.
#   - Optional independent eligible-target backgrounds (recommended).
#
# The GRN files contain all edges in each selected-method union; this script
# does not apply an additional Top 10% edge filter.
#
# R dependencies:
#   install.packages(c("data.table", "readxl", "ggplot2", "VennDiagram"))
#   BiocManager::install(c("AnnotationDbi", "org.Osativa.eg.db"))
#
# Run from the project root, or set GSGNA_INPUT_ROOT to the project input root:
#   source("code/Network validation/DEG_validation.R")
#
# IMPORTANT: TF perturbation DEGs suggest functional associations, not necessarily
# direct physical TF binding. A missing DEG is NOT a confirmed negative edge.
# =============================================================================

required <- c('data.table', 'readxl', 'ggplot2', 'VennDiagram')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop('Please install the required R packages with install.packages(c(',
       paste(shQuote(missing), collapse = ', '), '))')
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# ============================= SETTINGS ======================================
input_root <- Sys.getenv("GSGNA_INPUT_ROOT", unset = ".")
# Final XGBoost networks exported by gsGNA_Ensemble_Comparison.R.
# These files contain ALL edges in each stress-specific selected-method union;
network_base <- file.path(
  "gsGNA_FixedTop10",
  "FINAL_Locked_Test_And_Network"
)
network_files <- c(
  Drought = file.path(network_base, "Drought", "Final_Ranked_Network",
                      "XGBoost_Ranked_SelectedUnion_AllEdges.csv"),
  Alkaline = file.path(network_base, "Alkaline", "Final_Ranked_Network",
                       "XGBoost_Ranked_SelectedUnion_AllEdges.csv"),
  Cold = file.path(network_base, "Cold", "Final_Ranked_Network",
                   "XGBoost_Ranked_SelectedUnion_AllEdges.csv")
)
deg_files <- c(
  Drought = "DroughtDEG.xlsx",
  Alkaline = "AlkalineDEG.xlsx",
  Cold = "ColdDEG.xlsx"
)
# Display names for the three independent, single-TF DEG experiments.
# TF matching uses the actual TF IDs in XLSX/GRN, NOT these symbols.
tf_display_names <- c(Drought = "OsbZIP62", Alkaline = "OsWRKY1", Cold = "OsNAC5")

# =========================== COLOR SETTINGS =================================
stress_pal <- c(
  Drought = "#4DBBD5",
  Alkaline = "#E64B35",
  Cold = "#3C5488"
)
# Match the user-supplied VennDiagram figure style.
venn_grn_pal <- c(
  Drought = "#457B9D",
  Alkaline = "#E63946",
  Cold = "#1D3557"
)
venn_deg_fill <- "#F1FAEE"
# =============================================================================

# Optional: a file for each stress listing ALL genes that were eligible both
# for GRN prediction and DEG testing (including genes that were NOT DEGs).
# Accepted formats: .csv, .tsv, .txt or .xlsx; first column may be Gene / Target.
# Use NULL per stress to fall back to the UNION of Target genes in all 3 GRNs.
background_files <- list(Drought = NULL, Alkaline = NULL, Cold = NULL)
# e.g. background_files$Cold <- 'Cold_Eligible_Tested_Genes.txt'

# ID rule: ONLY Cold DEG Target is mapped RAP -> GID. All ID comparisons use
# uppercase (toupper) consistently, without deleting transcript suffixes or
# internal spaces. Cold RAP mapping happens BEFORE uppercase output normalization.
# Keep leading/trailing whitespace trim and the original DEG LOC-target filter.
bootstrap_reps <- 2000L  # exploratory TARGET-gene bootstrap for FE / overlap
random_seed <- 123L
plot_venns <- TRUE
output_dir <- file.path(input_root, 'gsGNA_DEG_Overlap_SingleTF')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# =============================================================================

# Case-normalized identifiers for equality matching across the GRNs and DEGs.
normalize_id <- function(x) {
  toupper(trimws(as.character(x)))
}

# Normalize column names for robust header recognition (e.g. TF ID -> tfid).
name_key <- function(x) {
  tolower(gsub('[^a-z0-9]', '', as.character(x)))
}

convert_cold_deg_target <- function(targets) {
  if (!requireNamespace('AnnotationDbi', quietly = TRUE) ||
      !requireNamespace('org.Osativa.eg.db', quietly = TRUE)) {
    stop('Cold DEG Target RAP-to-GID conversion requires Bioconductor packages. Run: ',
         'BiocManager::install(c("AnnotationDbi", "org.Osativa.eg.db"))')
  }
  raw <- trimws(as.character(targets))
  already_gid <- !is.na(raw) & grepl('^LOC', raw, ignore.case = TRUE)
  need_map <- !is.na(raw) & nzchar(raw) & !already_gid
  rap_keys <- unique(raw[need_map])
  mapped <- setNames(rep(NA_character_, length(rap_keys)), rap_keys)
  if (length(rap_keys)) {
    db <- org.Osativa.eg.db::org.Osativa.eg.db
    valid_keys <- AnnotationDbi::keys(db, keytype = 'RAP')
    # First obtain actual case-sensitive database keys without changing IDs.
    canonical <- valid_keys[match(toupper(rap_keys), toupper(valid_keys))]
    found <- unique(canonical[!is.na(canonical)])
    if (length(found)) {
      lookup <- AnnotationDbi::mapIds(
        db, keys = found, keytype = 'RAP', column = 'GID',
        multiVals = 'first')
      valid_idx <- which(!is.na(canonical))
      mapped[valid_idx] <- unname(lookup[canonical[valid_idx]])
    }
  }
  converted_raw <- raw
  mapped_for_rows <- unname(mapped[match(raw, names(mapped))])
  was_mapped <- need_map & !is.na(mapped_for_rows) & nzchar(mapped_for_rows)
  converted_raw[was_mapped] <- mapped_for_rows[was_mapped]
  result <- normalize_id(converted_raw)
  status <- rep('Unmapped_kept_original', length(raw))
  status[already_gid] <- 'Already_GID'
  status[was_mapped] <- 'Mapped_RAP_to_GID'
  status[is.na(raw) | !nzchar(raw)] <- 'Missing'
  audit <- unique(data.table(
    Original_Target = raw, Converted_Target = result, Status = status),
    by = c('Original_Target', 'Converted_Target', 'Status'))
  fwrite(audit, file.path(output_dir, 'Cold_DEG_Target_RAP_to_GID_Audit.csv'), na = 'NA')
  message(sprintf('[Cold DEG RAP->GID] unique raw IDs=%d; mapped=%d; already LOC=%d; unmapped=%d',
                  length(unique(raw[!is.na(raw)])),
                  sum(audit$Status == 'Mapped_RAP_to_GID'),
                  sum(audit$Status == 'Already_GID'),
                  sum(audit$Status == 'Unmapped_kept_original')))
  if (any(audit$Status == 'Unmapped_kept_original')) {
    warning('Cold stress has ', sum(audit$Status == 'Unmapped_kept_original'),
            ' Target IDs that could not be mapped. Original IDs were retained after uppercasing; the LOC-prefix filter is applied next. ',
            'Unmapped non-LOC IDs will be excluded; inspect the audit file.')
  }
  if (length(rap_keys) && !any(audit$Status == 'Mapped_RAP_to_GID') &&
      !any(audit$Status == 'Already_GID')) {
    stop('No Cold DEG Target IDs were successfully mapped from RAP or recognized as existing LOC IDs. ',
         'Check the RAP ID format and the RAP/GID mappings in org.Osativa.eg.db.')
  }
  result
}

# Distinct TF-target pairs; self loops are excluded consistently.
clean_pairs <- function(tf, target, deg_loc_filter = FALSE) {
  z <- data.table(TF = normalize_id(tf), Target = normalize_id(target))
  # All comparisons are case-normalized; original workflow dropped NA.
  z <- z[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target)]
  if (deg_loc_filter) z <- z[grepl('^LOC', Target)]
  unique(z, by = c('TF', 'Target'))
}

# Network files: first 3 columns = TF, Target, Weight. The code handles both
# genuine headerless CSV and common TF/Target/Weight header rows.
read_grn <- function(path) {
  if (!file.exists(path)) stop('GRN file does not exist: ', path)
  raw <- fread(path, header = FALSE, colClasses = 'character',
               fill = TRUE, showProgress = FALSE, check.names = FALSE)
  if (ncol(raw) < 2L || nrow(raw) == 0L) stop('GRN file is empty or contains fewer than two columns: ', path)
  is_header <- name_key(raw[[1L]][1L]) %in% c('tf','tfid','regulator','regulatorygene') &&
    name_key(raw[[2L]][1L]) %in% c('target','targetgene','targetid','regulatedgene')
  if (is_header) raw <- raw[-1L]
  pairs <- clean_pairs(raw[[1L]], raw[[2L]])
  if (!nrow(pairs)) stop('No valid TF-target pairs remain after cleaning the GRN file: ', path)
  message(sprintf('[GRN] %s: %s distinct pairs; header=%s',
                  basename(path), format(nrow(pairs), big.mark = ','), is_header))
  pairs
}

# DEG files: select TF and Target columns explicitly by recognizable names.
# If missing, use first 2 columns with a warning, as in the original workflow.
read_deg <- function(path, convert_cold_target = FALSE) {
  if (!file.exists(path)) stop('DEG file does not exist: ', path)
  raw <- as.data.table(readxl::read_excel(path, sheet = 1, .name_repair = 'unique'))
  if (ncol(raw) < 2L) stop('DEG file must contain at least two columns: ', path)
  nm <- name_key(names(raw))
  tf_candidate <- which(nm %in% c('tf','tfid','tfgene','regulator','transcriptionfactor'))
  tg_candidate <- which(nm %in% c('target','targetgene','targetid','targetgeneid',
                                  'gene','geneid','gid','rap','locus','locusid'))
  if (length(tf_candidate) == 1L && length(tg_candidate) == 1L &&
      tf_candidate != tg_candidate) {
    i <- tf_candidate; j <- tg_candidate
  } else {
    warning('Could not uniquely identify TF/Target columns in the DEG file; using the first two columns: ', path,
            '; verify that these columns represent TF and Target.')
    i <- 1L; j <- 2L
  }
  targets <- if (convert_cold_target) convert_cold_deg_target(raw[[j]]) else raw[[j]]
  pairs <- clean_pairs(raw[[i]], targets, deg_loc_filter = TRUE)
  if (!nrow(pairs)) stop('No valid TF-target pairs remain after cleaning the DEG file: ', path)
  message(sprintf('[DEG] %s: %s distinct pairs; %d TF(s)',
                  basename(path), format(nrow(pairs), big.mark = ','), uniqueN(pairs$TF)))
  pairs
}

read_background <- function(path) {
  if (!file.exists(path)) stop('Background file does not exist: ', path)
  if (grepl('\\.xlsx?$', path, ignore.case = TRUE)) {
    raw <- as.data.table(readxl::read_excel(path, sheet = 1))
  } else {
    raw <- fread(path, header = 'auto', colClasses = 'character',
                 showProgress = FALSE)
  }
  if (ncol(raw) < 1L) stop('Background file does not contain a gene-ID column: ', path)
  candidates <- which(name_key(names(raw)) %in%
                        c('gene','geneid','target','targetgene','targetid','gid',
                          'rap','locus','locusid'))
  i <- if (length(candidates) == 1L) candidates else 1L
  ids <- unique(normalize_id(raw[[i]]))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) stop('No valid gene IDs were found in the background file: ', path)
  ids
}

# Wilson intervals describe uncertainty under a binomial sampling assumption.
# Biological dependencies can make them too narrow, so we use exploratory gene
# bootstrap CIs for overlap ratios instead of implying TF-level independence.
quantile_ci <- function(x, n_valid_min = 50L) {
  x <- x[is.finite(x)]
  if (length(x) < n_valid_min) return(c(NA_real_, NA_real_))
  as.numeric(quantile(x, probs = c(.025, .975), names = FALSE, na.rm = TRUE))
}

# Gene resampling via multinomial draws of the 2x2 gene counts.
# This is the ordinary nonparametric bootstrap of gene-level binary statuses;
# it does NOT account for gene-gene dependence or selection of the universe.
boot_gene <- function(a, b, c, d, reps, seed) {
  empty <- list(Overlap_GRN_L95 = NA_real_, Overlap_GRN_U95 = NA_real_,
                Overlap_DEG_L95 = NA_real_, Overlap_DEG_U95 = NA_real_,
                FE_L95 = NA_real_, FE_U95 = NA_real_,
                N_Valid_FEBoot = 0L)
  N <- a + b + c + d
  if (N <= 1L || reps < 50L) return(empty)
  set.seed(seed)
  samp <- rmultinom(reps, size = N, prob = c(a,b,c,d) / N)
  aa <- samp[1L, ]; bb <- samp[2L, ]; cc <- samp[3L, ]
  pred_n <- aa + bb
  ref_n <- aa + cc
  support <- ifelse(pred_n > 0, aa / pred_n, NA_real_)
  recovery <- ifelse(ref_n > 0, aa / ref_n, NA_real_)
  expected <- pred_n * ref_n / N
  fe <- ifelse(expected > 0, aa / expected, NA_real_)
  support_ci <- quantile_ci(support)
  recovery_ci <- quantile_ci(recovery)
  fe_ci <- quantile_ci(fe)
  list(Overlap_GRN_L95 = support_ci[1L], Overlap_GRN_U95 = support_ci[2L],
       Overlap_DEG_L95 = recovery_ci[1L], Overlap_DEG_U95 = recovery_ci[2L],
       FE_L95 = fe_ci[1L], FE_U95 = fe_ci[2L],
       N_Valid_FEBoot = sum(is.finite(fe)))
}

# All counts use the SAME target universe for a particular TF and stress.
analyze_tf <- function(stress, tf, grn, deg, universe, bg_mode, index) {
  grn_raw <- unique(grn[TF == tf, Target])
  deg_raw <- unique(deg[TF == tf, Target])
  pred <- intersect(grn_raw, universe)
  ref <- intersect(deg_raw, universe)
  N <- length(universe)
  n <- length(pred)
  K <- length(ref)
  a <- length(intersect(pred, ref))
  b <- n - a
  c <- K - a
  d <- N - a - b - c
  if (any(c(a,b,c,d) < 0)) stop('Invalid 2x2 contingency table; check the background and TF: ', stress, '/', tf)
  
  expected <- if (N > 0) n * K / N else NA_real_
  FE <- if (is.finite(expected) && expected > 0) a / expected else NA_real_
  support <- if (n > 0) a / n else NA_real_
  recovery <- if (K > 0) a / K else NA_real_
  # OR = a*d/(b*c); apply Haldane-Anscombe correction ONLY if a zero cell.
  cells <- c(a,b,c,d)
  corrected <- any(cells == 0)
  orcells <- if (corrected) cells + .5 else cells
  OR <- (orcells[1L] * orcells[4L]) / (orcells[2L] * orcells[3L])
  
  # Hypergeometric one-sided test for OVER-representation, not independent
  # validation of direct regulation. Meaningful mainly with defensible universe.
  p_hyper <- if (N > 0 && n > 0 && K > 0 && a > 0)
    phyper(a - 1L, K, N - K, n, lower.tail = FALSE) else 1.0
  
  # Exact conditional OR CI (two-sided); standard independence assumptions
  # across genes apply. Fisher estimate may differ slightly from sample OR.
  fisher_ci <- tryCatch({
    tab <- matrix(c(a,c,b,d), nrow = 2L, byrow = FALSE)
    stats::fisher.test(tab, alternative = 'two.sided', conf.int = TRUE)$conf.int
  }, error = function(e) {
    warning('Fisher OR confidence-interval calculation failed for ', stress, '/', tf, ': ', conditionMessage(e))
    c(NA_real_, NA_real_)
  })
  gene_ci <- boot_gene(a,b,c,d, bootstrap_reps, random_seed + index * 37L)
  
  counts <- data.table(Stress = stress, TF = tf,
                       N_Universe = N,
                       N_GRN_Targets_Original = length(grn_raw),
                       N_DEG_Targets_Original = length(deg_raw),
                       N_GRN_Targets_Excluded = length(setdiff(grn_raw, universe)),
                       N_DEG_Targets_Excluded = length(setdiff(deg_raw, universe)),
                       N_GRN_Targets = n, N_DEG_Targets = K,
                       N_Overlap = a, N_GRN_Only = b,
                       N_DEG_Only = c, N_Neither = d,
                       GRN_Overlap_Fraction = support,
                       DEG_Recovery_Fraction = recovery,
                       Expected_Overlap = expected,
                       Fold_Enrichment = FE, Odds_Ratio = as.numeric(OR),
                       OR_Fisher_L95 = fisher_ci[1L],
                       OR_Fisher_U95 = fisher_ci[2L],
                       Used_Haldane_Correction = corrected,
                       P_Hyper_OneSided = p_hyper,
                       Background_Mode = bg_mode,
                       N_Gene_Bootstrap = bootstrap_reps)
  for (nm in names(gene_ci)) counts[, (nm) := gene_ci[[nm]]]
  list(summary = counts, grn_targets = pred, deg_targets = ref)
}

# ============================== MAIN =========================================
if (!identical(sort(names(network_files)), sort(names(deg_files))) ||
    !identical(sort(names(network_files)), sort(names(background_files))) ||
    !identical(sort(names(network_files)), sort(names(tf_display_names)))) {
  stop('network_files, deg_files, background_files, and tf_display_names must use identical stress names.')
}

stress_names <- names(network_files)

network_paths <- setNames(file.path(input_root, network_files), names(network_files))
deg_paths <- setNames(file.path(input_root, deg_files), names(deg_files))
if (any(!file.exists(network_paths))) {
  stop("Missing final XGBoost selected-union GRN files:\n",
       paste(network_paths[!file.exists(network_paths)], collapse = "\n"))
}
if (any(!file.exists(deg_paths))) {
  stop("Missing TF-perturbation DEG workbooks:\n",
       paste(deg_paths[!file.exists(deg_paths)], collapse = "\n"))
}
message("Using final XGBoost selected-union GRNs:")
for (s in stress_names) message("  ", s, ": ", network_paths[[s]])
message("Using single-TF DEG workbooks:")
for (s in stress_names) message("  ", s, ": ", deg_paths[[s]])
fwrite(
  data.table(
    Stress = stress_names,
    GRN_File = unname(network_paths[stress_names]),
    DEG_File = unname(deg_paths[stress_names])
  ),
  file.path(output_dir, "DEG_Validation_Input_Manifest.csv")
)

grns <- lapply(stress_names, function(s) read_grn(network_paths[[s]]))
names(grns) <- stress_names
degs <- lapply(stress_names, function(s) read_deg(
  deg_paths[[s]], convert_cold_target = identical(s, "Cold")))
names(degs) <- stress_names

fallback_universe <- sort(unique(unlist(lapply(grns, function(z) z$Target), use.names = FALSE)))
if (!length(fallback_universe)) stop('The global target-gene candidate background is empty.')
message('Global GRN-derived background (EXPLORATORY): ', length(fallback_universe), ' target genes')

results <- list()
venn_jobs <- list()
idx <- 0L
for (s in stress_names) {
  grn <- grns[[s]]
  deg <- degs[[s]]
  f <- background_files[[s]]
  if (is.null(f) || is.na(f) || !nzchar(f)) {
    universe <- fallback_universe
    bg_mode <- 'EXPLORATORY: union of Target genes from all three input GRNs'
  } else {
    universe <- read_background(file.path(input_root, f))
    bg_mode <- paste0('User-defined eligible tested genes: ', f)
  }
  if (length(universe) == 0L) stop('The evaluation background is empty for: ', s)
  # Each workbook contains DEG targets for exactly ONE TF.
  # Do not search/iterate over multiple TFs, and do not silently take TF[1]
  # if the data are unexpectedly mixed.
  deg_tfs <- unique(deg$TF)
  if (length(deg_tfs) != 1L) {
    stop('[', s, '] The DEG workbook contains ', length(deg_tfs),
         ' distinct TFs; exactly one TF is required per workbook: ',
         paste(head(deg_tfs, 8L), collapse=', '))
  }
  tf <- deg_tfs[[1L]]
  if (!(tf %in% grn$TF)) {
    stop('[', s, '] The unique DEG TF ', tf,
         ' was not found in the corresponding GRN. Verify that the original TF ID matches the GRN; ',
         ' only the Cold DEG Target column is converted from RAP to GID; all IDs are normalized to uppercase.')
  }
  message('[', s, '] Single DEG TF: ', tf, ' (', tf_display_names[[s]], ')')
  idx <- idx + 1L
  z <- analyze_tf(s, tf, grn, deg, universe, bg_mode, idx)
  z$summary[, TF_Display_Name := tf_display_names[[s]]]
  results[[idx]] <- z$summary
  venn_jobs[[idx]] <- list(stress=s,tf=tf,tf_display=tf_display_names[[s]],
                           pred=z$grn_targets,ref=z$deg_targets)
}

ans <- rbindlist(results, fill=TRUE)
# Exactly 3 planned tests: one unique DEG TF per stress.
# BH adjustment is applied across these three tests.
ans[, P_adjust_FDR := p.adjust(P_Hyper_OneSided, method = 'BH')]
ans[, Significance := fifelse(P_adjust_FDR < .001, '***',
                              fifelse(P_adjust_FDR < .01, '**',
                                      fifelse(P_adjust_FDR < .05, '*', 'ns')))]
setcolorder(ans, c('Stress','TF_Display_Name','TF','Background_Mode','N_Universe',
                   'N_GRN_Targets','N_DEG_Targets','N_Overlap',
                   'GRN_Overlap_Fraction','DEG_Recovery_Fraction',
                   'Expected_Overlap','Fold_Enrichment','FE_L95','FE_U95',
                   'Odds_Ratio','OR_Fisher_L95','OR_Fisher_U95',
                   'P_Hyper_OneSided','P_adjust_FDR','Significance'))

fwrite(ans, file.path(output_dir, 'DEG_Overlap_Full_Statistics.tsv'), sep='\t', na='NA')
fwrite(ans, file.path(output_dir, 'DEG_Overlap_Full_Statistics.csv'), na='NA')

contingency <- rbindlist(lapply(seq_len(nrow(ans)), function(i) {
  z <- ans[i]
  data.table(Stress=z$Stress, TF=z$TF,
             Cell=c('GRN+ / DEG+','GRN+ / DEG-', 'GRN- / DEG+','GRN- / DEG-'),
             N=c(z$N_Overlap,z$N_GRN_Only,z$N_DEG_Only,z$N_Neither))
}))
fwrite(contingency, file.path(output_dir, 'DEG_2x2_Contingency.csv'))

# =========================== PUBLICATION FIGURES ==============================
# Consistent three-stress labels and order.
panel_levels <- paste(stress_names, unname(tf_display_names[stress_names]), sep = " | ")
panel_levels_top_to_bottom <- rev(panel_levels)

# Figure 1: overlap proportions. The two panels deliberately retain their own
# denominators: GRN support = overlap / GRN targets; DEG recovery = overlap / DEGs.
plot_df <- rbindlist(list(
  ans[, .(Stress, TF, Metric = "GRN support",
          Estimate = GRN_Overlap_Fraction,
          Lower = Overlap_GRN_L95, Upper = Overlap_GRN_U95)],
  ans[, .(Stress, TF, Metric = "DEG recovery",
          Estimate = DEG_Recovery_Fraction,
          Lower = Overlap_DEG_L95, Upper = Overlap_DEG_U95)]
))
plot_df[, TF_Display_Name := unname(tf_display_names[as.character(Stress)])]
plot_df[, Panel := paste(Stress, TF_Display_Name, sep = " | ")]
plot_df[, Panel := factor(Panel, levels = panel_levels_top_to_bottom)]
plot_df[, Metric := factor(Metric, levels = c("GRN support", "DEG recovery"))]
plot_df[, Stress := factor(Stress, levels = stress_names)]

p_pct <- ggplot(plot_df, aes(x = Estimate * 100, y = Panel, color = Stress)) +
  geom_segment(
    aes(x = Lower * 100, xend = Upper * 100, y = Panel, yend = Panel),
    linewidth = 0.75, na.rm = TRUE
  ) +
  geom_point(size = 3.0) +
  geom_text(
    aes(label = sprintf("%.1f%%", Estimate * 100)),
    nudge_x = 3.0, hjust = 0, size = 3.0, show.legend = FALSE
  ) +
  facet_wrap(~Metric, nrow = 1) +
  scale_color_manual(values = stress_pal, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = c(0, 115), clip = "off") +
  labs(
    x = "Overlap proportion",
    y = NULL,
    title = "Concordance between gsGNA-predicted targets and TF perturbation DEGs",
    subtitle = "Point estimates with exploratory 95% target-gene bootstrap confidence intervals",
    caption = paste(
      "GRN support = overlap / predicted GRN targets; DEG recovery = overlap / perturbation DEGs.",
      "The two panels therefore have different denominators."
    )
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.45),
    strip.text = element_text(face = "bold", color = "black", size = 9.5),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 8.8),
    axis.title.x = element_text(face = "bold", margin = margin(t = 7)),
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 9, color = "grey25"),
    plot.caption = element_text(size = 7.5, color = "grey30", hjust = 0),
    plot.margin = margin(7, 18, 7, 7)
  )

ggsave(
  file.path(output_dir, "DEG_Overlap_Proportions.pdf"), p_pct,
  width = 7.6, height = 3.7, device = grDevices::pdf,
  useDingbats = FALSE, bg = "white"
)
ggsave(
  file.path(output_dir, "DEG_Overlap_Proportions.png"), p_pct,
  width = 7.6, height = 3.7, dpi = 600, bg = "white"
)

# Figure 2: forest plots for fold enrichment (FE) and odds ratio (OR).
# FE CI comes from the exploratory target-gene bootstrap; OR CI is Fisher exact.
effect_df <- rbindlist(list(
  ans[, .(Stress, TF, Metric = "Fold enrichment",
          Estimate = Fold_Enrichment, Lower = FE_L95, Upper = FE_U95)],
  ans[, .(Stress, TF, Metric = "Odds ratio",
          Estimate = Odds_Ratio, Lower = OR_Fisher_L95, Upper = OR_Fisher_U95)]
))
effect_df[, TF_Display_Name := unname(tf_display_names[as.character(Stress)])]
effect_df[, Panel := paste(Stress, TF_Display_Name, sep = " | ")]
effect_df[, Panel := factor(Panel, levels = panel_levels_top_to_bottom)]
effect_df[, Metric := factor(Metric, levels = c("Fold enrichment", "Odds ratio"))]
effect_df[, Stress := factor(Stress, levels = stress_names)]

# A logarithmic axis cannot display zero or non-finite values. Such estimates
# remain in DEG_Overlap_Full_Statistics.csv and are not silently changed.
effect_df <- effect_df[is.finite(Estimate) & Estimate > 0]
effect_df[!is.finite(Lower) | Lower <= 0, Lower := NA_real_]
effect_df[!is.finite(Upper) | Upper <= 0, Upper := NA_real_]

if (nrow(effect_df)) {
  p_effect <- ggplot(effect_df, aes(x = Estimate, y = Panel, color = Stress)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey40", linewidth = 0.55) +
    geom_segment(
      aes(x = Lower, xend = Upper, y = Panel, yend = Panel),
      linewidth = 0.75, na.rm = TRUE
    ) +
    geom_point(size = 3.0) +
    facet_wrap(~Metric, scales = "free_x", nrow = 1) +
    scale_color_manual(values = stress_pal, drop = FALSE) +
    scale_x_log10(
      breaks = scales::breaks_log(n = 4),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    labs(
      x = "Effect size (log10 scale)",
      y = NULL,
      title = "Enrichment of gsGNA-predicted targets among TF perturbation DEGs",
      subtitle = "Reference line indicates no enrichment (effect size = 1)",
      caption = paste(
        "FE: exploratory target-gene bootstrap 95% CI; OR: Fisher exact 95% CI.",
        "FE = 0 cannot be displayed on a logarithmic axis and is retained in the statistics table."
      )
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.45),
      strip.text = element_text(face = "bold", color = "black", size = 9.5),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 8.8),
      axis.title.x = element_text(face = "bold", margin = margin(t = 7)),
      plot.title = element_text(face = "bold", size = 11.5),
      plot.subtitle = element_text(size = 9, color = "grey25"),
      plot.caption = element_text(size = 7.5, color = "grey30", hjust = 0),
      plot.margin = margin(7, 10, 7, 7)
    )
  
  ggsave(
    file.path(output_dir, "DEG_Enrichment_FE_OR.pdf"), p_effect,
    width = 7.6, height = 3.7, device = grDevices::pdf,
    useDingbats = FALSE, bg = "white"
  )
  ggsave(
    file.path(output_dir, "DEG_Enrichment_FE_OR.png"), p_effect,
    width = 7.6, height = 3.7, dpi = 600, bg = "white"
  )
}

format_venn_q <- function(x) {
  if (!length(x) || !is.finite(x)) return("NA")
  if (x == 0) return("<1e-300")
  if (x < .001) return(formatC(x, format = "e", digits = 1))
  formatC(x, format = "f", digits = 3)
}

make_venn_grob <- function(job) {
  n_grn <- length(job$pred)
  n_deg <- length(job$ref)
  n_inter <- length(intersect(job$pred, job$ref))
  
  # With a zero-sized set, draw a textual placeholder instead of asking the
  # VennDiagram package to draw a degenerate pair of circles.
  if (n_grn == 0L || n_deg == 0L) return(NULL)
  
  VennDiagram::draw.pairwise.venn(
    area1 = n_grn,
    area2 = n_deg,
    cross.area = n_inter,
    category = c("Predicted targets", "Perturbation DEGs"),
    euler.d = TRUE,
    scaled = FALSE,
    offset = 0.15,
    fill = c(unname(venn_grn_pal[[job$stress]]), venn_deg_fill),
    alpha = c(0.75, 0.75),
    col = c("black", "black"),
    lwd = c(1.15, 1.15),
    cex = c(1.05, 1.05, 1.05),
    fontface = c("bold", "bold", "bold"),
    cat.cex = c(0.88, 0.88),
    cat.fontface = c("bold", "bold"),
    cat.pos = c(0, 180),
    cat.dist = c(0.035, 0.035),
    cat.col = c("black", "black"),
    ext.text = FALSE,
    ind = FALSE
  )
}

# Draw inside the CURRENT grid viewport, so the same function works for
# one figure per stress and a compact horizontal three-stress figure.
draw_venn_panel <- function(job, info, compact = FALSE) {
  panel_title <- paste0(job$stress, " | ", job$tf_display)
  fe <- info$Fold_Enrichment[[1L]]
  fdr <- info$P_adjust_FDR[[1L]]
  signif <- info$Significance[[1L]]
  fe_text <- if (is.finite(fe)) sprintf("%.2f", fe) else "NA"
  stat_label <- paste0(
    "FE = ", fe_text,
    "     BH q = ", format_venn_q(fdr),
    "     ", signif
  )
  
  grid::grid.rect(gp = grid::gpar(
    col = "grey55", fill = "white", lwd = 0.75
  ))
  grid::grid.text(
    panel_title, x = 0.50, y = 0.927,
    gp = grid::gpar(
      fontsize = if (compact) 9.5 else 14,
      fontface = "bold", col = "black"
    )
  )
  
  grob <- make_venn_grob(job)
  grid::pushViewport(grid::viewport(
    x = 0.50, y = 0.535,
    width = if (compact) 0.97 else 0.87,
    height = 0.675
  ))
  if (is.null(grob)) {
    grid::grid.text(
      "An empty target set; counts reported in the table",
      gp = grid::gpar(fontsize = if (compact) 8 else 11)
    )
  } else {
    grid::grid.draw(grob)
  }
  grid::popViewport()
  
  grid::grid.text(
    stat_label, x = 0.50, y = 0.125,
    gp = grid::gpar(
      fontsize = if (compact) 7.9 else 11,
      fontface = "bold", col = "black"
    )
  )
  grid::grid.text(
    "Schematic areas; exact counts",
    x = 0.50, y = 0.066,
    gp = grid::gpar(fontsize = if (compact) 7.1 else 9, col = "grey35")
  )
}

save_grid_venn <- function(path, width, height, draw_fun) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    grDevices::pdf(
      file = path, width = width, height = height,
      useDingbats = FALSE, bg = "white"
    )
  } else {
    grDevices::png(
      filename = path, width = width, height = height,
      units = "in", res = 600, bg = "white"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  draw_fun()
}

if (plot_venns) {
  # Verify consistency between the figures and the previously saved analysis.
  venn_stats <- rbindlist(lapply(venn_jobs, function(job) {
    data.table(
      Stress = job$stress,
      TF_Display_Name = job$tf_display,
      TF = job$tf,
      N_GRN_Targets = length(job$pred),
      N_DEG_Targets = length(job$ref),
      N_Overlap = length(intersect(job$pred, job$ref))
    )
  }))
  comparison <- merge(
    venn_stats,
    ans[, .(Stress, TF, N_GRN_Targets, N_DEG_Targets, N_Overlap)],
    by = c("Stress", "TF"), suffixes = c("_Venn", "_Analysis")
  )
  stopifnot(
    nrow(comparison) == length(venn_jobs),
    all(comparison$N_GRN_Targets_Venn == comparison$N_GRN_Targets_Analysis),
    all(comparison$N_DEG_Targets_Venn == comparison$N_DEG_Targets_Analysis),
    all(comparison$N_Overlap_Venn == comparison$N_Overlap_Analysis)
  )
  fwrite(venn_stats, file.path(output_dir, "DEG_Venn_ThreeStress_Counts.csv"))
  
  # 1) One combined figure: three panels, matching the reference drawing style.
  draw_combined_venns <- function() {
    grid::grid.text(
      "Overlap of predicted targets and TF perturbation DEGs",
      x = 0.5, y = 0.974,
      gp = grid::gpar(fontsize = 14.5, fontface = "bold")
    )
    grid::grid.text(
      "One TF-specific comparison per stress; areas are not proportional to set size",
      x = 0.5, y = 0.025,
      gp = grid::gpar(fontsize = 9, col = "grey30")
    )
    for (i in seq_along(venn_jobs)) {
      job <- venn_jobs[[i]]
      info <- ans[Stress == job$stress & TF == job$tf]
      stopifnot(nrow(info) == 1L)
      grid::pushViewport(grid::viewport(
        x = (i - 0.5) / length(venn_jobs),
        y = 0.51,
        width = (1 / length(venn_jobs)) - 0.008,
        height = 0.855
      ))
      draw_venn_panel(job, info, compact = TRUE)
      grid::popViewport()
    }
  }
  
  save_grid_venn(
    file.path(output_dir, "DEG_Venn_ThreeStress.pdf"),
    width = 10.5, height = 4.3, draw_fun = draw_combined_venns
  )
  save_grid_venn(
    file.path(output_dir, "DEG_Venn_ThreeStress.png"),
    width = 10.5, height = 4.3, draw_fun = draw_combined_venns
  )
  
  # 2) Three individual figures, as in the supplied reference code.
  for (job in venn_jobs) {
    info <- ans[Stress == job$stress & TF == job$tf]
    filename <- paste0(job$stress, "_", job$tf_display, "_DEG_Venn")
    draw_single_venn <- local({
      local_job <- job
      local_info <- info
      function() {
        grid::pushViewport(grid::viewport(
          x = 0.5, y = 0.5, width = 0.96, height = 0.96
        ))
        draw_venn_panel(local_job, local_info, compact = FALSE)
        grid::popViewport()
      }
    })
    save_grid_venn(
      file.path(output_dir, paste0(filename, ".pdf")),
      width = 7, height = 5.6, draw_fun = draw_single_venn
    )
    save_grid_venn(
      file.path(output_dir, paste0(filename, ".png")),
      width = 7, height = 5.6, draw_fun = draw_single_venn
    )
  }
}

stopifnot(nrow(ans) == 3L)
message('\n[DONE] Three single-TF DEG-overlap analyses and figures saved in: ', output_dir)
print(ans[, .(Stress,TF_Display_Name,TF,N_Universe,N_GRN_Targets,N_DEG_Targets,N_Overlap,
              GRN_Overlap_Fraction,DEG_Recovery_Fraction,Fold_Enrichment,
              FE_L95,FE_U95,Odds_Ratio,OR_Fisher_L95,OR_Fisher_U95,
              P_adjust_FDR)], digits=5)
