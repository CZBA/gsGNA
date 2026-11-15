# gsGNA TFBS overlap and enrichment for rice abiotic stress networks
#
# Run from the repository root:
#   Rscript "code/Network validation/TFBS_validation.R"
# Set GSGNA_INPUT_ROOT to the input root used by the ensemble comparison.
# Place TFTarget.txt, with a header row and TF/Target columns, in that root.
#
# Uses the top 10% of the ranked XGBoost selected-union edges per stress.
# Reports overlap, fold enrichment, odds ratios and TF-cluster bootstrap CIs.
# The target background should represent loci covered by the motif scan.
# TFBS motif predictions provide computational support for inferred edges.

required <- c('data.table', 'ggplot2')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install the required R packages with install.packages(c(', paste(shQuote(missing),collapse=','), '))')
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# =========================== SETTINGS ========================================
input_root <- Sys.getenv('GSGNA_INPUT_ROOT', unset='.')
ensemble_dir <- file.path(input_root, 'gsGNA_FixedTop10')
# This order determines the bootstrap seed assigned to each stress.
stress_names <- c('Drought', 'Alkaline', 'Cold')
network_files <- setNames(file.path(
  ensemble_dir, 'FINAL_Locked_Test_And_Network', stress_names, 'Final_Ranked_Network',
  'XGBoost_Ranked_SelectedUnion_AllEdges.csv'
), stress_names)
network_edge_fraction <- 0.10
tfbs_file <- 'TFTarget.txt'

# Optional file with a Gene or Target column of eligible motif-scanned genes.
# Paths for the TFBS resource and optional coverage files are relative to input_root.
background_targets_file <- NULL  # e.g. 'TFBS_Scanned_Eligible_Genes.txt'

# If the TFBS resource provides a list of TFs scanned including zero TFBS hits,
# use it here. Otherwise only TFs with at least one reference TFBS pair are used.
assayed_tfs_file <- NULL          # e.g. 'TFBS_Scanned_TFs.txt'

n_tf_boot <- 2000L               # TF-cluster percentile CI, NOT edge bootstrap
bootstrap_seed <- 123L
min_boot_tfs <- 3L
output_dir <- 'gsGNA_TFBS_Enrichment'
# ============================================================================

out <- file.path(input_root, output_dir)

clean_ids <- function(z) toupper(trimws(as.character(z)))
key <- function(s) gsub('[^a-z0-9]', '', tolower(s))
locate_col <- function(dt, possible, fallback, role) {
  ix <- which(key(names(dt)) %in% possible)
  if (length(ix) > 1L) stop(role, ' matches multiple columns: ', paste(names(dt)[ix],collapse=', '))
  if (length(ix) == 1L) ix else fallback
}
clean_pairs <- function(tf, target) {
  z <- data.table(TF=clean_ids(tf), Target=clean_ids(target))
  z <- z[!is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) &
           TF != 'NA' & Target != 'NA' & TF != Target]
  z <- unique(z, by=c('TF','Target'))
  setkey(z, TF, Target)
  z
}

read_pairs <- function(path, role) {
  if (!file.exists(path)) stop('Missing ', role, ' file: ', path)
  dt <- fread(path, header=TRUE, check.names=FALSE, showProgress=FALSE)
  if (ncol(dt) < 2L) stop('The ', role, ' file needs TF and Target columns: ', path)
  i <- locate_col(dt, c('tf','tfid','tfgene','regulator','regulatorygene'), 1L, 'TF')
  j <- locate_col(dt, c('target','targetgene','targetid','targetgeneid','regulatedgene'), 2L, 'Target')
  if (i == j) stop('TF and Target resolve to the same column in ', role, ': ', path)
  pairs <- clean_pairs(dt[[i]], dt[[j]])
  message(sprintf('[%s] %s | %s unique TF-target pairs', role, path, nrow(pairs)))
  pairs
}

read_ranked_network <- function(path, stress) {
  dt <- fread(path, header=TRUE, select=c('TF', 'Target', 'EdgeWeight'), showProgress=FALSE)
  if (!all(c('TF', 'Target', 'EdgeWeight') %in% names(dt)) || !nrow(dt)) {
    stop('Expected a nonempty ranked network with TF, Target, and EdgeWeight: ', path)
  }
  dt[, EdgeWeight := suppressWarnings(as.numeric(EdgeWeight))]
  if (any(!is.finite(dt$EdgeWeight))) stop('Non-finite ranked edge weights: ', path)
  # The ensemble script also writes a summary that identifies partial exports.
  export_summary <- file.path(dirname(path), 'Final_Network_Export_Summary.csv')
  if (file.exists(export_summary)) {
    audit <- fread(export_summary, showProgress=FALSE)
    required <- c('File', 'Exported_Edges', 'N_Union_Candidates')
    if (!all(required %in% names(audit)) || nrow(audit) != 1L ||
        !isTRUE(audit$File == basename(path)) || !isTRUE(audit$Exported_Edges == nrow(dt)) ||
        !isTRUE(audit$N_Union_Candidates == nrow(dt))) {
      stop('Expected the full ranked selected-union network: ', path)
    }
  }
  setorderv(dt, c('EdgeWeight', 'TF', 'Target'), c(-1L, 1L, 1L))
  n_keep <- ceiling(nrow(dt) * network_edge_fraction)
  retained <- dt[seq_len(n_keep)]
  pairs <- clean_pairs(retained$TF, retained$Target)
  if (nrow(pairs) != n_keep) stop('Duplicate or invalid edges in the retained network: ', path)
  message(sprintf('[%s] Retained %s of %s ranked GRN edges', stress, n_keep, nrow(dt)))
  pairs
}

read_ids <- function(path, role) {
  if (!file.exists(path)) stop(role,' file does not exist: ',path)
  dt <- fread(path, header=TRUE, showProgress=FALSE)
  if (!ncol(dt)) stop(role,' file has no columns: ',path)
  idx <- locate_col(dt,
                    if (role=='targets') c('target','gene','geneid','targetgene','locus')
                    else c('tf','tfid','regulator','tfgene'), 1L, role)
  ids <- unique(clean_ids(dt[[idx]]))
  ids[!is.na(ids) & nzchar(ids) & ids!='NA']
}

# Return 2x2 effect sizes from counts that already refer to the SAME universe.
calc_effect <- function(n_universe, n_pred, n_ref, n_overlap) {
  a <- as.numeric(n_overlap)
  b <- as.numeric(n_pred - a)
  c <- as.numeric(n_ref - a)
  d <- as.numeric(n_universe - a - b - c)
  if (any(!is.finite(c(a,b,c,d))) || any(c(a,b,c,d) < 0))
    stop('Invalid 2x2 table. Check the evaluable TF and target universes.')
  n <- a+b+c+d
  support <- if (a+b>0) a/(a+b) else NA_real_
  reference_coverage <- if (a+c>0) a/(a+c) else NA_real_
  prevalence <- if (n>0) (a+c)/n else NA_real_
  expected <- if (n>0) (a+b)*(a+c)/n else NA_real_
  fold <- if (is.finite(prevalence) && prevalence>0 && is.finite(support))
    support/prevalence else NA_real_

  # Edge-level approximate log(OR) CI; dependency across edges means these CIs
  # are potentially too narrow. TF-cluster bootstrapped CIs are preferred below.
  tab <- c(a,b,c,d)
  corrected <- any(tab == 0)
  if (corrected) tab <- tab + .5
  OR <- (tab[1]*tab[4])/(tab[2]*tab[3])
  se_logOR <- sqrt(sum(1/tab))
  OR_l <- exp(log(OR)-qnorm(.975)*se_logOR)
  OR_u <- exp(log(OR)+qnorm(.975)*se_logOR)
  # Risk ratio here compares TFBS fraction in predicted vs non-predicted pairs.
  RR <- if ((c+d)>0 && c>0 && a+b>0) (a/(a+b))/(c/(c+d)) else NA_real_

  list(a=a,b=b,c=c,d=d, support=support, tfbs_coverage=reference_coverage,
       prevalence=prevalence, expected=expected, fold=fold,
       OR=OR, OR_Wald_L95=OR_l, OR_Wald_U95=OR_u,
       RR=RR, Haldane_Correction=corrected)
}

boot_tf <- function(tf_stats, reps, seed) {
  names_of_metrics <- c('support','tfbs_coverage','fold','OR','RR')
  out <- as.list(setNames(rep(NA_real_,2L*length(names_of_metrics)),
                          as.vector(rbind(paste0(names_of_metrics,'_TFboot_L95'),
                                          paste0(names_of_metrics,'_TFboot_U95')))))
  out$N_Valid_TFBoot <- 0L
  n_tf <- nrow(tf_stats)
  if (n_tf < min_boot_tfs || reps < 20L) return(out)
  m <- as.matrix(tf_stats[,.(Universe_Pairs, GRN_Pairs, TFBS_Pairs, Overlap_Pairs)])
  storage.mode(m) <- 'double'
  set.seed(seed)
  vals <- matrix(NA_real_, nrow=reps, ncol=length(names_of_metrics),
                 dimnames=list(NULL,names_of_metrics))
  for (b in seq_len(reps)) {
    chosen <- sample.int(n_tf, size=n_tf, replace=TRUE)
    sums <- colSums(m[chosen,,drop=FALSE])
    z <- calc_effect(sums[1],sums[2],sums[3],sums[4])
    vals[b,] <- unlist(z[names_of_metrics],use.names=FALSE)
  }
  out$N_Valid_TFBoot <- sum(is.finite(vals[,'OR']))
  for (nm in names_of_metrics) {
    x <- vals[,nm]
    x <- x[is.finite(x)]
    if (length(x) >= max(20L,ceiling(reps*.75))) {
      out[[paste0(nm,'_TFboot_L95')]] <- as.numeric(quantile(x,.025,names=FALSE))
      out[[paste0(nm,'_TFboot_U95')]] <- as.numeric(quantile(x,.975,names=FALSE))
    }
  }
  out
}

# Minimal internal sanity check: known 2x2 counts must reproduce known metrics.
.selftest <- calc_effect(n_universe=100, n_pred=20, n_ref=10, n_overlap=5)
stopifnot(.selftest$a == 5, .selftest$b == 15, .selftest$c == 5,
          .selftest$d == 75, isTRUE(all.equal(.selftest$fold, 2.5)),
          isTRUE(all.equal(.selftest$OR, 5)))
rm(.selftest)

# ------------------------------- LOAD ---------------------------------------
missing_networks <- network_files[!file.exists(network_files)]
if (length(missing_networks)) {
  stop('Missing ranked XGBoost networks:\n', paste(missing_networks, collapse='\n'))
}
grns <- setNames(lapply(stress_names, function(stress) {
  read_ranked_network(network_files[[stress]], stress)
}), stress_names)
tfbs <- read_pairs(file.path(input_root,tfbs_file),'TFBS')

if (is.null(background_targets_file)) {
  target_universe <- unique(c(tfbs$Target, unlist(lapply(grns,`[[`,'Target'),
                                                  use.names=FALSE)))
  background_mode <- 'EXPLORATORY: data-derived union of TFBS/GRN targets'
  warning('Using the union-derived target background. Enrichment estimates are exploratory.')
} else {
  target_universe <- read_ids(file.path(input_root,background_targets_file),'targets')
  background_mode <- paste0('User-specified independent target universe: ',background_targets_file)
}
if (!length(target_universe)) stop('The evaluable target universe is empty.')
assayed_tfs <- if (is.null(assayed_tfs_file)) unique(tfbs$TF) else
  read_ids(file.path(input_root,assayed_tfs_file),'tfs')
if (!length(assayed_tfs)) stop('The TFBS reference contains no evaluable TFs.')

dir.create(out, recursive=TRUE, showWarnings=FALSE)
res <- list(); per_tf_list <- list(); contingency <- list(); plot_count <- list()
for (ix in seq_along(grns)) {
  ds <- names(grns)[ix]
  pred_all <- grns[[ds]]
  eval_tfs <- intersect(unique(pred_all$TF),assayed_tfs)
  if (!length(eval_tfs)) {
    warning(ds, ' has no evaluable TFs in common with the TFBS reference.')
    next
  }
  pred <- pred_all[TF %in% eval_tfs & Target %in% target_universe]
  ref <- tfbs[TF %in% eval_tfs & Target %in% target_universe]
  # Because TFBS and GRN use the same pair-key definition, no string
  # concatenation with separators can introduce spurious collisions.
  overlap <- pred[ref,on=.(TF,Target),nomatch=0L]
  overlap <- unique(overlap,by=c('TF','Target'))
  fwrite(overlap,file.path(out,paste0(ds,'_TFBS_Supported_Pairs.csv')))

  tf_summary <- data.table(TF=sort(eval_tfs))
  tf_summary[,Universe_Pairs:=as.numeric(length(target_universe) -
                                           as.integer(TF %in% target_universe))]
  pred_count <- pred[,.N,by=TF]
  ref_count <- ref[,.N,by=TF]
  over_count <- overlap[,.N,by=TF]
  tf_summary[,GRN_Pairs:=as.numeric(pred_count$N[match(TF,pred_count$TF)])]
  tf_summary[,TFBS_Pairs:=as.numeric(ref_count$N[match(TF,ref_count$TF)])]
  tf_summary[,Overlap_Pairs:=as.numeric(over_count$N[match(TF,over_count$TF)])]
  for (name in c('GRN_Pairs','TFBS_Pairs','Overlap_Pairs')) {
    set(tf_summary,which(is.na(tf_summary[[name]])),name,0)
  }
  tf_summary[,`:=`(Stress=ds, GRN_Supported_Fraction=ifelse(GRN_Pairs>0,
                                                            Overlap_Pairs/GRN_Pairs,NA_real_))]
  per_tf_list[[ds]] <- tf_summary
  counts <- tf_summary[,lapply(.SD,sum),.SDcols=c('Universe_Pairs','GRN_Pairs',
                                                  'TFBS_Pairs','Overlap_Pairs')]
  z <- calc_effect(counts$Universe_Pairs,counts$GRN_Pairs,
                   counts$TFBS_Pairs,counts$Overlap_Pairs)
  ci <- boot_tf(tf_summary,n_tf_boot,bootstrap_seed + ix*29L)

  summary <- data.table(
    Stress=ds, Reference='TFBS motif-supported candidate pairs',
    Prediction_Input=sprintf('Top %g%% of ranked XGBoost selected-union edges',
                             100*network_edge_fraction),
    Prediction_File=network_files[[ds]],
    Background_Mode=background_mode,
    TFBS_TF_Coverage_Mode=if(is.null(assayed_tfs_file))
      'Only TFs with >=1 TFBS database pair; zero-hit assayed TFs unknown'
    else paste('User-supplied assayed TF list:',assayed_tfs_file),
    N_All_GRN_TFs=uniqueN(pred_all$TF),
    N_Evaluable_TFs=length(eval_tfs),
    N_Universe_Targets=length(target_universe),
    N_GRN_Pairs_Original=nrow(pred_all),
    N_GRN_Pairs_Excluded=nrow(pred_all)-nrow(pred),
    N_TFBS_Pairs_Eligible=nrow(ref),
    N_Universe_Pairs=counts$Universe_Pairs,
    N_GRN_Pairs=counts$GRN_Pairs,
    N_TFBS_Pairs=counts$TFBS_Pairs,
    N_Overlap=counts$Overlap_Pairs,
    N_Predicted_Without_TFBS=z$b,
    N_TFBS_Without_Prediction=z$c,
    N_Neither=z$d,
    Predicted_TFBS_Support_Fraction=z$support,
    TFBS_Reference_Coverage_Fraction=z$tfbs_coverage,
    Background_TFBS_Fraction=z$prevalence,
    Expected_Overlap_Under_Uniform_Null=z$expected,
    Fold_Enrichment=z$fold,
    Odds_Ratio=z$OR,
    OR_Wald_L95=z$OR_Wald_L95,
    OR_Wald_U95=z$OR_Wald_U95,
    Risk_Ratio=z$RR,
    Used_Haldane_Correction=z$Haldane_Correction,
    P_Hypergeometric_Exploratory=if (z$a>0) {
      phyper(z$a-1,z$a+z$c,z$b+z$d,z$a+z$b,lower.tail=FALSE)
    } else 1.0,
    N_TF_Bootstrap=n_tf_boot
  )
  for (nm in names(ci)) summary[,(nm):=ci[[nm]]]
  res[[ds]] <- summary
  contingency[[ds]] <- data.table(
    Stress=ds, Cell=c('GRN+ / TFBS+','GRN+ / TFBS-',
                      'GRN- / TFBS+','GRN- / TFBS-'),
    N=c(z$a,z$b,z$c,z$d))
  plot_count[[ds]] <- rbindlist(list(
    data.table(Stress=ds,Denominator='Predicted GRN',
               Category=c('With TFBS support','Without TFBS support'),Count=c(z$a,z$b)),
    data.table(Stress=ds,Denominator='TFBS reference',
               Category=c('With GRN overlap','Without GRN overlap'),Count=c(z$a,z$c))))
  message(sprintf('[%s] overlap=%s / GRN=%s (%.2f%%); TFBS=%s; fold=%.3f; OR=%.3f',
                  ds,z$a,z$a+z$b,100*z$support,z$a+z$c,z$fold,z$OR))
}
if (!length(res)) stop('No evaluable TFBS comparisons for the three stresses. Check the gene IDs.')
ans <- rbindlist(res,fill=TRUE)
fwrite(ans,file.path(out,'TFBS_Enrichment_Summary.csv'))
fwrite(rbindlist(per_tf_list,fill=TRUE),file.path(out,'TFBS_Enrichment_PerTF.csv'))
fwrite(rbindlist(contingency,fill=TRUE),file.path(out,'TFBS_2x2_Contingency.csv'))
fwrite(rbindlist(plot_count,fill=TRUE),file.path(out,'TFBS_Overlap_Counts_Denominators.csv'))

# ------------------------------- FIGURES -------------------------------------
plot_df <- rbindlist(list(
  ans[,.(Stress,Metric='Predicted-edge TFBS support',Estimate=Predicted_TFBS_Support_Fraction,
         Lower=support_TFboot_L95,Upper=support_TFboot_U95)],
  ans[,.(Stress,Metric='Reference TFBS recovery',Estimate=TFBS_Reference_Coverage_Fraction,
         Lower=tfbs_coverage_TFboot_L95,Upper=tfbs_coverage_TFboot_U95)]))
plot_df[,Stress:=factor(Stress,levels=names(network_files))]
plot_df[,Metric:=factor(Metric,levels=c('Predicted-edge TFBS support',
                                        'Reference TFBS recovery'))]
p1 <- ggplot(plot_df,aes(Stress,Estimate,fill=Stress))+
  geom_col(width=.55)+
  geom_errorbar(aes(ymin=pmax(0,Lower),ymax=pmin(1,Upper)),width=.16,
                linewidth=.65,na.rm=TRUE)+
  geom_text(aes(label=sprintf('%.1f%%',100*Estimate)),vjust=-.5,size=3.8)+
  facet_wrap(~Metric,nrow=1)+
  scale_y_continuous(labels=function(x) sprintf('%.0f%%',100*x),
                     breaks=seq(0,1,.2))+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_manual(values=c(Drought='#2C967D',Alkaline='#D97C6D',Cold='#6097E7'))+
  labs(x=NULL,y='Fraction (denominator specified in each facet)',
       title='TFBS support and reference coverage of inferred GRN edges',
       subtitle='Error bars: 95% TF-cluster bootstrap intervals',
       caption='TFBS motif support is computational; it does not establish direct regulation.')+
  theme_bw(base_size=12)+theme(legend.position='none',panel.grid.minor=element_blank(),
                               strip.text=element_text(face='bold'))

# ------------------------- Figure 2: FE/OR forest plot ------------------------
# Read the just-exported TFBS summary, not a potentially stale global `ans`.
# ONLY plotting: point estimates and TF-cluster bootstrap CIs are not recalculated.
render_tfbs_fe_or_forest <- function(summary_file, output_dir=dirname(summary_file)) {
  if (!file.exists(summary_file)) stop('Missing TFBS summary file: ', summary_file)
  tfbs_summary <- data.table::fread(summary_file, check.names=FALSE, showProgress=FALSE,
                                    encoding='UTF-8')
  data.table::setnames(tfbs_summary, trimws(sub('^\ufeff', '', names(tfbs_summary))))
  needed <- c('Stress','Fold_Enrichment','fold_TFboot_L95','fold_TFboot_U95',
              'Odds_Ratio','OR_TFboot_L95','OR_TFboot_U95')
  absent <- setdiff(needed,names(tfbs_summary))
  if (length(absent)) {
    stop('Missing columns in the TFBS summary: ',paste(absent,collapse=', '),
         '\nSource: ',normalizePath(summary_file,winslash='/',mustWork=TRUE))
  }
  if (!nrow(tfbs_summary)) stop('The TFBS summary has no rows.')
  tfbs_summary[, Stress := trimws(as.character(Stress))]
  stress_order <- stress_names
  unexpected <- setdiff(unique(tfbs_summary$Stress),stress_order)
  if (length(unexpected)) stop('Unknown stress labels: ',paste(unexpected,collapse=', '))
  if (anyDuplicated(tfbs_summary,by='Stress')) stop('Expected one TFBS summary row per stress.')
  stress_colors <- c(Drought='#2C967D',Alkaline='#D97C6D',Cold='#6097E7')
  metric_order <- c('Fold enrichment','Odds ratio')
  effect_df <- data.table::rbindlist(list(
    tfbs_summary[,.(Stress,Metric='Fold enrichment',
                   Estimate=as.numeric(Fold_Enrichment),
                   Lower=as.numeric(fold_TFboot_L95),
                   Upper=as.numeric(fold_TFboot_U95))],
    tfbs_summary[,.(Stress,Metric='Odds ratio',
                   Estimate=as.numeric(Odds_Ratio),
                   Lower=as.numeric(OR_TFboot_L95),
                   Upper=as.numeric(OR_TFboot_U95))]
  ),use.names=TRUE)

  bad_estimate <- !is.finite(effect_df$Estimate) | effect_df$Estimate<=0
  if (any(bad_estimate)) warning(
    'Omitted nonpositive or nonfinite effects from the log-scale plot: ',
    paste(effect_df$Stress[bad_estimate],effect_df$Metric[bad_estimate],
          sep=' / ',collapse='; '))
  effect_df <- effect_df[!bad_estimate]
  if (!nrow(effect_df)) stop('No positive finite effect sizes to plot.')
  effect_df[,CI_Valid := is.finite(Lower)&is.finite(Upper)&
              Lower>0&Upper>0&Lower<=Estimate&Upper>=Estimate&Lower<=Upper]
  if (any(!effect_df$CI_Valid)) warning(
    'Showing point estimates for effects with unavailable confidence intervals: ',
    paste(effect_df$Stress[!effect_df$CI_Valid],
          effect_df$Metric[!effect_df$CI_Valid],sep=' / ',collapse='; '))
  effect_df[,Y := as.numeric(length(stress_order)-match(Stress,stress_order)+1L)]
  effect_df[,Stress := factor(Stress,levels=stress_order)]
  effect_df[,Metric := factor(Metric,levels=metric_order)]
  effect_df[,Label := sprintf('%.2f',Estimate)]
  valid_ci <- effect_df[CI_Valid==TRUE]
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  data.table::fwrite(effect_df,file.path(output_dir,'TFBS_FE_OR_Forest_Plot_Data.csv'))

  p2 <- ggplot2::ggplot(effect_df,ggplot2::aes(x=Estimate,y=Y))+
    ggplot2::geom_vline(xintercept=1,linetype='dashed',color='#6B7280',linewidth=0.55)+
    ggplot2::geom_segment(data=valid_ci,
       ggplot2::aes(x=Lower,xend=Upper,y=Y,yend=Y,color=Stress),
       linewidth=0.85,lineend='round',show.legend=FALSE)+
    ggplot2::geom_segment(data=valid_ci,
       ggplot2::aes(x=Lower,xend=Lower,y=Y-0.075,yend=Y+0.075,color=Stress),
       linewidth=0.75,show.legend=FALSE)+
    ggplot2::geom_segment(data=valid_ci,
       ggplot2::aes(x=Upper,xend=Upper,y=Y-0.075,yend=Y+0.075,color=Stress),
       linewidth=0.75,show.legend=FALSE)+
    ggplot2::geom_point(ggplot2::aes(fill=Stress),shape=21,size=3.5,
                        color='black',stroke=0.45,show.legend=FALSE)+
    ggplot2::geom_text(ggplot2::aes(y=Y+0.25,label=Label),
                       size=3.0,fontface='bold',color='black',show.legend=FALSE)+
    ggplot2::facet_wrap(~Metric,scales='free_x',nrow=1)+
    ggplot2::scale_fill_manual(values=stress_colors,drop=FALSE)+
    ggplot2::scale_color_manual(values=stress_colors,drop=FALSE)+
    ggplot2::scale_x_log10(labels=function(x)
      format(x,trim=TRUE,scientific=FALSE))+
    ggplot2::scale_y_continuous(breaks=seq_along(stress_order),
       labels=rev(stress_order),
       limits=c(0.6,length(stress_order)+0.60),expand=ggplot2::expansion(mult=0))+
    ggplot2::labs(
       x='Effect size (log10 scale)',y=NULL,
       title='TFBS motif-support enrichment of inferred GRN edges',
       subtitle='Points: effect estimates; horizontal lines: 95% TF-cluster bootstrap CIs',
       caption=paste('Dashed line: no enrichment (FE / OR = 1).',
                     'TFBS motif support is computational, not independent experimental validation.'))+
    ggplot2::theme_bw(base_size=10.8)+
    ggplot2::theme(
       legend.position='none',
       text=ggplot2::element_text(color='black'),
       panel.border=ggplot2::element_rect(color='black',fill=NA,linewidth=0.48),
       panel.grid.minor=ggplot2::element_blank(),
       panel.grid.major.y=ggplot2::element_blank(),
       panel.grid.major.x=ggplot2::element_line(color='#E9ECEF',linewidth=0.36),
       strip.background=ggplot2::element_rect(fill='#F2F5F8',color='black',linewidth=0.48),
       strip.text=ggplot2::element_text(face='bold',color='black',size=10.2),
       axis.text=ggplot2::element_text(color='black',size=9.1),
       axis.title.x=ggplot2::element_text(face='bold',color='black',size=10.2,
                                          margin=ggplot2::margin(t=8)),
       plot.title=ggplot2::element_text(face='bold',color='black',size=11.4),
       plot.subtitle=ggplot2::element_text(color='black',size=9),
       plot.caption=ggplot2::element_text(color='black',size=7.8,hjust=0,
                                          margin=ggplot2::margin(t=8)),
       plot.margin=ggplot2::margin(9,10,9,9))

  ggplot2::ggsave(file.path(output_dir,'TFBS_Enrichment_EffectSizes.pdf'),
      plot=p2,width=8.3,height=3.65,bg='white',device=grDevices::pdf,
      useDingbats=FALSE)
  ggplot2::ggsave(file.path(output_dir,'TFBS_Enrichment_EffectSizes.png'),
      plot=p2,width=8.3,height=3.65,dpi=600,bg='white')
  message('[TFBS forest] Source: ',normalizePath(summary_file,winslash='/',mustWork=TRUE))
  message('[TFBS forest] Output: ',normalizePath(output_dir,winslash='/',mustWork=TRUE))
  invisible(p2)
}

p2 <- render_tfbs_fe_or_forest(
  summary_file=file.path(out,'TFBS_Enrichment_Summary.csv'),output_dir=out)

# ------------------ Figure 3: TFBS reference-recovery stacked bars -----------
# Same denominator in every bar: eligible TFBS reference pairs N_TFBS_Pairs.
# Blue  = TFBS reference pairs NOT predicted by the GRN.
# Red   = reference pairs present in the GRN (N_Overlap).
# Percent above each bar = N_Overlap / N_TFBS_Pairs (reference recovery).
# This is NOT the predicted-edge TFBS support proportion (N_Overlap / N_GRN_Pairs).
plot_tfbs_reference_stack <- function(summary_dt, output_dir) {
  required_cols <- c('Stress','N_TFBS_Pairs','N_Overlap',
                     'TFBS_Reference_Coverage_Fraction')
  absent <- setdiff(required_cols,names(summary_dt))
  if (length(absent)) stop('Missing columns in the TFBS summary: ',paste(absent,collapse=', '))

  stress_order <- c('Alkaline','Drought','Cold')
  z <- data.table::as.data.table(data.table::copy(summary_dt))
  z[, Stress := trimws(as.character(Stress))]
  if (anyDuplicated(z,by='Stress') ||
      !setequal(z$Stress,stress_order) || nrow(z)!=3L)
    stop('The stacked bar plot requires one summary row for each of Alkaline, Drought, and Cold.')

  z[, `:=`(TFBS_Total=as.numeric(N_TFBS_Pairs),
           Overlap=as.numeric(N_Overlap),
           Recovery=as.numeric(TFBS_Reference_Coverage_Fraction))]
  if (any(!is.finite(z$TFBS_Total)) || any(!is.finite(z$Overlap)) ||
      any(z$TFBS_Total<=0) || any(z$Overlap<0) ||
      any(z$Overlap>z$TFBS_Total))
    stop('Invalid TFBS reference totals or overlap counts for the stacked bar plot.')
  if (any(!is.finite(z$Recovery)) ||
      any(abs(z$Recovery-z$Overlap/z$TFBS_Total)>1e-6))
    stop('The recorded TFBS reference coverage differs from overlap divided by total reference pairs.')

  z[, TFBS_Only := TFBS_Total-Overlap]
  z[, Stress := factor(Stress,levels=stress_order)]
  z[, TopLabel := sprintf('%.1f%% of TFBS',100*Recovery)]
  max_total <- max(z$TFBS_Total)
  z[, LabelY := TFBS_Total+0.065*max_total]

  stack <- data.table::rbindlist(list(
    z[,.(Stress,Category='TFBS only',Count=TFBS_Only)],
    z[,.(Stress,Category='GRN–TFBS overlap',Count=Overlap)]
  ))
  # ggplot2 stack order is reverse of fill factor: blue below, red above.
  stack[,Category := factor(Category,
                            levels=c('GRN–TFBS overlap','TFBS only'))]
  stack[,CountLabel := ifelse(Count>0,
    format(round(Count),big.mark=',',scientific=FALSE,trim=TRUE),'')]
  data.table::setorder(stack,Stress,Category)

  # Export the exact values underlying the figure; these are not new estimates.
  data.table::fwrite(z[,.(Stress,N_TFBS_Pairs=TFBS_Total,
                          N_Overlap=Overlap,N_TFBS_Only=TFBS_Only,
                          TFBS_Reference_Coverage_Fraction=Recovery)],
                     file.path(output_dir,'TFBS_Reference_Stacked_Overlap_Data.csv'))

  p3 <- ggplot2::ggplot(stack,
                        ggplot2::aes(x=Stress,y=Count,fill=Category)) +
    ggplot2::geom_col(width=0.62,color='white',linewidth=0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label=CountLabel),
      position=ggplot2::position_stack(vjust=0.5),
      size=3.5,fontface='bold',color='black',show.legend=FALSE) +
    ggplot2::geom_text(
      data=z,ggplot2::aes(x=Stress,y=LabelY,label=TopLabel),
      inherit.aes=FALSE,size=3.25,color='black') +
    ggplot2::scale_fill_manual(
      values=c('GRN–TFBS overlap'='#C93E55','TFBS only'='#3587AB'),
      breaks=c('GRN–TFBS overlap','TFBS only'),drop=FALSE) +
    ggplot2::scale_y_continuous(
      labels=function(x) format(x,big.mark=',',scientific=FALSE,trim=TRUE),
      expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::coord_cartesian(ylim=c(0,max(z$LabelY)+0.065*max_total),clip='off') +
    ggplot2::labs(tag='a',x='Stress condition',
                  y='Number of TFBS-reference pairs',
                  fill='Data category',
                  caption='Percentages indicate TFBS reference recovery (overlap / eligible TFBS pairs).') +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(
      text=ggplot2::element_text(color='black'),
      plot.tag=ggplot2::element_text(size=19,face='bold',color='black'),
      axis.text=ggplot2::element_text(color='black',size=10),
      axis.title=ggplot2::element_text(color='black',size=11),
      axis.line=ggplot2::element_line(color='#89939B',linewidth=0.4),
      axis.ticks=ggplot2::element_line(color='#89939B',linewidth=0.4),
      legend.position='bottom',legend.direction='horizontal',
      legend.title=ggplot2::element_text(color='black',size=9.5),
      legend.text=ggplot2::element_text(color='black',size=9.5),
      legend.key.width=grid::unit(4.5,'mm'),
      plot.caption=ggplot2::element_text(color='black',size=8,hjust=0,
                                         margin=ggplot2::margin(t=6)),
      plot.margin=ggplot2::margin(7,12,7,9))

  ggplot2::ggsave(file.path(output_dir,'TFBS_Reference_Overlap_Stacked_Bars.pdf'),
                  plot=p3,width=230,height=160,units='mm',
                  device=grDevices::pdf,useDingbats=FALSE,bg='white')
  ggplot2::ggsave(file.path(output_dir,'TFBS_Reference_Overlap_Stacked_Bars.png'),
                  plot=p3,width=230,height=160,units='mm',dpi=500,bg='white')
  message('[TFBS stacked] TFBS-reference stacked bars saved in: ',output_dir)
  invisible(p3)
}

p3 <- plot_tfbs_reference_stack(ans,out)

# Save the TFBS support and reference recovery proportions.
ggplot2::ggsave(file.path(out,'TFBS_Support_Proportions.pdf'),plot=p1,
                width=11,height=5,bg='white',device=grDevices::pdf,
                useDingbats=FALSE)
ggplot2::ggsave(file.path(out,'TFBS_Support_Proportions.png'),plot=p1,
                width=11,height=5,dpi=320,bg='white')

writeLines(notes,file.path(out,'TFBS_Interpretation_Notes.txt'),useBytes=TRUE)
message('[DONE] Summary, per-TF tables, 2x2 counts, overlap pairs, 3 figures (including stacked bars): ',out)
print(ans[,.(Stress,N_Evaluable_TFs,N_GRN_Pairs,N_TFBS_Pairs,N_Overlap,
             Predicted_TFBS_Support_Fraction,TFBS_Reference_Coverage_Fraction,
             Fold_Enrichment,Odds_Ratio)])
