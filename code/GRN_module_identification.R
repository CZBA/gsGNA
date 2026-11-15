# gsGNA module detection for rice abiotic-stress GRNs
#
# Run from the repository root:
#   Rscript code/GRN_module_identification.R
#
# Set GSGNA_METHOD to run one algorithm, or leave it unset to run all eight.
# Choices: FN, MCL, MCODE, Infomap, Prop, Eigen, Louvain, Walktrap.
# GSGNA_INPUT_ROOT points to the project directory (default: current directory).
# GSGNA_ENSEMBLE_RUN names the input directory under the project directory
# (default: gsGNA_FixedTop10).
# Inputs: precomputed Top 10% XGBoost edge lists with TF and Target columns
# and one score column (EdgeWeight, Weight, XGB_Score, or Score).
# Outputs: results/GRN_Modules/<method>/<stress>/.
#
# Each stress network is clustered independently as an undirected, unweighted
# graph. Retained modules contain 10–300 distinct genes. All modules are saved.

required_packages <- c("data.table", "miRspongeR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Install data.table from CRAN and miRspongeR from Bioconductor.")
}
suppressPackageStartupMessages(library(data.table))

input_root <- Sys.getenv("GSGNA_INPUT_ROOT", ".")
ensemble_run <- Sys.getenv("GSGNA_ENSEMBLE_RUN", "gsGNA_FixedTop10")
output_root <- file.path(input_root, "results", "GRN_Modules")
stresses <- c("Drought", "Alkaline", "Cold")
algorithms <- c(
  FN = "FN",
  MCL = "MCL",
  MCODE = "MCODE",
  Infomap = "infomap",
  Prop = "prop",
  Eigen = "eigen",
  Louvain = "louvain",
  Walktrap = "walktrap"
)
requested_method <- Sys.getenv("GSGNA_METHOD")
if (nzchar(requested_method) && !requested_method %in% names(algorithms)) {
  stop("GSGNA_METHOD must be one of: ", paste(names(algorithms), collapse = ", "))
}
selected_methods <- if (nzchar(requested_method)) requested_method else names(algorithms)
thread_setting <- Sys.getenv("GSGNA_DT_THREADS", "16")
if (!grepl("^[1-9][0-9]*$", thread_setting)) {
  stop("GSGNA_DT_THREADS must be a positive integer.")
}
dt_threads <- suppressWarnings(as.integer(thread_setting))
if (is.na(dt_threads)) stop("GSGNA_DT_THREADS exceeds the integer range.")
setDTthreads(dt_threads)

seed_base <- 123L
min_module_size <- 10L
max_module_size <- 300L

network_files <- setNames(
  file.path(
    input_root, ensemble_run, "FINAL_Locked_Test_And_Network",
    stresses, "Final_Ranked_Network",
    "XGBoost_SelectedUnion_Top10Percent_Edges.csv"
  ),
  stresses
)
missing_files <- network_files[!file.exists(network_files)]
if (length(missing_files)) {
  stop("Missing final Top 10% edge list(s):\n", paste(missing_files, collapse = "\n"))
}

read_network <- function(path) {
  columns <- names(fread(path, nrows = 0L, showProgress = FALSE))
  score_columns <- intersect(c("Weight", "EdgeWeight", "XGB_Score", "Score"), columns)
  if (!all(c("TF", "Target") %in% columns) || !length(score_columns)) {
    stop("Expected TF, Target, and a score column in ", path)
  }

  score_column <- score_columns[[1L]]
  edges <- fread(path, select = c("TF", "Target", score_column), showProgress = FALSE)
  setnames(edges, score_column, "Score")
  edges[, `:=`(
    TF = toupper(trimws(as.character(TF))),
    Target = toupper(trimws(as.character(Target))),
    Score = suppressWarnings(as.numeric(Score))
  )]
  edges <- edges[
    !is.na(TF) & !is.na(Target) & nzchar(TF) & nzchar(Target) &
      !TF %chin% c("NA", "NAN") & !Target %chin% c("NA", "NAN") &
      TF != Target & is.finite(Score)
  ]
  setorder(edges, TF, Target, -Score)
  pairs <- unique(edges, by = c("TF", "Target"))[, .(TF, Target)]
  if (nrow(pairs) < 2L) stop("Fewer than two valid TF–target pairs in ", path)
  pairs
}

normalize_modules <- function(raw, graph_genes) {
  if (!is.list(raw)) stop("netModule() must return a list of gene vectors.")
  modules <- lapply(raw, function(values) {
    ids <- unique(toupper(trimws(as.character(unlist(values, use.names = FALSE)))))
    ids[ids %chin% graph_genes]
  })
  modules <- modules[lengths(modules) > 0L]
  if (!length(modules)) stop("netModule() returned no valid modules.")
  setNames(modules, sprintf("M%05d", seq_along(modules)))
}

module_memberships <- function(modules) {
  data.table(
    Module = rep(as.character(names(modules)), lengths(modules)),
    Gene = as.character(unlist(modules, use.names = FALSE))
  )
}

write_modules <- function(modules, method, stress, dest) {
  sizes <- lengths(modules)
  keep <- sizes >= min_module_size & sizes <= max_module_size
  retained <- modules[keep]

  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(Module = names(modules), Size = as.integer(sizes),
                    Method = method, Stress = stress, Retained_10_to_300 = keep),
         file.path(dest, "Module_Sizes_All.csv"))
  fwrite(module_memberships(modules), file.path(dest, "Module_Memberships_All.csv"))
  fwrite(module_memberships(retained),
         file.path(dest, "Module_Memberships_Retained_10_300.csv"))
  saveRDS(modules, file.path(dest, "All_Modules.rds"))
  saveRDS(retained, file.path(dest, "Retained_Modules_10_300.rds"))

  txt_dir <- file.path(dest, "Retained_Module_TXT")
  if (length(retained)) dir.create(txt_dir, showWarnings = FALSE)
  module_files <- sprintf("%s.txt", names(retained))
  for (i in seq_along(retained)) {
    writeLines(retained[[i]], file.path(txt_dir, module_files[[i]]))
  }
  existing <- list.files(txt_dir, pattern = "^M[0-9]+\\.txt$")
  unlink(file.path(txt_dir, setdiff(existing, module_files)))
  message("  Modules: ", length(modules), "; retained (10–300 genes): ", length(retained))
}

failures <- character()
for (stress in stresses) {
  edges <- read_network(network_files[[stress]])
  graph_genes <- unique(c(edges$TF, edges$Target))
  graph_pairs <- as.data.frame(edges)
  for (method in selected_methods) {
    message(stress, " / ", method)
    set.seed(seed_base + 100L * match(method, names(algorithms)) + match(stress, stresses))
    error <- tryCatch({
      raw <- miRspongeR::netModule(
        spongenetwork = graph_pairs,
        method = algorithms[[method]],
        directed = FALSE,
        modulesize = 1,
        save = FALSE
      )
      modules <- normalize_modules(raw, graph_genes)
      write_modules(modules, method, stress, file.path(output_root, method, stress))
      NULL
    }, error = conditionMessage)
    if (!is.null(error)) {
      failure <- paste0(stress, " / ", method, ": ", error)
      message("  ", failure)
      failures <- c(failures, failure)
    }
    raw <- modules <- NULL
  }
  edges <- graph_pairs <- graph_genes <- NULL
  invisible(gc())
}
if (length(failures)) stop("Clustering failed for:\n", paste(failures, collapse = "\n"))
message("Module files: ", normalizePath(output_root, winslash = "/"))
