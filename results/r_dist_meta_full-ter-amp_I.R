# ==============================================================================
# DDM FITTING PIPELINE 
# ==============================================================================
# Author: Luc Vermeylen
# Date:   2026-01-12
# 
# Description:
# This script fits custom variants of the DDM to observed data.
# It supports multi-condition designs via parameter expansion. Parameters 
# defined in VARYING_PARAMS will be estimated separately for each condition 
# level, while others are shared across conditions.
#
# Usage:
# 1. Ensure 'helper_functions.R' and 'models.cpp' are in the working directory.
# 2. Adjust "USER SETTINGS" section below.
# 3. Run the script:
#    - Single fit: Rscript Fitting_Pipeline.R SUBJECT_IDX or just Run/Source in RStudio
#    - HPC Batch: Use the provided batch.slurm file to loop SUBJECT_IDX from 1:N
#
# for each new fit, make a new version of this file with different file name
# this avoids overwriting/duplicating previous results

# Note 1: SUBJECT_IDX refers to the index in sorted(unique(subject_id)).
# Set SUBJECT_IDX to NULL for a group-level fit (all data combined).

# Note 2: Column expectations:
# - rt (s for DDM, ms for DMC), acc (0=error, 1=correct)
# - congruency (-1=incongruent, 1=congruent) for conflict models
# - condition_id / meta_bin etc. for varying parameters
# - rtconf, cj (1:6, low to high confidence)
# - s = 1 for DDM, s = 4 for DMC (= sigma = diffusion noise)
# ==============================================================================

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================
Sys.setenv(LC_ALL = "C")
Sys.setenv(LANG = "C")
suppressPackageStartupMessages({
  library(Rcpp); library(DEoptim); library(dplyr); library(RcppZiggurat)
  source("helper_functions.R")
  sourceCpp("models.cpp")
})

# ==============================================================================
# 1. USER SETTINGS (EDIT ONLY THIS SECTION)
# ==============================================================================

# --- A. SET FOLDER ---
OUTPUT_DIR <- "r_dist_meta_full-ter-amp_I" # new folder for each fit to avoid overwriting. Will be created if it doesn't exist.

# --- B. DATA INPUT ---
DATA_FILE <- "observed_data_eren_metabinned2.csv"
SUBJECT_COL   <- "sub_id"

# Local Mode: Change this manually when running directly in RStudio (If running from CLI, this is overridden by the first argument)
SUBJECT_IDX <- NULL  # e.g., 1, or NULL for Group Fit

# --- C. MODEL SELECTION ---
MODEL_NAME <- "DMC" # C++ function name from models.cpp, and a key in 'model_params' (helper_functions.R) where parameters bounds are specified

# --- D. VARYING PARAMETERS ---
# Example: list(v = "condA", a = "condB", ter = c("condA", "condB"))
VARYING_PARAMS <- list(
  v_c = ~ as.factor(meta_bin),
  a = ~ distance * as.factor(meta_bin),
  #ter_mean = ~ as.factor(meta_bin),
  amp = ~ distance
)

# --- E. PARAMETERS TO FIX ---
# Constant values (overrides bounds). Use "ParamName" = value.
FIXED_PARAMS <- list(
  tau = 64, beta = 2.5, ter_sd = 39 
)

# --- F. COST FUNCTION TARGETS ---
# split_cols must include all variables used in VARYING_PARAMS.
FIT_TARGETS <- list(
  list(rt_col = "rt", split_cols = c("acc", "congruency", "distance", "meta_bin"), weight = 1)
  #list(rt_col = "rtconf", split_cols = c("cj","acc"), weight = 0.25,
  #list(rt_col = "rtconf", split_cols = "acc", weight = 0.25)
)
COST_METHOD = "gsquare"

# --- G. SIMULATION CONSTANTS ---
sim_constants <- list(
  ntrials = 5000, # trials used for simulating the predictions in each fitting iterations
  s       = 4, # diffusion constant, note: Typically 1/.1 for standard DDM, 4 for DMC!
  dt      = 1 # time step, .001 or 1 depending on model scale (seconds vs ms)
)

# --- H. OPTIMIZER SETTINGS ---
ITER_MAX  <- 500   # number of optimization iterations (min 500)
USE_CORES <- 1     # 0 = Single Core, 1 = Parallel
N_PRED    <- 10000 # Number of trials for final predictions (high-res)

# ==============================================================================
# 2. DATA PREPARATION (STOP EDITING FROM HERE)
# ==============================================================================

# 1. Resolve Subject Index from CLI
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (!is.na(args[1]) && args[1] != "NULL") SUBJECT_IDX <- as.numeric(args[1])
}
# Print configuration to confirm
print(paste("Active Configuration -> Subject:", ifelse(is.null(SUBJECT_IDX), "Group", SUBJECT_IDX)))

# 2. Loading and Filtering Data
if (!file.exists(DATA_FILE)) stop(paste("File not found:", DATA_FILE))
raw_data <- read.csv(DATA_FILE)

# Handle subject filtering
if (!is.null(SUBJECT_IDX)) {
  all_subjects <- sort(unique(raw_data[[SUBJECT_COL]]))
  target_sub_id <- all_subjects[SUBJECT_IDX]
  SUBJECT_ID_VAL <- target_sub_id
  observations <- raw_data[raw_data[[SUBJECT_COL]] == target_sub_id, ]
  print(paste("Filtering Subject Index:", SUBJECT_IDX, "-> ID:", target_sub_id))
} else {
  observations <- raw_data
  SUBJECT_ID_VAL <- "GROUP"
  print("Fitting Group/All Data.")
}
if (nrow(observations) == 0) stop("No data found for the selected subject...")

# 3. Configure parameter space
config <- get_param_config(MODEL_NAME, VARYING_PARAMS, observations)
config <- fix_parameters(config, FIXED_PARAMS)
fit_constants <- c(config$fixed, sim_constants)

# 5. STANDARDIZATION OF FACTORS
observations$acc <- factor(observations$acc, levels = c(0, 1))
if("congruency" %in% names(observations)) observations$congruency <- factor(observations$congruency, levels = c(-1, 1))

# Gebruik all.vars om variabelen uit formules te halen, of gebruik de string direct
cond_col_names <- if(length(VARYING_PARAMS) > 0) {
  unique(unlist(lapply(VARYING_PARAMS, function(x) if(inherits(x, "formula")) all.vars(x) else x)))
} else { NULL }

if (!is.null(cond_col_names)) {
  for (col_name in cond_col_names) {
    # Check of de variabele wel echt een kolom is in de data (voorkomt errors bij transformaties)
    if (col_name %in% names(raw_data)) {
      master_lvls <- sort(unique(raw_data[[col_name]]))
      observations[[col_name]] <- factor(observations[[col_name]], levels = master_lvls)
      cat(sprintf("Factorized column: %s | Intercept baseline: %s\n", col_name, master_lvls[1]))
    }
  }
}

print(paste("Observations ready:", nrow(observations), "trials."))

# ==============================================================================
# 3. RUN OPTIMIZATION
# ==============================================================================

TRACE <- ITER_MAX/10 # how often to print progress (every TRACE iterations)
POP_SCALE <- 10    # NP = POP_SCALE * n_parameters (10 is standard)

ctrl <- DEoptim.control(
  itermax = ITER_MAX, 
  trace   = TRACE,
  NP      = POP_SCALE * length(config$lower),
  parallelType = USE_CORES, 
  packages     = c("Rcpp", "RcppZiggurat", "dplyr"),
  parVar       = c("FIT_TARGETS", "VARYING_PARAMS", "MODEL_NAME") 
)

print(paste("Starting DEoptim on", length(config$lower), "parameters..."))
start_time <- Sys.time()

# --- PRE-FLIGHT VALIDATION ---
cat("\nValidating configuration...")
is_regression <- any(sapply(VARYING_PARAMS, function(x) inherits(x, "formula")))
all_design_vars <- unique(unlist(lapply(VARYING_PARAMS, all.vars)))

# 1. Check if all design variables are in FIT_TARGETS split_cols
for (target in FIT_TARGETS) {
  missing_in_target <- setdiff(all_design_vars, target$split_cols)
  if (length(missing_in_target) > 0) {
    stop(paste0("\n[STRICT ERROR] Variables '", paste(missing_in_target, collapse=", "), 
                "' are used in VARYING_PARAMS but missing in FIT_TARGETS split_cols.",
                "\nThis is required for valid likelihood blocks."))
  }
}

# 2. Print mapping overview
cat("\nParameter mapping overview:")
for (p in names(get_model_params(MODEL_NAME))) {
  if (p %in% names(VARYING_PARAMS)) {
    cat(sprintf("\n  [VARYING] %-10s -> Formula: %s", p, deparse(VARYING_PARAMS[[p]])))
  } else if (p %in% names(fit_constants)) {
    cat(sprintf("\n  [FIXED]   %-10s -> Value: %s", p, fit_constants[[p]]))
  } else {
    cat(sprintf("\n  [SHARED]  %-10s -> Estimated globally", p))
  }
}
cat("\n\nAll checks passed. Starting optimization...\n")

cat("\n==================================================\n")
cat("FREE PARAMETERS (Optimizing):\n")
cat(sprintf("  %-4s %-25s %-20s\n", "Idx", "Parameter Name", "Bounds [Min, Max]"))
cat("  --------------------------------------------------\n")
for (i in seq_along(config$lower)) {
  cat(sprintf("  [%02d] %-25s [%.4f, %.4f]\n", i, config$names[i], config$lower[i], config$upper[i]))
}
cat("==================================================\n")

result <- DEoptim(
  fn = objective_function, 
  lower = config$lower, 
  upper = config$upper,
  observations   = observations,
  param_names    = config$names,
  constants      = fit_constants,
  model_fun      = MODEL_NAME,
  targets        = FIT_TARGETS, 
  cost_method    = COST_METHOD,  
  varying_params = VARYING_PARAMS,
  control        = ctrl
)

duration <- Sys.time() - start_time
print(paste("Optimization finished in:", round(duration, 2), units(duration)))

# ==============================================================================
# 4. GENERATE HIGH RESOLUTION PREDICTIONS
# ==============================================================================

print(paste("Generating best-fit predictions with", N_PRED, "trials..."))

best_params <- result$optim$bestmem
best_betas  <- NULL

# convert the regressors to cell means if the model was specified with regression formulas in VARYING_PARAMS
is_regression <- any(sapply(VARYING_PARAMS, function(x) inherits(x, "formula")))

if (is_regression) {
  cat("Translating Regression Betas to Marginal Cell Means for downstream compatibility...\n")
  best_betas <- best_params 
  reconstructed_params <- c()
  
  # identify all design columns used across all parameters
  all_cond_cols <- unique(unlist(lapply(VARYING_PARAMS, all.vars)))
  cond_grid <- unique(observations[, all_cond_cols, drop = FALSE])
  p_base_names <- names(get_model_params(MODEL_NAME))
  
  for (i in 1:nrow(cond_grid)) {
    current_cell <- cond_grid[i, , drop = FALSE]
    
    # pluck the values for this specific cell (handles X * Beta)
    p_lvl <- pluck_params(best_betas, fit_constants, current_cell, VARYING_PARAMS, MODEL_NAME)
    
    for (p in p_base_names) {
      if (p %in% names(VARYING_PARAMS)) {
        # --- SMART NAMING: use only variables relevant to this specific parameter ---
        p_formula <- VARYING_PARAMS[[p]]
        p_vars <- all.vars(p_formula)
        
        # construct label based ONLY on the columns in this parameter's formula
        if (length(p_vars) > 0) {
          p_lvl_val <- paste(unlist(current_cell[, p_vars, drop = FALSE]), collapse = ".")
          reconstructed_params[paste0(p, ":", p_lvl_val)] <- p_lvl[[p]]
        } else {
          # intercept only (~ 1)
          reconstructed_params[p] <- p_lvl[[p]]
        }
      } else if (!p %in% names(fit_constants)) {
        # shared free parameter
        reconstructed_params[p] <- p_lvl[[p]]
      }
    }
  }
  # cleanup duplicates and overwrite best_params
  best_params <- reconstructed_params[!duplicated(names(reconstructed_params))]
}

# Override trial count for high-res predictions
pred_constants <- fit_constants
pred_constants$ntrials <- N_PRED

# we use the new varying_params format here if it was regression, 
# so the objective_function knows to expect classic cell-means (strings) now!
sim_varying_params <- list()
for (p in names(VARYING_PARAMS)) {
  # extract the specific variables for THIS parameter
  sim_varying_params[[p]] <- all.vars(VARYING_PARAMS[[p]])
}

final_predictions <- objective_function(
  params         = best_params,
  observations   = observations,
  param_names    = names(best_params),
  constants      = pred_constants,
  model_fun      = MODEL_NAME,
  targets        = FIT_TARGETS, 
  varying_params = sim_varying_params,
  returnFit      = 0 
)

# Align Factor Levels for Plotting
for (col in names(observations)) {
  if (col %in% names(final_predictions) && is.factor(observations[[col]])) {
    final_predictions[[col]] <- factor(final_predictions[[col]], levels = levels(observations[[col]]))
  }
}

# returnFit = 2 returns the actual proportions (Obs & Pred) used in G-Square
final_proportions <- objective_function(
  params         = best_params,
  observations   = observations,
  param_names    = names(best_params),
  constants      = pred_constants,
  model_fun      = MODEL_NAME,
  targets        = FIT_TARGETS, 
  varying_params = sim_varying_params, 
  returnFit      = 2
)

# Extract correct number of bins from attribute
n_bins_total <- attr(result$optim$bestval, "n_bins")
if (is.null(n_bins_total)) { # Fallback
  n_bins_total <- sum(unlist(lapply(final_proportions, function(target_list) {
    sum(sapply(target_list, function(level) length(level$obs)))
  })))
}
# ==============================================================================
# 5. REPORTING & SAVING
# ==============================================================================

# A. Calculate Metrics
n_free_params  <- length(config$lower) 
fit_metrics <- list(
  n_observations = nrow(observations),
  n_free_params  = n_free_params,
  n_bins         = n_bins_total,
  best_cost      = result$optim$bestval,
  bic = result$optim$bestval + n_free_params * log(nrow(observations)),
  aic = result$optim$bestval + 2 * n_free_params
)

# B. Bundle Output
fit_output <- list(
  best_params  = best_params,
  best_betas   = best_betas,  
  fit_metrics  = fit_metrics,
  optim_full   = result,
  param_info   = config,
  constants    = fit_constants,
  targets      = FIT_TARGETS,
  cost_method  = COST_METHOD,
  observations = observations,
  predictions  = final_predictions,
  final_proportions = final_proportions,
  info = list(
    model          = MODEL_NAME,
    timestamp      = Sys.time(),
    subject        = SUBJECT_ID_VAL,
    varying_params = VARYING_PARAMS
  )
)

# C. Save RDS
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
file_base <- paste0(MODEL_NAME, "_FIT_sub-", SUBJECT_ID_VAL, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
saveRDS(fit_output, file.path(OUTPUT_DIR, paste0(file_base, ".rds")))
print(paste("Best Cost:", round(fit_metrics$best_cost, 4)))

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================
cat("\nGenerating diagnostic reports...\n")

# determine columns for faceting based on the design
is_regression <- any(sapply(VARYING_PARAMS, function(x) inherits(x, "formula")))
cond_cols_plotting <- if(length(VARYING_PARAMS) > 0) {
  if(is_regression) unique(unlist(lapply(VARYING_PARAMS, all.vars))) else unique(unlist(VARYING_PARAMS))
} else { NULL }

# setup arguments for the standard plot_fit
plot_args <- list(
  obs            = fit_output$observations,
  pred           = fit_output$predictions,
  varying_params = VARYING_PARAMS,
  model_name     = MODEL_NAME,
  best_params    = fit_output$best_params,
  constants      = fit_output$constants,
  types          = c("dist", "delta", "caf", "mechanism") 
)

# --- A. SCREEN OUTPUT ---
# 1. standard plots
do.call(plot_fit, plot_args)

# 2. new bin-diagnostics
p_cdf <- plot_defective_cdf_mirror(fit_output$observations, fit_output$predictions, 
                                   cond_cols = cond_cols_plotting, 
                                   has_conflict = ("congruency" %in% names(fit_output$observations)))
p_bin_mass <- plot_bin_mass_mirror(fit_output$final_proportions)

print(p_cdf)
print(p_bin_mass)

# --- B. PDF OUTPUT ---
pdf_path <- file.path(OUTPUT_DIR, paste0(file_base, ".pdf"))
pdf(file = pdf_path, width = 10, height = 7)

# page 1-N: standard plots
do.call(plot_fit, plot_args)

# page N+1: defective cdf mirror
print(p_cdf)

# page N+2: bin mass mirror
print(p_bin_mass)

dev.off() 

cat(sprintf("\nPDF Report saved to: %s\n", pdf_path))