# ==============================================================================
# V-RATIO DEMO
# ==============================================================================
# Demo script for estimating metacognitive efficiency (v-ratio)
# using the Flexible Confidence Boundary (FCB) models from
# Herregods, S., Le Denmat, P., Vermeylen, L., & Desender, K. (2025). 
# Modeling speed–accuracy trade-offs in the stopping rule for confidence judgments. 
# Psychological Review.
# ==============================================================================
# Author: Luc Vermeylen
# Date:   05/06/2026
#
# Usage:
# 1. Adjust "USER SETTINGS" section below.
# 2. Run the script:
#    - Single fit: Rscript Fitting_Pipeline.R SUBJECT_IDX or just Run/Source in RStudio
#    - HPC Batch: Use the provided batch_fit.slurm file to loop SUBJECT_IDX from 1:N
#
# for each new fit, make a new version of this file with different file name
# this avoids overwriting/duplicating previous results

# Note 1: SUBJECT_IDX refers to the index in sorted(unique(subject_id)).
# Set SUBJECT_IDX to NULL for a group-level fit (all data combined).

# Note 2: Your data needs to contain the following columns:
# - subject column, name can be changed below
# - rt (reaction time; in seconds)
# - acc (accuracy; 0=error, 1=correct)
# - rtconf (confidence RT; in seconds), 
# - cj (confidence judgement; 0-1 or 1-6, from low to high confidence)
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

# --- A. FOLDER & DATA INPUT ---
OUTPUT_FOLDER <- "test_vratio" # output folders will appear within /results/
DATA_NAME     <- "Exp2_cj6_Herregods2025.csv" # .csv or .rds supported
SUBJECT_COL   <- "sub_id"
SUBJECT_IDX   <- 1  # Change manually (e.g., 1), or NULL for Group Fit

# --- B. MODEL SELECTION ---
# "FCB_cj2" (binary confidence) or "FCB_cj6" (6-point confidence scale)
MODEL_NAME    <- "FCB_cj6" 
print(get_model_params(MODEL_NAME)) # print parameters and their bounds

# --- C. VARYING PARAMETERS ---
# Define which parameters vary by experimental condition.
# Example: list(v = ~ as.factor(Difficulty), a = ~ as.factor(Difficulty))
VARYING_PARAMS <- list() # Leave empty if you are not fitting conditions

# --- D. PARAMETERS TO FIX ---
# Fix parameters to a constant value (i.e., not estimated).
# Example: list(starting_point_confidence = 0.5)
FIXED_PARAMS <- list() 

# --- E. OPTIMIZER SETTINGS ---
ITER_MAX  <- 500   # number of optimization iterations (1000 recommended)
USE_CORES <- 1     # 0 = Single Core, 1 = Parallel

# ==============================================================================
# 2. DATA PREPARATION (STOP EDITING FROM HERE)
# ==============================================================================

# 0. Pipeline internal settings
COST_METHOD   <- "gsquare"
sim_constants <- list(ntrials = 5000, s = 1, dt = .001)
TRACE         <- ITER_MAX / 10 
POP_SCALE     <- 10
N_PRED        <- 10000

# Auto-detect experimental conditions to split G-Square likelihood blocks safely
cond_col_names <- if(length(VARYING_PARAMS) > 0) {
  unique(unlist(lapply(VARYING_PARAMS, function(x) if(inherits(x, "formula")) all.vars(x) else x)))
} else { NULL }

# define how the cost function should be built
FIT_TARGETS <- list(
  # Eq 4: Primary RT shape, split by Accuracy
  list(rt_col = "rt", split_cols = c("acc", cond_col_names), weight = 1),
  
  # Eq 5: Confidence RT shape, split by Accuracy
  list(rt_col = "rtconf", split_cols = c("acc", cond_col_names), weight = 1),
  
  # Eq 6: Pure confidence proportions
  list(rt_col = "cj", split_cols = c("acc", "cj", cond_col_names), weight = 1) 
)

# 1. Resolve Subject Index from CLI, if given
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (!is.na(args[1]) && args[1] != "NULL") SUBJECT_IDX <- as.numeric(args[1])
}
# Print configuration to confirm
print(paste("Active Configuration -> Subject:", ifelse(is.null(SUBJECT_IDX), "Group", SUBJECT_IDX)))

# 2. Loading and Filtering Data
DATA_FILE <- file.path("data", DATA_NAME) 
if (!file.exists(DATA_FILE)) stop(paste("File not found:", DATA_FILE))
# Handle RDS or CSV
if (grepl("\\.rds$", DATA_FILE, ignore.case = TRUE)) {
  raw_data <- readRDS(DATA_FILE)
} else {
  raw_data <- read.csv(DATA_FILE)
}
OUTPUT_DIR <- file.path("results", OUTPUT_FOLDER)

# 3. Handle subject filtering
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

# 4. Configure parameter space
config <- get_param_config(MODEL_NAME, VARYING_PARAMS, observations)
config <- fix_parameters(config, FIXED_PARAMS)
fit_constants <- c(config$fixed, sim_constants)

# 5. Standardization of the factors
observations$acc <- factor(observations$acc, levels = c(0, 1))
if("cj" %in% names(observations)) {
  cj_levels <- if(MODEL_NAME == "FCB_cj6") 1:6 else 0:1
  observations$cj <- factor(observations$cj, levels = cj_levels)
}
# Ensure condition columns are factors
if (!is.null(cond_col_names)) {
  for (col_name in cond_col_names) {
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

# --- VALIDATION ---
cat("\nValidating configuration...")
is_regression <- any(sapply(VARYING_PARAMS, function(x) inherits(x, "formula")))
all_design_vars <- unique(unlist(lapply(VARYING_PARAMS, all.vars)))

# Print mapping overview
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
cat("\n\nStarting optimization...\n")

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
  cat("Translating Regression Betas to Marginal Cell Means...\n")
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

# --- A. SCREEN OUTPUT ---
# 1. Primary RT Distribution
plot_dist(obs = fit_output$observations, pred = fit_output$predictions, 
          val_col = "rt", split_by_acc = TRUE, main_title = "Decision RT")

# 2. Confidence RT Distribution
plot_dist(obs = fit_output$observations, pred = fit_output$predictions, 
          val_col = "rtconf", split_by_acc = TRUE, main_title = "Confidence RT")

# 3. Confidence Rating (CJ) Distribution
p_cj <- plot_cj_distribution(fit_output$observations, fit_output$predictions)
print(p_cj)


# --- B. PDF OUTPUT ---
pdf_path <- file.path(OUTPUT_DIR, paste0(file_base, ".pdf"))
pdf(file = pdf_path, width = 10, height = 6)

plot_dist(obs = fit_output$observations, pred = fit_output$predictions, 
          val_col = "rt", split_by_acc = TRUE, main_title = "Decision RT")

plot_dist(obs = fit_output$observations, pred = fit_output$predictions, 
          val_col = "rtconf", split_by_acc = TRUE, main_title = "Confidence RT")

print(p_cj)

dev.off() 

cat(sprintf("\nPDF Report saved to: %s\n", pdf_path))