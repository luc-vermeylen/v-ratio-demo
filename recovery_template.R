# ==============================================================================
# DDM PARAMETER RECOVERY PIPELINE
# ==============================================================================
# Author: Luc Vermeylen
# Description: 
# This script performs parameter recovery to test identifiability.
# 1. Samples random "Ground Truth" parameters from the model bounds.
# 2. Generates synthetic data based on a specified design.
# 3. Fits the data using DEoptim.
# 4. Compares "True" vs "Recovered" parameters in a diagnostic report.
# ==============================================================================

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================
Sys.setenv(LC_ALL = "C")
Sys.setenv(LANG = "C")
suppressPackageStartupMessages({
  library(Rcpp); library(DEoptim); library(dplyr); library(RcppZiggurat); library(data.table)
  source("helper_functions.R")
  sourceCpp("models.cpp")
})

# ==============================================================================
# 1. USER SETTINGS (EDIT ONLY THIS SECTION)
# ==============================================================================

# --- A. SET FOLDER ---
OUTPUT_FOLDER <- "recovery_test" # Files will appear in /results/recovery_test_v2/

# --- B. DESIGN SPECIFICATION ---
# Define the experimental factors and levels you want to simulate
RECOVERY_DESIGN <- list(
  meta_bin = c("0", "1")
)

# --- C. MODEL SELECTION ---
MODEL_NAME <- "DMC" 

# --- D. VARYING PARAMETERS ---
VARYING_PARAMS <- list(
  v_c      = ~ meta_bin,
  a        = ~ meta_bin
)

# --- E. PARAMETERS TO FIX ---
# Parameters held constant during both simulation and fitting.
FIXED_PARAMS <- list(
  tau = 64, beta = 2.5, ter_sd = 39 
)

# --- F. GROUND TRUTH SAMPLING ---
# How do we pick the 'True' parameters?
# If empty, parameters are sampled uniformly from the model bounds.
# You can provide specific values or restricted distributions here.
CUSTOM_TRUTH <- list(
  # "v_c:Int" = 0.5,
  # "v_c:distant" = 20
)

# --- G. COST FUNCTION TARGETS ---
FIT_TARGETS <- list(
  list(rt_col = "rt", split_cols = c("acc", "congruency", "meta_bin"), weight = 1)
)
COST_METHOD = "gsquare"

# --- H. SIMULATION CONSTANTS ---
sim_constants <- list(
  ntrials = 5000, # N trials per cell (e.g., 500 * 6 = 3000 trials total)
  s       = 4,
  dt      = 1
)

# --- I. OPTIMIZER SETTINGS ---
ITER_MAX  <- 500  
TRACE     <- ITER_MAX/10
POP_SCALE <- 10
USE_CORES <- 1     
N_PRED    <- 10000 

# ==============================================================================
# 2. GENERATE GROUND TRUTH & DATA (STOP EDITING FROM HERE)
# ==============================================================================

# 1. Create the design matrix
dummy_obs <- expand.grid(RECOVERY_DESIGN)

# 2. Build the parameter configuration
config <- get_param_config(MODEL_NAME, VARYING_PARAMS, dummy_obs)
config <- fix_parameters(config, FIXED_PARAMS)

# 3 & 4. SAMPLE GROUND TRUTH & GENERATE DATA
# We loop until we find a set of parameters that are physically valid
# as the objective function does not allow out of bound parameters
# which can occur in regression designs due to the additive nature of parameters
observations <- 1e9
attempts     <- 0
cat("\nSampling valid Ground Truth...")

while(!is.data.frame(observations)) {
  attempts <- attempts + 1
  
  # 1. draw_random_params returns ALL params (Free + Fixed)
  full_truth_draw <- draw_random_params(config)
  
  # 2. Apply any custom overrides if defined
  if (length(CUSTOM_TRUTH) > 0) {
    for (p in names(CUSTOM_TRUTH)) {
      if (p %in% names(full_truth_draw)) full_truth_draw[p] <- CUSTOM_TRUTH[[p]]
    }
  }
  
  # 3. Split the vector before passing to objective_function
  # We subset using config$names to get ONLY the free parameters.
  # This ensures the length matches the param_names argument.
  true_free_params <- full_truth_draw[config$names]
  
  # 4. Attempt to generate data
  # params gets ONLY the free ones. constants gets the fixed ones from config.
  observations <- objective_function(
    params         = true_free_params,  
    observations   = dummy_obs,
    param_names    = config$names,
    constants      = c(config$fixed, sim_constants),
    model_fun      = MODEL_NAME,
    targets        = FIT_TARGETS,
    varying_params = VARYING_PARAMS,
    returnFit      = 0 
  )
  
  if (attempts > 50) stop("Could not find a valid Ground Truth after 50 attempts.")
}

# SAVE the full set for later comparison (so we can compare both betas and cell means)
# We use the full draw from step 1
true_full_params <- full_truth_draw

cat(paste(" Success after", attempts, "attempt(s).\n"))
cat(true_free_params)

# 5. STANDARDIZATION OF FACTORS (Now safe to run)
observations$acc <- factor(observations$acc, levels = c(0, 1))
if("congruency" %in% names(observations)) observations$congruency <- factor(observations$congruency, levels = c(-1, 1))
for (col in names(RECOVERY_DESIGN)) {
  observations[[col]] <- factor(observations[[col]], levels = RECOVERY_DESIGN[[col]])
}

SUBJECT_ID_VAL <- paste0("REC_", paste0(sample(c(0:9, letters), 4), collapse = ""))
print(paste("Ground Truth ready. Synthetic N =", nrow(observations), "trials."))

# ==============================================================================
# 3. RUN OPTIMIZATION
# ==============================================================================
fit_constants <- c(config$fixed, sim_constants)

ctrl <- DEoptim.control(
  itermax = ITER_MAX, trace = TRACE,
  NP = POP_SCALE * length(config$lower),
  parallelType = USE_CORES, 
  packages = c("Rcpp", "RcppZiggurat", "dplyr", "data.table"),
  parVar = c("FIT_TARGETS", "VARYING_PARAMS", "MODEL_NAME") 
)

cat("\nStarting Parameter Recovery Fit...")
cat("\n==================================================")
cat("\nGROUND TRUTH TO RECOVER (BETAS):")
cat(sprintf("\n  %-25s %-10s %-20s", "Parameter", "Truth", "Bounds [Min, Max]"))
cat("\n  --------------------------------------------------")
for (i in seq_along(config$names)) {
  p_name <- config$names[i]
  cat(sprintf("\n  %-25s %-10.4f [%.2f, %.2f]", 
              p_name, true_free_params[p_name], config$lower[i], config$upper[i]))
}
cat("\n==================================================\n")

start_time <- Sys.time()

result <- DEoptim(
  fn             = objective_function, 
  lower          = config$lower, 
  upper          = config$upper,
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
# 4. PARAMETER RECONSTRUCTION & PREDICTIONS
# ==============================================================================
cat("\nGenerating best-fit predictions and reconstructing cell means...")

# 1. Identify the raw Beta vectors
best_betas <- result$optim$bestmem
true_betas <- true_free_params[config$names]

# 2. Translate Betas into readable Cell Means (e.g. amp:close.1)
# We use the new helper function for both Truth and Recovered sets
best_params_cell_means <- reconstruct_cell_means(best_betas, fit_constants, VARYING_PARAMS, observations, MODEL_NAME)
true_params_cell_means <- reconstruct_cell_means(true_betas, fit_constants, VARYING_PARAMS, observations, MODEL_NAME)

# 3. Generate high-resolution predictions 
# We pass the BETAS to the objective function because VARYING_PARAMS contains formulas
pred_constants <- fit_constants
pred_constants$ntrials <- N_PRED

final_predictions <- objective_function(
  params         = best_betas,
  observations   = observations,
  param_names    = config$names,
  constants      = pred_constants,
  model_fun      = MODEL_NAME,
  targets        = FIT_TARGETS, 
  varying_params = VARYING_PARAMS, 
  returnFit      = 0 
)

# Align factor levels for plotting consistency
for (col in names(observations)) {
  if (col %in% names(final_predictions) && is.factor(observations[[col]])) {
    final_predictions[[col]] <- factor(final_predictions[[col]], levels = levels(observations[[col]]))
  }
}

# ==============================================================================
# 5. RECOVERY REPORTING & SAVING
# ==============================================================================

# Comparison Table for the Coefficients (Betas)
comparison_betas <- data.frame(
  Parameter = config$names,
  True      = as.numeric(true_betas),
  Recovered = as.numeric(best_betas),
  stringsAsFactors = FALSE
) %>%
  mutate(Error = Recovered - True,
         Pct_Error = round(abs(Error / (config$upper - config$lower)) * 100, 2))

# Save the comprehensive bundle
fit_output <- list(
  # Optimization Results (Betas)
  best_betas        = best_betas,
  true_betas        = true_betas,
  
  # Readable Results (Cell Means)
  best_params       = best_params_cell_means,
  true_params       = true_params_cell_means,
  
  # Core fit data
  comparison        = comparison_betas,
  fit_metrics       = list(best_cost = result$optim$bestval, n_obs = nrow(observations)),
  param_info        = config,
  constants         = fit_constants,
  observations      = observations,
  predictions       = final_predictions,
  info = list(model = MODEL_NAME, timestamp = Sys.time(), mode = "RECOVERY", varying_params = VARYING_PARAMS)
)

OUTPUT_DIR <- file.path("results", OUTPUT_FOLDER) 
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
saveRDS(fit_output, file.path(OUTPUT_DIR, paste0(MODEL_NAME, "_RECOVERY_", SUBJECT_ID_VAL, ".rds")))

cat("\n==================================================\n")
cat("RECOVERY SUMMARY (Betas):\n")
print(comparison_betas)
cat("==================================================\n")

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================
# We pass the reconstructed Cell Means to plot_fit so it can find individual conditions
plot_fit(obs            = observations, 
         pred           = final_predictions, 
         varying_params = VARYING_PARAMS, 
         model_name     = MODEL_NAME, 
         best_params    = fit_output$best_params, 
         constants      = fit_constants)