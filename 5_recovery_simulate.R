# ==============================================================================
# 5_recovery_simulate.R
# PARAMETER RECOVERY: Synthetic Data Simulation
# Author: Luc Vermeylen
# Description:
# Generates synthetic subject data with known ground-truth parameters.
# Supports 3 parameter generation scenarios:
#   1. "deoptim":   Master optimization search bounds.
#   2. "custom":    User-defined parameter intervals (with fallback to DEoptim).
#   3. "empirical": Empirical envelope from a previous fit report with outlier trimming.
# ==============================================================================

# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

# --- A. GENERAL & FILE OUTPUT ---
OUTPUT_FOLDER <- "recovery_fcb6_demo"          # Folder in /results where ground_truth.csv is saved
DATA_NAME     <- "synthetic_recovery_demo.csv" # CSV file created in /data (matches 1_run_pipeline.R)
MODEL_NAME    <- "FCB_cj6"                    # "FCB_cj2" (binary) or "FCB_cj6" (6-point scale)
N_SUBJECTS    <- 2                            # Number of synthetic subjects (2 for smoke test; 20-50 for paper)
RANDOM_SEED   <- 42                           # Set seed for reproducible synthetic subjects

# --- B. COGNITIVE MODEL SPECIFICATION ---
CONDITIONS     <- c("Difficulty")             # e.g., c("Difficulty") or empty c() if no conditions
VARYING_PARAMS <- list(v = ~ Difficulty)      # Formulas: list(v = ~ Difficulty) or empty list()
FIXED_PARAMS   <- list(a_slope = 0, sv = 0, st = 0, sz = 0) # Parameters held constant

# ------------------------------------------------------------------------------
# --- C. PARAMETER GENERATION SCENARIO (Choose 1 of 3) ---
# ------------------------------------------------------------------------------
# "deoptim"   -> SCENARIO 1: Sample across the full DEoptim search bounds (stress test).
# "custom"    -> SCENARIO 2: Sample across user-specified intervals below (targeted test).
#                            Any parameter omitted falls back to its default DEoptim bound.
# "empirical" -> SCENARIO 3: Sample inside the observed empirical envelope (Approach 3A)
#                            derived from a previous fit's master_parameters_report.csv.

PARAM_SOURCE <- "deoptim" # "deoptim", "custom", or "empirical"

# >>> SCENARIO 2: CUSTOM BOUNDS (Only used if PARAM_SOURCE == "custom") <<<
CUSTOM_BOUNDS <- list(
  a      = c(1.0, 2.5),
  v      = c(0.8, 2.5),
  vratio = c(0.5, 3.0),
  ter    = c(0.2, 0.8)
)

# >>> SCENARIO 3: EMPIRICAL FIT SETTINGS (Only used if PARAM_SOURCE == "empirical") <<<
EMPIRICAL_FOLDER        <- "vratio_exp2_fit" # Subfolder in /results containing master_parameters_report.csv
EMPIRICAL_QUANTILE_TRIM <- c(0.05, 0.95)     # Percentile trimming (5% to 95%) to discard extreme human tails
FILTER_COST_OUTLIERS    <- TRUE              # Exclude degenerate fits (cost > median + 3*IQR or cost >= 1e6)

# ------------------------------------------------------------------------------
# --- D. TRIAL DESIGN: CHOOSE OPTION A OR OPTION B ---
# ------------------------------------------------------------------------------
# OPTION A: Set to CSV/RDS in /data to clone an empirical subject's trial count/balance.
TEMPLATE_DATA    <- NULL 

# OPTION B: Used only if TEMPLATE_DATA is NULL. Fixed trials per condition combination.
CONDITION_LEVELS <- list(Difficulty = c("Easy", "Hard"))
TRIALS_PER_CELL  <- 100

# ==============================================================================
# 2. SETUP & ENGINE INITIALIZATION
# ==============================================================================
library(here)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppZiggurat)

source(here("R", "helper_functions.R"))
sourceCpp(here("R", "models.cpp"))

set.seed(RANDOM_SEED)

OUTPUT_DIR <- here("results", OUTPUT_FOLDER)
DATA_FILE  <- here("data", DATA_NAME)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Start with default master boundaries
master_bounds <- get_bounds(MODEL_NAME)
active_lower  <- master_bounds$lower
active_upper  <- master_bounds$upper

# ==============================================================================
# 3. CONFIGURE PARAMETER BOUNDS BY SCENARIO
# ==============================================================================
cat("\n=== CONFIGURING PARAMETER SAMPLING BOUNDS ===\n")

if (PARAM_SOURCE == "deoptim") {
  cat("[SCENARIO 1: DEoptim Search Bounds] Sampling across full optimization search space.\n")
  
} else if (PARAM_SOURCE == "custom") {
  cat("[SCENARIO 2: Custom Bounds with Fallback] Overriding user-specified parameters:\n")
  for (p in names(CUSTOM_BOUNDS)) {
    if (p %in% names(active_lower)) {
      active_lower[p] <- CUSTOM_BOUNDS[[p]][1]
      active_upper[p] <- CUSTOM_BOUNDS[[p]][2]
      cat(sprintf("  %-15s -> [%.3f, %.3f]\n", p, active_lower[p], active_upper[p]))
    } else {
      warning("Custom parameter '", p, "' not recognized in model '", MODEL_NAME, "'. Ignored.")
    }
  }
  unspecified <- setdiff(names(active_lower), names(CUSTOM_BOUNDS))
  cat(sprintf("  (Fallback to DEoptim bounds for: %s)\n", paste(unspecified, collapse = ", ")))
  
} else if (PARAM_SOURCE == "empirical") {
  cat(sprintf("[SCENARIO 3: Empirical Envelope] Extracting human ranges from folder: '%s'\n", EMPIRICAL_FOLDER))
  emp_path <- here("results", EMPIRICAL_FOLDER, "master_parameters_report.csv")
  if (!file.exists(emp_path)) stop("Empirical report not found: ", emp_path)
  
  emp_df <- fread(emp_path)
  
  # A. Filter Out Degenerate Fits / Cost Outliers
  if (FILTER_COST_OUTLIERS && "cost" %in% names(emp_df)) {
    initial_n <- nrow(emp_df)
    valid_costs <- emp_df$cost[is.finite(emp_df$cost) & emp_df$cost < 1e6]
    cost_threshold <- median(valid_costs, na.rm = TRUE) + 3 * IQR(valid_costs, na.rm = TRUE)
    emp_df <- emp_df[cost < cost_threshold & cost < 1e6]
    cat(sprintf("  Filtered degenerate fits by cost: retained %d of %d fits (threshold: %.1f)\n", 
                nrow(emp_df), initial_n, cost_threshold))
  }
  
  # B. Calculate Robust Envelope (Quantile Trimming)
  for (p in names(active_lower)) {
    # Match parameter column (e.g., 'a', 'ter', 'v', 'v.Easy', 'v:Easy')
    matched_cols <- grep(paste0("^", p, "($|[:._])"), names(emp_df), value = TRUE)
    
    if (length(matched_cols) > 0) {
      all_vals <- unlist(emp_df[, matched_cols, with = FALSE])
      all_vals <- all_vals[is.finite(all_vals)]
      
      if (length(all_vals) >= 3 && sd(all_vals) > 1e-6) {
        emp_low  <- as.numeric(quantile(all_vals, probs = EMPIRICAL_QUANTILE_TRIM[1], na.rm = TRUE))
        emp_high <- as.numeric(quantile(all_vals, probs = EMPIRICAL_QUANTILE_TRIM[2], na.rm = TRUE))
        
        # Clamp strictly within physical DEoptim bounds
        active_lower[p] <- max(master_bounds$lower[p], emp_low)
        active_upper[p] <- min(master_bounds$upper[p], emp_high)
        cat(sprintf("  %-15s -> [%.3f, %.3f] (from empirical %s)\n", 
                    p, active_lower[p], active_upper[p], paste(matched_cols, collapse = "/")))
      } else {
        cat(sprintf("  %-15s -> Fallback to DEoptim [%.3f, %.3f] (low variance/N in empirical data)\n", 
                    p, active_lower[p], active_upper[p]))
      }
    } else {
      cat(sprintf("  %-15s -> Fallback to DEoptim [%.3f, %.3f] (not found in empirical report)\n", 
                  p, active_lower[p], active_upper[p]))
    }
  }
} else {
  stop("Invalid PARAM_SOURCE! Choose 'deoptim', 'custom', or 'empirical'.")
}

# ==============================================================================
# 4. BUILD EXPERIMENTAL DESIGN PER SUBJECT
# ==============================================================================
cat("\n=== CONFIGURING EXPERIMENTAL DESIGN ===\n")

if (!is.null(TEMPLATE_DATA)) {
  cat(sprintf("[MODE: OPTION A] Mirroring trial design from template: %s\n", TEMPLATE_DATA))
  tmpl_path <- here("data", TEMPLATE_DATA)
  if (!file.exists(tmpl_path)) stop("Template data file not found in /data: ", TEMPLATE_DATA)
  
  tmpl_df <- if (grepl("\\.rds$", tmpl_path, ignore.case = TRUE)) readRDS(tmpl_path) else as.data.frame(fread(tmpl_path))
  
  missing_conds <- setdiff(CONDITIONS, names(tmpl_df))
  if (length(missing_conds) > 0) {
    stop("Conditions specified in CONDITIONS are missing from template: ", paste(missing_conds, collapse = ", "))
  }
  
  sub_col <- if ("sub_id" %in% names(tmpl_df)) "sub_id" else names(tmpl_df)[1]
  first_sub_id <- unique(tmpl_df[[sub_col]])[1]
  subject_design <- tmpl_df[tmpl_df[[sub_col]] == first_sub_id, CONDITIONS, drop = FALSE]
  rownames(subject_design) <- NULL
  
  for (cd in CONDITIONS) subject_design[[cd]] <- factor(subject_design[[cd]])
  cat(sprintf("Extracted %d empirical trials per subject.\n", nrow(subject_design)))
  
} else {
  cat("[MODE: OPTION B] Building balanced design from scratch with fixed trials per cell.\n")
  if (length(CONDITIONS) > 0) {
    missing_levels <- setdiff(CONDITIONS, names(CONDITION_LEVELS))
    if (length(missing_levels) > 0) {
      stop("Define factor levels in CONDITION_LEVELS for: ", paste(missing_levels, collapse = ", "))
    }
    cell_grid <- expand.grid(CONDITION_LEVELS[CONDITIONS])
    subject_design <- cell_grid[rep(seq_len(nrow(cell_grid)), each = TRIALS_PER_CELL), , drop = FALSE]
    rownames(subject_design) <- NULL
    
    for (cd in CONDITIONS) {
      subject_design[[cd]] <- factor(subject_design[[cd]], levels = CONDITION_LEVELS[[cd]])
    }
    cat(sprintf("Generated %d cells x %d trials = %d trials per subject.\n",
                nrow(cell_grid), TRIALS_PER_CELL, nrow(subject_design)))
  } else {
    subject_design <- data.frame(matrix(ncol = 0, nrow = TRIALS_PER_CELL))
    cat(sprintf("No conditions specified. Generating %d trials per subject.\n", TRIALS_PER_CELL))
  }
}

# ==============================================================================
# 5. SIMULATION LOOP (GROUND TRUTH & SYNTHETIC DATA)
# ==============================================================================
all_sim_trials <- list()
ground_truth_records <- list()

cat(sprintf("\nSimulating %d synthetic subjects for model '%s'...\n", N_SUBJECTS, MODEL_NAME))

base_param_names <- names(master_bounds$lower)

for (s_idx in 1:N_SUBJECTS) {
  
  # A. Setup Config with Active Bounds
  config <- get_param_config(MODEL_NAME, VARYING_PARAMS, subject_design)
  config <- fix_parameters(config, FIXED_PARAMS)
  
  # Inject active boundaries into config
  for (pn in names(config$lower)) {
    # Match root parameter name
    root_p <- sub("(:|\\.).*$", "", pn)
    if (root_p %in% names(active_lower)) {
      if (grepl("(:Int$|^[a-zA-Z0-9_]+$)", pn)) {
        config$lower[pn] <- active_lower[root_p]
        config$upper[pn] <- active_upper[root_p]
      } else {
        # Slope range based on active bounds
        p_rng <- active_upper[root_p] - active_lower[root_p]
        config$lower[pn] <- -p_rng
        config$upper[pn] <- p_rng
      }
    }
  }
  
  # B. Draw Valid Parameters via Rejection Sampling
  valid_draw <- FALSE
  attempts <- 0
  
  while (!valid_draw && attempts < 500) {
    attempts <- attempts + 1
    raw_params <- draw_random_params(config)
    
    if (length(CONDITIONS) > 0) {
      cond_grid <- unique(subject_design[, CONDITIONS, drop = FALSE])
      cell_checks <- logical(nrow(cond_grid))
      
      for (cg_idx in 1:nrow(cond_grid)) {
        cur_cell <- cond_grid[cg_idx, , drop = FALSE]
        p_cell <- pluck_params(raw_params, config$fixed, cur_cell, VARYING_PARAMS, MODEL_NAME)
        
        in_bounds <- TRUE
        for (pn in base_param_names) {
          val <- p_cell[[pn]]
          if (is.null(val) || val < active_lower[pn] || val > active_upper[pn]) {
            in_bounds <- FALSE; break
          }
        }
        cell_checks[cg_idx] <- in_bounds
      }
      if (all(cell_checks)) valid_draw <- TRUE
    } else {
      valid_draw <- TRUE
    }
  }
  
  if (!valid_draw) stop("Failed to draw valid parameters within active bounds. Check constraints.")
  
  # C. Record Ground Truth (Deduplicated)
  cell_means <- reconstruct_cell_means(raw_params, config$fixed, VARYING_PARAMS, subject_design, MODEL_NAME)
  combined_p <- c(cell_means, unlist(config$fixed))
  gt_row <- as.data.frame(as.list(combined_p[!duplicated(names(combined_p))]))
  gt_row$subject_id <- s_idx
  ground_truth_records[[s_idx]] <- gt_row
  
  # D. Simulate Behavioral Trials
  sub_trials <- list()
  
  if (length(CONDITIONS) > 0) {
    cond_grid <- unique(subject_design[, CONDITIONS, drop = FALSE])
    
    for (cg_idx in 1:nrow(cond_grid)) {
      cur_cell <- cond_grid[cg_idx, , drop = FALSE]
      
      cell_match <- rep(TRUE, nrow(subject_design))
      for (cd in CONDITIONS) {
        cell_match <- cell_match & (subject_design[[cd]] == cur_cell[[cd]][1])
      }
      n_cell_trials <- sum(cell_match)
      if (n_cell_trials == 0) next
      
      p_cell <- pluck_params(raw_params, config$fixed, cur_cell, VARYING_PARAMS, MODEL_NAME)
      
      sim_args <- p_cell
      sim_args$ntrials <- n_cell_trials
      sim_args$s       <- 1
      sim_args$dt      <- 0.001
      sim_args$seed    <- sample.int(2e9, 1)
      
      cell_sim <- do.call(MODEL_NAME, sim_args)
      for (cd in CONDITIONS) {
        cell_sim[[cd]] <- cur_cell[[cd]][1]
      }
      sub_trials[[cg_idx]] <- cell_sim
    }
    sub_df <- bind_rows(sub_trials)
  } else {
    sim_args <- as.list(c(raw_params, unlist(config$fixed)))
    sim_args$ntrials <- nrow(subject_design)
    sim_args$s       <- 1
    sim_args$dt      <- 0.001
    sim_args$seed    <- sample.int(2e9, 1)
    sub_df <- do.call(MODEL_NAME, sim_args)
  }
  
  sub_df$sub_id <- s_idx
  all_sim_trials[[s_idx]] <- sub_df
}

final_dataset   <- bind_rows(all_sim_trials) %>% select(sub_id, rt, acc, rtconf, cj, everything())
ground_truth_df <- bind_rows(ground_truth_records) %>% select(subject_id, everything())

# ==============================================================================
# 6. VERIFICATION & EXPORT
# ==============================================================================
fwrite(final_dataset, DATA_FILE)
gt_file <- file.path(OUTPUT_DIR, "ground_truth.csv")
fwrite(ground_truth_df, gt_file)

cat("\n====================================================================\n")
cat("VERIFICATION: TRIAL COUNTS PER CONDITION\n")
cat("====================================================================\n")
if (length(CONDITIONS) > 0) {
  print(table(final_dataset$sub_id, final_dataset[[CONDITIONS[1]]], dnn = c("Subject", CONDITIONS[1])))
} else {
  print(table(final_dataset$sub_id, dnn = "Subject"))
}
cat("====================================================================\n")

cat(sprintf("[SUCCESS] Dataset saved to:      %s (%d rows)\n", DATA_FILE, nrow(final_dataset)))
cat(sprintf("[SUCCESS] Ground truth saved to: %s (%d subjects)\n", gt_file, nrow(ground_truth_df)))
cat("\nNext step: Run '1_run_pipeline.R' using this synthetic dataset.\n")