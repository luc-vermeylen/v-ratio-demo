# helper_functions.R
# author: Luc Vermeylen

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

#### PARAMETER FUNCTIONS ####

# This list defines the master boundaries for all models in the pipeline.
# Each entry contains a lower and upper bound c(min, max).
# For parameters you wish to keep constant across all studies, set min == max.
model_params <- list(
  "FCB_cj2" = list(
    a                         = c(0.3, 4),
    v                         = c(0, 3),
    a_slope                   = c(0, 0),     # fixed to 0 (no decisional boundary collapse)
    ter                       = c(0, 1.5),
    a2                        = c(0, 20),
    vratio                    = c(0.2, 5),
    a2_slope_upper            = c(0, 20),
    a2_slope_lower            = c(0, 20),
    ter2                      = c(-1.5, 1.5),
    starting_point_confidence = c(0, 1)
  ),
  "FCB_cj6" = list(
    a                         = c(0.3, 4),
    v                         = c(0, 3),
    a_slope                   = c(0, 0),     # fixed to 0 (no decisional boundary collapse)
    ter                       = c(0, 1.5),
    a2                        = c(0, 20),
    vratio                    = c(0.2, 5),
    a2_slope_upper            = c(0, 20),
    a2_slope_lower            = c(0, 20),
    ter2                      = c(-1.5, 1.5),
    starting_point_confidence = c(0, 1)
  )
)

# Aliases: Allow different names to point to the same parameter set
# model_params[["DDM_metaconflict_bounds_fast"]] <- model_params[["DDM_metaconflict_bounds"]]

# Helper to retrieve the parameter list for a specific model name.
get_model_params <- function(model_name) {
  if (model_name %in% names(model_params)) {
    return(model_params[[model_name]])
  } else {
    stop(paste("Model", model_name, "not found in model_params list."))
  }
}

# Extracts lower and upper bounds separately into a list of vectors.
get_bounds <- function(model_name) {
  param_list <- get_model_params(model_name)
  lower <- sapply(param_list, `[`, 1)
  upper <- sapply(param_list, `[`, 2)
  return(list(lower = lower, upper = upper))
}

# Generates a single set of parameters by drawing from a uniform distribution
# between the lower and upper bounds defined in the provided config.
# This function is design-aware: if the config contains expanded parameters
# (e.g. amp_near, amp_far), it will draw independent values for each.
draw_random_params <- function(config, rounding = 4) {
  
  # 1. Draw random values for all parameters currently marked as FREE
  vals <- runif(length(config$lower), min = config$lower, max = config$upper)
  names(vals) <- names(config$lower)
  
  # 2. Combine with the parameters currently marked as FIXED (constants)
  all_params <- c(vals, unlist(config$fixed))
  
  return(round(all_params, rounding))
}

# This is the core function for setting up an optimization run.
# 1. It expands parameters that should vary by condition (e.g. amp becomes amp_level1, amp_level2).
# 2. It identifies which parameters are fixed (where lower == upper).
# 3. It returns a structured object containing only the "Free" parameters for DEoptim.
# INPUTS:
# - model_name (String): The key for the model_params list (e.g., "DMC").
# - varying_params (List): A named list defining the experimental design. 
#   Format: list(parameter_name = "column_name"). 
#   Example: list(amp = "condition_id") tells the function to create a 
#   unique 'amp' parameter for every level in the "condition_id" column.
# - observations (Dataframe): The dataset to be fit. Required if varying_params 
#   is not empty, as it is used to detect the factor levels of the design.
#
# OUTPUTS:
# A structured list used by the objective function and the optimizer:
# - lower (Named Vector): Lower bounds for all parameters DEoptim should search.
# - upper (Named Vector): Upper bounds for all parameters DEoptim should search.
# - names (Character Vector): The names of the free parameters (e.g., "amp_near").
# - fixed (Named List): Parameters held constant (where lower == upper).
# - varying (List): The design list passed in (stored for reference).
# - model_name (String): The name of the model (stored for reference).
get_param_config <- function(model_name, varying_params = list(), observations = NULL) {
  
  bounds <- get_bounds(model_name)
  model_param_names <- names(bounds$lower)
  
  # ensure all user-defined varying params actually exist in the model
  user_params <- names(varying_params)
  invalid_params <- setdiff(user_params, model_param_names)
  
  if (length(invalid_params) > 0) {
    stop(paste0("\n[CONFIG ERROR] invalid parameter names in VARYING_PARAMS: '", 
                paste(invalid_params, collapse = "', '"), "'",
                "\nlegal parameters for model '", model_name, "' are: ", 
                paste(model_param_names, collapse = ", ")))
  }
  new_lower <- c()
  new_upper <- c()
  
  for (p in names(bounds$lower)) {
    if (p %in% names(varying_params)) {
      mapping <- varying_params[[p]]
      
      # branch a: regression (formulas)
      if (inherits(mapping, "formula")) {
        if (is.null(observations)) stop("observations required for formula designs.")
        
        # create design matrix to extract beta names
        X <- model.matrix(mapping, observations)
        beta_names <- clean_beta_names(colnames(X), p, mapping)
        
        # calculate the maximum possible shift based on the base bounds.
        # this is a principled way to bound betas without hardcoding values.
        p_range <- bounds$upper[p] - bounds$lower[p]
        
        for (i in seq_along(beta_names)) {
          b_name <- beta_names[i]
          
          if (grepl(":Int$", b_name)) {
            # intercept gets the strict physical boundaries of the parameter
            new_lower[b_name] <- bounds$lower[p]
            new_upper[b_name] <- bounds$upper[p]
          } else {
            # slopes can shift the intercept at most from one end of the range to the other
            new_lower[b_name] <- -p_range
            new_upper[b_name] <- p_range
          }
        }
        
        # branch b: cell-mean (strings)
      } else {
        if (is.null(observations)) stop("observations required for cell-mean designs.")
        lvls <- levels(interaction(observations[mapping], drop = TRUE, sep = "."))
        
        for (l in lvls) {
          full_name <- paste0(p, ":", l)
          new_lower[full_name] <- bounds$lower[p]
          new_upper[full_name] <- bounds$upper[p]
        }
      }
      
    } else {
      # branch c: shared parameters
      new_lower[p] <- bounds$lower[p]
      new_upper[p] <- bounds$upper[p]
    }
  }
  
  # identify fixed parameters where lower equals upper
  is_fixed <- abs(new_upper - new_lower) < 1e-9
  
  return(list(
    lower   = new_lower[!is_fixed],
    upper   = new_upper[!is_fixed],
    names   = names(new_lower[!is_fixed]),
    fixed   = as.list(new_lower[is_fixed]),
    varying = varying_params
  ))
}

# Moves parameters from the "Free" list in the config to the "Fixed" list.
# This function is design-aware:
# - If you fix "amp", it will fix "amp" and all expanded versions (amp_near, amp_far).
# - If you fix "amp_near", it will only fix that specific level.
fix_parameters <- function(config, overrides) {
  if (length(overrides) == 0) return(config)
  
  for (p in names(overrides)) {
    val <- as.numeric(overrides[[p]])
    
    # Identify which parameters in the current config match this override.
    # Regex "^p(:|$)" ensures we match "v_c" exactly OR "v_c:level" 
    # It will NOT match "v_c" if the search term was just "v", preventing accidental partial matches.
    targets <- names(config$lower)[grepl(paste0("^", p, "(:|$)"), names(config$lower))]
    
    
    if (length(targets) > 0) {
      for (t in targets) {
        # Move from Free to Fixed
        config$fixed[[t]] <- val
        
        keep_mask <- names(config$lower) != t
        config$lower <- config$lower[keep_mask]
        config$upper <- config$upper[keep_mask]
        config$names <- config$names[keep_mask]
        message(paste("Fixed:", t, "=", val))
      }
    } else if (p %in% names(config$fixed)) {
      config$fixed[[p]] <- val
      message(paste("Updated Fixed:", p, "=", val))
    } else {
      warning(paste("Parameter", p, "not found. Cannot fix."))
    }
  }
  return(config)
}

# Helper to clean R beta names for consistent lookup
clean_beta_names <- function(raw_names, p_prefix, formula) {
  clean <- gsub("\\(Intercept\\)", "Int", raw_names)
  clean <- gsub("as.factor\\(", "", clean) 
  clean <- gsub("\\)", "", clean)
  
  vars <- all.vars(formula)
  for (v in vars) {
    stripped <- gsub(paste0("^", v), "", clean)
    clean <- ifelse(stripped == "" | stripped == ":", clean, stripped)
    clean <- gsub(paste0(":", v), ":", clean)
  }
  return(paste0(p_prefix, ":", clean))
}

# Extracts cell-specific parameters for a given experimental condition.
# Acts as a 'translator' between design-aware expanded parameters (e.g., amp:close.high) 
# and model-aware base parameters (e.g., amp). Supports multi-factor interactions.
#
# INPUTS:
# - best_params (Named Vector): The free parameters found by the optimizer.
# - constants (List): The parameters held constant during the fit (including Fixed model params).
# - current_cell (Dataframe/List): A single row representing the specific condition to pluck for.
# - varying_params (List): The design definition (e.g., list(amp = c("dist", "meta"))).
# - model_name (String): The name of the model to determine required base names.
#
# OUTPUTS:
# - A named list of parameters ready to be passed to a simulation function or plotting utility.
# hybrid parameter plucker
# calculates regression math or does string lookup depending on the formula
pluck_params <- function(best_params, constants, current_cell, varying_params, model_name) {
  
  p_base_names <- names(get_model_params(model_name))
  all_p <- as.list(c(best_params, unlist(constants)))
  out <- list()
  
  for (p in p_base_names) {
    if (p %in% names(varying_params)) {
      mapping <- varying_params[[p]]
      
      if (inherits(mapping, "formula")) {
        # --- REGRESSION ---
        X <- model.matrix(mapping, current_cell)
        clean_b_names <- clean_beta_names(colnames(X), p, mapping)
        
        # Check of we de betas wel echt hebben in de lijst
        if (all(clean_b_names %in% names(all_p))) {
          p_betas <- unlist(all_p[clean_b_names])
          out[[p]] <- sum(X[1, ] * p_betas)
        } else {
          # FALLBACK: Als betas missen, probeer de cel-naam direct (voor geconverteerde data)
          lvl_val <- paste(as.character(unlist(current_cell[, all.vars(mapping), drop=FALSE])), collapse=".")
          p_ext <- paste0(p, ":", lvl_val)
          out[[p]] <- if (!is.null(all_p[[p_ext]])) all_p[[p_ext]] else all_p[[p]]
        }
      } else {
        # --- CELL MEAN ---
        lvl_val <- paste(as.character(unlist(current_cell[, mapping, drop=FALSE])), collapse=".")
        p_ext <- paste0(p, ":", lvl_val)
        out[[p]] <- if (!is.null(all_p[[p_ext]])) all_p[[p_ext]] else all_p[[p]]
      }
    } else {
      # --- SHARED / FIXED ---
      out[[p]] <- if (!is.null(all_p[[p]])) all_p[[p]] else constants[[p]]
    }
    
    # Noodrem
    if (is.null(out[[p]])) {
      if (p == "alpha") { out[[p]] <- 2 } 
      else { stop(paste("Parameter lookup failed for:", p)) }
    }
  }
  return(out)
}

# Translates regression coefficients (Betas) into actual parameter values for every 
# unique experimental condition (Cell Means). 
# This is used to make the model results readable and compatible with plotting.
#
# INPUTS:
# - betas (Named Vector): The raw coefficients (Intercepts and Slopes).
# - constants (List): Fixed parameters (Fixed model params + Env constants).
# - varying_params (List): The design definition (e.g. list(amp = ~ distance)).
# - observations (Dataframe): The dataset containing the design columns.
# - model_name (String): The name of the model (e.g. "DMC").
#
# OUTPUTS:
# - A named vector of parameters in "Param:Level" format (e.g. "amp:close.1").
#   If a parameter was shared, it returns the base name (e.g. "tau").
reconstruct_cell_means <- function(betas, constants, varying_params, observations, model_name) {
  
  is_regression <- any(sapply(varying_params, function(x) inherits(x, "formula")))
  # If no formulas are used, the parameters are already cell-means
  if (!is_regression) return(betas) 
  
  reconstructed <- c()
  
  # Identify all design columns (e.g. distance, meta_bin)
  all_cond_cols <- unique(unlist(lapply(varying_params, all.vars)))
  # Get every unique combination of these factors found in the data
  cond_grid <- unique(observations[, all_cond_cols, drop = FALSE])
  p_base_names <- names(get_model_params(model_name))
  
  for (i in 1:nrow(cond_grid)) {
    current_cell <- cond_grid[i, , drop = FALSE]
    
    # Use pluck_params to calculate the specific value for this cell (X * Beta)
    p_lvl <- pluck_params(betas, constants, current_cell, varying_params, model_name)
    
    for (p in p_base_names) {
      if (p %in% names(varying_params)) {
        p_formula <- varying_params[[p]]
        p_vars <- all.vars(p_formula)
        
        # Create a label specific to the variables used for THIS parameter
        if (length(p_vars) > 0) {
          p_lvl_val <- paste(unlist(current_cell[, p_vars, drop = FALSE]), collapse = ".")
          reconstructed[paste0(p, ":", p_lvl_val)] <- p_lvl[[p]]
        } else {
          # Handle intercept-only varyings (~ 1)
          reconstructed[p] <- p_lvl[[p]]
        }
      } else if (p %in% names(betas)) {
        # Catch shared free parameters
        reconstructed[p] <- p_lvl[[p]]
      }
    }
  }
  # Remove duplicate shared parameters created during the grid loop
  return(reconstructed[!duplicated(names(reconstructed))])
}

#### FITTING FUNCTIONS ####

# Calculates dynamic cut-points (quantiles) to create RT bins based on sample size.
# High trial counts receive higher resolution bins to capture subtle shape changes.
#
# INPUTS:
# - observed_rt_vector (Numeric): A vector of raw reaction times from observed data.
# - min_trials (Integer): Minimum required trials to attempt binning.
# - med_threshold (Integer): N trials required for 4 bins.
# - high_threshold (Integer): N trials required for 6 bins.
# - xhigh_threshold (Integer): N trials required for 10 bins.
#
# OUTPUTS:
# - A numeric vector of quantile values, or NULL if data is insufficient/degenerate.
determine_quantiles <- function(observed_rt_vector, 
                                min_trials = 10,
                                med_threshold = 30,
                                high_threshold = 60,
                                xhigh_threshhold = 200) {
  
  # 1. Clean data
  if (is.null(observed_rt_vector)) return(NULL)
  
  # If the target is a factor (e.g., cj), we cannot make RT quantiles. 
  # Returning NULL will triggers pure proportion fitting in the pipeline
  if (!is.numeric(observed_rt_vector)) return(NULL) 
  
  observed_rt_vector <- observed_rt_vector[!is.na(observed_rt_vector)]
  n <- length(observed_rt_vector)
  
  # 2. Data quality check 
  # If N is too small or variance is near zero, we cannot estimate a distribution shape.
  if (n <= min_trials || (max(observed_rt_vector) - min(observed_rt_vector) < 1e-6)) {
    return(NULL)
  }
  
  # 3. Select resolution based on sample size
  if (n >= xhigh_threshhold) {
    probs <- seq(0.1, 0.9, by = 0.1) # X high resolution: 10 bins
  } else if (n > high_threshold) {
    probs <- c(.1,.3, .5, .7, .9)  # High resolution: 6 bins
  } else if (n > med_threshold) {
    probs <- c(.3, .5, .7)         # Medium resolution: 4 bins
  } else {
    probs <- c(.5)                 # Low resolution: 2 bins (Median split)
  }
  
  # 4. Compute quantiles
  qs <- quantile(observed_rt_vector, probs = probs, names = FALSE)
  
  # 5. Robustness: check for tied quantiles
  # If data is heavily discretized, q.1 might equal q.3. 
  # This creates 0-width bins which break the cost function
  # If detected, fallback to Median or NULL.
  if (any(duplicated(qs))) {
    # Try falling back to just the median
    qs <- quantile(observed_rt_vector, probs = c(.5), names = FALSE)
    # If even the median is duplicated (should be caught by step 2, but another check)
    if (length(qs) > 1 && any(duplicated(qs))) return(NULL) 
  }
  
  return(qs)
}

# Calculates the proportion of RTs falling into specific quantile bins.
# Optimized using fast C-level R functions (findInterval and tabulate).
#
# INPUTS:
# - rts (Numeric): Vector of reaction times (usually from simulation).
# - qs (Numeric): Vector of quantile cut-points defining the bins.
#
# OUTPUTS:
# - A numeric vector of proportions (length = length(qs) + 1). Sums to 1.0.
get_bin_proportions <- function(rts, qs) {
  
  # 1. Clean data
  if (is.null(rts)) return(numeric(length(qs) + 1))
  rts <- rts[!is.na(rts)]
  n_rts <- length(rts)
  
  # 2. Handle Empty/NULL Data
  if (n_rts == 0) return(numeric(length(qs) + 1))
  if (is.null(qs)) return(1.0) # If qs is NULL, we treat the entire range as 1 bin of proportion 1.0. 
  
  # 3. Binning
  # Boundaries are set from -Inf to Inf to ensure every RT is caught.
  if (is.unsorted(qs)) qs <- sort(qs)
  bounds <- c(-Inf, qs, Inf)
  
  # findInterval identifies which bin each RT belongs to
  bin_indices <- findInterval(rts, bounds, rightmost.closed = TRUE, all.inside = TRUE)
  
  # tabulate counts occurrences in each bin
  counts <- tabulate(bin_indices, nbins = length(qs) + 1)
  
  return(counts / n_rts) # Convert counts to proportions
}

# Get Weighted Proportions for a specific data split
# This function is the heart of the cost calculation. 
# It transforms raw RTs and discrete conditions into a probability distribution that can be
# compared between observed and predicted data, taking into account the requested 
# ways to split the data (e.g., accuracy or experimental cells).
# The result represents the joint probability mass across all cells
# (the whole vector sums to 1.0).
#
# INPUTS:
# - obs (Dataframe): The observed subject data.
# - pred (Dataframe): The simulated model data.
# - split_cols (Character Vector): Columns defining the cells (e.g. c("acc", "congruency")).
# - rt_col (String): The RT column to use (e.g. "rt").
#
# OUTPUTS:
# A list containing:
# - obs: Proportions of observed data across all cells.
# - pred: Proportions of predicted data across the same cells.
get_split_proportions <- function(obs, pred, split_cols, rt_col) {
  
  # 1. Validation
  
  # A. Check Observed Data
  if (!all(split_cols %in% names(obs))) {
    stop(paste("split_cols missing from Observed data:", 
               paste(setdiff(split_cols, names(obs)), collapse=", ")))
  }
  
  # B. Check Predicted Data
  if (!all(split_cols %in% names(pred))) {
    stop(paste("The simulation failed to produce columns needed for fitting:", 
               paste(setdiff(split_cols, names(pred)), collapse=", ")))
  }
  
  # C. Check RT column
  if (!rt_col %in% names(obs)) stop(paste("RT column", rt_col, "missing from Obs"))
  if (!rt_col %in% names(pred)) stop(paste("RT column", rt_col, "missing from Pred"))
  
  # 2. Create Grouping Factors
  # interaction() creates a factor "LevelA.LevelB". 
  # drop = FALSE ensures we respect all possible factor levels 
  obs_factor  <- interaction(obs[split_cols], drop = FALSE, sep = ".")
  pred_factor <- interaction(pred[split_cols], drop = FALSE, sep = ".")
  
  n_obs_total  <- nrow(obs)
  n_pred_total <- nrow(pred)
  
  # --- SAFETY CHECK ---
  # Ensure both dataframes operate in the exact same "Universe" of levels.
  # This catches upstream errors where the user forgot to set factor levels identically.
  if (!identical(levels(obs_factor), levels(pred_factor))) {
    stop(paste("Factor Mismatch! Observed and Predicted data have different factor definitions for:", 
               paste(split_cols, collapse=", "), 
               "\nHint: Ensure you used factor(..., levels=...) on both dataframes before calling this function."))
  }
  
  # We loop over ALL theoretical levels defined in the factors
  levels_to_check <- levels(obs_factor)
  
  all_obs_props  <- numeric()
  all_pred_props <- numeric()
  
  # 3. Loop through conditions
  for (lvl in levels_to_check) {
    
    # Vector Indexing (Fast)
    idx_obs  <- which(obs_factor == lvl)
    idx_pred <- which(pred_factor == lvl)
    
    # --- A. Marginal Probability ---
    # If idx_obs is empty (0 trials), w_obs becomes 0.
    w_obs  <- length(idx_obs) / n_obs_total
    w_pred <- length(idx_pred) / n_pred_total
    
    # Extract RTs
    rts_obs  <- obs[[rt_col]][idx_obs]
    rts_pred <- pred[[rt_col]][idx_pred]
    
    # --- B. Quantiles & Binning ---
    qs <- determine_quantiles(rts_obs)
    
    # Logic: If Observed data is missing or sparse (qs is NULL),
    # we cannot match shapes. We fallback to matching Marginal Probabilities only.
    if (is.null(qs)) {
      
      # Case 1: Real Sparse Data (e.g. 3 trials) -> w_obs is small positive.
      # Case 2: Missing Level (0 trials) -> w_obs is 0.
      
      # In both cases, we compare w_obs vs w_pred directly.
      # If Obs=0 and Pred=0.2, the cost function will see: (0 - 0.2)^2. Correct Penalty.
      
      all_obs_props  <- c(all_obs_props, w_obs)
      all_pred_props <- c(all_pred_props, w_pred)
      next
    }
    
    # Calculate bin proportions
    bin_props_obs  <- get_bin_proportions(rts_obs, qs)
    bin_props_pred <- get_bin_proportions(rts_pred, qs)
    
    names(bin_props_obs) <- rep(lvl, length(bin_props_obs))
    names(bin_props_pred) <- rep(lvl, length(bin_props_pred))
    
    # --- C. Apply Weighting ---
    weighted_obs  <- bin_props_obs * w_obs
    weighted_pred <- bin_props_pred * w_pred
    
    all_obs_props  <- c(all_obs_props, weighted_obs)
    all_pred_props <- c(all_pred_props, weighted_pred)
  }
  
  # 4. FINAL CHECKS
  
  # Check 1: Alignment
  if (length(all_obs_props) != length(all_pred_props)) {
    stop("Error: Observed and Predicted output vectors have different lengths.")
  }
  
  # Check 2: Observed Sum (Must be ~1.0)
  # Real data shouldn't have NAs in accuracy columns
  if (abs(sum(all_obs_props) - 1.0) > 1e-4) {
    warning(paste("Observed proportions sum to", round(sum(all_obs_props), 4), 
                  ". Check input data for missing factor levels."))
  }
  
  # Check 3: Predicted Sum
  # If trials time out (NA), they are excluded from the bins.
  # Sum will be < 1.0. This is expected behavior.
  sum_pred <- sum(all_pred_props)
  
  # Error if > 1 
  if (sum_pred > 1.001) {
    warning(paste("Predicted proportions > 1.0 (", sum_pred, "). Check simulation logic."))
  } 
  
  # Warn if predicted sum is too low
  # This implies too many trials are timing out
  if (sum_pred < 0.8) {
    # no stop() just warn. DEoptim should see a huge cost and move away naturally.
    warning(paste("High rate of undecided trials. Valid Prediction Mass =", round(sum_pred, 4)))
  }
  
  return(list(obs = all_obs_props, pred = all_pred_props))
}

# Calculate Model Fit Cost
#
# INPUTS:
# - observed (Numeric Vector): Proportions from data.
# - predicted (Numeric Vector): Proportions from model.
# - method (String): Optimization metric ("gsquare", "chisquare", "mse", etc.).
# - N (Integer): Scaling factor 
# - epsilon (Numeric): Small constant to prevent log(0).
#
# OUTPUTS:
# - A single numeric cost value. Lower is better.
calculate_cost <- function(observed, predicted, method = "chisquare", N = 1, epsilon = 1e-9) {
  
  # Ensure no negative predictions (Safety check)
  predicted <- pmax(predicted, 0) 
  
  # Select method and compute cost
  cost <- switch(method,
                 
                 # --- Pearson Chi-Square ---
                 # Formula: N * Sum( (Obs - Pred)^2 / Pred )
                 "chisquare" = {
                   # Epsilon prevents division by zero if Pred is 0
                   sum( (observed - predicted)^2 / (predicted + epsilon) ) * N
                 },
                 
                 # --- G-Square (Likelihood Ratio Chi-Square) ---
                 # Formula: 2 * N * Sum( Obs * log(Obs / Pred) )
                 "gsquare" = {
                   # 1. We only care about terms where Observed > 0. 
                   #    (If Obs=0, then 0 * log(...) is 0).
                   idx <- observed > 0
                   
                   # 2. Calculate log term only for valid indices
                   #    (predicted + epsilon) handles the P=0 case (infinite penalty)
                   sum_term <- sum(observed[idx] * log(observed[idx] / (predicted[idx] + epsilon)))
                   
                   2 * N * sum_term
                 },
                 
                 # --- Mean Squared Error (MSE) ---
                 "mse" = {
                   mean( (observed - predicted)^2 )
                 },
                 
                 # --- Root Mean Squared Error (RMSE) ---
                 "rmse" = {
                   sqrt(mean( (observed - predicted)^2 ))
                 },
                 
                 # --- Sum of Squared Errors (SSE) ---
                 "sse" = {
                   sum( (observed - predicted)^2 )
                 },
                 
                 stop(paste("Unknown cost method:", method))
  )
  
  return(cost)
}

# Generalized Objective Function
# 
# This is the objective function called by the optimizer (DEoptim). 
# It can run in parallel, manages the C++ simulation, aligns data structures,
# and calculates the total weighted cost across multiple to be specified
# features or "targets" of the data.
#
# params Numeric vector. Parameters generated by the optimizer (unnamed).
# observations Dataframe. The observed data. 
# Note: Must have columns converted to Factors to define the "Universe" of possible outcomes.
# param_names Character vector. Used to map the unnamed 'params' vector to named arguments for the model.
# constants List. Fixed simulation constants (ntrials, dt).
# model_fun Character string or Function. The name of the C++ model function. 
# targets List of lists. Defines the specific distributions to fit.
#        Each element represents a "Feature" of the data to minimize error on.
#        Example Structure:
#        targets = list(
#          # Target 1: RT (weighted 0.5)
#          list(rt_col = "rt", split_cols = c("acc", "congruency"), weight = 0.5),
#          
#          # Target 2: Confidence RT (weighted 0.5)
#          list(rt_col = "rtconf", split_cols = "cj", weight = 0.5)
#        )
# Each target can have an optional 'weight' element.
# cost_method String. The statistic to minimize ("gsquare", "chisquare", "mse", ...).
# returnFit Integer. 1 = return Cost (scalar, for optimization), 
# 0 = return Predictions (dataframe, for plotting), 
# 2 = return the proportions list (list, for plotting/debugging).
# 
# Scalar cost (for optimizer) or Dataframe of predictions (e.g., for plotting).
# generalized multi-factor objective function.
# orchestrates simulation batches, runs the linear predictors for regression,
# and calculates independent likelihood blocks to prevent parameter mimicry.
objective_function <- function(params, observations, param_names, constants, 
                               model_fun, targets, cost_method = "gsquare", 
                               varying_params = list(), returnFit = 1) {
  
  # load c++ environment if running on parallel workers
  if (!exists(".cpp_initialized", envir = .GlobalEnv)) {
    Sys.setenv(LC_ALL = "C") 
    library(Rcpp); library(RcppZiggurat); library(dplyr)
    source(file.path("R", "helper_functions.R")) 
    sourceCpp(file.path("R", "models.cpp"), rebuild = FALSE) 
  }
  
  # merge parameters into one searchable vector
  names(params) <- param_names
  all_p_vec <- c(params, unlist(constants))
  
  # figure out which columns drive the experimental design
  is_regression <- any(sapply(varying_params, function(x) inherits(x, "formula")))
  cond_cols <- if(length(varying_params) > 0) {
    if(is_regression) unique(unlist(lapply(varying_params, all.vars))) else unique(unlist(varying_params))
  } else {
    NULL
  }
  
  all_preds <- list() 
  
  # --- 1. SIMULATION BATCHES ---
  if (!is.null(cond_cols)) {
    
    # build the grid of unique experimental cells
    cond_grid <- unique(observations[, cond_cols, drop = FALSE])
    
    # pre-calculate design matrices to save time inside the loop
    X_list <- list()
    if (is_regression) {
      for (p in names(varying_params)) {
        if (inherits(varying_params[[p]], "formula")) {
          X_list[[p]] <- model.matrix(varying_params[[p]], cond_grid)
        }
      }
    }
    
    b_lims <- get_bounds(model_fun)
    base_names <- names(b_lims$lower)
    
    # simulate each cell
    for (i in seq_len(nrow(cond_grid))) {
      current_cell <- cond_grid[i, , drop = FALSE]
      current_sim_args <- list()
      
      for (p in base_names) {
        if (p %in% names(varying_params)) {
          mapping <- varying_params[[p]]
          
          # branch a: regression (linear predictor)
          if (inherits(mapping, "formula")) {
            clean_b_names <- clean_beta_names(colnames(X_list[[p]]), p, mapping)
            p_betas <- all_p_vec[clean_b_names]
            
            # calculate the actual parameter value (intercept + slopes)
            val <- sum(X_list[[p]][i, ] * p_betas)
            
            # strict boundary penalty: if the sum of betas results in an impossible 
            # physical value (e.g. negative drift), kill this iteration instantly.
            if (is.na(val) || val < b_lims$lower[p] || val > b_lims$upper[p]) {
              return(1e9)
            }
            current_sim_args[[p]] <- val
            
            # branch b: cell-mean (direct lookup)
          } else {
            lvls_in_cell <- sapply(current_cell[, mapping, drop = FALSE], as.character)
            lvl_val <- paste(lvls_in_cell, collapse = ".")
            p_ext <- paste0(p, ":", lvl_val)
            current_sim_args[[p]] <- all_p_vec[[p_ext]]
          }
        } else {
          # shared / fixed parameter
          current_sim_args[[p]] <- all_p_vec[[p]]
        }
      }
      
      # inject simulation constants
      current_sim_args$ntrials <- constants$ntrials
      current_sim_args$s       <- constants$s
      current_sim_args$dt      <- constants$dt
      current_sim_args$seed    <- sample.int(2e9, 1)
      
      batch_df <- do.call(model_fun, current_sim_args)
      
      # tag the simulated data with the condition labels
      for (col in cond_cols) { 
        batch_df[[col]] <- current_cell[[col]] 
      }
      all_preds[[i]] <- batch_df
    }
    preds <- do.call(rbind, all_preds)
    
  } else {
    # simple simulation without conditions
    sim_args <- c(as.list(params), constants)
    sim_args$seed <- sample.int(2e9, 1)
    preds <- do.call(model_fun, sim_args)
  }
  row.names(preds) <- NULL 
  
  if (returnFit == 0) return(preds)
  
  # --- 2. FACTOR ALIGNMENT ---
  # align factor levels to prevent splitting errors during likelihood computation
  cols_to_check <- unique(unlist(lapply(targets, `[[`, "split_cols")))
  for (col in cols_to_check) {
    if (col %in% names(preds)) {
      preds[[col]] <- factor(preds[[col]], levels = levels(observations[[col]]))
    }
  }
  
  # --- 3. COST CALCULATION (Independent Likelihood Blocks) ---
  total_cost <- 0
  total_bins_used <- 0
  all_proportions <- list() 
  
  for (target in targets) {
    target_id <- target$rt_col 
    
    # we use the actual design columns to construct independent blocks
    overlap_cols <- cond_cols 
    
    if (length(overlap_cols) > 0) {
      sub_target_cols <- setdiff(target$split_cols, overlap_cols)
      
      obs_split_factor   <- interaction(observations[overlap_cols], drop = TRUE)
      preds_split_factor <- interaction(preds[overlap_cols], drop = TRUE)
      
      lvls <- levels(obs_split_factor)
      condition_sum_cost <- 0
      
      for (lvl in lvls) {
        obs_sub   <- observations[obs_split_factor == lvl, , drop = FALSE]
        preds_sub <- preds[preds_split_factor == lvl, , drop = FALSE]
        
        if (nrow(obs_sub) == 0) next 
        
        props <- get_split_proportions(obs_sub, preds_sub, sub_target_cols, target$rt_col)
        if (returnFit == 2) all_proportions[[target_id]][[lvl]] <- props
        
        # carefully track exactly how many bins we use for this block
        total_bins_used <- total_bins_used + length(props$obs)
        
        sub_cost <- calculate_cost(props$obs, props$pred, cost_method, nrow(obs_sub))
        condition_sum_cost <- condition_sum_cost + sub_cost
      }
      current_target_cost <- condition_sum_cost
      
    } else {
      props <- get_split_proportions(observations, preds, target$split_cols, target$rt_col)
      if (returnFit == 2) all_proportions[[target_id]][["Joint"]] <- props
      
      total_bins_used <- total_bins_used + length(props$obs)
      current_target_cost <- calculate_cost(props$obs, props$pred, cost_method, nrow(observations))
    }
    
    w <- if(!is.null(target$weight)) target$weight else 1.0
    total_cost <- total_cost + (current_target_cost * w)
  }
  
  if (returnFit == 2) return(all_proportions)
  
  # return a heavy penalty if the final cost exploded
  if (!is.finite(total_cost)) return(1e9)
  
  # attach the exact number of bins used so the template can compute a valid bic
  attr(total_cost, "n_bins") <- total_bins_used
  
  return(total_cost)
}

#### PLOTTING FUNCTIONS ####

# Master plotting controller for DMC fit results.
# Automatically handles multi-factor designs by creating interaction-based 
# subsets and generating a sequence of diagnostic plots for each experimental cell.
#
# INPUTS:
# - obs (Dataframe): Observed trial data.
# - pred (Dataframe): Predicted trial data from the model.
# - varying_params (List): The design definition used during fitting.
# - model_name (String): Name of the fitted model (e.g., "DMC").
# - best_params (Named Vector): The optimal parameters found by the optimizer.
# - constants (List): Fixed parameters and simulation constants.
# - types (Character Vector): Types of plots to generate (dist, caf, delta, mechanism).
#
# OUTPUTS:
# - Sequence of plots sent to the current device (Screen or PDF).
# Master plotting controller - Hybrid Version (Supports Strings & Formulas)
plot_fit <- function(obs, pred, 
                     varying_params = list(), 
                     model_name     = "DMC", 
                     best_params    = NULL, 
                     constants      = NULL,
                     types          = c("dist", "caf", "delta", "mechanism")) {
  
  # 1. Identify unique conditions based on varying_params.
  is_regression <- any(sapply(varying_params, function(x) inherits(x, "formula")))
  cond_cols <- if(length(varying_params) > 0) {
    if(is_regression) unique(unlist(lapply(varying_params, all.vars))) else unique(unlist(varying_params))
  } else { NULL }
  
  obs_df <- as.data.frame(obs); pred_df <- as.data.frame(pred)
  
  if (!is.null(cond_cols)) {
    # Only use columns that are actually present in the data to avoid errors.
    valid_cols <- intersect(cond_cols, names(obs_df))
    obs_df$plot_grp  <- interaction(obs_df[, valid_cols, drop=FALSE], drop = TRUE, sep = ".")
    pred_df$plot_grp <- interaction(pred_df[, valid_cols, drop=FALSE], drop = TRUE, sep = ".")
    lvls <- levels(obs_df$plot_grp)
  } else {
    obs_df$plot_grp <- "Overall"; pred_df$plot_grp <- "Overall"; lvls <- "Overall"
  }
  
  has_conflict <- "congruency" %in% names(obs_df)
  
  # 2. Loop through each unique condition and generate plots
  for (l in lvls) {
    d_o <- obs_df[which(obs_df$plot_grp == l), ]
    d_p <- pred_df[which(pred_df$plot_grp == l), ]
    if (nrow(d_o) == 0) next
    
    if ("dist" %in% types) plot_dist(d_o, d_p, facet_col = if(has_conflict) 'congruency' else NULL, main_title = l)
    if ("caf" %in% types) plot_conditional_function(d_o, d_p, group_col = if(has_conflict) 'congruency' else NULL, main_title = l)
    if ("delta" %in% types && has_conflict) plot_delta(d_o, d_p, main_title = l)
    
    if ("mechanism" %in% types && !is.null(best_params) && !is.null(constants)) {
      if (model_name %in% c("DMC", "DDM_metaconflict_bounds")) {
        # To plot the mechanism, we need to extract the relevant parameters for this specific condition 'l'.
        grid_full <- unique(obs_df[, c(cond_cols, "plot_grp"), drop = FALSE])
        row_idx <- which(grid_full$plot_grp == l)[1]
        
        if (!is.na(row_idx)) {
          # pluck_params will handle the logic of matching the current cell's condition to the appropriate parameters, whether they are global or cell-specific.
          current_cell <- grid_full[row_idx, cond_cols, drop = FALSE]
          p_lvl <- pluck_params(best_params, constants, current_cell, varying_params, model_name)
          plot_dmc_mechanism(p_lvl, constants, main_title = paste("Mechanism:", l))
        }
      }
    }
  }
}

# Calculate Binned Means
# 
# Splits data into equal nbins, calculate Mean X and Mean Y.
# df Dataframe.
# x_col String. Name of X variable (e.g. "rt").
# y_col String. Name of Y variable (e.g. "acc", "cj").
# nbins Integer. Number of bins.
get_binned_means <- function(df, x_col, y_col, nbins) {
  
  # Safety checks
  if (nrow(df) < nbins * 2) return(NULL)
  
  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]
  
  # Ensure Y is numeric (handle factors like Accuracy/CJ)
  if (is.factor(y_vals)) y_vals <- as.numeric(as.character(y_vals))
  
  # Calculate quantiles on X
  probs  <- seq(0, 1, length.out = nbins + 1)
  breaks <- unique(quantile(x_vals, probs = probs, names = FALSE, na.rm = TRUE))
  
  if (length(breaks) < nbins + 1) return(NULL)
  
  # Binning
  bin_idx <- .bincode(x_vals, breaks, include.lowest = TRUE)
  
  # Aggregate
  mean_x <- tapply(x_vals, bin_idx, mean, na.rm = TRUE)
  mean_y <- tapply(y_vals, bin_idx, mean, na.rm = TRUE)
  
  return(data.frame(x = as.numeric(mean_x), y = as.numeric(mean_y)))
}

# Plot Conditional Function
# 
# Plots the mean of Y (e.g., Accuracy, Confidence) against bins of X (e.g., RT).
# Separates lines possible by Group (e.g., Congruency).
# Solid = Observed, Dashed = Predicted.
# 
# obs Dataframe of observed data.
# pred Dataframe of predicted data (Optional).
# x_col String. Variable for X-axis bins (default "rt").
# y_col String. Variable for Y-axis means (default "acc").
# group_col String. Variable to split lines by (eg, "congruency" default NULL).
# nbins Integer. Number of bins.
# main_title String. Plot title.
# color_palette Vector. Colors for the group column.
plot_conditional_function <- function(obs, pred = NULL, 
                                      x_col = "rt", 
                                      y_col = "acc", 
                                      group_col = NULL, 
                                      nbins = 5, 
                                      main_title = NULL,
                                      color_palette = c("black", "red3", "blue3", "orange2")) {
  
  # 1. Checks
  if (!x_col %in% names(obs)) stop(paste(x_col, "not found in Obs"))
  if (!y_col %in% names(obs)) stop(paste(y_col, "not found in Obs"))
  if (!is.null(pred)) {
    if (!x_col %in% names(pred)) stop(paste(x_col, "not found in Pred"))
    if (!y_col %in% names(pred)) stop(paste(y_col, "not found in Pred"))
  }
  
  # 2. Identify Groups
  if (is.null(group_col) || !group_col %in% names(obs)) {
    obs$dummy_grp <- "Overall"
    if(!is.null(pred)) pred$dummy_grp <- "Overall"
    group_col <- "dummy_grp"
  }
  
  groups <- sort(unique(obs[[group_col]]))
  if (length(groups) > length(color_palette)) warning("More groups than colors in palette.")
  
  # 3. Calculate Coordinates
  plot_data <- list()
  all_x <- c()
  all_y <- c()
  
  for (i in seq_along(groups)) {
    grp <- groups[i]
    col <- color_palette[((i-1) %% length(color_palette)) + 1] # Cycle colors
    
    # Obs
    dat_o <- obs[obs[[group_col]] == grp, ]
    stats_o <- get_binned_means(dat_o, x_col, y_col, nbins)
    
    # Pred
    stats_p <- NULL
    if (!is.null(pred)) {
      dat_p <- pred[pred[[group_col]] == grp, ]
      stats_p <- get_binned_means(dat_p, x_col, y_col, nbins)
    }
    
    # Store
    plot_data[[paste0("g", i)]] <- list(stats_o=stats_o, stats_p=stats_p, col=col, name=grp)
    
    # Collect limits
    if(!is.null(stats_o)) { all_x <- c(all_x, stats_o$x); all_y <- c(all_y, stats_o$y) }
    if(!is.null(stats_p)) { all_x <- c(all_x, stats_p$x); all_y <- c(all_y, stats_p$y) }
  }
  
  # 4. Setup Canvas
  if (length(all_x) == 0) {
    plot(1, type="n", axes=F, xlab="", ylab="", main="No Data")
    return()
  }
  
  # Smart Y-Limits:
  # If Y is Accuracy (0-1), usually clamp 0.5-1.0. 
  # If Y is Confidence (1-6), show full range or data range.
  ylim_val <- range(all_y, na.rm = TRUE)
  
  # Heuristic for Accuracy plots to look standard
  if (y_col %in% c("acc", "cor", "accuracy") && max(ylim_val) <= 1) {
    ylim_val <- c(min(0.5, min(ylim_val)), 1.0) 
  } else {
    # Add small padding for other variables
    ylim_val <- ylim_val + c(-0.05, 0.05) * diff(ylim_val)
  }
  
  if (is.null(main_title)) main_title <- paste("Conditional Function")
  
  plot(1, type="n", xlim=range(all_x, na.rm=T), ylim=ylim_val,
       xlab = paste("Mean", x_col), ylab = paste("Mean", y_col),
       main = main_title, las = 1)
  
  grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
  
  # 5. Draw Lines
  legend_txt <- c()
  legend_col <- c()
  legend_pch <- c()
  legend_lty <- c()
  
  for (item in plot_data) {
    
    # Draw Obs (Solid)
    if (!is.null(item$stats_o)) {
      lines(item$stats_o$x, item$stats_o$y, col=item$col, lwd=2, type="b", pch=19)
      
      # Add to legend
      legend_txt <- c(legend_txt, paste("Obs", item$name))
      legend_col <- c(legend_col, item$col)
      legend_pch <- c(legend_pch, 19)
      legend_lty <- c(legend_lty, 1)
    }
    
    # Draw Pred (Dashed)
    if (!is.null(item$stats_p)) {
      lines(item$stats_p$x, item$stats_p$y, col=item$col, lwd=2, type="b", pch=1, lty=2)
      
      legend_txt <- c(legend_txt, paste("Pred", item$name))
      legend_col <- c(legend_col, item$col)
      legend_pch <- c(legend_pch, 1)
      legend_lty <- c(legend_lty, 2)
    }
  }
  
  # 6. Legend
  if (length(legend_txt) > 0) {
    legend("bottomright", legend = legend_txt, col = legend_col, 
           pch = legend_pch, lty = legend_lty, lwd = 2, bty = "n", cex = 0.7)
  }
}

# Distribution Plot
# 
# Visualizes distributions of any variable (RT or ConfRT), optionally split by 
# conditions (facets) and/or accuracy (with errors displayed as negative RTs).
# 
# obs Dataframe of observed data.
# pred Dataframe of predicted data (Optional).
# val_col String. The dependent variable (x-axis). E.g. "rt", "rtconf".
# facet_col String. The variable to create subplots for. E.g. "congruency", "cj". If NULL, one plot.
# split_by_acc Boolean. If TRUE, splits distribution by Accuracy (Errors as negative RTs). If FALSE, plots overall density.
# acc_col String. Name of accuracy column (default "acc").
# main_title String.
# global_weighting Boolean. 
#  If TRUE, area sums to 1 across ALL panels (good for comparing relative frequency of panels).
#  If FALSE (default), area sums to 1 within EACH panel (good for comparing shape fit when obs/pred frequencies differ).
plot_dist <- function(obs, pred = NULL, 
                      val_col = "rt", 
                      facet_col = NULL, 
                      split_by_acc = TRUE,
                      acc_col = "acc",
                      main_title = NULL,
                      global_weighting = FALSE) {
  
  # 1. Checks
  if (!val_col %in% names(obs)) stop(paste("Column", val_col, "not found."))
  
  # 2. Panels
  if (is.null(facet_col)) {
    obs$dummy_facet <- "All Data"
    if(!is.null(pred)) pred$dummy_facet <- "All Data"
    facet_col <- "dummy_facet"
    panels <- "All Data"
  } else {
    # Add safety check if column is missing
    if (!facet_col %in% names(obs)) {
      obs$dummy_facet <- "All Data"; facet_col <- "dummy_facet"; panels <- "All Data"
    } else {
      panels <- sort(unique(obs[[facet_col]]))
    }
  }
  
  # 3. Global Limits (X-axis)
  pool <- obs[[val_col]]
  if (!is.null(pred)) pool <- c(pool, pred[[val_col]])
  pool <- pool[!is.na(pool)]
  if(length(pool) == 0) { warning("No valid data."); return() }
  
  max_val <- quantile(pool, 0.995, na.rm = TRUE)
  if (split_by_acc) xlim_val <- c(-max_val, max_val) else xlim_val <- c(0, max_val)
  
  # 4. Helper: Density Calculation
  get_dens <- function(vals, base_n, flip=FALSE) {
    vals <- vals[!is.na(vals)]
    if (length(vals) < 2) return(list(x=numeric(0), y=numeric(0)))
    d <- density(vals, from = 0, to = max_val, n = 512)
    d$y <- d$y * (length(vals) / base_n)
    if (flip) d$x <- -d$x
    return(d)
  }
  
  # 5. PRE-CALCULATION LOOP (Crucial for Y-Axis Scaling)
  plot_data <- list()
  global_max_y <- 0
  
  n_obs_total  <- nrow(obs)
  n_pred_total <- if(!is.null(pred)) nrow(pred) else 0
  
  for (lvl in panels) {
    # Subset
    d_o <- obs[obs[[facet_col]] == lvl, ]
    d_p <- if (!is.null(pred)) pred[pred[[facet_col]] == lvl, ] else NULL
    
    # Determine Denominator
    if (global_weighting) {
      n_base_o <- n_obs_total
      n_base_p <- n_pred_total
    } else {
      n_base_o <- nrow(d_o)
      n_base_p <- if(!is.null(d_p)) nrow(d_p) else 0
    }
    
    curves <- list()
    
    if (split_by_acc) {
      curves$o_corr <- get_dens(d_o[[val_col]][d_o[[acc_col]]==1], n_base_o, flip=F)
      curves$o_err  <- get_dens(d_o[[val_col]][d_o[[acc_col]]==0], n_base_o, flip=T)
      
      if (!is.null(d_p)) {
        curves$p_corr <- get_dens(d_p[[val_col]][d_p[[acc_col]]==1], n_base_p, flip=F)
        curves$p_err  <- get_dens(d_p[[val_col]][d_p[[acc_col]]==0], n_base_p, flip=T)
      }
    } else {
      curves$o_all <- get_dens(d_o[[val_col]], n_base_o, flip=F)
      if (!is.null(d_p)) {
        curves$p_all <- get_dens(d_p[[val_col]], n_base_p, flip=F)
      }
    }
    
    plot_data[[as.character(lvl)]] <- curves
    
    # Global maximum
    all_ys <- unlist(lapply(curves, `[[`, "y"))
    if (length(all_ys) > 0) {
      curr_max <- max(all_ys, na.rm=TRUE)
      if (curr_max > global_max_y) global_max_y <- curr_max
    }
  }
  
  if (global_max_y == 0) global_max_y <- 1 # Safety
  
  
  # 6. PLOTTING LOOP
  
  # layout
  n_panels <- length(panels)
  if (n_panels == 1) par(mfrow=c(1,1))
  else if (n_panels == 2) par(mfrow=c(1,2))
  else if (n_panels <= 4) par(mfrow=c(2,2))
  else par(mfrow=c(2,3)) 
  
  par(oma = c(0, 0, 3, 0), mar = c(5, 4, 4, 1))
  
  for (lvl in panels) {
    curves <- plot_data[[as.character(lvl)]]
    
    # --- DETERMINE Y-LIMIT ---
    if (global_weighting) {
      # Use the global maximum we found earlier. 
      # This ensures rare conditions look small compared to frequent ones.
      ylim_val <- c(0, global_max_y * 1.1)
    } else {
      # Use local max to zoom in to fit the shape, regardless of frequency.
      all_ys <- unlist(lapply(curves, `[[`, "y"))
      local_max <- if(length(all_ys)>0) max(all_ys, na.rm=TRUE) else 0
      ylim_val <- c(0, local_max * 1.1)
    }
    
    # Plot canvas
    title_txt <- if(facet_col == "dummy_facet") "All Data" else paste(facet_col, "=", lvl)
    xlabel <- if(split_by_acc) paste(val_col, "(Error <--> Correct)") else val_col
    
    plot(1, type="n", xlim=xlim_val, ylim=ylim_val, 
         main=title_txt, xlab=xlabel, ylab="Density", las=1)
    
    if (split_by_acc) abline(v=0, lwd=1, lty=2, col="gray")
    
    draw_poly <- function(d, col) {
      if(length(d$x) > 1) {
        polygon(c(d$x[1], d$x, d$x[length(d$x)]), c(0, d$y, 0), 
                col=adjustcolor(col, 0.4), border=NA)
      }
    }
    
    # Draw observations
    if (split_by_acc) {
      draw_poly(curves$o_corr, "gray")
      draw_poly(curves$o_err,  "gray")
    } else {
      draw_poly(curves$o_all, "lightblue")
    }
    
    # Draw predictions
    if (!is.null(d_p)) {
      if (split_by_acc) {
        if(length(curves$p_corr$x)>1) lines(curves$p_corr, col="darkgreen", lwd=2.5)
        if(length(curves$p_err$x)>1)  lines(curves$p_err,  col="red",       lwd=2.5)
      } else {
        if(length(curves$p_all$x)>1)  lines(curves$p_all,  col="blue",      lwd=2.5, lty=1)
      }
    }
  }
  
  if(is.null(main_title)) main_title <- ""
  mtext(main_title, outer = TRUE, cex = 1.2, line = 1)
  par(mfrow=c(1,1))
}

# Calculate Delta Plot Coordinates
# 
# Splits Correct trials into nbins, calculates Mean RT
# for each bin per condition, and computes the delta (Incong - Cong).
#
# df Dataframe with columns: rt, acc, congruency
# nbins Integer. Number of bins (e.g. 5 for quintiles).
# Returns: Dataframe with columns x (Mean RT) and y (Delta). Or NULL if invalid.
get_delta_coords <- function(df, nbins) {
  
  # 1. Filter Correct Trials
  # Assumes acc=1 is correct; congruency 1=Cong, -1=Incong
  rt_c <- df$rt[df$acc == 1 & df$congruency == 1]
  rt_i <- df$rt[df$acc == 1 & df$congruency == -1]
  
  # Safety Check: Need enough trials (2 per bin minimum)
  if (length(rt_c) < nbins * 2 || length(rt_i) < nbins * 2) return(NULL)
  
  # 2. Internal Helper: Bin and Average
  calc_bin_means <- function(x) {
    # Calculate quantile breaks (0%, 20%, 40%...)
    breaks <- quantile(x, probs = seq(0, 1, length.out = nbins + 1), names = FALSE)
    # Handle discrete data (tied quantiles)
    if (length(unique(breaks)) < length(breaks)) return(NULL)
    
    # Assign bins
    indices <- .bincode(x, breaks, include.lowest = TRUE)
    # Calculate means
    tapply(x, indices, mean)
  }
  
  # 3. Compute Means
  means_c <- calc_bin_means(rt_c)
  means_i <- calc_bin_means(rt_i)
  
  if (is.null(means_c) || is.null(means_i)) return(NULL)
  
  # 4. Compute Coordinates
  # X = Average RT of the two conditions
  # Y = Conflict Effect (Incongruent - Congruent)
  data.frame(
    x = (means_c + means_i) / 2,
    y = means_i - means_c
  )
}

# Plot Delta Function
# 
# Plots the conflict effect (Incong - Cong) as a function of time.
# Solid = Observed, Dashed = Predicted.
# 
# obs Dataframe. Observed data.
# pred Dataframe. Predicted data (optional).
# nbins Integer. Number of bins.
# main_title String.
plot_delta <- function(obs, 
                       pred = NULL, 
                       nbins = 5, 
                       main_title = "Delta Plot") {
  
  # 1. STANDARDIZATION
  # Helper to clean column names and types safely
  clean_df <- function(df) {
    if (is.null(df)) return(NULL)
    
    # Force Factors to Numeric (-1/1 and 0/1)
    # We use as.character first to avoid factor index issues (e.g. factor level 1 becoming '1')
    if (is.factor(df$acc))        df$acc        <- as.numeric(as.character(df$acc))
    if (is.factor(df$congruency)) df$congruency <- as.numeric(as.character(df$congruency))
    
    return(df)
  }
  
  obs  <- clean_df(obs)
  pred <- clean_df(pred)
  
  # 2. CALCULATE COORDINATES
  coords_o <- get_delta_coords(obs, nbins)
  coords_p <- if(!is.null(pred)) get_delta_coords(pred, nbins) else NULL
  
  if (is.null(coords_o)) {
    plot(1, type="n", axes=F, xlab="", ylab="", main="Insufficient Data")
    return()
  }
  
  # 3. SETUP CANVAS
  # Combine data to find full range
  all_x <- c(coords_o$x, coords_p$x)
  all_y <- c(coords_o$y, coords_p$y, 0) # Ensure 0 is included
  
  # Add 5% padding
  x_lim <- range(all_x, na.rm = TRUE)
  y_lim <- range(all_y, na.rm = TRUE) + c(-0.02, 0.02)
  
  plot(1, type = "n", xlim = x_lim, ylim = y_lim,
       xlab = "Mean Reaction Time (s)", 
       ylab = expression(Delta ~ "RT (Incongruent - Congruent)"),
       main = main_title, las = 1)
  
  grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
  abline(h = 0, lty = 2, col = "gray50")
  
  # 4. DRAW LINES
  # Observed (Solid, Filled)
  lines(coords_o$x, coords_o$y, type = "b", pch = 19, lwd = 2, col = "black")
  
  # Predicted (Dashed, Open)
  if (!is.null(coords_p)) {
    lines(coords_p$x, coords_p$y, type = "b", pch = 1, lty = 2, lwd = 2, col = "red")
  }
  
  # 5. LEGEND
  leg_txt <- "Observed"
  leg_col <- "black"
  leg_pch <- 19
  leg_lty <- 1
  
  if (!is.null(coords_p)) {
    leg_txt <- c(leg_txt, "Predicted")
    leg_col <- c(leg_col, "red")
    leg_pch <- c(leg_pch, 1)
    leg_lty <- c(leg_lty, 2)
  }
  
  legend("topleft", legend = leg_txt, col = leg_col, 
         pch = leg_pch, lty = leg_lty, lwd = 2, bty = "n", cex = 0.8)
}


plot_dmc_mechanism <- function(params, constants = NULL, 
                               t_max = 1000, alpha = 2, main_title = NULL) {
  
  # 1. Merge parameters
  p <- c(as.list(params), as.list(constants))
  
  # Check for required parameters
  req <- c("v_c", "amp", "tau", "a")
  if (!all(req %in% names(p))) {
    stop(paste("Missing params:", paste(setdiff(req, names(p)), collapse=", ")))
  }
  
  if (!"alpha" %in% names(p)) p$alpha <- alpha
  
  # 2. Time Vector
  t <- seq(0, t_max, length.out = 500)
  
  # 3. Calculate Gamma Pulse
  t_peak <- p$tau * (p$alpha - 1)
  if (t_peak <= 1e-6) {
    pulse <- p$amp * exp(-t / p$tau)
  } else {
    term_power <- (t / t_peak)^(p$alpha - 1)
    term_exp   <- exp((t_peak - t) / p$tau)
    pulse      <- p$amp * term_power * term_exp
  }
  
  # 4. Calculate Trajectories
  start <- if ("z" %in% names(p)) (p$z * 2 * p$a) - p$a else 0 
  drift <- p$v_c * t
  ev_cong   <- start + drift + pulse
  ev_incong <- start + drift - pulse
  
  # 5. AUTO-SCALE
  hit_incong <- which(abs(ev_incong) >= p$a)[1]
  min_hit <- min(c(hit_incong, length(t)), na.rm = TRUE)
  cutoff_idx <- min(length(t), floor(min_hit * 1.5))
  min_show <- max(200, t_peak * 1.5)
  if (t[cutoff_idx] < min_show) {
    target_idx <- which(t >= min_show)[1]
    cutoff_idx <- if(is.na(target_idx)) length(t) else target_idx
  }
  
  plot_t <- t[1:cutoff_idx]; plot_pulse <- pulse[1:cutoff_idx]
  plot_cong  <- ev_cong[1:cutoff_idx]; plot_incong <- ev_incong[1:cutoff_idx]
  
  # --- PLOTTING ---
  # We use oma (outer margin area) to leave space for the main_title at the top
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), mar = c(5, 4, 2, 1))
  
  # A. The Conflict Pulse
  plot(plot_t, plot_pulse, type = "l", lwd = 3, col = "purple",
       xlab = "Time (ms)", ylab = "Evidence",
       main = "Conflict Pulse", las = 1,
       ylim = c(min(0, min(plot_pulse)), max(plot_pulse) * 1.1))
  grid()
  
  # B. Trajectories
  y_min <- min(-p$a, min(c(plot_cong, plot_incong))) * 1.1
  y_max <- max(p$a,  max(c(plot_cong, plot_incong))) * 1.1
  plot(plot_t, plot_cong, type = "l", lwd = 2, col = "darkgreen",
       ylim = c(y_min, y_max), xlab = "Time (ms)", ylab = "Evidence",
       main = "Expected Trajectories", las = 1)
  grid()
  lines(plot_t, plot_incong, lwd = 2, col = "red")
  abline(h = p$a, lwd = 2); abline(h = -p$a, lwd = 2)
  abline(h = start, lty = 3, col = "blue")
  
  # --- ADD THE MAIN TITLE ---
  if (!is.null(main_title)) {
    mtext(main_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
  }
  
  par(mfrow = c(1, 1))
}

# defective cdf plot (quantile probability plot)
# plots observed quantiles (circles) vs predicted quantiles (crosses + lines)
# mirrored defective cdf plot (quantile probability plot)
# x-axis: errors (negative rt) <--- 0 ---> corrects (positive rt)
# y-axis: cumulative probability x accuracy
plot_defective_cdf_mirror <- function(obs, pred, cond_cols = NULL, has_conflict = TRUE) {
  
  get_cdf_data <- function(df, source_name) {
    if(nrow(df) == 0) return(NULL)
    
    actual_has_conflict <- has_conflict && ("congruency" %in% names(df))
    base_groups <- if(is.null(cond_cols)) character(0) else cond_cols
    if(actual_has_conflict) base_groups <- c(base_groups, "congruency")
    
    cell_totals <- df %>% group_by(across(all_of(base_groups))) %>% summarise(n_cell = n(), .groups = "drop")
    probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    
    cdf_df <- df %>%
      group_by(across(all_of(c(base_groups, "acc")))) %>%
      summarise(
        n_acc = n(),
        q_rt = list(quantile(rt, probs = probs, names = FALSE, na.rm = TRUE)),
        q_prob = list(probs),
        .groups = "drop"
      ) %>%
      left_join(cell_totals, by = base_groups) %>%
      mutate(p_acc = n_acc / n_cell) %>%
      unnest(cols = c(q_rt, q_prob)) %>%
      mutate(
        # CRUCIAAL: Spiegel de Y-as massa, niet de X-as tijd
        defective_p = if_else(as.numeric(as.character(acc)) == 1, q_prob * p_acc, -(q_prob * p_acc)),
        source = source_name,
        accuracy_label = if_else(as.numeric(as.character(acc)) == 1, "Correct", "Error")
      )
    return(cdf_df)
  }
  
  all_cdf <- bind_rows(
    get_cdf_data(as.data.frame(obs), "Observed"),
    get_cdf_data(as.data.frame(pred), "Predicted")
  )
  
  # Sorteer op RT voor vloeiende lijnen
  all_cdf <- all_cdf %>% arrange(source, accuracy_label, q_rt)
  
  group_vars <- c("accuracy_label", "source")
  if(has_conflict && "congruency" %in% names(all_cdf)) group_vars <- c(group_vars, "congruency")
  
  p <- ggplot(all_cdf, aes(x = q_rt, y = defective_p, group = interaction(!!!syms(group_vars)))) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray60", size = 0.5) # De nul-lijn
  
  if (has_conflict && "congruency" %in% names(all_cdf)) {
    all_cdf$congruency <- as.factor(all_cdf$congruency)
    p <- p + geom_line(data = all_cdf %>% filter(source == "Predicted"), aes(color = congruency), size = 0.8, alpha = 0.7) +
      geom_point(data = all_cdf %>% filter(source == "Predicted"), aes(color = congruency, shape = accuracy_label), shape = 4, size = 3, stroke = 1.1) +
      geom_point(data = all_cdf %>% filter(source == "Observed"), aes(color = congruency, shape = accuracy_label), size = 2.5, alpha = 0.6) +
      scale_color_manual(values = c("1" = "dodgerblue3", "-1" = "firebrick3"), 
                         labels = c("1" = "Congruent", "-1" = "Incongruent"), name = "Condition")
  } else {
    p <- p + geom_line(data = all_cdf %>% filter(source == "Predicted"), color = "black", size = 0.8, alpha = 0.7) +
      geom_point(data = all_cdf %>% filter(source == "Predicted"), aes(shape = accuracy_label), shape = 4, size = 3, stroke = 1.1) +
      geom_point(data = all_cdf %>% filter(source == "Observed"), aes(shape = accuracy_label), size = 2.5, alpha = 0.6, color = "black")
  }
  
  p <- p + scale_shape_manual(values = c("Correct" = 16, "Error" = 17), name = "Outcome") +
    theme_bw(base_size = 12) +
    # Zoom in op het relevante RT gebied (bijv. vanaf 250ms tot max)
    coord_cartesian(xlim = c(250, max(all_cdf$q_rt) * 1.05)) +
    labs(title = "Quantile Probability Plot",
         subtitle = "Corrects (>0) <--- 0 ---> Errors (<0)",
         x = "Reaction Time (ms)", 
         y = "Defective Cumulative Probability") +
    theme(strip.background = element_rect(fill = "gray95", color = NA),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  if (!is.null(cond_cols)) p <- p + facet_wrap(vars(!!!syms(cond_cols)))
  
  return(p)
}

# g-square mirror plot
# facets by experimental condition and highlights which bins drive the cost
plot_gsquare_mirror <- function(final_proportions) {
  bin_list <- list()
  
  for (t_name in names(final_proportions)) {
    for (b_name in names(final_proportions[[t_name]])) {
      obs_vec  <- final_proportions[[t_name]][[b_name]]$obs
      pred_vec <- final_proportions[[t_name]][[b_name]]$pred
      
      # calculate local g-square contribution per bin
      eps <- 1e-9
      g2_part <- obs_vec * log(obs_vec / (pred_vec + eps))
      
      bin_list[[paste0(t_name, b_name)]] <- data.frame(
        Obs      = obs_vec,
        Pred     = pred_vec,
        G2_Pijn  = abs(g2_part), 
        Bin_Idx  = 1:length(obs_vec),
        Block    = b_name,
        Target   = t_name
      )
    }
  }
  
  df <- bind_rows(bin_list)
  
  # calculate shared limit for both axes to keep it a perfect square
  max_val <- max(c(df$Obs, df$Pred), na.rm = TRUE) * 1.1
  
  p <- ggplot(df, aes(x = Obs, y = Pred, color = Bin_Idx)) +
    # perfection line
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
    # the bins
    geom_point(aes(size = G2_Pijn), alpha = 0.6) +
    # scaling & colors
    scale_color_viridis_c(option = "plasma", name = "Bin Order\n(Fast -> Slow)") +
    scale_size_continuous(name = "G-Square\nContribution") +
    # FIX: Remove scales = "free" to allow coord_fixed()
    facet_wrap(~Block) + 
    coord_fixed(xlim = c(0, max_val), ylim = c(0, max_val)) +
    theme_minimal(base_size = 12) +
    labs(title = "G-Square Diagnostic: Observed vs Predicted Mass",
         subtitle = "Faceted by condition. Larger points indicate areas of poor fit.",
         x = "Observed Probability Mass",
         y = "Predicted Probability Mass") +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "gray95", color = NA),
          strip.text = element_text(face = "bold"))
  
  return(p)
}

# plots the probability mass per bin, split by accuracy and congruency
# universal bin mass mirror plot
# handles both conflict (congruency) and standard (accuracy only) models
plot_bin_mass_mirror <- function(final_proportions) {
  all_data <- list()
  
  for (t_name in names(final_proportions)) {
    for (b_name in names(final_proportions[[t_name]])) {
      
      block_data <- final_proportions[[t_name]][[b_name]]
      obs_vec  <- as.numeric(block_data$obs)
      pred_vec <- as.numeric(block_data$pred)
      bin_names <- names(block_data$obs)
      
      # recovery logic if names are stripped during aggregation
      if (is.null(bin_names)) {
        n_bins <- length(obs_vec)
        # assume standard 4-way split if length matches, else 2-way
        if (n_bins %% 4 == 0) {
          bin_names <- rep(c("0.-1", "1.-1", "0.1", "1.1"), each = n_bins/4)
        } else {
          bin_names <- rep(c("0", "1"), each = n_bins/2)
        }
      }
      
      tmp_df <- data.frame(raw_label = bin_names, obs = obs_vec, pred = pred_vec, block = b_name)
      
      # --- smart parsing of labels ---
      # check if the labels contain a dot (indicating acc.cong interaction)
      has_dot <- any(grepl("\\.", tmp_df$raw_label))
      
      if (has_dot) {
        tmp_df <- tmp_df %>%
          separate(raw_label, into = c("acc_val", "cong_val"), sep = "\\.", extra = "drop") %>%
          mutate(congruency_label = case_when(cong_val == "1" ~ "Congruent", 
                                              cong_val == "-1" ~ "Incongruent", 
                                              TRUE ~ as.character(cong_val)))
      } else {
        tmp_df <- tmp_df %>%
          mutate(acc_val = raw_label, congruency_label = "Standard")
      }
      
      tmp_df <- tmp_df %>%
        group_by(acc_val, congruency_label) %>%
        mutate(
          bin_idx = 1:n(),
          accuracy_label = ifelse(acc_val == "1", "Correct", "Error"),
          plot_x = if_else(acc_val == "1", as.numeric(bin_idx), -as.numeric(bin_idx))
        ) %>%
        ungroup()
      
      all_data[[paste0(t_name, b_name)]] <- tmp_df
    }
  }
  
  df <- bind_rows(all_data) %>% arrange(block, accuracy_label, congruency_label, plot_x)
  
  # dynamic color palette: if only one 'congruency', use black/gray
  n_congs <- n_distinct(df$congruency_label)
  color_vals <- if(n_congs > 1) c("Congruent" = "dodgerblue3", "Incongruent" = "firebrick3") else c("Standard" = "black")
  
  ggplot(df, aes(x = plot_x, color = congruency_label, group = interaction(accuracy_label, congruency_label, block))) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80") +
    geom_point(aes(y = obs, shape = accuracy_label), size = 3, alpha = 0.5) +
    geom_line(aes(y = pred), size = 1) +
    geom_point(aes(y = pred, shape = accuracy_label), shape = 4, size = 3, stroke = 1.2) +
    facet_wrap(~block) +
    scale_color_manual(values = color_vals, name = "Condition") +
    scale_shape_manual(values = c("Correct" = 16, "Error" = 17), name = "Outcome") +
    theme_minimal(base_size = 12) +
    labs(title = "G-Square Fit: Probability Mass per Bin",
         subtitle = "Errors (Left) <--- 0 ---> Corrects (Right) | Symbols = Obs, Crosses = Pred",
         x = "Bin Index (Fastest near 0)", y = "Joint Probability Mass") +
    theme(strip.background = element_rect(fill = "gray95", color = NA),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

# Plot Confidence Judgment (CJ) Distribution
# Generates a barplot of observed vs predicted confidence rating frequencies,
# separated for Correct and Error trials.
plot_cj_distribution <- function(obs, pred, main_title = "Confidence Rating Distribution") {
  
  # Helper to calculate proportions
  calc_props <- function(df, source_name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      filter(!is.na(cj), !is.na(acc)) %>%
      group_by(acc, cj) %>%
      summarise(count = n(), .groups = "drop_last") %>%
      mutate(prop = count / sum(count),
             source = source_name,
             Accuracy = ifelse(as.numeric(as.character(acc)) == 1, "Correct", "Error")) %>%
      ungroup()
  }
  
  obs_df  <- calc_props(as.data.frame(obs), "Observed")
  pred_df <- calc_props(as.data.frame(pred), "Predicted")
  plot_df <- bind_rows(obs_df, pred_df)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = as.factor(cj), y = prop, fill = source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
    facet_wrap(~Accuracy) +
    scale_fill_manual(values = c("Observed" = "gray50", "Predicted" = "firebrick3")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 14) +
    labs(title = main_title, x = "Confidence Rating", y = "Relative Frequency", fill = "Data") +
    theme(strip.background = element_rect(fill = "gray90"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          panel.grid.major.x = element_blank())
  
  return(p)
}
