# ==============================================================================
# fit_stats.R: Group Statistics & Parameter Inspection
# ==============================================================================
# Author: Luc Vermeylen
# Description: 
# 1. Reconstructs marginal cell means from parameter betas.
# 2. Plots Varying Parameters by conditions if present (Spaghetti & Group Means).
# 3. Plots Global Parameters (i.e. not varying by condition).
# 4. Performs T-Test/RM-ANOVA for all parameters and saves reports.
# ==============================================================================

# ==============================================================================
# 0. DIRECTORY SETUP
# ==============================================================================
# --- 1. SETTINGS ---
OUTPUT_FOLDER   <- "" 

# Toggle to TRUE to see the red DEoptim bounds, or FALSE to zoom in on the data
SHOW_BOUNDARIES <- FALSE     

# Define derived parameters to calculate for every condition cell
# Example: list(v2 = function(p) p$v1 * p$v_ratio)
DERIVED_PARAMS <- list()

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)
source(here("R", "helper_functions.R"))

# ------------------------------------------------------------------------------
# 2. DATA LOADING & RECONSTRUCTION
# ------------------------------------------------------------------------------
RESULTS_DIR <- here("results", OUTPUT_FOLDER)
files <- list.files(RESULTS_DIR, pattern = "\\.rds$", full.names = TRUE)
valid_files <- c(); varying_info <- NULL; meta <- NULL

# First pass: identify valid files and the experimental design
for (f in files) {
  tmp <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(tmp) && is.list(tmp) && "info" %in% names(tmp)) {
    valid_files <- c(valid_files, f)
    if (is.null(varying_info)) { 
      meta <- tmp
      varying_info <- tmp$info$varying_params 
    }
  }
}

if (length(valid_files) == 0) stop("No valid fit result files found in: ", RESULTS_DIR)
cat("Found", length(valid_files), "valid fits. Reconstructing cell-means...\n")

# Detect all design variables from formulas or strings
all_design_vars <- unique(unlist(lapply(varying_info, function(x) if(inherits(x, "formula")) all.vars(x) else x)))

if (length(all_design_vars) == 0) {
  # No conditions: create a 1-row, 0-column dataframe so the loop runs exactly once
  full_design_grid <- data.frame(row.names = 1) 
} else {
  # Conditions exist: build the full interaction grid
  full_design_grid <- expand.grid(lapply(meta$observations[all_design_vars], function(x) levels(as.factor(x))))
  colnames(full_design_grid) <- all_design_vars
}

marginal_list <- list()
for (f in valid_files) {
  fit <- readRDS(f)
  
  # If this is a regression fit, use the 'best_betas'.
  # If it's a legacy cell-mean fit, use 'best_params'.
  params_to_use <- if(!is.null(fit$best_betas)) fit$best_betas else fit$best_params
  
  for (i in 1:nrow(full_design_grid)) {
    current_cell <- full_design_grid[i, , drop = FALSE]
    
    # Pluck parameters for this specific cell using the backend plucker
    p_lvl <- pluck_params(params_to_use, unlist(fit$constants), current_cell, 
                          fit$info$varying_params, fit$info$model)
    
    # Calculate Derived Parameters
    for (new_p in names(DERIVED_PARAMS)) { p_lvl[[new_p]] <- DERIVED_PARAMS[[new_p]](p_lvl) }
    
    res_row <- as.data.frame(p_lvl)
    res_row$subject_id <- as.character(fit$info$subject)
    res_row <- cbind(res_row, current_cell)
    marginal_list[[paste0(fit$info$subject, "_", i)]] <- res_row
  }
}

all_data_long <- bind_rows(marginal_list) %>%
  pivot_longer(cols = -c(subject_id, all_of(all_design_vars)), names_to = "Parameter", values_to = "Value")

# ------------------------------------------------------------------------------
# 3. ANALYSIS & PLOTTING ENGINE
# ------------------------------------------------------------------------------
model_bounds <- get_bounds(meta$info$model)
fixed_names  <- names(meta$param_info$fixed)

# Parameters to analyze (Exclude fixed constants, keep estimated ones)
params_to_analyze <- all_data_long %>%
  filter(!Parameter %in% fixed_names) %>%
  group_by(Parameter) %>% filter(sd(Value, na.rm = TRUE) > 1e-12) %>%
  pull(Parameter) %>% unique()

plots_spaghetti <- list(); plots_average <- list(); plots_global <- list(); stats_results <- list()
n_subs <- n_distinct(all_data_long$subject_id)

for (p in params_to_analyze) {
  p_data <- all_data_long %>% filter(Parameter == p)
  p_data$subject_id <- as.factor(p_data$subject_id) # Ensure factor for ANOVA
  
  # Detect design factors for this specific parameter
  p_mapping <- if(p %in% names(varying_info)) varying_info[[p]] else NULL
  active_vars <- if(!is.null(p_mapping)) {
    if(inherits(p_mapping, "formula")) all.vars(p_mapping) else p_mapping
  } else if (p %in% names(DERIVED_PARAMS)) {
    vars_found <- c()
    for (v in all_design_vars) {
      if (n_distinct(p_data %>% group_by(.data[[v]]) %>% summarise(m=mean(Value)) %>% pull(m)) > 1) vars_found <- c(vars_found, v)
    }
    vars_found
  } else { c() }
  
  # --- A. STATISTICS ---
  if (n_subs > 1 && length(active_vars) > 0) {
    # Reconstruct the formula string to match the fit
    rhs <- if(!is.null(p_mapping) && inherits(p_mapping, "formula")) {
      gsub("^~\\s*", "", deparse(p_mapping))
    } else {
      paste(active_vars, collapse = " * ")
    }
    
    anova_form <- as.formula(paste("Value ~", rhs, "+ Error(subject_id / (", rhs, "))"))
    res_aov <- tryCatch({ summary(aov(anova_form, data = p_data)) }, error = function(e) NULL)
    
    if (!is.null(res_aov)) {
      for (stratum in names(res_aov)) {
        tab <- res_aov[[stratum]][[1]]
        if (!is.null(tab) && "F value" %in% colnames(tab)) {
          for (row_idx in 1:nrow(tab)) {
            term <- trimws(rownames(tab)[row_idx])
            if (term %in% c("Residuals", "(Intercept)")) next
            f_val <- tab[row_idx, "F value"]; p_val <- tab[row_idx, "Pr(>F)"]
            if (!is.null(f_val) && !is.na(f_val)) {
              pes <- tab[row_idx, "Sum Sq"] / (tab[row_idx, "Sum Sq"] + tab["Residuals", "Sum Sq"])
              stats_results[[paste0(p, term)]] <- data.frame(Parameter = p, Term = term, F = round(as.numeric(f_val), 2),
                                                             df = paste0(tab[row_idx, "Df"], ",", tab["Residuals", "Df"]), p = round(as.numeric(p_val), 4), PES = round(as.numeric(pes), 3))
            }
          }
        }
      }
    }
  }
  
  # --- B. PLOTTING ---
  if (length(active_vars) > 0) {
    # --- VARYING PARAMETER PLOTS ---
    x_var <- active_vars[1]
    col_var <- if(length(active_vars) > 1) setdiff(active_vars, x_var)[1] else NULL
    
    p_data_summarised <- p_data %>% group_by(subject_id, across(all_of(active_vars))) %>% summarise(Value = mean(Value, na.rm=T), .groups = "drop")
    p_base <- ggplot(p_data_summarised, aes(x = as.factor(.data[[x_var]]), y = Value)) + theme_classic() + labs(title = p, x = x_var, y = NULL)
    g_aes  <- if(is.null(col_var)) aes(group = 1) else aes(group = .data[[col_var]], color = .data[[col_var]])
    
    # 1. Spaghetti Plot
    p_spag <- p_base
    if (!is.null(col_var)) {
      p_spag <- p_spag + geom_line(aes(group = interaction(subject_id, .data[[col_var]]), color = .data[[col_var]]), alpha = 0.08) + scale_color_brewer(palette = "Set1")
    } else {
      p_spag <- p_spag + geom_line(aes(group = subject_id), alpha = 0.1, color="gray80")
    }
    
    # 2. Average Plot
    p_cln <- p_base
    if (!is.null(col_var)) p_cln <- p_cln + scale_color_brewer(palette = "Set1")
    
    # Add Summary layers
    add_summary <- function(plt, mapping) {
      plt + stat_summary(mapping, fun.data = mean_se, geom = "errorbar", width = 0.1, linewidth = 1, color = "black", position = position_dodge(0.1)) +
        stat_summary(mapping, fun = mean, geom = "line", linewidth = 1.2, position = position_dodge(0.1)) +
        stat_summary(mapping, fun = mean, geom = "point", size = 4, position = position_dodge(0.1))
    }
    plots_spaghetti[[p]] <- add_summary(p_spag, g_aes)
    plots_average[[p]]   <- add_summary(p_cln, g_aes)
    
  } else {
    # --- GLOBAL PARAMETER PLOTS ---
    p_data_unique <- p_data %>% distinct(subject_id, Value)
    
    p_bounds <- NULL
    if (SHOW_BOUNDARIES && p %in% names(model_bounds$lower)) {
      p_bounds <- c(model_bounds$lower[p], model_bounds$upper[p])
    }
    
    plots_global[[p]] <- ggplot(p_data_unique, aes(x = "", y = Value)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.4) +
      geom_jitter(width = 0.15, alpha = 0.5, color = "firebrick") +
      {if(!is.null(p_bounds)) geom_hline(yintercept = p_bounds, linetype = "dashed", color = "red", alpha = 0.4)} +
      theme_classic() + 
      labs(title = p, x = NULL, y = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
}

# ------------------------------------------------------------------------------
# 4. OUTPUT GENERATION
# ------------------------------------------------------------------------------
n_v <- length(plots_spaghetti)
grid_cols <- if(n_v == 4) 2 else ceiling(sqrt(n_v))

# Compile the plots into patchwork objects
p_global <- NULL; p_spag <- NULL; p_avg <- NULL

if(length(plots_global) > 0) {
  title_text <- if(SHOW_BOUNDARIES) "Global Parameters (and DEoptim Boundaries)" else "Global Parameters"
  p_global <- wrap_plots(plots_global, ncol = 3) + plot_annotation(title = title_text)
}

if(n_v > 0) {
  p_spag <- wrap_plots(plots_spaghetti, ncol = grid_cols, guides = "collect") + 
    plot_annotation(title = "Individual Trends") & theme(legend.position = 'bottom')
  
  p_avg <- wrap_plots(plots_average, ncol = grid_cols, guides = "collect") + 
    plot_annotation(title = "Group Means") & theme(legend.position = 'bottom')
}

# --- A. SCREEN OUTPUT ---
if (!is.null(p_global)) print(p_global)
if (!is.null(p_spag)) print(p_spag)
if (!is.null(p_avg)) print(p_avg)

# --- B. PDF & CSV EXPORT ---
pdf_path <- file.path(RESULTS_DIR, "parameter_statistics_report.pdf")
pdf(file = pdf_path, width = 10, height = 7)
if (!is.null(p_global)) print(p_global)
if (!is.null(p_spag)) print(p_spag)
if (!is.null(p_avg)) print(p_avg)
dev.off()
cat(sprintf("\nPlots saved to: %s\n", pdf_path))

if (length(stats_results) > 0) {
  cat("\n=================================================================================\n")
  cat("STATISTICAL SUMMARY\n")
  cat("=================================================================================\n")
  
  stats_df <- bind_rows(stats_results) %>% 
    mutate(sig = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", p < .1 ~ ".", TRUE ~ "ns"))
  
  print(as.data.frame(stats_df), row.names = FALSE)
  
  csv_path <- file.path(RESULTS_DIR, "parameter_anova_results.csv")
  write.csv(stats_df, csv_path, row.names = FALSE)
  cat(sprintf("Stats saved to: %s\n", csv_path))
}