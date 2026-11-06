# ==============================================================================
# 2_fit_assess.R
# Visual Fit Assessment & Data Aggregation (Model-Agnostic)
# Author: Luc Vermeylen 
# Description: 
# 1. Aggregates fitted model results (.rds) from a results directory.
# 2. Generates behavioral diagnostic plots adapted to the model family:
#    - Conflict models: RT distribution (by congruency), CAF, Delta plot, Mechanism.
#    - Confidence models: Decision RT, Confidence RT, Confidence Rating Mass.
# 3. Produces summary CSVs necessary for fit_stats.R and fit_compare.R.
# ==============================================================================

# ==============================================================================
# 0. DIRECTORY SETUP
# ==============================================================================

OUTPUT_FOLDER <- "fit_eeg46_ro_coupled_L2_fixSr7" # <--- Edit your folder name directly here!

# Optional CLI override (only applies if run from terminal with arguments)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && nzchar(args[1])) {
  OUTPUT_FOLDER <- args[1]
}

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)
source(here("R", "helper_functions.R"))

# ------------------------------------------------------------------------------
# 1. DATA AGGREGATION 
# ------------------------------------------------------------------------------
SAVE_CSV    <- TRUE    
RESULTS_DIR <- here("results", OUTPUT_FOLDER)
files <- list.files(RESULTS_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No RDS files found in: ", RESULTS_DIR)

obs_list <- list(); pred_list <- list(); param_list <- list(); metrics_list <- list()

cat("Aggregating data from", length(files), "files...\n")

for (f in files) {
  fit <- tryCatch(readRDS(f), error = function(e) NULL)
  
  if (is.null(fit) || !is.list(fit) || !"info" %in% names(fit)) {
    next
  }
  
  sub_id   <- as.character(fit$info$subject)
  fit_cond <- if(!is.null(fit$info$conditions)) paste(fit$info$conditions, collapse = "_") else "All"
  fname    <- basename(f)
  
  try({
    # A. Trial data
    obs_list[[f]]  <- as.data.table(fit$observations)[, `:=`(subject_id = sub_id, fit_condition = fit_cond, file_id = fname)]
    pred_list[[f]] <- as.data.table(fit$predictions)[, `:=`(subject_id = sub_id, fit_condition = fit_cond, file_id = fname)]
    
    # B. Parameters
    full_p <- c(fit$best_params, unlist(fit$constants))
    param_list[[f]] <- as.data.frame(t(full_p)) %>% 
      mutate(subject_id = sub_id, fit_condition = fit_cond, file_id = fname, cost = fit$fit_metrics$best_cost)
    
    # C. Fit Metrics
    metrics_list[[f]] <- data.frame(
      subject_id = sub_id, fit_condition = fit_cond, file_id = fname,
      model = fit$info$model, 
      n_free_params = fit$fit_metrics$n_free_params,
      n_bins = if(!is.null(fit$fit_metrics$n_bins)) fit$fit_metrics$n_bins else NA,
      n_obs = fit$fit_metrics$n_observations, 
      cost = fit$fit_metrics$best_cost,
      bic = fit$fit_metrics$bic, aic = fit$fit_metrics$aic, 
      stringsAsFactors = FALSE
    )
  }, silent = FALSE)
}

# Combine and select the best fit attempt per subject
master_table_full <- bind_rows(param_list) %>% 
  left_join(bind_rows(metrics_list), by = c("subject_id", "fit_condition", "file_id", "cost"))

unique_subs <- unique(master_table_full$subject_id)
if ("GROUP" %in% unique_subs && length(unique_subs) > 1) {
  warning(
    "\n[SAFETY TRIGGERED] Detected a 'GROUP' fit mixed with individual fits in the same folder.\n",
    "Excluding 'GROUP' fit from this individual aggregation.\n",
    call. = FALSE
  )
  master_table_full <- master_table_full %>% filter(subject_id != "GROUP")
}

master_table <- master_table_full %>% 
  group_by(subject_id, fit_condition) %>% 
  filter(cost == min(cost)) %>% 
  slice(1) %>% 
  ungroup()

all_obs  <- rbindlist(obs_list)[file_id %in% master_table$file_id]
all_pred <- rbindlist(pred_list)[file_id %in% master_table$file_id]

# ------------------------------------------------------------------------------
# 2. DESIGN & MODEL CAPABILITIES DETECTION
# ------------------------------------------------------------------------------
meta <- readRDS(here("results", OUTPUT_FOLDER, master_table$file_id[1])) 
varying_info <- meta$info$varying_params
model_name   <- meta$info$model
model_info   <- get_model_info(model_name)

is_regression <- any(sapply(varying_info, function(x) inherits(x, "formula")))
cond_cols <- if(length(varying_info) > 0) {
  if(is_regression) unique(unlist(lapply(varying_info, all.vars))) else unique(unlist(varying_info))
} else { NULL }

# Create interaction column for plotting across conditions
if (!is.null(cond_cols)) {
  existing_cond_cols <- intersect(cond_cols, names(all_obs))
  all_obs[,  interaction_col := interaction(.SD, drop = TRUE, sep = "."), .SDcols = existing_cond_cols]
  all_pred[, interaction_col := interaction(.SD, drop = TRUE, sep = "."), .SDcols = existing_cond_cols]
  lvls <- levels(all_obs$interaction_col)
} else {
  all_obs[, interaction_col := "Overall"]; all_pred[, interaction_col := "Overall"]; lvls <- "Overall"
}

plot_lvls <- unique(c("Overall", lvls))

# Group means for latent trajectory / mechanism plotting
non_param_cols <- c("subject_id", "fit_condition", "file_id", "model", "cost", 
                    "bic", "aic", "n_obs", "n_free_params", "n_bins", "folder")

group_means_all <- master_table %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(non_param_cols)) %>%
  colMeans(na.rm = TRUE)

# ------------------------------------------------------------------------------
# 3. MASTER FIT PLOTS (Adaptive by Model Family)
# ------------------------------------------------------------------------------
draw_all_diagnostics <- function() {
  for (l in plot_lvls) {
    d_o <- if(l=="Overall") all_obs else all_obs[interaction_col == l]
    d_p <- if(l=="Overall") all_pred else all_pred[interaction_col == l]
    if (nrow(d_o) == 0) next
    
    title_prefix <- if(l=="Overall") "Overall Data" else paste("Condition:", l)
    
    # --- A. PRIMARY DECISION RT ---
    facet_col <- if(model_info$has_conflict) "congruency" else NULL
    plot_dist(obs = as.data.frame(d_o), pred = as.data.frame(d_p), 
              val_col = "rt", facet_col = facet_col, split_by_acc = TRUE, 
              main_title = paste(title_prefix, "- Decision RT"))
    
    # --- B. CONFLICT DIAGNOSTICS ---
    if (model_info$has_conflict && "congruency" %in% names(d_o)) {
      # 1. Conditional Accuracy Function (CAF)
      plot_conditional_function(obs = as.data.frame(d_o), pred = as.data.frame(d_p), 
                                x_col = "rt", y_col = "acc", group_col = "congruency", 
                                main_title = paste(title_prefix, "- CAF (Accuracy vs RT)"))
      
      # 2. Delta Plot
      plot_delta(obs = as.data.frame(d_o), pred = as.data.frame(d_p), 
                 main_title = paste(title_prefix, "- Delta Plot"))
      
      # 3. DMC Mechanism Trajectory
      if (model_name == "DMC") {
        if (!is.null(cond_cols) && l != "Overall") {
          grid_full <- unique(all_obs[, c(cond_cols, "interaction_col"), with = FALSE])
          row_idx <- which(grid_full$interaction_col == l)[1]
          if (!is.na(row_idx)) {
            current_cell <- as.data.frame(grid_full[row_idx, cond_cols, with = FALSE])
            p_cell <- pluck_params(group_means_all, meta$constants, current_cell, varying_info, model_name)
            plot_dmc_mechanism(p_cell, meta$constants, main_title = paste("DMC Mechanism:", l))
          }
        } else {
          plot_dmc_mechanism(group_means_all, meta$constants, main_title = paste("DMC Mechanism:", title_prefix))
        }
      }
    }
    
    # --- C. CONFIDENCE DIAGNOSTICS ---
    if (model_info$has_confidence) {
      if ("rtconf" %in% names(d_o)) {
        plot_dist(obs = as.data.frame(d_o), pred = as.data.frame(d_p), 
                  val_col = "rtconf", split_by_acc = TRUE, 
                  main_title = paste(title_prefix, "- Confidence RT"))
      }
      if ("cj" %in% names(d_o)) {
        p_cj <- plot_cj_distribution(as.data.frame(d_o), as.data.frame(d_p), 
                                     main_title = paste(title_prefix, "- Confidence Rating Mass"))
        if (!is.null(p_cj)) print(p_cj)
      }
    }
  }
}

cat("Generating grand-average diagnostic plots...\n")

# --- SCREEN OUTPUT ---
draw_all_diagnostics()

# --- PDF OUTPUT ---
pdf_path <- file.path(RESULTS_DIR, "grand_average_fit_assessment.pdf")
pdf(file = pdf_path, width = 10, height = 6)
draw_all_diagnostics()
dev.off()
cat(sprintf("PDF Report saved to: %s\n", pdf_path))

# ------------------------------------------------------------------------------
# 4. EXPORT SUMMARY REPORTS
# ------------------------------------------------------------------------------
if(SAVE_CSV) {
  metrics_report <- master_table %>%
    select(subject_id, model, cost, bic, aic, n_obs, n_bins, n_free_params) %>%
    mutate(folder = RESULTS_DIR)
  
  write.csv(master_table, file.path(RESULTS_DIR, "master_parameters_report.csv"), row.names = FALSE)
  write.csv(metrics_report, file.path(RESULTS_DIR, "fit_metrics_summary.csv"), row.names = FALSE)
  
  cat("\nResults aggregated and reports saved to:", RESULTS_DIR, "\n")
}

print("Assessment Complete.")
