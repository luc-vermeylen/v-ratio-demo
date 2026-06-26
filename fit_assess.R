# ==============================================================================
# fit_assess.R
# Visual Fit Assessment & Data Aggregation
# Author: Luc Vermeylen 
# Description: 
# 1. Aggregates fitted model results (.rds) from a given directory.
# 2. Generates behavioral diagnostic plots (RT Distributions, Confidence Rating Mass).
# 3. Produces summary CSVs
# ==============================================================================

# ==============================================================================
# 0. DIRECTORY SETUP
# ==============================================================================

OUTPUT_FOLDER <- "vratio" # specify the directory with your fit results here

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)
source(here("helper_functions.R"))

# ------------------------------------------------------------------------------
# 1. DATA AGGREGATION 
# ------------------------------------------------------------------------------
SAVE_CSV <- TRUE    
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
  fit_cond <- if(!is.null(fit$info$condition)) as.character(fit$info$condition) else "All"
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

# Combine and select the best fit attempt per subject (if multiple exist)
master_table_full <- bind_rows(param_list) %>% 
  left_join(bind_rows(metrics_list), by = c("subject_id", "fit_condition", "file_id", "cost"))

master_table <- master_table_full %>% 
  group_by(subject_id, fit_condition) %>% 
  filter(cost == min(cost)) %>% 
  slice(1) %>% 
  ungroup()

all_obs  <- rbindlist(obs_list)[file_id %in% master_table$file_id]
all_pred <- rbindlist(pred_list)[file_id %in% master_table$file_id]

# ------------------------------------------------------------------------------
# 2. DESIGN DETECTION
# ------------------------------------------------------------------------------
meta <- readRDS(here("results", OUTPUT_FOLDER, master_table$file_id[1])) 
varying_info <- meta$info$varying_params

is_regression <- any(sapply(varying_info, function(x) inherits(x, "formula")))
cond_cols <- if(length(varying_info) > 0) {
  if(is_regression) unique(unlist(lapply(varying_info, all.vars))) else unique(unlist(varying_info))
} else { NULL }

# Create interaction column for plotting and effects
if (!is.null(cond_cols)) {
  existing_cond_cols <- intersect(cond_cols, names(all_obs))
  all_obs[,  interaction_col := interaction(.SD, drop = TRUE, sep = "."), .SDcols = existing_cond_cols]
  all_pred[, interaction_col := interaction(.SD, drop = TRUE, sep = "."), .SDcols = existing_cond_cols]
  lvls <- levels(all_obs$interaction_col)
} else {
  all_obs[, interaction_col := "Overall"]; all_pred[, interaction_col := "Overall"]; lvls <- "Overall"
}

plot_lvls <- unique(c("Overall", lvls))

# ------------------------------------------------------------------------------
# 3. MASTER FIT PLOTS
# ------------------------------------------------------------------------------
draw_all_distributions <- function() {
  for (l in plot_lvls) {
    d_o <- if(l=="Overall") all_obs else all_obs[interaction_col == l]
    d_p <- if(l=="Overall") all_pred else all_pred[interaction_col == l]
    
    title_prefix <- if(l=="Overall") "Overall Data" else paste("Condition:", l)
    
    # 1. Primary RT Distribution
    plot_dist(obs = as.data.frame(d_o), pred = as.data.frame(d_p), 
              val_col = "rt", split_by_acc = TRUE, 
              main_title = paste(title_prefix, "- Decision RT"))
    
    # 2. Confidence RT Distribution
    if ("rtconf" %in% names(d_o)) {
      plot_dist(obs = as.data.frame(d_o), pred = as.data.frame(d_p), 
                val_col = "rtconf", split_by_acc = TRUE, 
                main_title = paste(title_prefix, "- Confidence RT"))
    }
    
    # 3. Confidence Rating Mass
    if ("cj" %in% names(d_o)) {
      p_cj <- plot_cj_distribution(as.data.frame(d_o), as.data.frame(d_p), 
                                   main_title = paste(title_prefix, "- Confidence Rating Mass"))
      if (!is.null(p_cj)) print(p_cj)
    }
  }
}

cat("Generating grand-average fit plots...\n")

# --- A. SCREEN OUTPUT ---
draw_all_distributions()

# --- B. PDF OUTPUT ---
pdf_path <- file.path(RESULTS_DIR, "grand_average_fit_assessment.pdf")
pdf(file = pdf_path, width = 10, height = 6)

draw_all_distributions()     # Put all distributions in PDF

dev.off()
cat(sprintf("PDF Report saved to: %s\n", pdf_path))

# ------------------------------------------------------------------------------
# 4. EXPORT
# ------------------------------------------------------------------------------
if(SAVE_CSV) {
  metrics_report <- master_table %>%
    select(subject_id, model, cost, bic, aic, n_obs, n_bins, n_free_params) %>%
    mutate(folder = RESULTS_DIR)
  
  write.csv(master_table, file.path(RESULTS_DIR, "master_parameters_report.csv"), row.names=F)
  write.csv(metrics_report, file.path(RESULTS_DIR, "fit_metrics_summary.csv"), row.names=F)

  cat("\nResults aggregated and reports saved to:", RESULTS_DIR, "\n")
}

print("Assessment Complete.")