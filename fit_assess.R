# ==============================================================================
# fit_assess.R
# Visual Fit Assessment & Data Aggregation
# Author: Luc Vermeylen 
# Description: 
# 1. Aggregates fitted model results (.rds) from a directory.
# 2. Generates behavioral diagnostic plots (Distributions, CAF, Delta).
# 3. Assesses the "Mechanism" (Latent Trajectories) for the group.
# 4. Produces summary CSVs required for fit_compare.R and fit_stats.R.
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)
source("helper_functions.R")

# --- SETTINGS ---
RESULTS_DIR <- "r_dist_meta_full-ter-amp_I" 
SAVE_CSV    <- TRUE    

# ------------------------------------------------------------------------------
# 1. DATA AGGREGATION 
# ------------------------------------------------------------------------------
RESULTS_DIR <- file.path("results", RESULTS_DIR)
files <- list.files(RESULTS_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No RDS files found in: ", RESULTS_DIR)

obs_list <- list(); pred_list <- list(); param_list <- list(); metrics_list <- list()

cat("Aggregating data from", length(files), "files...\n")

for (f in files) {
  fit <- tryCatch(readRDS(f), error = function(e) NULL)
  
  # Structural integrity check (skip summary files or corrupt files)
  if (is.null(fit) || !is.list(fit) || !"info" %in% names(fit)) {
    next
  }
  
  sub_id  <- fit$info$subject
  fit_cond <- if(!is.null(fit$info$condition)) as.character(fit$info$condition) else "All"
  fname   <- basename(f)
  
  try({
    # A. Trial data
    obs_list[[f]]  <- as.data.table(fit$observations)[, `:=`(subject_id = sub_id, fit_condition = fit_cond, file_id = fname)]
    pred_list[[f]] <- as.data.table(fit$predictions)[, `:=`(subject_id = sub_id, fit_condition = fit_cond, file_id = fname)]
    
    # B. Parameters (using best_params which are already translated to cell-means in template)
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
# 2. DESIGN DETECTION (Hybrid: handles strings and formulas)
# ------------------------------------------------------------------------------
meta <- readRDS(file.path(RESULTS_DIR, master_table$file_id[1])) 
varying_info <- meta$info$varying_params

is_regression <- any(sapply(varying_info, function(x) inherits(x, "formula")))
cond_cols <- if(length(varying_info) > 0) {
  if(is_regression) unique(unlist(lapply(varying_info, all.vars))) else unique(unlist(varying_info))
} else { NULL }

has_conflict <- "congruency" %in% names(all_obs)

# Create interaction column for plotting and effects
if (!is.null(cond_cols)) {
  existing_cond_cols <- intersect(cond_cols, names(all_obs))
  all_obs[,  interaction_col := interaction(.SD, drop = TRUE, sep = "."), .SDcols = existing_cond_cols]
  all_pred[, interaction_col := interaction(.SD, drop = TRUE, sep = "."), .SDcols = existing_cond_cols]
  lvls <- levels(all_obs$interaction_col)
} else {
  all_obs[, interaction_col := "Overall"]; all_pred[, interaction_col := "Overall"]; lvls <- "Overall"
}

# ------------------------------------------------------------------------------
# 3. BEHAVIORAL EFFECT SUMMARY (CONGRUENCY EFFECTS)
# ------------------------------------------------------------------------------
if (has_conflict) {
  cat("Generating behavioral effect barplots...\n")
  
  calc_eff <- function(d, type="RT") {
    acc_n <- as.numeric(as.character(d$acc))
    cong  <- as.numeric(as.character(d$congruency))
    if(type=="RT") {
      return(mean(d$rt[cong == -1 & acc_n == 1], na.rm=T) - mean(d$rt[cong == 1 & acc_n == 1], na.rm=T))
    }
    if(type=="ER") {
      er_i <- (1 - mean(acc_n[cong == -1], na.rm=T)) * 100
      er_c <- (1 - mean(acc_n[cong == 1], na.rm=T)) * 100
      return(er_i - er_c)
    }
  }
  
  plot_lvls <- c("Overall", lvls)
  grid_dims <- c(ceiling(length(plot_lvls)/4), min(length(plot_lvls), 4))
  
  # RT EFFECT
  par(mfrow = grid_dims, mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  for(l in plot_lvls) {
    d_o <- if(l=="Overall") all_obs else all_obs[interaction_col == l]
    d_p <- if(l=="Overall") all_pred else all_pred[interaction_col == l]
    barplot(c(calc_eff(d_o, "RT"), calc_eff(d_p, "RT")), names.arg=c("Obs","Pred"), 
            main=l, col=c("gray40", "firebrick3"), ylab="ms", border=NA)
    abline(h=0)
  }
  mtext("RT Congruency Effect (Incong - Cong)", outer=T, font=2, line=1)
  
  # ERROR EFFECT
  par(mfrow = grid_dims, mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  for(l in plot_lvls) {
    d_o <- if(l=="Overall") all_obs else all_obs[interaction_col == l]
    d_p <- if(l=="Overall") all_pred else all_pred[interaction_col == l]
    barplot(c(calc_eff(d_o, "ER"), calc_eff(d_p, "ER")), names.arg=c("Obs","Pred"), 
            main=l, col=c("gray60", "dodgerblue3"), ylab="%", border=NA)
    abline(h=0)
  }
  mtext("Error Rate Effect (% Incong - % Cong)", outer=T, font=2, line=1)
  par(mfrow=c(1,1)) 
}

# ------------------------------------------------------------------------------
# 4. MASTER FIT PLOTS
# ------------------------------------------------------------------------------
cat("Generating grand-average fit plots...\n")

# filter out non-parameter columns for group mean calculation
non_param_cols <- c("subject_id", "fit_condition", "file_id", "model", "cost", 
                    "bic", "aic", "n_obs", "n_free_params", "n_bins", "seed", "ntrials", "dt", "s")

group_means_all <- master_table %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(non_param_cols)) %>%
  colMeans(na.rm = TRUE)

# 1. Collapsed Grand Average (Across all conditions)
plot_fit(
  obs            = all_obs, 
  pred           = all_pred, 
  varying_params = list(), # forces overall mode
  model_name     = meta$info$model,
  best_params    = group_means_all, 
  constants      = meta$constants,
  types          = c("dist", "caf", "delta") 
)

# 2. Detailed Grid (Split by condition)
plot_fit(
  obs            = all_obs, 
  pred           = all_pred, 
  varying_params = varying_info, 
  model_name     = meta$info$model,
  best_params    = group_means_all, 
  constants      = meta$constants,
  types          = c("dist", "caf", "delta", "mechanism") 
)

# ------------------------------------------------------------------------------
# 5. EXPORT
# ------------------------------------------------------------------------------
if(SAVE_CSV) {
  metrics_report <- master_table %>%
    select(subject_id, model, cost, bic, aic, n_obs, n_bins, n_free_params) %>%
    mutate(folder = RESULTS_DIR)
  
  write.csv(master_table, file.path(RESULTS_DIR, "master_parameters_report.csv"), row.names=F)
  write.csv(metrics_report, file.path(RESULTS_DIR, "fit_metrics_summary.csv"), row.names=F)
  saveRDS(master_table, file.path(RESULTS_DIR, "aggregated_results.rds"))
  
  cat("\nResults aggregated and reports saved to:", RESULTS_DIR, "\n")
}

print("Assessment Complete.")