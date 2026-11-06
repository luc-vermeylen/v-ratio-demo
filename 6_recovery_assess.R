# ==============================================================================
# 6_recovery_assess.R
# PARAMETER RECOVERY: Assessment, Statistics & Diagnostic Plots (Robust & Clean)
# Author: Luc Vermeylen
# ==============================================================================

# ==============================================================================
# 0. DIRECTORY SETUP & SETTINGS
# ==============================================================================
# Unconditional assignment: Change your folder name directly here!
OUTPUT_FOLDER <- "recov_ro_L0_full9" 

# Optional CLI override
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && nzchar(args[1])) {
  OUTPUT_FOLDER <- args[1]
}

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(data.table)

RESULTS_DIR <- here("results", OUTPUT_FOLDER)
cat(sprintf("\n=== RUNNING RECOVERY ASSESSMENT FOR: %s ===\n", OUTPUT_FOLDER))

GT_FILE  <- file.path(RESULTS_DIR, "ground_truth.csv")
FIT_FILE <- file.path(RESULTS_DIR, "master_parameters_report.csv")

if (!file.exists(GT_FILE)) stop("ground_truth.csv not found in: ", RESULTS_DIR)
if (!file.exists(FIT_FILE)) stop("master_parameters_report.csv not found in: ", RESULTS_DIR)

# ==============================================================================
# 1. LOAD & CANONICAL PARAMETER MATCHING
# ==============================================================================
cat("Loading ground truth and fitted parameters...\n")

df_gt  <- fread(GT_FILE)
df_fit <- fread(FIT_FILE)

df_gt[, subject_id := as.character(subject_id)]
df_fit[, subject_id := as.character(subject_id)]

ignore_cols <- c("subject_id", "fit_condition", "file_id", "cost", "model", 
                 "n_free_params", "n_bins", "n_obs", "bic", "aic", "folder",
                 "ntrials", "s", "dt", "seed", "st", "sz", "z", "gate", "prop_congruent")

raw_gt_names  <- setdiff(names(df_gt), ignore_cols)
raw_fit_names <- setdiff(names(df_fit), ignore_cols)

# Decode mangled names cleanly:
# Restores '..' back to negative sign '-' so v..1 becomes v:-1 and never collides with v.1
decode_name <- function(nm) {
  # If name has "..", it came from make.names turning ':-' into '..'
  nm <- sub("^(v|a|ter|sv)\\.\\.([0-9]+)$", "\\1:-\\2", nm) # v..1 -> v:-1
  nm <- sub("^(v|a|ter|sv)\\.([0-9]+)$",   "\\1:\\2",  nm) # v.1  -> v:1
  # Normalize separators
  nm <- gsub("[.]", "_", nm)
  return(nm)
}

# Canonical key extractor (extracts base parameter and level suffix)
get_canonical_key <- function(nm) {
  d <- decode_name(nm)
  if (grepl(":", d)) {
    parts <- strsplit(d, ":")[[1]]
    p_base <- parts[1]
    # Extract condition level (e.g. from 'congruency_-1' -> '-1', from 'Difficulty_Hard' -> 'Hard')
    lvl <- sub(".*_", "", parts[2])
    return(paste0(p_base, "@", lvl))
  } else {
    return(d)
  }
}

key_gt  <- setNames(sapply(raw_gt_names, get_canonical_key), raw_gt_names)
key_fit <- setNames(sapply(raw_fit_names, get_canonical_key), raw_fit_names)

# Match on canonical keys (BaseParam @ Level)
shared_keys <- intersect(unname(key_gt), unname(key_fit))

param_pairs <- list()
for (k in shared_keys) {
  g_orig <- names(key_gt)[which(key_gt == k)[1]]
  f_orig <- names(key_fit)[which(key_fit == k)[1]]
  param_pairs[[g_orig]] <- f_orig
}

if (length(param_pairs) == 0) stop("Could not find matching parameters between ground truth and fit report!")

# Filter out fixed constants (parameters with zero variance)
active_params <- c()
for (gp in names(param_pairs)) {
  vals_gt <- as.numeric(df_gt[[gp]])
  if (sd(vals_gt, na.rm = TRUE) > 1e-6) {
    active_params <- c(active_params, gp)
  }
}

cat(sprintf("Identified %d estimated parameters to evaluate: %s\n", 
            length(active_params), paste(active_params, collapse = ", ")))

# ==============================================================================
# 2. MERGE & RESTRUCTURE DATA (Distinct Readable Labels)
# ==============================================================================
format_label <- function(raw_nm) {
  d <- decode_name(raw_nm)
  if (grepl(":", d)) {
    parts <- strsplit(d, ":")[[1]]
    p_base <- parts[1]
    lvl <- sub(".*_", "", parts[2])
    return(paste0(p_base, " [", lvl, "]"))
  } else {
    return(raw_nm)
  }
}

long_records <- list()
subs_in_both <- intersect(df_gt$subject_id, df_fit$subject_id)

if (length(subs_in_both) == 0) stop("No matching subject IDs between ground truth and fit files!")
cat(sprintf("Evaluating recovery across %d subjects...\n", length(subs_in_both)))

for (gp in active_params) {
  fp <- param_pairs[[gp]]
  
  sub_gt  <- df_gt[subject_id %in% subs_in_both, .(subject_id, True = as.numeric(get(gp)))]
  sub_fit <- df_fit[subject_id %in% subs_in_both, .(subject_id, Estimated = as.numeric(get(fp)))]
  
  merged_p <- merge(sub_gt, sub_fit, by = "subject_id")
  merged_p[, Parameter := format_label(fp)]
  long_records[[gp]] <- merged_p
}

recovery_df <- rbindlist(long_records)

# ==============================================================================
# 3. COMPUTE RECOVERY METRICS
# ==============================================================================
metrics_list <- list()
unique_labels <- unique(recovery_df$Parameter)

for (p_lbl in unique_labels) {
  sub_d <- recovery_df[Parameter == p_lbl]
  n_pts <- nrow(sub_d)
  
  r_val   <- if (n_pts > 2) cor(sub_d$True, sub_d$Estimated, method = "pearson", use = "complete.obs") else NA
  rho_val <- if (n_pts > 2) cor(sub_d$True, sub_d$Estimated, method = "spearman", use = "complete.obs") else NA
  
  rmse_val <- sqrt(mean((sub_d$Estimated - sub_d$True)^2, na.rm = TRUE))
  bias_val <- mean(sub_d$Estimated - sub_d$True, na.rm = TRUE)
  
  lm_fit    <- if (n_pts > 2) tryCatch(lm(Estimated ~ True, data = sub_d), error = function(e) NULL) else NULL
  slope_val <- if (!is.null(lm_fit)) coef(lm_fit)[2] else NA
  r2_val    <- if (!is.null(lm_fit)) summary(lm_fit)$r.squared else NA
  
  metrics_list[[p_lbl]] <- data.frame(
    Parameter    = p_lbl,
    N_Subjects   = n_pts,
    Pearson_r    = round(r_val, 3),
    Spearman_rho = round(rho_val, 3),
    RMSE         = round(rmse_val, 3),
    Bias         = round(bias_val, 3),
    Slope        = round(slope_val, 3),
    R2           = round(r2_val, 3),
    stringsAsFactors = FALSE
  )
}

metrics_table <- bind_rows(metrics_list)

# ==============================================================================
# 4. SCATTER PLOTS
# ==============================================================================
scatter_plots <- list()

for (p_lbl in unique_labels) {
  p_data   <- recovery_df[Parameter == p_lbl]
  p_metric <- metrics_table %>% filter(Parameter == p_lbl)
  
  all_vals <- c(p_data$True, p_data$Estimated)
  lim_min  <- min(all_vals, na.rm = TRUE)
  lim_max  <- max(all_vals, na.rm = TRUE)
  pad      <- (lim_max - lim_min) * 0.08
  lims     <- c(lim_min - pad, lim_max + pad)
  
  ann_text <- if (!is.na(p_metric$Pearson_r)) {
    sprintf("r = %.2f\nRMSE = %.2f\nBias = %.2f", p_metric$Pearson_r, p_metric$RMSE, p_metric$Bias)
  } else {
    sprintf("RMSE = %.2f\nBias = %.2f", p_metric$RMSE, p_metric$Bias)
  }
  
  plt <- ggplot(p_data, aes(x = True, y = Estimated)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    {if(nrow(p_data) > 2) geom_smooth(method = "lm", formula = y ~ x, color = "dodgerblue3", fill = "dodgerblue1", alpha = 0.15, linewidth = 1)} +
    geom_point(color = "black", fill = "gray30", shape = 21, size = 3, alpha = 0.8) +
    annotate("text", x = lims[1] + (lims[2] - lims[1]) * 0.05, y = lims[2] - (lims[2] - lims[1]) * 0.05,
             label = ann_text, hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
    coord_fixed(xlim = lims, ylim = lims) +
    theme_classic(base_size = 11) +
    labs(title = p_lbl, x = "Ground Truth", y = "Estimated") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text  = element_text(color = "black"))
  
  scatter_plots[[p_lbl]] <- plt
}

n_p <- length(scatter_plots)
n_cols <- min(3, ceiling(sqrt(n_p)))
p_grid <- wrap_plots(scatter_plots, ncol = n_cols) + 
  plot_annotation(title = paste("Parameter Recovery:", OUTPUT_FOLDER),
                  theme = theme(plot.title = element_text(size = 14, face = "bold"),
                                plot.subtitle = element_text(size = 11, color = "gray30")))

# ==============================================================================
# 5. PARAMETER CONFUSION MATRIX (Safeguarded against duplicates)
# ==============================================================================
conf_matrix_plot <- NULL

if (length(subs_in_both) >= 5 && length(active_params) > 1) {
  
  df_gt_aligned  <- df_gt[match(subs_in_both, df_gt$subject_id), ]
  df_fit_aligned <- df_fit[match(subs_in_both, df_fit$subject_id), ]
  
  mat_true <- as.matrix(df_gt_aligned[, active_params, with = FALSE])
  mat_est  <- as.matrix(df_fit_aligned[, unlist(param_pairs[active_params]), with = FALSE])
  
  cross_cor <- cor(mat_true, mat_est, use = "pairwise.complete.obs")
  
  clean_labels <- sapply(unlist(param_pairs[active_params]), format_label)
  
  rownames(cross_cor) <- paste0("True_", clean_labels)
  colnames(cross_cor) <- paste0("Est_", clean_labels)
  
  cor_long <- as.data.frame(as.table(cross_cor)) %>%
    rename(True_Param = Var1, Est_Param = Var2, Correlation = Freq)
  
  # Defensive: unique() ensures no duplicated levels can crash ggplot
  u_labels_rev <- rev(unique(paste0("True_", clean_labels)))
  u_labels     <- unique(paste0("Est_", clean_labels))
  
  cor_long$True_Param <- factor(cor_long$True_Param, levels = u_labels_rev)
  cor_long$Est_Param  <- factor(cor_long$Est_Param,  levels = u_labels)
  
  conf_matrix_plot <- ggplot(cor_long, aes(x = Est_Param, y = True_Param, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "firebrick3", mid = "white", high = "dodgerblue3", 
                         midpoint = 0, limits = c(-1, 1), name = "Correlation") +
    theme_minimal(base_size = 11) +
    labs(title = "Parameter Recovery Confusion Matrix",
         x = "Estimated Parameters", y = "True Parameters") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank())
}

# ==============================================================================
# 6. REPORTING & EXPORT
# ==============================================================================
cat("\n====================================================================\n")
cat("PARAMETER RECOVERY SUMMARY METRICS\n")
cat("====================================================================\n")
print(as.data.frame(metrics_table), row.names = FALSE)
cat("====================================================================\n")

metrics_csv <- file.path(RESULTS_DIR, "parameter_recovery_metrics.csv")
write.csv(metrics_table, metrics_csv, row.names = FALSE)

pdf_path <- file.path(RESULTS_DIR, "parameter_recovery_report.pdf")
pdf(file = pdf_path, width = 11, height = max(7, ceiling(n_p / n_cols) * 3.5))
print(p_grid)
if (!is.null(conf_matrix_plot)) print(conf_matrix_plot)
dev.off()

cat(sprintf("Reports saved in: %s\n", RESULTS_DIR))

# Print both to screen
print(p_grid)
if (!is.null(conf_matrix_plot)) print(conf_matrix_plot)