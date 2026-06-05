# ==============================================================================
# fit_compare.R: Robust Model Comparison (Hybrid Compatible)
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# --- 1. SETTINGS ---
model_folders <- c(
  #"r_dist_meta_full",
  #"r_dist_meta_full-v",
  #"r_dist_meta_full-a",
  #"r_dist_meta_full-ter",
  #"r_dist_meta_full-amp"
  #"r_dist_meta_full-ter-v",
  #"r_dist_meta_full-ter-a",
  "r_dist_meta_full-ter-amp",
  #"r_dist_meta_full-ter-amp-v",
  #"r_dist_meta_full-ter-amp-a"
  "r_dist_meta_full-ter-amp_I"
  # Add other folders here
)

model_labels <- NULL 
CSV_NAME <- "fit_metrics_summary.csv" 

# ------------------------------------------------------------------------------
# 2. DATA LOADING & CLEANING
# ------------------------------------------------------------------------------
cat("Loading model data...\n")
all_model_data <- list()
model_folders <- list.dirs("results", full.names = TRUE, recursive = FALSE)

if (is.null(model_labels)) model_labels <- model_folders

for (i in seq_along(model_folders)) {
  folder <- model_folders[i]
  path <- file.path(folder, CSV_NAME)
  if (!file.exists(path)) path <- file.path(folder, "master_parameters_report.csv")
  
  if (!file.exists(path)) {
    warning(paste("Folder skipped (no CSV found):", folder))
    next
  }
  
  df <- read.csv(path) %>% 
    select(any_of(c("subject_id", "cost", "bic", "aic", "n_obs", "n_free_params", "n_bins"))) %>%
    mutate(model_label = model_labels[i], folder_id = folder)
  
  all_model_data[[folder]] <- df
}

comparison_df <- bind_rows(all_model_data) %>% mutate(subject_id = as.character(subject_id))

# --- Common Subjects Check ---
n_models_loaded <- length(unique(comparison_df$model_label))
common_subs <- comparison_df %>% 
  group_by(subject_id) %>% 
  filter(n() == n_models_loaded) %>% 
  pull(subject_id) %>% unique()

if (length(common_subs) < length(unique(comparison_df$subject_id))) {
  cat(sprintf("Note: %d subjects excluded (incomplete data across models).\n", 
              length(unique(comparison_df$subject_id)) - length(common_subs)))
  comparison_df <- comparison_df %>% filter(subject_id %in% common_subs)
}

# ------------------------------------------------------------------------------
# 3. MODEL SELECTION CALCULATION
# ------------------------------------------------------------------------------

# Subject level winner
subject_selection <- comparison_df %>%
  group_by(subject_id) %>%
  mutate(delta_bic = bic - min(bic), is_winner = (delta_bic == 0)) %>%
  ungroup()

# Group level summary
group_summary <- comparison_df %>%
  group_by(model_label) %>%
  summarise(
    N         = n(),
    # We take the mean of n_bins because dynamic binning can cause slight variations per sub
    Avg_Bins  = round(mean(n_bins, na.rm = TRUE), 1),
    Pars      = first(n_free_params),
    Cost      = mean(cost),                        
    Pen       = mean(bic - cost),                  
    mBIC      = mean(bic),                         
    BIC       = sum(bic),                          
    Wins      = sum(subject_selection$is_winner[subject_selection$model_label == model_label]),
    .groups = 'drop'
  ) %>%
  mutate(dBIC = BIC - min(BIC)) %>%
  mutate(Weight = exp(-0.5 * dBIC) / sum(exp(-0.5 * dBIC)),
         Raw_Ratio = max(Weight) / Weight) %>%
  mutate(Ratio = case_when(dBIC == 0 ~ "1.0", dBIC > 20 ~ ">1000", TRUE ~ as.character(round(Raw_Ratio, 1)))) %>%
  select(Model = model_label, BIC, dBIC, Weight, Ratio, Wins, Cost, Pen, mBIC, Pars, Bins = Avg_Bins, N) %>%
  arrange(BIC)

# ------------------------------------------------------------------------------
# 4. VISUALIZATION (Robust to Outliers)
# ------------------------------------------------------------------------------
# Clean data for plotting (remove penalty infinities if any)
plot_df <- comparison_df %>% filter(bic < 1e6)

p1 <- ggplot(plot_df, aes(x = reorder(model_label, bic), y = bic, fill = model_label)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.8) +
  coord_cartesian(ylim = c(quantile(plot_df$bic, 0.05, na.rm=T)*0.8, quantile(plot_df$bic, 0.95, na.rm=T)*1.2)) +
  theme_minimal() + labs(title = "Fit Quality (BIC)", x = "Model", y = "BIC") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(group_summary, aes(x = reorder(Model, -Wins), y = Wins, fill = Model)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() + labs(title = "Subject Preferences", x = "Model", y = "N Wins") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(subject_selection %>% filter(delta_bic < 1e6), aes(x = reorder(model_label, delta_bic), y = delta_bic, fill = model_label)) +
  geom_violin(alpha = 0.5, draw_quantiles = c(0.5)) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "red") + 
  coord_cartesian(ylim = c(0, 50)) +
  theme_minimal() + labs(title = "Evidence Gap (dBIC)", x = "Model", y = "dBIC from Winner") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print((p1 | p2) / p3)

# ------------------------------------------------------------------------------
# 5. FINAL REPORTING
# ------------------------------------------------------------------------------
cat("\n============================================================================================\n")
cat("MODEL COMPARISON SUMMARY\n")
cat("============================================================================================\n")

report_table <- group_summary %>%
  mutate(across(c(BIC, dBIC, Cost, Pen, mBIC), ~round(., 1)),
         Weight = sprintf("%.3f", Weight))

print(as.data.frame(report_table), row.names = FALSE)

# Robust Validity Check
unique_bin_avgs <- unique(group_summary$Bins)
if (length(unique_bin_avgs) > 1) {
  # We allow a small margin for dynamic binning (e.g., difference < 2 is OK)
  if (max(unique_bin_avgs) - min(unique_bin_avgs) > 2) {
    cat("\n!!! WARNING: Bins differ significantly between models. BIC comparison may be invalid. !!!\n")
  } else {
    cat("\nNote: Small variations in bin counts detected due to dynamic binning. Comparison remains valid.\n")
  }
}

cat("============================================================================================\n")
write.csv(group_summary, "summary_model_comparison_detailed.csv", row.names = FALSE)