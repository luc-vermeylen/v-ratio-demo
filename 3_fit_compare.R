# ==============================================================================
# fit_compare.R: Model Comparison
# Author: Luc Vermeylen
# ==============================================================================

# ==============================================================================
# 0. DIRECTORY SETUP
# ==============================================================================

# Define which model folders you want to compare. 
# They must be located inside the "results/" directory.
# Set to "ALL" to automatically compare every folder present in the results directory.

FOLDERS_TO_COMPARE <- c("", "") # Example for specific folders: c("test_vratio_1", "test_vratio_2")

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# 2. DATA LOADING & CLEANING
# ------------------------------------------------------------------------------
cat("Loading model data...\n")

# Auto-detect or use specified folders
if (length(FOLDERS_TO_COMPARE) == 1 && FOLDERS_TO_COMPARE == "ALL") {
  model_folders <- list.dirs(here("results"), full.names = TRUE, recursive = FALSE)
} else {
  model_folders <- here("results", FOLDERS_TO_COMPARE)
}

if (length(model_folders) < 2) {
  stop("You need at least 2 model folders to run a comparison!")
}

# Clean up labels for plotting (remove path details from the name)
model_labels <- basename(model_folders)

CSV_NAME <- "fit_metrics_summary.csv" 

all_model_data <- list()

for (i in seq_along(model_folders)) {
  folder <- model_folders[i]
  path <- file.path(folder, CSV_NAME)
  
  # Fallback for older versions
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

if (length(all_model_data) == 0) stop("No valid data found to compare.")

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
  theme_minimal() + labs(title = "Fit Quality (Sum BIC)", x = "Model", y = "BIC") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(group_summary, aes(x = reorder(Model, -Wins), y = Wins, fill = Model)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() + labs(title = "Subject Preferences", x = "Model", y = "N Wins") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

safe_dbic <- subject_selection$delta_bic[subject_selection$delta_bic < 1e6]
dbic_upper_limit <- max(50, quantile(safe_dbic, 0.95, na.rm = TRUE) * 1.1)

p3 <- ggplot(subject_selection %>% filter(delta_bic < 1e6), aes(x = reorder(model_label, delta_bic), y = delta_bic, fill = model_label)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "red", linewidth = 0.8) + 
  coord_cartesian(ylim = c(0, dbic_upper_limit)) +
  theme_minimal() + labs(title = "Evidence Gap (delta BIC)", x = "Model", y = "dBIC from Winner") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


# Compile plot
final_plot <- (p1 | p2) / p3

# Print to screen
print(final_plot)

# Save to PDF
pdf_path <- here("results", "model_comparison_plots.pdf")
pdf(file = pdf_path, width = 10, height = 7)
print(final_plot)
dev.off()
cat(sprintf("\nPlots saved to: %s\n", pdf_path))

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

cat("============================================================================================\n")

csv_path <- here("results", "summary_model_comparison_detailed.csv")
write.csv(group_summary, csv_path, row.names = FALSE)
cat(sprintf("Summary table saved to: %s\n", csv_path))

# Robust Validity Check
unique_bin_avgs <- unique(group_summary$Bins)
if (length(unique_bin_avgs) > 1) {
  # We allow a small margin for dynamic binning (e.g., difference < 2 is OK)
  if (max(unique_bin_avgs) - min(unique_bin_avgs) > 2) {
    cat("\n!!! WARNING: Bins differ significantly between models. BIC comparison may be invalid. !!!\n")
  } else {
    cat("\nNote: Small variations in bin counts detected. Comparison remains valid.\n")
  }
}