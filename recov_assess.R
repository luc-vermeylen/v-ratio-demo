# ==============================================================================
# parameter recovery analysis
# ==============================================================================
# author: luc vermeylen
#
# description: 
# 1. aggregates all recovery rds files from a target directory.
# 2. extracts true and recovered coefficients (betas).
# 3. calculates pearson correlation and normalized root mean square error.
# 4. generates a faceted scatterplot for performance assessment.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. configuration ---
# folder name defined in the recovery_template.r
session_name <- "prereg_recov_n70single"
results_dir  <- file.path("results", session_name)

# --- 2. data collection ---
files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop(paste("no rds files found in:", results_dir))

cat(paste("aggregating", length(files), "recovery iterations\n"))

recovery_list <- list()

for (f in files) {
  # safety check for corrupted files
  fit <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(fit)) next
  
  # extract the comparison table saved in each fit_output bundle
  df <- fit$comparison
  
  # attach parameter bounds for nrmse calculation
  df$min <- as.numeric(fit$param_info$lower[df$Parameter])
  df$max <- as.numeric(fit$param_info$upper[df$Parameter])
  
  # track file source
  df$file_id <- basename(f)
  recovery_list[[f]] <- df
}

all_data <- bind_rows(recovery_list)

# --- 3. statistical summaries ---
# calculate descriptive metrics per parameter/coefficient
stats_summary <- all_data %>%
  group_by(Parameter) %>%
  summarise(
    r = cor(True, Recovered, use = "complete.obs"),
    rmse = sqrt(mean((Recovered - True)^2, na.rm = TRUE)),
    range_width = mean(max - min, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    nrmse_pct = (rmse / range_width) * 100,
    label = paste0("r = ", round(r, 2))
                   #"\nnrmse = ", round(nrmse_pct, 1), "%")
  )

# --- 4. visualization ---
# generate faceted scatterplot comparing true vs recovered values
p <- ggplot(all_data, aes(x = True, y = Recovered)) +
  # identity line representing perfect recovery
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  # individual iteration results
  geom_point(alpha = 0.4, color = "midnightblue") +
  # linear regression trend
  geom_smooth(method = "lm", color = "firebrick3", se = FALSE, size = 0.7) +
  # independent axes per facet due to varying parameter scales
  facet_wrap(~Parameter, scales = "free") +
  # statistical annotations in the top-left of each facet
  geom_text(data = stats_summary, aes(label = label), 
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, 
            size = 3.2, fontface = "bold", inherit.aes = FALSE) +
  theme_minimal(base_size = 12) +
  labs(
    title = "parameter recovery assessment",
    subtitle = paste("directory:", session_name, "| n =", length(unique(all_data$file_id))),
    x = "generative value (ground truth)",
    y = "recovered value (estimate)"
  ) +
  theme(
    strip.background = element_rect(fill = "gray96", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# --- 5. export results ---
# write summary statistics to csv
summary_csv_path <- file.path(results_dir, "_recovery_stats_summary.csv")
write.csv(stats_summary, summary_csv_path, row.names = FALSE)

# save high-resolution plot
plot_pdf_path <- file.path(results_dir, "_recovery_performance_plot.pdf")
ggsave(plot_pdf_path, p, width = 11, height = 8.5)

# print summary table to console
print(p)
cat("\nrecovery analysis complete\n")
print(stats_summary %>% select(Parameter, r, nrmse_pct))
