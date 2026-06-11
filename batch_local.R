# ==============================================================================
# LOCAL BATCH RUNNER
# ==============================================================================
# This script automatically determines the number of unique subjects in your 
# dataset and runs the fitting script for each subject sequentially on your computer.
# make sure to also adjust the "USER SETTINGS" section in the fit_vratio_demo.R script
# ==============================================================================

# --- 1. SETTINGS ---
DATA_FILE     <- "data/Exp2_cj6_Herregods2025.csv" # Ensure this matches your data file
SUBJECT_COL   <- "sub_id"
SCRIPT_TO_RUN <- "fit_vratio_demo.R"

# --- 2. COUNT UNIQUE SUBJECTS ---
cat("Checking dataset:", DATA_FILE, "...\n")

if (!file.exists(DATA_FILE)) {
  stop("File not found! Please check the path and name.")
}

# Load the data to count subjects
if (grepl("\\.rds$", DATA_FILE, ignore.case = TRUE)) {
  raw_data <- readRDS(DATA_FILE)
} else {
  raw_data <- read.csv(DATA_FILE)
}

if (!SUBJECT_COL %in% names(raw_data)) {
  stop(paste("Column", SUBJECT_COL, "not found in the dataset!"))
}

n_subjects <- length(unique(raw_data[[SUBJECT_COL]]))
cat("Found", n_subjects, "unique subjects.\n")

# --- 3. RUN FITS SEQUENTIALLY ---
cat("Starting local batch fit...\n")
start_time_global <- Sys.time()

for (i in 1:n_subjects) {
  cat("\n=================================================================\n")
  cat(sprintf(" FITTING SUBJECT %d OF %d\n", i, n_subjects))
  cat("=================================================================\n")
  
  exit_status <- system2("Rscript", args = c("--vanilla", SCRIPT_TO_RUN, i))
  
  if (exit_status != 0) {
    warning(sprintf("Subject %d encountered an error and did not finish properly.", i))
  }
}

duration_global <- Sys.time() - start_time_global
cat("\n=================================================================\n")
cat("ALL FITS COMPLETED!\n")
cat("Total time taken:", round(duration_global, 2), units(duration_global), "\n")
cat("=================================================================\n")