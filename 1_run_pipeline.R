# ==============================================================================
# V-RATIO FITTING PIPELINE 
# ==============================================================================
# Directions:
# 1. Edit your settings in Section 1.
# 2. Choose your RUN_MODE in Section 2.
# 3. Run this file
# ==============================================================================

# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================
OUTPUT_FOLDER <- "vratio" # the folder where fits will be stored (the model identifier)
DATA_NAME     <- "Exp2_cj6_Herregods2025.csv" # Must be in the /data folder
SUBJECT_COL   <- "sub_id"

MODEL_NAME    <- "FCB_cj6" # FCB_cj2 or FCB_cj6 depending on if you have binary or 6 levels of confidence

VARYING_PARAMS <- list() # e.g., list(v = ~ as.factor(Difficulty))
FIXED_PARAMS   <- list() # e.g., list(starting_point_confidence = 0.5)

ITER_MAX  <- 10 # 1000 recommended
USE_CORES <- 1  # 1 for parallel processing within DEoptim, 0 for turning this off

# ==============================================================================
# 2. EXECUTION MODE
# ==============================================================================
# Choose one of the following:
# "single"      -> Fits ONLY the specific subject index defined below.
# "group"       -> Fits all data merged at the group level.
# "local_batch" -> Fits ALL subjects sequentially on this computer.
# "hpc"         -> Saves settings and prints batch command for running on HPC.

RUN_MODE <- "local_batch"
SUBJECT_IDX_TO_FIT <- 1  # Only used if RUN_MODE is "single"

# ==============================================================================
# 3. BACKGROUND EXECUTION (Do not edit below)
# ==============================================================================
library(here)

# Define absolute paths dynamically
OUTPUT_DIR  <- here("results", OUTPUT_FOLDER)
CONFIG_PATH <- here("results", OUTPUT_FOLDER, "config.RData")
DATA_FILE   <- here("data", DATA_NAME)
ENGINE_FILE <- here("R","fit.R")

# Create results folder and save config
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat(sprintf("\nSaving configuration to %s...\n", CONFIG_PATH))
save(list = c("OUTPUT_FOLDER", "DATA_NAME", "SUBJECT_COL", "MODEL_NAME", 
              "VARYING_PARAMS", "FIXED_PARAMS", "ITER_MAX", "USE_CORES"), 
     file = CONFIG_PATH)

# Load data to count subjects
if (!file.exists(DATA_FILE)) stop("File not found! Did you put your data in the 'data' folder?")

if (grepl("\\.rds$", DATA_FILE, ignore.case = TRUE)) {
  raw_data <- readRDS(DATA_FILE)
} else {
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  raw_data <- as.data.frame(data.table::fread(DATA_FILE))
}
n_subjects <- length(unique(raw_data[[SUBJECT_COL]]))

# Execute based on mode
if (RUN_MODE == "single") {
  cat(sprintf("\n=== FITTING SINGLE SUBJECT (Index %d) ===\n", SUBJECT_IDX_TO_FIT))
  system2("Rscript", args = c("--vanilla", shQuote(ENGINE_FILE), SUBJECT_IDX_TO_FIT, shQuote(CONFIG_PATH)))
  
} else if (RUN_MODE == "local_batch") {
  cat("\n=== STARTING LOCAL BATCH (", n_subjects, "SUBJECTS ) ===\n")
  start_time <- Sys.time()
  for (i in 1:n_subjects) {
    cat(sprintf("\n--- FITTING SUBJECT %d OF %d ---\n", i, n_subjects))
    system2("Rscript", args = c("--vanilla", shQuote(ENGINE_FILE), i, shQuote(CONFIG_PATH)))
  }
  cat("\nALL FITS COMPLETED in", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes.\n")
  
} else if (RUN_MODE == "group") {
  cat("\n=== FITTING GROUP (All subjects pooled) ===\n")
  system2("Rscript", args = c("--vanilla", shQuote(ENGINE_FILE), "GROUP", shQuote(CONFIG_PATH)))
  
} else if (RUN_MODE == "hpc") {
  # For the HPC printout, use a relative path so it works when copy-pasted onto Linux!
  CONFIG_REL <- file.path("results", OUTPUT_FOLDER, "config.RData")
  
  cat("\n=== HPC PREPARATION COMPLETE ===\n")
  cat(sprintf("1. Your settings are safely saved in:\n   %s\n", CONFIG_PATH))
  cat("2. To run this on the HPC cluster, simply copy and paste this command into your terminal:\n\n")
  cat(sprintf("sbatch --job-name=vratio --array=1-%d --export=NONE,R_SCRIPT=R/fit.R,CONFIG_PATH=%s hpc/batch_fit.slurm\n\n", n_subjects, CONFIG_REL))
  
} else {
  stop("Invalid RUN_MODE! Choose 'single', 'local_batch', 'group', or 'hpc'.")
}