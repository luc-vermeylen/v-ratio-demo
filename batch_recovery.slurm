#!/bin/bash -l
#SBATCH --job-name=DMC_Recovery
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2g
#SBATCH --output=slurm_oe/recovery-%A_%a.o
#SBATCH --error=slurm_oe/recovery-%A_%a.e
#SBATCH --account=lp_desenderlab
#SBATCH --cluster=wice

# ==============================================================================
# RECOVERY BATCH RUNNER (v2.0)
# ==============================================================================
: <<'HELP'
1. PURPOSE:
   This script runs multiple independent parameter recovery iterations.
   Each run samples a new random Ground Truth and attempts to recover it.

2. HOW TO RUN:
   To run 100 independent recovery iterations:
   sbatch --array=1-100 --export=NONE,R_SCRIPT=recovery_template.R batch_recovery.slurm

3. PARAMETERS:
   - --array=1-100: Determines the number of data points in your final 
     recovery scatterplot. 50-100 is recommended for a stable R-value.
   - R_SCRIPT: Your recovery script (e.g., recovery_template.R).

4. OUTPUT:
   Results will be saved in the results/ folder specified in your R script.
   Filenames will include the Task ID and a random ID to prevent collisions.
HELP

# ==============================================================================
# CONFIGURATION
# ==============================================================================
CONDA_ENV_NAME="r_env"

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================

# 1. Ensure an R script was specified
if [ -z "$R_SCRIPT" ]; then
    echo "ERROR: No R script provided. Use --export=NONE,R_SCRIPT=your_script.R"
    exit 1
fi

# 2. Setup Directories
cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p slurm_oe

# 3. Initialize Conda
module --force purge
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

conda activate "$CONDA_ENV_NAME" || { echo "ERROR: Could not activate $CONDA_ENV_NAME"; exit 1; }

# ==============================================================================
# 1. JOB INFO / LOGGING
# ==============================================================================
echo "=============================================="
echo "RECOVERY Iteration started at: $(date)"
echo "Iteration ID:   ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_COUNT}"
echo "R script:       $R_SCRIPT"
echo "Working dir:    $SLURM_SUBMIT_DIR"
echo "Cores used:     $SLURM_CPUS_PER_TASK"
echo "=============================================="

# ==============================================================================
# 2. RUN THE RECOVERY
# ==============================================================================
# We pass the Task ID to R. In the recovery script, this serves as 
# the SUBJECT_ID_VAL for unique filename generation.
Rscript --vanilla "$R_SCRIPT" ${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "RECOVERY Iteration finished at: $(date)"
echo "=============================================="