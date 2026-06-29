#!/bin/bash -l
#SBATCH --job-name=DDM
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2g
#SBATCH --output=slurm_oe/slurm-%A_%a.o
#SBATCH --error=slurm_oe/slurm-%A_%a.e
#SBATCH --account=lp_desenderlab
#SBATCH --cluster=wice

# ==============================================================================
# HOW TO USE THIS BATCH SCRIPT ON THE KU LEUVEN HPC
# ==============================================================================
: <<'HELP'
1. SETUP YOUR ENVIRONMENT (Run once on the login node):
   conda create -n r_env r-base=4.3.1 r-essentials -c conda-forge
   conda activate r_env
   R -e "install.packages(c('here','Rcpp','RcppZiggurat','DEoptim','dplyr','tidyr','data.table','ggplot2','patchwork'), repos='https://cloud.r-project.org')"

2. LOCATIONS:
- You must submit this script from the ROOT folder of this repository.
- Slurm logs (errors/output) will appear in a subfolder: slurm_oe/
- Model results (.rds/.pdf) will appear in the results/ folder.

3. THE "EXPORT=NONE" TRICK:
We use --export=NONE to prevent "user env retrieval" errors. When running 
large arrays, Slurm can crash trying to copy your terminal settings. 
NONE ensures a clean, stable startup for every job.

4. SUBMIT JOBS:
You do not need to type the sbatch command manually! 
Instead, configure your model in `1_run_pipeline.R` and set RUN_MODE <- "hpc". 
Running `1_run_pipeline.R` will automatically print the exact `sbatch` command 
you need to copy and paste into the terminal.

Example of an auto-generated command:
sbatch --job-name=vratio --array=1-51 --export=NONE,R_SCRIPT=R/fit.R,CONFIG_PATH=results/vratio/config.RData hpc/batch_fit.slurm
HELP

# ==============================================================================
# CONFIGURATION
# ==============================================================================
CONDA_ENV_NAME="r_env"

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================

# 1. Ensure an R script and Config Path were specified
if [ -z "$R_SCRIPT" ] || [ -z "$CONFIG_PATH" ]; then
    echo "ERROR: Missing arguments. Use --export=NONE,R_SCRIPT=R/fit.R,CONFIG_PATH=results/.../config.RData"
    exit 1
fi

# 2. Set the working directory to where you submitted the job (Root of Repo)
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
echo "SLURM job started at: $(date)"
echo "User:           $USER"
echo "Job ID:         ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "R script:       $R_SCRIPT"
echo "Config path:    $CONFIG_PATH"
echo "Working dir:    $SLURM_SUBMIT_DIR"
echo "Using R from:   $(which R)"
echo "=============================================="

# ==============================================================================
# 2. RUN THE R SCRIPT
# ==============================================================================
Rscript --vanilla "$R_SCRIPT" ${SLURM_ARRAY_TASK_ID} "$CONFIG_PATH"

echo "=============================================="
echo "SLURM job finished at: $(date)"
echo "=============================================="