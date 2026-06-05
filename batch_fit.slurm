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
# HOW TO USE THIS PIPELINE
# ==============================================================================
: <<'HELP'
1. SETUP YOUR ENVIRONMENT (Run once on the login node):
   conda create -n r_env r-base=4.3.1 r-essentials -c conda-forge
   conda activate r_env
   R -e "install.packages(c('Rcpp','RcppZiggurat','DEoptim','dplyr','tidyr','data.table','ggplot2','patchwork'), repos='https://cloud.r-project.org')"
2. LOCATIONS:
- Navigate to your project folder, where you have the R scripts, data and this .slurm file
- In that project folder, run the sbatch commands explained below.
- Slurm logs (errors/output) will appear in a subfolder: slurm_oe/
- Model results (.rds/.pdf) will appear in a subfolder you named INSIDE your R script.
- Tip: name results folder the same as the fit script

3. THE "EXPORT=NONE" TRICK:
We use --export=NONE to prevent "user env retrieval" errors. When running 
large arrays, Slurm can crash trying to copy your terminal settings. 
NONE ensures a clean, stable startup for every job.

4. SUBMIT JOBS:
Copy and paste these into your terminal from in the project folder
Change your array to the subject index range and R_SCRIPT to your fit file.

--- To run a single model ---
sbatch --job-name=my_fit --array=1-51 --export=NONE,R_SCRIPT=model1.R batch_fit.slurm

--- Run multiple models ---
for script in model1.R model2.R model3.R; do
  sbatch --job-name=$script --array=1-51 --export=NONE,R_SCRIPT=$script batch_fit.slurm
done

--- To run all scripts in your folder matching a pattern ---
for script in model-*.R; do
  sbatch --job-name=$script --array=1-51 --export=NONE,R_SCRIPT=$script batch_fit.slurm
done
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

# 2. Set the working directory to where you submitted the job
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
echo "Working dir:    $SLURM_SUBMIT_DIR"
echo "Using R from:   $(which R)"
echo "=============================================="

# ==============================================================================
# 2. RUN THE R SCRIPT
# ==============================================================================
Rscript --vanilla "$R_SCRIPT" ${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "SLURM job finished at: $(date)"
echo "=============================================="