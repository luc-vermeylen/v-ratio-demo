***

# DMC Fitting Pipeline

A pipeline for fitting and analyzing DMC models. It allows you to specify your own models in C++ and handles multi-condition designs, custom cost functions and is optimized for the HPC cluster.

## Required Files
*   **`models.cpp`**: The C++ simulation engine. Where you can define new models.
*   **`helper_functions.R`**: Core logic for parameter mapping, fitting and plotting. Don't touch, except for defining the parameters of a newly added model!
*   **`fit_template.R`**: The main script where you define the model and fit settings. Take a copy and give it a name that relates to your model, don't use the original.
*   **`batch_fit.slurm`**: Submission script for the HPC.
*   **Post-Fit Analysis scripts**: `fit_assess.R`, `fit_stats.R`, and `fit_compare.R`.

## 1. Configuring a Fit (`fit_template.R`)
Open `fit_template.R` and adjust the **USER SETTINGS**:
*   **`OUTPUT_DIR`**: Set a unique folder name for each model variant (e.g., `"full_model"`, `"fixed_amp"`,  `"v_by_condition"`). Best to also take a copy of fit_template.R and name it like the **`OUTPUT_DIR`** folder for reference later and keeping the template itself clean.
*   **`VARYING_PARAMS`**: Define which parameters vary by which column(s) in your CSV. Use the format `parameter = ~ factor1 * factor2`.
    *   *Example:* `list(amp = ~ distance * meta_bin` fits a full interaction for the amp parameter.
*   **`FIT_TARGETS`**: Ensure all columns used in `VARYING_PARAMS` are also included in `split_cols`. This is critical for correct Likelihood calculation.
*   **`FIXED_PARAMS`**: Use this to hold specific parameters constant (e.g., `tau = 60`).

## 2. Running on the HPC (`batch_fit.slurm`)

Navigate to the folder on the HPC where you store these scripts, including the batch_fit.slurm.

**Single model submission:**
```bash
sbatch --job-name=my_model --array=1-51 --export=NONE,R_SCRIPT=fit_template.R,WORKDIR=$VSC_DATA/my_project batch_fit.slurm
```

**Multiple models submission at once:**
```bash
for script in fit_variant_*.R; do
  sbatch --job-name=$script --array=1-51 --export=NONE,R_SCRIPT=$script,WORKDIR=$VSC_DATA/my_project batch_fit.slurm
done
```

Don't forget to set up your conda environment with the necessary packages. See batch_fit.slurm for more information on how to. 

## 3. Analysis Workflow

### Step 1: Visual Assessment (`fit_assess.R`)
Run this first to check if the model actually describes the behavior.
*   Aggregates `.rds` files into a master table.
*   Generates "Mirror" RT distributions, CAFs, and Delta plots.
*   Saves `fit_metrics_summary.csv` (required for further stats and model comparison).

### Step 2: Parameter Statistics (`fit_stats.R`)
Analyze the estimated parameters across the group.
*   **Varying Parameters**: Generates Spaghetti plots and clean Group-Mean plots. Automatically runs RM-ANOVAs and calculates Effect Sizes (PES).
*   **Global Parameters**: Generates boxplots of shared parameters against model boundaries to check for "boundary swarming."
*   **Derived Parameters**: Add custom formulas (e.g., `amp / a`) in the `DERIVED_PARAMS` list.

### Step 3: Model Comparison (`fit_compare.R`)
Compare different model folders to find the best fit.
*   Ranks models using **BIC** and **BIC Weights**.
*   Includes a **Validity Check** to ensure all models were fitted on the same bin structures, which is necessary to validly compare G2/BIC.

## Important Notes
*   **BIC Validity**: Only compare BICs if the models used the same `FIT_TARGETS` splits.
*   **Subject IDs**: Ensure IDs are unique in your data for the RM-ANOVA to work.

  ## Requirements
The following R packages are required:
* `Rcpp`, `RcppZiggurat` (for simulation)
* `DEoptim` (for fitting)
* `dplyr`, `tidyr`, `data.table` (for data handling)
* `ggplot2`, `patchwork` (for plotting)
