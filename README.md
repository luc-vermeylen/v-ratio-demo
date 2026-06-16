
# V-Ratio Fitting Pipeline (FCB Model)

A plug-and-play computational modeling pipeline for estimating metacognitive efficiency (**v-ratio**). 

This pipeline fits **v-ratio** using the **Flexible Confidence Boundary (FCB)** models described in:

 *Desender, K., Vermeylen, L., & Verguts, T. (2022). Dynamic influences on static measures of metacognition. Nature communications, 13(1), 4208.*
 
 *Herregods, S., Le Denmat, P., Vermeylen, L., & Desender, K. (2025). Modeling speed–accuracy trade-offs in the stopping rule for confidence judgments. Psychological Review.*

The pipeline uses differential evolution (`DEoptim`) to fit the joint distributions of primary choices, reaction times, confidence judgments, and confidence reaction times using a G-Square cost function. It natively supports multi-condition experimental designs and is optimized to run in batch mode on both local laptops and High-Performance Computing (HPC) clusters.

---

## Required Data Format
Before you begin, your dataset (CSV or RDS) **must** contain the following columns exactly as named:

*   **`sub_id`**: The subject identifier (or change `SUBJECT_COL` in the settings).
*   **`rt`**: Primary decision reaction time (**in seconds**).
*   **`acc`**: Primary decision accuracy (`0` = Error, `1` = Correct).
*   **`rtconf`**: Confidence reaction time (**in seconds**).
*   **`cj`**: Confidence judgment rating. 
    * If using binary confidence, this must be coded as `0` (Low) and `1` (High).
    * If using a 6-point scale, this must be coded as `1` through `6`.

---

## Core Files

*   **`fit_vratio_demo.R`**: The main configuration script. **This is the only file you need to edit.**
*   **`run_batch_local.R`**: A handy script to automatically loop through all subjects and fit them on your local computer.
*   **`batch_fit.slurm`**: Submission script for running massive parallel fits on an HPC cluster (optimized for the KU Leuven HPC)
*   **Post-Fit Analysis scripts**: `fit_assess.R`, `fit_compare.R`, and `fit_stats.R`.
*   *Backend Files (Do not touch)*: `models.cpp` (C++ engine) and `helper_functions.R` (core R logic).

---

## Step-by-Step Workflow

### Step 1: Configure & Fit

Open **`fit_vratio_demo.R`** and adjust the **USER SETTINGS**:

*   **`OUTPUT_FOLDER`**: Give your model a unique folder name (e.g., `"vratio_baseline"` or `"vratio_by_emotion"`).
*   **`DATA_NAME`**: The name of your dataset inside the `data/` folder.
*   **`MODEL_NAME`**: Choose `"FCB_cj2"` (binary confidence) or `"FCB_cj6"` (6-point scale).
*   **`VARYING_PARAMS`**: Define if parameters change between experimental conditions using R formulas. 
    * *Example:* `list(v = ~ as.factor(Difficulty))` will estimate a separate drift rate for each difficulty level. The pipeline will automatically split your data and likelihood blocks safely!

**To Run the Fits:**

*   **Locally (Laptop):** Open `run_batch_local.R` and source it. It will safely loop through every subject in your dataset.
*   **On the HPC:** Navigate to the folder in your terminal and submit the `batch_fit.slurm` using sbatch:

```bash
sbatch --job-name=vratio --array=1-50 --export=NONE,R_SCRIPT=fit_vratio_demo.R batch_fit.slurm
```

More information on how to use this on the HPC is included in the `batch_fit.slurm` file.

### Step 2: Visual Assessment (`fit_assess.R`)

**⚠️ You MUST run this script before comparing models or running statistics!**

Open `fit_assess.R`, set `RESULTS_DIR` to your output folder, and run it.

*   It aggregates all the heavy individual `.rds` fit files into a lightweight master dataset.
*   It generates a PDF report containing:
    1. Decision RT distributions.
    2. Confidence RT distributions.
    3. Confidence mass distributions.
*   It saves `fit_metrics_summary.csv` (which is strictly required for Step 3).

### Step 3: Model Comparison (`fit_compare.R`)

If you fit multiple models (e.g., one where `v` varies by condition, and one where `vratio` varies by condition), this script will tell you which one fits the data best.

*   Open `fit_compare.R` and list the folders you want to compare in `FOLDERS_TO_COMPARE`.
*   It reads the CSVs generated in Step 2 and ranks the models using **BIC** and **BIC Weights**.
*   It outputs a visual comparison plot and a detailed summary table.

### Step 4: Parameter Statistics (`fit_stats.R`)

Once you have selected your winning model, use `fit_stats.R` to analyze the actual parameters.

*   **Global Parameters**: Generates boxplots of shared parameters against model boundaries to check for "boundary swarming" (e.g., ensuring parameters didn't crash into the parameter limits).
*   **Varying Parameters**: Generates Spaghetti plots and Group-Mean plots for any parameters you manipulated. Automatically runs RM-ANOVAs or T-Tests and calculates Effect Sizes.

---

## Important Notes

*   **Time Scaling:** Ensure your `rt` and `rtconf` columns are in **SECONDS** (e.g., `0.450`, not `450`). The integration step is set to `dt = 0.001` to match RTs in seconds.
*   **Subject IDs:** Ensure your `sub_id` column contains unique identifiers for the Repeated-Measures ANOVAs to work correctly in `fit_stats.R`.

## Requirements

The following R packages are required. *(If running on the HPC, ensure your conda environment contains these).*

*   `Rcpp`, `RcppZiggurat` (for the C++ simulation engine)
*   `DEoptim` (for global parameter optimization)
*   `dplyr`, `tidyr`, `data.table` (for data wrangling)
*   `ggplot2`, `patchwork` (for plotting)

## Details of the Cost Function ($G^2$)

The pipeline minimizes the **Likelihood Ratio Chi-Square statistic ($G^2$)**. For categorical and binned multinomial data (such as discretized RT distributions), minimizing $G^2$ is equivalent to Maximum Likelihood Estimation (MLE).

The standard formula for $G^2$ is:

$$
G^2 = 2 \times N \sum_{i=1}^{k} o_i \times \ln\left(\frac{o_i}{p_i}\right)
$$

(Where $N$ is the number of trials, $o_i$ is the observed probability mass in a specific bin, and $p_i$ is the predicted probability mass simulated by the model).

### Evaluating the FCB Model (The 3 Targets)
To capture both the primary decision and the confidence judgement, the pipeline mirrors the fitting architecture introduced in **Herregods et al. (2025)**. It calculates a joint cost function across three specific targets, which are weighted equally:

1. **Target 1: Primary Decision RTs**
   The pipeline calculates dynamic RT quantiles separately for Correct and Error trials. To ensure robust parameter estimation even when certain conditions (or error rates) have very few trials, the pipeline automatically adjusts the binning resolution based on the exact sample size available in that specific cell:
   * **$\ge$ 200 trials:** 10 bins (quantiles: 0.1, 0.2, ..., 0.9)
   * **> 60 trials:** 6 bins (quantiles: 0.1, 0.3, 0.5, 0.7, 0.9)
   * **> 30 trials:** 4 bins (quantiles: 0.3, 0.5, 0.7)
   * **11 to 30 trials:** 2 bins (median split: 0.5)
   * **$\le$ 10 trials (or zero variance):** The pipeline safely abandons RT shape fitting and evaluates the pure marginal probability mass of the cell (e.g., fitting the overall error rate without penalizing RT shape).
   
   It then calculates the $G^2$ difference between the observed and predicted probability mass over these primary decision bins.

2. **Target 2: Confidence RTs**
   The pipeline performs the same dynamic RT quantile binning for the Confidence Reaction Times (`rtconf`), again evaluated separately for correct and error trials.
3. **Target 3: Pure Confidence Proportions**
   Instead of looking at time, this target evaluates the raw choices. It calculates the proportion of trials falling into each confidence rating bin (e.g., ratings 1 through 6), mapped separately for Corrects and Errors.

### Independent Likelihood Blocks for Multiple Conditions
If a user maps parameters to experimental conditions (e.g., a `Difficulty` factor), the pipeline splits the data into **Independent Likelihood Blocks**. It calculates the $G^2$ cost for "Hard" trials independently from "Easy" trials (where the probability mass sums to 1.0 within each specific condition), and then sums the costs together. 

---
