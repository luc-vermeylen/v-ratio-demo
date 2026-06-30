
------------------------------------------------------------------------

# V-Ratio Fitting Pipeline (FCB Model)

A computational modeling pipeline for estimating metacognitive efficiency (**v-ratio**).

This pipeline fits the **Flexible Confidence Boundary (FCB)** models described in: 

> *Herregods, S., Le Denmat, P., Vermeylen, L., & Desender, K. (2025). Modeling speed–accuracy trade-offs in the stopping rule for confidence judgments. Psychological Review.*

> *Desender, K., Vermeylen, L., & Verguts, T. (2022). Dynamic influences on static measures of metacognition. Nature communications, 13(1), 4208.*

To help you get started, **we have included the original datasets from Herregods et al. (2025)** in the `data/` folder.

------------------------------------------------------------------------

## Repository Overview

- 📁 **`data/`**: Place your dataset (.csv or .RDS) here.
- 📁 **`results/`**: Each fit will get a new folder here where configuration (config.Rdata), fit information (.rds), plots (.pdf) and summaries (.csv) will appear.
- 📁 **`R/`**: Backend C++ simulation code and R helper functions. *(Do not edit these).*
- 📁 **`hpc/`**: Slurm scripts (made for the KU Leuven HPC infrastructure). Detailed information on how to run this on the hpc can be found in below and within `hpc/batch_fit.slurm`.
- 📄 **`1_run_pipeline.R`**: The main file for running a fit. Configure settings and launch fits here.
- 📄 **`2_fit_assess.R`**: Aggregates fits and generates behavioral diagnostic plots.
- 📄 **`3_fit_compare.R`**: Compares different models to find the best fit via BIC.
- 📄 **`4_fit_stats.R`**: Runs statistics (T-Tests, ANOVAs) and generates parameter plots.

------------------------------------------------------------------------

## Required Data Format

Your dataset **must** contain the following columns exactly as named: 

* **`sub_id`**: The subject identifier (or change `SUBJECT_COL` in the settings).
* **`rt`**: Primary decision reaction time (**strictly in seconds**, e.g., `0.450`). 
* **`acc`**: Primary decision accuracy (`0` = Error, `1` = Correct). 
* **`rtconf`**: Confidence reaction time (**strictly in seconds**). 
* **`cj`**: Confidence judgment rating. If using binary confidence (`FCB_cj2`), code this as `0` (Low) and `1` (High). If using a 6-point scale (`FCB_cj6`), code this as `1` through `6`.

------------------------------------------------------------------------

## Step-by-Step Workflow

### Step 1: Configure & Fit (`1_run_pipeline.R`)

Open Rstudio using the `.Rproj` (to set the paths correctly), open `1_run_pipeline.R` and adjust the User Settings. You will specify your dataset, the `OUTPUT_FOLDER` for the results, your model (`FCB_cj2` or `FCB_cj6`, depending on your number of confidence response options), which parameters to vary by experimental conditions and which parameters to fix. When the file is ready, save and simply run the file by using ctrl+a and the run button or press the source button in RStudio or type ```source('1_run_pipeline.R')``` in the console.

> **The Importance of the `OUTPUT_FOLDER`**
> The `OUTPUT_FOLDER` (which will show up within the `/results/` folder) is the most critical organizational setting in the entire pipeline. It acts as the permanent "home" for your specific model configuration. When the pipeline runs, every `.rds` file, `.pdf` report, and `.csv` summary is saved inside `results/OUTPUT_FOLDER/`. Therefore, the approach is to use one folder for one fit. You can create different folders with different types of fit (e.g., conditions varying or different parameters fixed) and then later compare those using the `3_fit_compare.R` script.
> Crucially, **all subsequent analysis scripts (Steps 2, 3, and 4) rely entirely on this exact folder name** to locate your data. The idea is to give it a descriptive, unique name representing the current model's hypothesis (e.g., `"vratio_by_emotion"` or `"static_boundaries"`). 

**Specify the `CONDITIONS` and `VARYING_PARAMS`:** If you wish to include experimental conditions, you can do so here. In order to guarantee fair model comparison (see more details below), simply follow this rule:

1. **`CONDITIONS`**: Look at the *most complex* model you plan to test. What experimental factors does it use? List them here as a vector of strings (e.g., `c("Difficulty", "Emotion")`). **You must use this exact same `CONDITIONS` list for every single model you fit in your comparison set**, even for the Null model. This controls the splitting of the Likelihood.
2. **`VARYING_PARAMS`**: Which parameters are actually allowed to change across those conditions in *this specific fit*? Use R formulas here (e.g., `list(v = ~ as.factor(Difficulty))` or `list(vratio = ~ as.factor(Emotion) * as.factor(Validity))`). To fit a Null model, simply leave this empty as `list()`. This controls the splitting of the Parameters.

> **Quick Rationale:** We separate *splitting the likelihood* (`CONDITIONS`) from *splitting the parameters* (`VARYING_PARAMS`) in order to be able to compare complex models to complete Null models. If your Null model doesn't use the exact same `CONDITIONS` as your Complex model, they are "graded" on different ways of splitting the data, and therefore, their BICs should not be compared! See the **Good Model Comparison Practices** section below for a detailed explanation.

**Specify the `FIXED_PARAMS`:** Provide a list with parameter = value to fix the parameter to that value.

* `list()`: (Empty list). All parameters are free to vary.
* `list(starting_point_confidence = 0.5)`. starting point for confidence is fixed to 0.5.

**Select a `RUN_MODE`'s:** Change the `RUN_MODE` variable to choose how to run the pipeline: 

* `"single"`: Runs a test fit on 1 subject so you can check for errors. 
* `"group"`: Pools all data together to fit one massive "mega-subject" (Group Fit).
* `"local_batch"`: Automatically loops through all subjects and fits them on your local computer.
* `"hpc"`: Saves your configuration file and created the fit folder (`config.RData`) and prints the exact `sbatch` command you need to copy/paste into the cluster terminal to run a massive parallel array job.

> **What is saved?** For each subject, a `.pdf` file with diagnostic plots and a `.rds` file is saved in `results/YOUR_FOLDER_NAME`. Here is a detailed breakdown of its contents:
>*   **`best_params`**: A named vector of the winning parameter values. If you used experimental conditions (formulas), the pipeline automatically translates the regression coefficients back into readable "marginal cell means" (e.g., `v:Hard = 1.2`, `v:Easy = 2.1`) so they are ready for easy plotting, interpretation and ANOVAs.
>*   **`best_betas`**: The raw regression coefficients (Intercepts and Slopes) that the optimizer actually found. If you did not use formulas in `VARYING_PARAMS`, this will be `NULL`.
>*   **`best_cost`**: The final raw $G^2$ value minimized by the optimizer.
>*   **`bic` & `aic`**: The Bayesian and Akaike Information Criterions, strictly penalized for the number of free parameters and the amount of data.
>*   **`n_bins`**: The exact number of probability mass bins the data was split into for the cost calculation.
>*   **`n_free_params` & `n_observations`**: The complexity of the model and the sample size.
>*   **`observations`**: A clean, subsetted dataframe containing the exact empirical trial data used to fit this specific subject.
>*   **`predictions`**: A high-resolution simulated dataset (by default, 10,000 trials per condition) generated by the C++ engine using the `best_params`. This is extremely useful if you want to create your own custom behavioral plots!
>*   **`final_proportions`**: A nested list containing the exact probability mass distributions (observed vs. predicted) across all dynamic RT and Confidence bins. These are the raw binning vectors used for the $G^2$ calculation.
>*   **`info`**: A metadata list recording the `model` name, the exact `timestamp` the fit completed, the `subject` ID, and the exact `VARYING_PARAMS` formulas requested by the user.
>*   **`param_info`**: The boundary configurations passed to DEoptim (the absolute lower and upper search limits for each parameter).
>*   **`constants`**: All fixed mathematical constants used during the simulation (e.g., `dt`, `s`, `ntrials`, and any values set in `FIXED_PARAMS`).
>*   **`targets` & `cost_method`**: A record of how the likelihood blocks and variables were mapped (the `FIT_TARGETS` list).
>*   **`optim_full`**: The massive, raw output object returned directly by the `DEoptim` package, which contains the trace of all generations (useful for advanced optimization diagnostics).

> **Running Multiple Fit Attempts**
> * If you suspect a model didn't fit well, **you can simply run the exact same fit again in the same folder**. 
> * The pipeline uses precise timestamps in the `.rds` filenames, meaning old fits are never accidentally overwritten.
> * When you run Step 2 (`2_fit_assess.R`), the script automatically groups all files by subject, **identifies the fit attempt with the lowest $G^2$ cost**, and discards the rest. The pipeline automatically handles the cleanup and guarantees you only analyze the absolute best fits!

### Step 2: Visual Assessment (`2_fit_assess.R`)

**⚠️ You MUST run this script before comparing models or running statistics!** Make sure to change OUTPUT_FOLDER to match the exact folder you created in Step 1! This script will then read all the `.rds` files in that folder, aggregate them into a lightweight CSV, and provide plots that help you evaluate whether the model actually captured human behavior adequately.

**Outputs:** 
* **`fit_metrics_summary.csv`**: A spreadsheet of BICs.
* **`master_parameter_report.csv`**: A spreadsheet with all parameters.
* **`grand_average_fit_assessment.pdf`**: A visual report containing: 
    * *Decision RT Distributions:* A mirrored density plot of observed vs predicted RTs.
    * *Confidence RT Distributions:* The speed of the confidence judgment.
    * *Confidence Rating Mass:* A barplot showing if the model predicted the exact empirical frequencies of ratings (e.g., 1-6).

### Step 3: Model Comparison (`3_fit_compare.R`)

If you fit multiple models to test competing hypotheses (e.g., one folder where `v` varies, and another where `vratio` varies), this script ranks them. Open it and list your result folders in `FOLDERS_TO_COMPARE`. Note: You can only validly compare models that were fitted on the exact same observations, the same number of subjects, and the same binning structure (for this reason, n_bins used for the the $G^2$ calculation is stored and printed in the model comparison table).

**Outputs:** 
* `model_comparison_plots.pdf`: Shows boxplots of sum BICs, Subject "Wins", and the Evidence Gap (delta-BIC) by model. 
* `summary_model_comparison_detailed.csv`: A table ranking models by sum BIC, dBIC, Akaike Weights and Ratio, Raw cost, model complexity penalty, number of parameters, observations and bins used for the $G^2$ calculation.

### Step 4: Parameter Statistics (`4_fit_stats.R`)

Once you have selected your winning model, this script automatically extracts and analyzes the estimated parameters.

**Outputs:** 
* **Global Parameters Plot**: Generates boxplots of shared ("global") parameters. You can optionally set `SHOW_BOUNDARIES <- TRUE` to plot them against the DEoptim optimization boundaries. This is useful for checking for "boundary swarming", ensuring the optimizer didn't get trapped against the parameter boundaries (otherwise this suggests that a larger parameter space is warranted). 
* **Individual Trends & Group Means**: Generates Spaghetti plots and summary plots for any `VARYING_PARAMS`.
* **`parameter_anova_results.csv`**: Automatically runs T-Tests or Repeated-Measures ANOVAs on your varying parameters and saves the T/F-values, p-values, and effect sizes.

------------------------------------------------------------------------

## Running on an HPC Cluster (KU Leuven / VSC)

If you are running this pipeline on a supercomputer using Slurm (e.g., the WICE cluster at KU Leuven), the repository includes a dedicated `hpc/batch_fit.slurm` script. 

### 1. Initial Environment Setup (Run Once)
Before running the pipeline for the first time, you must log into the login node and set up an R environment containing all the necessary packages. Copy and paste these commands into your terminal:

```bash
# Create a fresh conda environment with base R
conda create -n r_env r-base=4.3.1 r-essentials -c conda-forge

# Activate the environment
conda activate r_env

# Install the required R packages
R -e "install.packages(c('here','Rcpp','RcppZiggurat','DEoptim','dplyr','tidyr','data.table','ggplot2','patchwork'), repos='https://cloud.r-project.org')"
```

### 2. Submitting Jobs
You **do not** need to type out the `sbatch` command manually. The pipeline will write it for you! 
1. Go to your data folder on the cluster, e.g., ```cd $VSC_DATA```
2. Clone the repository here: ```git clone https://github.com/luc-vermeylen/v-ratio-demo.git```
3. Move into the newly downloaded folder: ```cd v-ratio-demo``` (*note: If you are using your own dataset, make sure to upload your CSV file into the data/ folder now*)
4. Activate your Conda environment: ```conda activate r_env```
5. open and configure your model: ```nano 1_run_pipeline.R``` *(note: change your settings, ensure RUN_MODE <- "hpc", then press Ctrl+O and Enter to save, and Ctrl+X to exit)*
6. Run the script on the login node: ```Rscript 1_run_pipeline.R```
7. The script will create the new results folder, save your configuration there (`config.Rdata`) and print the exact `sbatch` command you need. Simply copy and paste that command into your terminal.

### 3. Understanding the `sbatch` Command (for more flexibility)
If you want to construct the command yourself, you must execute it from the **root folder** of the repository. Here is the anatomy of the submission command:

```bash
sbatch --job-name=vratio --array=1-51 --export=NONE,R_SCRIPT=R/fit.R,CONFIG_PATH=results/vratio/config.RData hpc/batch_fit.slurm
```

**What these arguments do:**
*   `--job-name=vratio`: The name of your job as it will appear in the Slurm queue (`squeue`).
*   `--array=1-51`: Tells Slurm to launch 51 parallel nodes (one for each subject). The script will automatically assign Subject 1 to Node 1, Subject 2 to Node 2, etc.
*   `--export=NONE`: This prevents Slurm from copying your current terminal's environment variables to the compute nodes. When running large arrays, omitting this can cause the cluster to crash. `NONE` ensures a clean, stable startup for every job.
*   `R_SCRIPT=R/fit.R`: Points the Slurm file to the fitting engine file.
*   `CONFIG_PATH=...`: Points the engine to the specific settings file you created for this model fit.
*   `hpc/batch_fit.slurm`: The actual Slurm script containing the hardware requests (memory, CPUs, time limits).

### 4. Advanced: Submitting Multiple Models at Once
You can easily queue up multiple different models without duplicating any actions. 

**Step 1:** Open `1_run_pipeline.R`, set your first model's settings (e.g., `OUTPUT_FOLDER <- "model_A"`), and run it in `"hpc"` mode to generate its `config.RData`. Repeat this for `"model_B"`, `"model_C"`, etc.

**Step 2:** Use a bash `for` loop in your terminal to submit them all simultaneously.

**Example A: Submit a specific list of models**
```bash
for folder in results/model_A results/model_B results/model_C; do
  # Extracts the folder name to use as the job name
  job_name=$(basename $folder)
  
  sbatch --job-name=$job_name --array=1-51 --export=NONE,R_SCRIPT=R/fit.R,CONFIG_PATH=$folder/config.RData hpc/batch_fit.slurm
done
```

**Example B: Submit ALL models that match a pattern**
If you created 10 models and named them all starting with `vratio_`, you can queue all 10 arrays (510 jobs) with one loop:
```bash
for folder in results/vratio_*; do
  job_name=$(basename $folder)
  sbatch --job-name=$job_name --array=1-51 --export=NONE,R_SCRIPT=R/fit.R,CONFIG_PATH=$folder/config.RData hpc/batch_fit.slurm
done
```

### 5. Tracking HPC Output
*   **Results:** As jobs finish, `.rds` files will quietly appear inside your designated `results/OUTPUT_FOLDER/`.
*   **Error Logs:** Slurm will generate text logs for every subject (e.g., `slurm-12345_1.o` and `slurm-12345_1.e`) in a folder named **`slurm_oe/`**. If a subject crashes, check the `.e` (error) file in this directory to see what went wrong.
  
------------------------------------------------------------------------

## The FCB Model Parameters Dictionary

The models rely on the following parameters. (Note: `a_slope` is fixed to 0 by default in the helper functions). Parameter boundaries can be changed in R/helper_functions.R, in the model_params list.

| Parameter | Stage | Description |
|:-----------------------|:-----------------------|:-----------------------|
| **`a`** | Decision | **Decision Boundary**: The total evidence required to make the primary choice. |
| **`v`** | Decision | **Drift Rate**: The speed of information accumulation for the primary choice. |
| **`ter`** | Decision | **Non-Decision Time**: Encoding and motor execution time (seconds). |
| **`a_slope`** | Decision | **Decision Boundary Slope**: How fast the decision boundaries collapse. |
| **`a2`** | Confidence | **Confidence Boundary**: The total evidence space for the confidence judgment. |
| **`vratio`** | Confidence | **V-Ratio**: The ratio of post-decisional drift to primary drift ($v_{post} = v \times vratio$). |
| **`ter2`** | Confidence | **Confidence Non-Decision Time**: Delay specific to confidence reporting. Can be negative, implying metacognitive processing starts during the primary motor execution. |
| **`a2_slope_upper`** | Confidence | **Upper Bound Collapse**: How fast the upper confidence boundary collapses. |
| **`a2_slope_lower`** | Confidence | **Lower Bound Collapse**: How fast the lower confidence boundary collapses. |
| **`starting_point_confidence`** | Confidence | **Starting Point**: The fraction of `a2` where confidence accumulation begins (between 0 and 1). |

------------------------------------------------------------------------

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

------------------------------------------------------------------------

## Good Model Comparison Practices: Decoupling Likelihood from Parameters

When conducting model comparisons, we typically evaluate a **Null Model** (where a cognitive process is invariant across experimental conditions) against a more **Complex Model** (where a parameter, such as drift rate, varies across conditions). 

Within this pipeline, specifying these models requires two distinct, independent arguments in `1_run_pipeline.R`:
*   **`CONDITIONS`**: Specifies the experimental factors used to partition the empirical data into independent likelihood blocks prior to the cost function calculation.
*   **`VARYING_PARAMS`**: Specifies which computational parameters are permitted to freely vary across those predefined blocks.

This decoupling is a statistical requirement to ensure that the BIC remain comparable across models.

### The Problem: Incomparable Likelihood Spaces
To validly compare two models using BIC or the Likelihood Ratio Chi-Square ($G^2$), both models must be evaluated over the exact same data space and binning architecture. 

If data partitioning is only applied when a parameter varies, the Null model is evaluated against the *marginal* distribution (pooled data), while the Complex model is evaluated against the *conditional* distributions (data split by experimental manipulation). This causes two fatal issues:
1. **Probability Normalization:** Matching a single parameter to a pooled distribution is a mathematically distinct objective compared to matching parameters to multiple, separately-normalized conditional distributions (where each block's probability mass sums to 1.0 locally). 
2. **Degrees of Freedom in Binning:** Because the pipeline utilizes dynamic, sample-size-dependent adaptive binning, splitting the data inherently reduces the trial count per cell. This alters the discrete binning resolution (e.g., a pooled dataset might use 10 bins, while split conditions might drop to 4 bins). 

Comparing a model evaluated on pooled data to a model evaluated on split data violates the foundational assumptions of model selection. The resulting $G^2$ and BIC values would be evaluated on entirely different scales, making them incomparable.

### The Solution: Constant Likelihood Spaces, Flexible Parameters
By decoupling data partitioning from parameter flexibility, the pipeline ensures valid model selection:
*   `CONDITIONS` defines the fixed evaluation architecture. It partitions the data into independent conditional likelihood blocks regardless of the model's complexity.
*   `VARYING_PARAMS` dictates the flexibility the model has to accommodate those blocks.

Under this framework, **the Null model is evaluated condition-by-condition**. However, because its parameters are tied globally, it predicts identical probability distributions for every block. If the empirical data structurally differ across conditions, the invariant parameters of the Null model will fail to capture the local variance, appropriately yielding a high $G^2$ deviance penalty. This allows the Complex model to rightfully win the BIC comparison based purely on explained variance, rather than differences in likelihood architecture.

### Best Practices for Model Comparison
> **To ensure valid model comparisons, the `CONDITIONS` vector must be identical across all candidate models.** You should always set `CONDITIONS` to reflect the most fine-grained experimental partition used by the most complex model in your comparison set, even for models that do not utilize those factors in `VARYING_PARAMS`. Perhaps you don't know this yet at the start of your fitting, but once you know the most complex model, you can go back to more simple models and make sure you partition them in the same manner as your most complex model.

### Built-in Safety Check: Tracking n_bins
Because identical likelihood spaces are important for comparison, the pipeline automatically keeps track of the exact number of probability mass bins (n_bins) used during the $G^2$ calculation for every single fit. When you run 3_fit_compare.R, the script displays the average n_bins for each model in the summary table. If the script detects that models were evaluated on a different number of bins (which happens if CONDITIONS were not identical), it will print a warning that your BIC comparisons are potentially invalid.

### Implementation Example

If you wanted to test whether Drift Rate (`v`) is affected by Task Difficulty, you would fit two models in separate folders and then compare them using `3_fit_compare.R`. You must set up `1_run_pipeline.R` like this for each fit:

```r
# =========================================================
# FIT 1: THE NULL MODEL (No parameters vary)
# =========================================================
OUTPUT_FOLDER <- "vratio_null_model"

# The data is partitioned into independent likelihood blocks by the experimental manipulation.
CONDITIONS <- c("Difficulty") 

# Parameters are restricted. The model must utilize a single global drift rate across both Difficulty blocks.
VARYING_PARAMS <- list()


# =========================================================
# FIT 2: THE COMPLEX MODEL (Drift rate varies)
# =========================================================
OUTPUT_FOLDER <- "vratio_v_varying"

# The data partitioning structure remains strictly identical. This ensures the G^2 bins match the Null model exactly!
CONDITIONS <- c("Difficulty")

# The model is permitted to estimate a conditional Drift Rate (v) for Easy and Hard trials respectively.
VARYING_PARAMS <- list(v = ~ Difficulty)

# =========================================================
# Because both models were evaluated across identical likelihood 
# blocks and identical adaptive bins, the resulting Delta-BIC 
# accurately reflects the statistical evidence for the parameter shift.
```

------------------------------------------------------------------------

## Requirements
The following R packages are required. (If running on the HPC, ensure your conda environment contains these).
- `here` (for robust directory management)
- `Rcpp`, `RcppZiggurat` (for the C++ simulation engine)
- `DEoptim` (for global parameter optimization)
- `dplyr`, `tidyr`, `data.table` (for data wrangling)
- `ggplot2`, `patchwork` (for plotting)
