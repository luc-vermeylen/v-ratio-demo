---

editor_options: 
  markdown: 
    wrap: sentence
---

------------------------------------------------------------------------

# V-Ratio Fitting Pipeline (FCB Model)

A computational modeling pipeline for estimating metacognitive efficiency (**v-ratio**).

This pipeline fits the **Flexible Confidence Boundary (FCB)** models described in: \> *Herregods, S., Le Denmat, P., Vermeylen, L., & Desender, K. (2025). Modeling speed–accuracy trade-offs in the stopping rule for confidence judgments. Psychological Review.*

> *Desender, K., Vermeylen, L., & Verguts, T. (2022). Dynamic influences on static measures of metacognition. Nature communications, 13(1), 4208.*

To help users get started, **we have included the original datasets from Herregods et al. (2025)** in the `data/` folder.

------------------------------------------------------------------------

## 📂 Repository Overview

- 📁 **`data/`**: Place your dataset (CSV or RDS) here.
- 📁 **`results/`**: Model fits, CSV summaries, and PDF plots will appear here.
- 📁 **`R/`** and 📁 **`hpc/`**: Backend C++ engine, R helper functions, and Slurm scripts. *(Do not edit these).*
- 📄 **`1_run_pipeline.R`**: The main file. Configure settings and launch fits here.
- 📄 **`2_fit_assess.R`**: Aggregates fits and generates behavioral diagnostic plots.
- 📄 **`3_fit_compare.R`**: Compares different models to find the best fit via BIC.
- 📄 **`4_fit_stats.R`**: Runs statistics (T-Tests, ANOVAs) and generates parameter plots.

------------------------------------------------------------------------

## ⚠️ Required Data Format

Your dataset **must** contain the following columns exactly as named: \* **`sub_id`**: The subject identifier (or change `SUBJECT_COL` in the settings). \* **`rt`**: Primary decision reaction time (**strictly in seconds**, e.g., `0.450`). \* **`acc`**: Primary decision accuracy (`0` = Error, `1` = Correct). \* **`rtconf`**: Confidence reaction time (**strictly in seconds**). \* **`cj`**: Confidence judgment rating. \* If using binary confidence (`FCB_cj2`), code this as `0` (Low) and `1` (High). \* If using a 6-point scale (`FCB_cj6`), code this as `1` through `6`.

------------------------------------------------------------------------

## 🛠️ Step-by-Step Workflow

### Step 1: Configure & Fit (`1_run_pipeline.R`)

Open this file and adjust the User Settings. You will specify your dataset, your model (`FCB_cj2` or `FCB_cj6`, depending on your number of confidence response options), and which parameters vary by experimental conditions.

**How to use `VARYING_PARAMS`:** The pipeline uses standard R formula syntax to map parameters to your data columns. \* `list()`: (Empty list). Fits one global parameter per subject. \* `list(v = ~ as.factor(Difficulty))`: Estimates a separate Drift Rate for each level of "Difficulty". \* `list(vratio = ~ as.factor(Emotion) * as.factor(Validity))`: Fits a full interaction for the v-ratio parameter.

**How to execute:** Change the `RUN_MODE` variable to choose how to run the pipeline: \* `"single"`: Runs a test fit on 1 subject so you can check for errors. \* `"local_batch"`: Automatically loops through all subjects and fits them on your computer. \* `"hpc"`: Saves your configuration and prints the exact `sbatch` command you need to copy/paste into the cluster terminal to run a massive parallel array job.

> **What is saved?** For each subject, a `.rds` file is saved in `results/`. This file contains the `$best_params`, `$fit_metrics` (BIC/Cost), the original `$observations`, the C++ simulated `$predictions` (10,000 trials of the winning model), and the `$final_proportions` used for the likelihood calculation.

### Step 2: Visual Assessment (`2_fit_assess.R`)

**⚠️ You MUST run this script before comparing models or running statistics!** This script reads all the `.rds` files, aggregates them into a lightweight CSV, and evaluates if the model actually captured human behavior.

**Outputs:** 1. **`fit_metrics_summary.csv`**: A spreadsheet of BICs and parameters. 2. **`grand_average_fit_assessment.pdf`**: A visual report containing: \* *Decision RT Distributions:* A mirrored density plot of observed vs predicted RTs. \* *Confidence RT Distributions:* The speed of the confidence judgment. \* *Confidence Rating Mass:* A barplot showing if the model predicted the exact empirical frequencies of ratings (e.g., 1-6).

### Step 3: Model Comparison (`3_fit_compare.R`)

If you fit multiple models to test competing hypotheses (e.g., one folder where `v` varies, and another where `vratio` varies), this script ranks them. Open it and list your result folders in `FOLDERS_TO_COMPARE`.

**Outputs:** \* `model_comparison_plots.pdf`: Shows boxplots of BICs, Subject "Wins", and the Evidence Gap (delta-BIC). \* `summary_model_comparison_detailed.csv`: A table ranking models by overall BIC, Akaike Weights, and penalizations.

### Step 4: Parameter Statistics (`4_fit_stats.R`)

Once you have selected your winning model, this script automatically extracts and analyzes the estimated parameters.

**Outputs:** \* **Global Parameters Plot**: Generates boxplots of shared parameters against the theoretical optimization boundaries. This is crucial for checking for "boundary swarming" (ensuring the optimizer didn't get trapped against a limit). \* **Individual Trends & Group Means**: Generates Spaghetti plots and summary plots for any parameters you manipulated experimentally. \* **`parameter_anova_results.csv`**: Automatically runs T-Tests or Repeated-Measures ANOVAs on your varying parameters and saves the F-values, p-values, and Partial Eta Squared (PES) effect sizes.

------------------------------------------------------------------------

## 📖 The FCB Model Parameters Dictionary

The models rely on the following parameters. (Note: `a_slope` is fixed to 0 by default in the helper functions).

| Parameter | Stage | Description |
|:-----------------------|:-----------------------|:-----------------------|
| **`a`** | Primary | **Decision Boundary**: The total evidence required to make the primary choice. |
| **`v`** | Primary | **Drift Rate**: The speed of information accumulation for the primary choice. |
| **`ter`** | Primary | **Non-Decision Time**: Encoding and motor execution time (seconds). |
| **`a2`** | Confidence | **Confidence Boundary**: The total evidence space for the confidence judgment. |
| **`vratio`** | Confidence | **V-Ratio**: The ratio of post-decisional drift to primary drift ($v_{post} = v \times vratio$). If an error is made, $v_{post}$ is automatically reversed (multiplied by -1) in the C++ code to drive confidence downwards. |
| **`a2_slope_upper`** | Confidence | **Upper Bound Collapse**: How fast the upper confidence boundary collapses. |
| **`a2_slope_lower`** | Confidence | **Lower Bound Collapse**: How fast the lower confidence boundary collapses. |
| **`ter2`** | Confidence | **Confidence Non-Decision Time**: Delay specific to confidence reporting. Can be negative, implying metacognitive processing starts during the primary motor execution. |
| **`starting_point_confidence`** | Confidence | **Starting Point**: The fraction of `a2` where confidence accumulation begins (between 0 and 1). |

------------------------------------------------------------------------

## 🧮 Details of the Cost Function ($G^2$)

The pipeline uses R's `DEoptim` package to optimize parameters, calling a fast C++ simulation engine to generate predictions. It minimizes the **Likelihood Ratio Chi-Square statistic (**$G^2$). For binned multinomial data, minimizing $G^2$ is mathematically equivalent to Maximum Likelihood Estimation (MLE).

$$G^2 = 2 \times N \sum_{i=1}^{k} o_i \times \ln\left(\frac{o_i}{p_i}\right)$$
*(Where* $N$ is the number of trials, $o_i$ is the observed probability mass, and $p_i$ is the predicted probability mass).

### Evaluating the FCB Model (The 3 Targets)

To capture both stages, the pipeline calculates a joint cost function across three empirical targets. Following Herregods et al. (2025), all three targets are weighted equally (Weight = 1.0):

1.  **Target 1: Primary Decision RTs** The pipeline calculates dynamic RT quantiles separately for Correct and Error trials. To ensure robust parameter estimation even when conditions have very few trials, the pipeline automatically adjusts the binning resolution based on the exact sample size in that cell (ranging from 10 bins for \>200 trials, down to 2 bins for $\le$ 30 trials). *If a cell has* $\le$ 10 trials or zero variance, the pipeline safely abandons RT shape fitting and fits the pure marginal probability mass instead.
2.  **Target 2: Confidence RTs** The pipeline performs the exact same dynamic RT quantile binning for the Confidence Reaction Times (`rtconf`), again evaluated separately for correct and error trials.
3.  **Target 3: Pure Confidence Proportions** This target evaluates raw choices instead of time. It calculates the proportion of trials falling into each specific confidence rating bin (e.g., ratings 1 through 6), mapped separately for Corrects and Errors.

### Independent Likelihood Blocks

If you map parameters to experimental conditions, the pipeline does not lump all data into one giant distribution. Instead, it splits the data into **Independent Likelihood Blocks**. The $G^2$ cost is calculated for each condition independently (where probability mass sums to 1.0 locally), and then summed. This prevents **Parameter Mimicry**, ensuring that massive global RT shifts do not wash out subtle, condition-specific parameter effects.
