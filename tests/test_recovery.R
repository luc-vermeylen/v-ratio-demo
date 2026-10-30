# ==============================================================================
# test_recovery.R
# Unit Tests for Parameter Recovery Suite
# ==============================================================================

library(here)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppZiggurat)

source(here("R", "helper_functions.R"))
sourceCpp(here("R", "models.cpp"))

cat("\n===================================================\n")
cat("RUNNING PARAMETER RECOVERY UNIT TESTS\n")
cat("===================================================\n\n")

n_passed <- 0
n_failed <- 0

assert <- function(desc, condition) {
  if (isTRUE(condition)) {
    cat(sprintf("  [PASS] %s\n", desc))
    n_passed <<- n_passed + 1
  } else {
    cat(sprintf("  [FAIL] %s\n", desc))
    n_failed <<- n_failed + 1
  }
}

# ------------------------------------------------------------------------------
# TEST 1: Custom Bounds Override and Fallback Logic
# ------------------------------------------------------------------------------
cat("--- Test Suite 1: Custom Bounds Fallback Logic ---\n")

master_b <- get_bounds("FCB_cj6")
custom_b <- list(a = c(1.2, 2.0), v = c(1.0, 1.8))

test_lower <- master_b$lower
test_upper <- master_b$upper

for (p in names(custom_b)) {
  test_lower[p] <- custom_b[[p]][1]
  test_upper[p] <- custom_b[[p]][2]
}

assert("Custom 'a' lower bound updated to 1.2", test_lower["a"] == 1.2)
assert("Custom 'a' upper bound updated to 2.0", test_upper["a"] == 2.0)
assert("Custom 'v' lower bound updated to 1.0", test_lower["v"] == 1.0)
assert("Unspecified 'vratio' retained master lower bound (0.2)", test_lower["vratio"] == master_b$lower["vratio"])
assert("Unspecified 'ter' retained master upper bound (1.5)", test_upper["ter"] == master_b$upper["ter"])

# ------------------------------------------------------------------------------
# TEST 2: Empirical Envelope Outlier Filtering (Approach 3A)
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 2: Empirical Outlier Filtering ---\n")

# Mock empirical results table with 20 good fits and 2 bad outliers
set.seed(123)
mock_emp <- data.frame(
  subject_id = 1:22,
  a = c(rnorm(20, mean = 1.5, sd = 0.2), 999.0, 1.4),     # Outlier 1: insane boundary
  v = c(rnorm(20, mean = 1.8, sd = 0.3), 1.8, 1.7),
  cost = c(rnorm(20, mean = 500, sd = 50), 450, 1e9)       # Outlier 2: exploded cost penalty
)

# Step A: Filter by cost
valid_costs <- mock_emp$cost[is.finite(mock_emp$cost) & mock_emp$cost < 1e6]
cost_thresh <- median(valid_costs) + 3 * IQR(valid_costs)
filtered_emp <- mock_emp[mock_emp$cost < cost_thresh & mock_emp$cost < 1e6, ]

assert("Cost outlier (1e9) successfully filtered out", nrow(filtered_emp) == 21)

# Step B: Percentile trimming (5% to 95%)
q_trim <- quantile(filtered_emp$a, probs = c(0.05, 0.95))
clamped_low  <- max(master_b$lower["a"], q_trim[1])
clamped_high <- min(master_b$upper["a"], q_trim[2])

assert("Extreme 'a' value (999.0) trimmed from upper envelope", clamped_high < 3.0)
assert("Empirical upper bound is within reasonable range (< 2.5)", clamped_high < 2.5)
assert("Empirical lower bound is above master bound (> 0.3)", clamped_low > 0.3)

# ------------------------------------------------------------------------------
# TEST 3: Ground Truth Boundedness Under Rejection Sampling
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 3: Parameter Sampling Invariant Checks ---\n")

dummy_design <- data.frame(Difficulty = factor(c("Easy", "Hard")))
cfg <- get_param_config("FCB_cj6", list(v = ~ Difficulty), dummy_design)
cfg <- fix_parameters(cfg, list(a_slope = 0, sv = 0, st = 0, sz = 0))

draws <- replicate(50, {
  raw <- draw_random_params(cfg)
  p_cell <- pluck_params(raw, cfg$fixed, dummy_design[1, , drop = FALSE], list(v = ~ Difficulty), "FCB_cj6")
  p_cell$a >= master_b$lower["a"] && p_cell$a <= master_b$upper["a"]
})

assert("All 50 sampled draws obey master parameter boundaries for 'a'", all(draws))

# ------------------------------------------------------------------------------
# TEST 4: Trial Data Dimensions and Factor Structure
# ------------------------------------------------------------------------------
cat("\n--- Test Suite 4: Simulation Data Generation Output ---\n")

sim_res <- FCB_cj6(
  ntrials = 50, a = 1.5, v = 2.0, a_slope = 0, ter = 0.3,
  a2 = 10, vratio = 1.2, a2_slope_upper = 5, a2_slope_lower = 5,
  ter2 = 0.2, starting_point_confidence = 0.5,
  sv = 0, st = 0, sz = 0, s = 1, dt = 0.001, tmax = 10.0, seed = 42
)

assert("C++ simulator produces exactly 50 trials", nrow(sim_res) == 50)
assert("All required columns present (rt, acc, rtconf, cj)", all(c("rt", "acc", "rtconf", "cj") %in% names(sim_res)))
assert("Accuracy strictly binary (0 or 1)", all(sim_res$acc %in% c(0, 1)))
assert("Confidence ratings in 1..6 range for FCB_cj6", all(sim_res$cj %in% 1:6))

# ------------------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------------------
cat("\n===================================================\n")
cat(sprintf("TEST RESULTS: %d PASSED, %d FAILED\n", n_passed, n_failed))
cat("===================================================\n")

if (n_failed > 0) stop("Unit test suite failed!")