#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;

static Ziggurat::Ziggurat::Ziggurat zigg;

// ============================================================================
// MODEL 1: Flexible Confidence Bounds Model (BINARY CONFIDENCE)
// ============================================================================
// [[Rcpp::export]]
DataFrame FCB_cj2(
    int ntrials, 
    double a, double v, double a_slope, double ter, 
    double a2, double vratio, double a2_slope_upper, double a2_slope_lower, 
    double ter2, double starting_point_confidence,
    double sv = 0, double st = 0, double sz = 0,
    double s = 1, double dt = 0.001, double tmax = 10.0, unsigned long seed = 0) {
  
  RNGScope scope; 
  if (seed != 0) zigg.setSeed(seed);
  
  NumericMatrix DATA(ntrials, 4);
  double z_base = 0.5;
  
  for (int i = 0; i < ntrials; i++) {
    
    // --- TRIAL-BY-TRIAL VARIABILITY ---
    double trial_v = v;
    if (sv > 0) trial_v = v + sv * zigg.norm();
    
    double trial_start = a * z_base;
    if (sz > 0) {
      trial_start += sz * (R::runif(0.0, 1.0) - 0.5);
      if (trial_start > a) trial_start = a;
      if (trial_start < 0) trial_start = 0;
    }
    
    double trial_ter = ter;
    if (st > 0) {
      trial_ter += st * (R::runif(0.0, 1.0) - 0.5);
      if (trial_ter < 0) trial_ter = 0;
    }
    
    // --- DECISIONAL PROCESSING ---
    double evidence = trial_start;
    double t = 0;
    int cor = -1;
    
    while (evidence <= a && evidence >= 0) {
      t += dt;
      evidence += trial_v * dt + s * std::sqrt(dt) * zigg.norm();
      
      if (evidence >= a - t * a_slope) {
        cor = 1; break;
      } else if (evidence <= 0 + t * a_slope) {
        cor = 0; break;
      }
      if (t >= tmax) break;
    }
    
    DATA(i, 0) = t + trial_ter; // rt
    DATA(i, 1) = cor;           // acc
    
    // --- POST-DECISIONAL PROCESSING (Confidence) ---
    double t2 = 0;
    double v_post = v * vratio;
    evidence = a2 * starting_point_confidence;
    
    if (cor == 0) { v_post = -1 * v_post; }
    
    while ((evidence < a2 - t2 * a2_slope_upper) && (evidence > t2 * a2_slope_lower)) {
      t2 += dt;
      evidence += v_post * dt + s * std::sqrt(dt) * zigg.norm();
      if (t2 >= tmax) break;
    }
    
    DATA(i, 2) = t2 + ter2; // rtconf
    
    if (evidence >= a2 / 2) {
      DATA(i, 3) = 1; // High Confidence
    } else { 
      DATA(i, 3) = 0; // Low Confidence
    }
  }
  
  return DataFrame::create(
    Named("rt") = DATA(_, 0),
    Named("acc") = DATA(_, 1),
    Named("rtconf") = DATA(_, 2),
    Named("cj") = DATA(_, 3)
  );
} 

// ============================================================================
// MODEL 2: Flexible Confidence Bounds Model (6-POINT SCALE CONFIDENCE)
// ============================================================================
// [[Rcpp::export]]
DataFrame FCB_cj6(
    int ntrials, 
    double a, double v, double a_slope, double ter, 
    double a2, double vratio, double a2_slope_upper, double a2_slope_lower, 
    double ter2, double starting_point_confidence,
    double sv = 0, double st = 0, double sz = 0,
    double s = 1, double dt = 0.001, double tmax = 10.0, unsigned long seed = 0) {
  
  RNGScope scope; 
  if (seed != 0) zigg.setSeed(seed);
  
  NumericMatrix DATA(ntrials, 4);
  double z_base = 0.5;
  
  for (int i = 0; i < ntrials; i++) {
    
    // --- TRIAL-BY-TRIAL VARIABILITY ---
    double trial_v = v;
    if (sv > 0) trial_v = v + sv * zigg.norm();
    
    double trial_start = a * z_base;
    if (sz > 0) {
      trial_start += sz * (R::runif(0.0, 1.0) - 0.5);
      if (trial_start > a) trial_start = a;
      if (trial_start < 0) trial_start = 0;
    }
    
    double trial_ter = ter;
    if (st > 0) {
      trial_ter += st * (R::runif(0.0, 1.0) - 0.5);
      if (trial_ter < 0) trial_ter = 0;
    }
    
    // --- DECISIONAL PROCESSING ---
    double evidence = trial_start;
    double t = 0;
    int cor = -1;
    
    while (evidence <= a && evidence >= 0) {
      t += dt;
      evidence += trial_v * dt + s * std::sqrt(dt) * zigg.norm();
      
      if (evidence >= a - t * a_slope) {
        cor = 1; break;
      } else if (evidence <= 0 + t * a_slope) {
        cor = 0; break;
      }
      if (t >= tmax) break;
    }
    
    DATA(i, 0) = t + trial_ter; 
    DATA(i, 1) = cor;    
    
    // --- POST-DECISIONAL PROCESSING (Confidence) ---
    double t2 = 0;
    double v_post = v * vratio;
    evidence = a2 * starting_point_confidence;
    
    if (cor == 0) { v_post = -1 * v_post; }
    
    while ((evidence < a2 - t2 * a2_slope_upper) && (evidence > t2 * a2_slope_lower)) {
      t2 += dt;
      evidence += v_post * dt + s * std::sqrt(dt) * zigg.norm();
      if (t2 >= tmax) break;
    }
    
    DATA(i, 2) = t2 + ter2; // rtconf
    
    if (evidence < a2 / 6) {
      DATA(i, 3) = 1;
    } else if (evidence < 2 * a2 / 6) { 
      DATA(i, 3) = 2;
    } else if (evidence < 3 * a2 / 6) { 
      DATA(i, 3) = 3;
    } else if (evidence < 4 * a2 / 6) { 
      DATA(i, 3) = 4;
    } else if (evidence < 5 * a2 / 6) { 
      DATA(i, 3) = 5;
    } else { 
      DATA(i, 3) = 6;
    }
  } 
  
  return DataFrame::create(
    Named("rt") = DATA(_, 0),
    Named("acc") = DATA(_, 1),
    Named("rtconf") = DATA(_, 2),
    Named("cj") = DATA(_, 3)
  );
} 

// ============================================================================
// MODEL 3: Diffusion Model for Conflict Tasks (DMC)
// ============================================================================
// [[Rcpp::export]]
DataFrame DMC(
    double v_c = .5,            // Controlled drift rate (task-relevant)
    double a = 50,              // Bounds are [ -a , +a ], total separation = 2a
    double ter_mean = 300,      // Mean non-decision time (ms)
    double amp = 20,            // Amplitude of gamma pulse (automatic process)
    double tau = 130,           // Time-to-peak of gamma pulse (scale, ms)
    double beta = 3,            // Starting point variability shape (0 = no variability)
    double ter_sd = 20,         // Non-decision time variability (SD, ms)
    double alpha = 2,           // Shape of gamma pulse (fixed at 2)
    double s = 4,               // Within-trial noise (standard DMC scale: 4)
    double dt = 1,              // Time step size (1 ms)
    double tmax = 10000,        // Max trial time (ms)
    int ntrials = 5000,         // Number of trials to simulate
    unsigned long seed = 0) {       
  
  RNGScope scope; 
  if (seed != 0) zigg.setSeed(seed);
  
  NumericMatrix DATA(ntrials, 3);
  
  // 1. SETUP GAMMA LOOKUP TABLE (Pre-calculated for speed)
  int max_steps = (int)std::ceil(tmax / dt);
  std::vector<double> gamma_lookup(max_steps);
  
  double ae = alpha - 1.0;
  double const_term = std::pow(std::exp(1.0) / (ae * tau), ae);
  
  for (int i = 0; i < max_steps; ++i) {
    double t_val = (i + 1) * dt; 
    double shape = std::exp(-t_val / tau) * 
      std::pow(t_val, ae) * const_term * 
      (ae / t_val - 1.0 / tau);
    gamma_lookup[i] = shape; 
  } 
  
  double s_sqrt_dt = s * std::sqrt(dt);
  double b_upper = a;
  double b_lower = -a;
  double width = b_upper - b_lower;
  
  // 2. TRIAL LOOP (50% Congruent (+1), 50% Incongruent (-1))
  std::vector<double> congruency_vec(ntrials);
  for (int i = 0; i < ntrials; i++) {
    congruency_vec[i] = (i % 2 == 0) ? 1.0 : -1.0;
  }
  
  for (int i = 0; i < ntrials; i++) {
    
    // A. Starting Point Variability (Beta Distribution)
    double sp = 0.0;
    if (beta > 0) {
      sp = R::rbeta(beta, beta) * width + b_lower;
    }  
    
    // B. Non-Decision Time Variability (Normal Distribution)
    double ter = R::rnorm(ter_mean, ter_sd);
    if (ter < 0) ter = 0;
    
    // C. Evidence Accumulation
    double evidence = sp;
    double t = 0.0;
    int step_idx = 0;
    int accuracy = NA_REAL;
    double congruency = congruency_vec[i];
    double current_amp = amp * congruency;
    
    while (t < tmax) {
      t += dt;
      
      double v_automatic = 0.0;
      if (step_idx < max_steps) {
        v_automatic = current_amp * gamma_lookup[step_idx];
      } 
      step_idx++;
      
      double v_total = v_c + v_automatic;
      evidence += v_total * dt + s_sqrt_dt * zigg.norm();
      
      if (evidence >= b_upper) {
        accuracy = 1;
        break;
      } 
      if (evidence <= b_lower) {
        accuracy = 0;
        break;
      }
    } 
    
    DATA(i, 0) = t + ter;
    DATA(i, 1) = accuracy;
    DATA(i, 2) = congruency;
  }
  
  return DataFrame::create(
    Named("rt") = DATA(_, 0), 
    Named("acc") = DATA(_, 1),
    Named("congruency") = DATA(_, 2)
  );
} 

// ============================================================================
// MODEL 4: Standard Drift Diffusion Model (DDM)
// ============================================================================
// [[Rcpp::export]]
DataFrame DDM(
    int ntrials,
    double a, double v, double ter, double z = 0.5,
    double sv = 0, double st = 0, double sz = 0,
    double s = 1, double dt = 0.001, double tmax = 10.0, unsigned long seed = 0) {
  
  RNGScope scope; 
  if (seed != 0) zigg.setSeed(seed);
  
  NumericMatrix DATA(ntrials, 2);
  
  for (int i = 0; i < ntrials; i++) {
    
    double trial_v = v;
    if (sv > 0) trial_v = v + sv * zigg.norm();
    
    double trial_start = a * z;
    if (sz > 0) {
      trial_start += sz * (R::runif(0.0, 1.0) - 0.5);
      if (trial_start > a) trial_start = a;
      if (trial_start < 0) trial_start = 0;
    }
    
    double trial_ter = ter;
    if (st > 0) {
      trial_ter += st * (R::runif(0.0, 1.0) - 0.5);
      if (trial_ter < 0) trial_ter = 0;
    }
    
    double evidence = trial_start;
    double t = 0;
    int cor = -1;
    
    while (evidence <= a && evidence >= 0) {
      t += dt;
      evidence += trial_v * dt + s * std::sqrt(dt) * zigg.norm();
      
      if (evidence >= a) {
        cor = 1; break;
      } else if (evidence <= 0) {
        cor = 0; break;
      }
      if (t >= tmax) break;
    }
    
    DATA(i, 0) = t + trial_ter;
    DATA(i, 1) = cor;
  }
  
  return DataFrame::create(
    Named("rt") = DATA(_, 0),
    Named("acc") = DATA(_, 1)
  );
}

// ============================================================================
// MODEL 5: RuleOnset -- Late, variable rule onset (Lognormal)
// ============================================================================
static Ziggurat::Ziggurat::Ziggurat zigg_rule;

// [[Rcpp::export]]
DataFrame RuleOnset(
    int ntrials,
    double v, 
    double sv, 
    double thr, 
    double ter, 
    double st,
    double alpha, 
    double mR, 
    double sR, 
    double tau, 
    double gamma,
    double gate = 1.0, 
    double prop_congruent = 0.5,
    double s = 1.0, double dt = 0.001, double tmax = 5.0, unsigned long seed = 0) {
  
  RNGScope scope;
  if (seed != 0) zigg_rule.setSeed(seed);
  
  // Convert mean and SD of R to lognormal parameters
  double sdlog   = std::sqrt(std::log(1.0 + (sR * sR) / (mR * mR)));
  double meanlog = std::log(mR) - 0.5 * sdlog * sdlog;
  double decay   = std::exp(-dt / tau);   // Per-step fading of location activation
  double noise   = s * std::sqrt(dt);
  int n_congruent = (int) std::round(ntrials * prop_congruent);
  
  NumericVector RT(ntrials), ACC(ntrials), CONG(ntrials), RULE(ntrials), HWRESP(ntrials);
  
  for (int i = 0; i < ntrials; i++) {
    
    // --- TRIAL-LEVEL DRAWS ---
    double r = (i < n_congruent) ? 1.0 : -1.0;                  
    double trial_v = v;
    if (sv > 0) trial_v = v + sv * zigg_rule.norm();             // drift, N(v, sv)
    double rule_t = std::exp(meanlog + sdlog * zigg_rule.norm()); // rule onset, lognormal
    double trial_ter = ter;
    if (st > 0) {
      trial_ter += st * (R::runif(0.0, 1.0) - 0.5);              // uniform jitter on ter
      if (trial_ter < 0) trial_ter = 0;
    }
    
    // --- ACCUMULATION ---
    double y_loc = 0.0, y_rule = 0.0, t = 0.0, hw_resp = NA_REAL;
    int acc = -1;
    
    while (t < tmax) {
      t += dt;
      double dE = trial_v * dt + noise * zigg_rule.norm();       // colour evidence
      bool before_rule = (t <= rule_t);
      
      if (before_rule) {
        y_loc += alpha * r * dE;                                 // location route
      } else {
        y_rule += dE;                                            // colour route
        y_loc  *= decay;                                         // location fading
      }
      
      double H_c = std::max(y_rule, 0.0) + std::max(y_loc, 0.0);   // correct hand
      double H_w = std::max(-y_rule, 0.0) + std::max(-y_loc, 0.0); // wrong hand
      
      if (!before_rule || gate < 0.5) {                          // the gate
        if (H_c - gamma * H_w >= thr) { acc = 1; hw_resp = H_w; break; }
        if (H_w - gamma * H_c >= thr) { acc = 0; hw_resp = H_w; break; }
      }
    }
    
    RT[i]     = t + trial_ter;
    ACC[i]    = acc;
    CONG[i]   = r;
    RULE[i]   = rule_t;
    HWRESP[i] = hw_resp / thr;
  }
  
  return DataFrame::create(
    Named("rt")                 = RT,
    Named("acc")                = ACC,
    Named("congruency")         = CONG,
    Named("rule_onset")         = RULE,
    Named("wrong_hand_at_resp") = HWRESP
  );
}