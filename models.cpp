#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
DataFrame DMC(
    double v_c = .5,            // response-relevant drift
    double a = 50,              // Bounds are [ -a , +a ], hence boundary separation is 2a
    double ter_mean = 300,      // mean non-decision time
    double amp = 20,            // amplitude of gamma pulse (strength of irrelevant evidence)
    double tau = 130,           // peak time of gamma pulse (scale)
    double beta = 3,            // Starting point variability shape, beta = 0 means no variability
    double ter_sd = 20,         // non-decision time variability (normal distribution)
    double alpha = 2,           // shape of gamma pulse
    double s = 4,               // within-trial noise
    double dt = 1,              // time step size
    double tmax = 10000,        // maximum time per trial (trials at max are NA, i.e. "undecided")
    int ntrials = 5000,         // number of trials to simulate
    unsigned long seed = 0) {       
  
  // Set seeds (important for parallel workers)
  RNGScope scope; // synchronizes RNG with R
  if (seed != 0) zigg.setSeed(seed);
  
  // Output container
  NumericMatrix DATA(ntrials, 3);
   
  // 1. SETUP GAMMA LOOKUP TABLE (for speedup)
  // ---------------------------------------------------------
  // Pre-calculate the gamma pulse.
  // Formula: exp(-t/tau) * (t*e / ((a-1)*tau))^(a-1) * ((a-1)/t - 1/tau)
  
  int max_steps = (int)std::ceil(tmax / dt);
  std::vector<double> gamma_lookup(max_steps);
  
  // Constants for Gamma
  double ae = alpha - 1.0;
  // Constant part: (e / ( (alpha-1)*tau )) ^ (alpha-1)
  double const_term = std::pow(std::exp(1.0) / (ae * tau), ae);
  
  for(int i = 0; i < max_steps; ++i) {
    double t_val = (i + 1) * dt; 
    double shape = std::exp(-t_val / tau) * 
      std::pow(t_val, ae) * const_term * 
      (ae / t_val - 1.0 / tau);
    gamma_lookup[i] = shape; 
  } 
  
  double s_sqrt_dt = s * std::sqrt(dt);
  double b_upper = a;
  double b_lower = -a;
  double width = b_upper - b_lower; // Total separation (2*a)
   
  // 2. TRIAL LOOP
  // ---------------------------------------------------------
  
  // half congruent (1) and half incongruent (-1) trials
  std::vector<double> congruency_vec;
  int per_cond = ntrials / 2;
  for (int i = 0; i < per_cond; i++) {
    congruency_vec.push_back(1);
    congruency_vec.push_back(-1);
  }
   
  // Loop through trials
  for (int i = 0; i < ntrials; i++) {
    
    // A. Starting Point Variability (Beta Distribution)
    // Map Beta(0-1) to Range (-a to +a)
    double sp = 0.0;
    if (beta > 0) { // If beta=0, sp=0 (Unbiased)
      sp = R::rbeta(beta, beta) * width + b_lower; // Scale to [b_lower, b_upper]
    }  
    
    // B. Non-Decision Time Variability (Normal Distribution)
    double ter = R::rnorm(ter_mean, ter_sd);
    if (ter < 0) ter = 0; // Safety check to avoid negative NDT
    
    // C. Accumulate Evidence
    double evidence = sp;
    double t = 0.0;
    int step_idx = 0;
    int accuracy = NA_REAL; // NA = undecided/max_time reached
    double congruency = congruency_vec[i];
    double current_amp = amp * congruency; // Signed amplitude
     
    while (t < tmax) {
      t += dt;
      
      // Lookup Drift
      double v_automatic = 0.0;
      if (step_idx < max_steps) {
        v_automatic = current_amp * gamma_lookup[step_idx];
      } 
      step_idx++;
      
      // Total Drift = Controlled (v_c) + Automatic (Pulse)
      double v_total = v_c + v_automatic;
      
      evidence += v_total * dt + s_sqrt_dt * zigg.norm();
      
      // Check Bounds
      if (evidence >= b_upper) {
        accuracy = 1; // Correct
        break;
      } 
      if (evidence <= b_lower) {
        accuracy = 0; // Error
        break;
      }
    } 
    
    // Save Results
    DATA(i, 0) = t + ter;
    DATA(i, 1) = accuracy;
    DATA(i, 2) = congruency;
  }
   
  DataFrame df = DataFrame::create(Named("rt") = DATA(_, 0), 
                                   Named("acc") = DATA(_, 1),
                                   Named("congruency") = DATA(_, 2));
  return df;
} 

