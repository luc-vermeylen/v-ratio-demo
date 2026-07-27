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
    // Drift Rate Variability (Normal distribution: N(v, sv))
    double trial_v = v;
    if (sv > 0) trial_v = v + sv * zigg.norm();
    
    // Starting Point Variability (Uniform distribution width sz around z_base)
    double trial_start = a * z_base;
    if (sz > 0) {
      trial_start += sz * (R::runif(0.0, 1.0) - 0.5);
      // Clamp to boundaries just in case sz is very large
      if (trial_start > a) trial_start = a;
      if (trial_start < 0) trial_start = 0;
    }
    
    // Non-Decision Time Variability (Uniform distribution width st around ter)
    double trial_ter = ter;
    if (st > 0) {
      trial_ter += st * (R::runif(0.0, 1.0) - 0.5);
      if (trial_ter < 0) trial_ter = 0; // Prevent negative physical time
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
      if (t >= tmax) break; // Safety break
    }
    
    DATA(i, 0) = t + trial_ter; // rt
    DATA(i, 1) = cor;           // acc
     
     // --- POST-DECISIONAL PROCESSING (Confidence) ---
    double t2 = 0;
    double v_post = v * vratio;
    evidence = a2 * starting_point_confidence;
    
    if (cor == 0) { v_post = -1 * v_post; } // reverse drift for errors
     
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