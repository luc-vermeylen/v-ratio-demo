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
    double s = 1, double dt = 0.001, double tmax = 10.0, unsigned long seed = 0) {
  
  // Parallel RNG Setup
  RNGScope scope; 
  if (seed != 0) zigg.setSeed(seed);
  
  NumericMatrix DATA(ntrials, 4);
  double z = 0.5;
  
  for (int i = 0; i < ntrials; i++) {
    // 1. Decisional processing
    double evidence = a * z;
    double t = 0;
    int cor = -1;
     
    while (evidence <= a && evidence >= 0) {
      t += dt;
      evidence += v * dt + s * std::sqrt(dt) * zigg.norm();
      
      if (evidence >= a - t * a_slope) {
        cor = 1; break;
      } else if (evidence <= 0 + t * a_slope) { 
        cor = 0; break;
      } 
      if (t >= tmax) break;
    }
     
    DATA(i, 0) = t + ter; // rt
    DATA(i, 1) = cor;     // acc
     
    // 2. Post-decisional processing (Confidence)
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
    double s = 1, double dt = 0.001, double tmax = 10.0, unsigned long seed = 0) {
  
  // Parallel RNG Setup
  RNGScope scope; 
  if (seed != 0) zigg.setSeed(seed);
  
  NumericMatrix DATA(ntrials, 4);
  double z = 0.5;
  
  for (int i = 0; i < ntrials; i++) {
    // 1. Decisional processing
    double evidence = a * z;
    double t = 0;
    int cor = -1;
     
    while (evidence <= a && evidence >= 0) {
      t += dt;
      evidence += v * dt + s * std::sqrt(dt) * zigg.norm();
      
      if (evidence >= a - t * a_slope) {
        cor = 1; break;
      } else if (evidence <= 0 + t * a_slope) { 
        cor = 0; break;
      } 
      if (t >= tmax) break;
    }
     
    DATA(i, 0) = t + ter; // rt
    DATA(i, 1) = cor;     // acc
     
    // 2. Post-decisional processing (Confidence)
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