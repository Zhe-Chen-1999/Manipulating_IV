# IMPLEMENT SIMULATION STUDY 3
################################################################################
# Import libraries
library(Rlab)
library(dplyr)
library(tidyr)
library(xtable)
library(locfit)
library(nbpMatching)
################################################################################
set.seed(333)
################################################################################
# Write a function that generates covariates
gen_covariates = function(n) {
  x = rbinom(n, size = 1, p = 0.5)
  #  x = rnorm(n, mean = 0, sd = 5)
  x = data.frame(x)
  return(x)
}

# Write a function that generates observed IV dose
gen_Z_obs = function(n) {
  Z.obs = runif(n, min = 5, max = 50)
  return(Z.obs)
}

# Write a function that generates potential outcomes for all the units
gen_r = function(n) {
  r_d0 = rnorm(n, mean = 0, sd = 1)
  psi = runif(n, min = 4, max = 6)
  r_d1 = r_d0 + psi
  return(cbind(r_d0, r_d1))
}

# Write a function that adds a dose caliper to a distance matrix
add_dose_caliper = function(Z.obs, penalty, caliper, distmat) {
  distmat.Z.obs = abs(outer(Z.obs, Z.obs, "-"))
  penalty.dose.caliper = penalty * (distmat.Z.obs <= caliper)
  distmat = distmat + penalty.dose.caliper
  return(distmat)
}

# Write a function that performs nbp matching
nbp_match = function(x, Z.obs, penalty, caliper.1, caliper.2, caliper.3, caliper.4) {
  # Construct a distance matrix
  dist.mat.list = gendistance(x)
  dist.mat = dist.mat.list$dist
  # Match 1: Add a small IV dose caliper to the distance matrix
  dist.mat.caliper.1 = add_dose_caliper(Z.obs, penalty, caliper.1, dist.mat)
  # Match 2: Add a large IV dose caliper to the distance matrix
  dist.mat.caliper.2 = add_dose_caliper(Z.obs, penalty, caliper.2, dist.mat)
  # Match 3: Add an even larger IV dose caliper to the distance matrix
  dist.mat.caliper.3 = add_dose_caliper(Z.obs, penalty, caliper.3, dist.mat)
  # Match 4: Add the largest IV dose caliper to the distance matrix
  dist.mat.caliper.4 = add_dose_caliper(Z.obs, penalty, caliper.4, dist.mat)
  
  # Match 0: Do matching without IV dose caliper
  dist.mat.no.caliper = distancematrix(dist.mat.list)
  matching.no.caliper = nonbimatch(dist.mat.no.caliper)$halves
  # Match 1: Do matching with the small IV dose caliper
  dist.mat.caliper.1 = distancematrix(dist.mat.caliper.1)
  matching.caliper.1 = nonbimatch(dist.mat.caliper.1)$halves
  # Match 2: Do matching with the large IV dose caliper
  dist.mat.caliper.2 = distancematrix(dist.mat.caliper.2)
  matching.caliper.2 = nonbimatch(dist.mat.caliper.2)$halves
  # Match 3: Do matching with even larger IV dose caliper
  dist.mat.caliper.3 = distancematrix(dist.mat.caliper.3)
  matching.caliper.3 = nonbimatch(dist.mat.caliper.3)$halves
  # Match 4: Do matching with the largest IV dose caliper
  dist.mat.caliper.4 = distancematrix(dist.mat.caliper.4)
  matching.caliper.4 = nonbimatch(dist.mat.caliper.4)$halves
  
  return(list("matching.caliper.0" = matching.no.caliper, 
              "matching.caliper.1" = matching.caliper.1,
              "matching.caliper.2" = matching.caliper.2,
              "matching.caliper.3" = matching.caliper.3,
              "matching.caliper.4" = matching.caliper.4))
}

# Write a function that extracts the indices of matched pairs
extract_match = function(matches, caliper.idx) {
  matched.pairs = matches[[caliper.idx]] # Caliper.idx takes values btw {1,2,3,4}
  first.idx = matched.pairs$Group1.Row
  second.idx = matched.pairs$Group2.Row
  return(cbind(first.idx, second.idx))
}

# Write a function that extracts IV doses based on matched pairs
extract_IV_doses = function(matched.idx, Z.obs) {
  first.idx = matched.idx[, 1]
  second.idx = matched.idx[, 2]
  Z.obs.first = Z.obs[first.idx]
  Z.obs.second = Z.obs[second.idx]
  Z.obs.matched = cbind(Z.obs.first, Z.obs.second)
  Z.obs.min = apply(Z.obs.matched, 1, FUN = min)
  Z.obs.max = apply(Z.obs.matched, 1, FUN = max)
  Z.obs.matched = cbind(Z.obs.min, Z.obs.max)
  return(Z.obs.matched) # indices of matched pairs are from the original
}

# Generate IV dose assignment probability for each matched pair
gen_pi = function(Z, gamma) { # If gamma == 0, randomization i.e., Prop. 1
  I = nrow(Z)
  Z.min = Z[, 1]
  Z.max = Z[, 2]
  Z.diff = Z.max - Z.min
  exp.term = exp(gamma * Z.diff)
  LB = 0.5
  UB = exp.term / (exp.term + 1)
  pi = runif(n = I, min = LB, max = UB)
  
  return(pi)
}

# Write a function that computes the average difference in the IV dose  
compute_IV_dose_avg_diff = function(Z.obs.matched) {
  return(mean(Z.obs.matched[,2] - Z.obs.matched[,1]))
}

# Write a function that generates indicator of which matched unit gets the larger IV dose
gen_obs_IV_dose_indi = function(Z.min, Z.max, pi) {
  J = rbern(n = 1, prob = pi)
  if (J == 1) {
    Z.obs = c(Z.max, Z.min)
  } else {
    Z.obs = c(Z.min, Z.max)
  }
  return(Z.obs)
}

# Write a function that generates observed IV dose (pi is a list)
gen_obs_IV_dose_unstacked = function(Z, pi) {
  Z.obs.raw = mapply(FUN = gen_obs_IV_dose_indi, Z[, 1], Z[, 2], pi)
  Z.obs.raw = t(Z.obs.raw)
  colnames(Z.obs.raw) = c("Z_i1", "Z_i2")
  return(Z.obs.raw)
}

# Write a function that stacks matched IV dose s.t. the first I rows are IV doses for the first unit in each matched pair and the next I rows are IV doses for the second unit
gen_obs_IV_dose_stacked = function(Z.obs.unstacked) {
  Z.obs.stacked = data.frame(Z = c(Z.obs.unstacked[,1], Z.obs.unstacked[,2]))
  return(Z.obs.stacked)
}

# Write a function that creates binary Z for each unit
create_binary_Z_indi = function(Z_i1, Z_i2) {
  if (Z_i1 < Z_i2) {
    return(c(0,1)) 
  } else {
    return(c(1,0))
  }
}

# Write a function that creates binary Z for all the units
create_binary_Z = function(Z.obs.unstacked) {
  Z.ind = mapply(FUN = create_binary_Z_indi, Z.obs.unstacked[,1], Z.obs.unstacked[,2])
  Z.ind.unstacked = t(Z.ind)
  Z.ind.stacked = c(Z.ind.unstacked[,1], Z.ind.unstacked[,2])
  return(data.frame(Z.ind.stacked))
}

# Generate compliance status for each unit
gen_compliance_indi = function(IV.diff, a) {
  if (IV.diff >= a) {
    comp.status = 1 # Complier
  } else {
    comp.status = sample(x = c(2,3), size = 1, prob = c(0.5, 0.5)) # Flip a fair coin for 
  }
  return(comp.status)
}

# Generate compliance status for all the units
gen_compliance = function(Z.obs.matched.caliper) {
  I = nrow(Z.obs.matched.caliper)
  a = runif(I, 15, 18)
  IV.diff = rep(Z.obs.matched.caliper[, 2] - Z.obs.matched.caliper[, 1], 2)
  compliance = mapply(gen_compliance_indi, IV.diff, a)
  compliance = data.frame(S = compliance)
  return(compliance)
}

# Write a function that computes compliance rate
calc_compl_rate = function(S) {
  n = nrow(S)
  n.compliers = length(S[S$S == 1, ])
  iota.C = n.compliers / n
  return(iota.C)
}

# Write a function that generates potential treatments for each unit
gen_poten_trt_indi = function(S_ij) {
  # Compliance status: 1 = complier, 2 = always-taker, 3 = never-taker
  if (S_ij == 1) {
    d_T = 1
    d_C = 0
  } else if (S_ij == 2) {
    d_T = 1
    d_C = 1
  } else {
    d_T = 0
    d_C = 0
  }
  return(c(d_T, d_C))
}

# Write a function that generates potential treatments for all the units
gen_poten_trt = function(S) {
  df = t(apply(S, MARGIN=1, FUN = gen_poten_trt_indi))
  colnames(df) = c("d_T", "d_C")
  return(df)
}

# Generate potential outcome based on encouragement of treatment or control for each individual
gen_poten_out_enc_indi = function(Sij, r_d1, r_d0) {
  r_T = case_when(Sij == 1 ~ r_d1, # Complier
                  Sij == 2 ~ r_d1, # Always-taker
                  Sij == 3 ~ r_d0) # Never-taker
  r_C = case_when(Sij == 1 ~ r_d0,
                  Sij == 2 ~ r_d1,
                  Sij == 3 ~ r_d0)
  return(c(r_T, r_C))
}

# Generate potential outcome based on encouragement of treatment or control
gen_poten_out_enc = function(S, r_d) {
  df = mapply(FUN = gen_poten_out_enc_indi, S$S, r_d[ ,1], r_d[ ,2])
  df = t(df)
  colnames(df) = c("r_T", "r_C")
  return(df)
}

# Write a function that generates observed treatment
gen_obs_trt = function(Z.ind.stacked, d) {
  d_T = d[ ,1]
  d_C = d[ ,2]
  Z.ind.stacked*d_T + (1-Z.ind.stacked)*d_C
}

# Write a function that generates observed outcome 
gen_obs_outcome = function(Z.ind.stacked, r) {
  r_T = r[ ,1]
  r_C = r[ ,2]
  Z.ind.stacked*r_T + (1-Z.ind.stacked)*r_C 
}

# Write a function that computes IV dose-dependent odds
compute_IV_odds = function(gamma, Z.obs.unstacked) {
  Z.obs.max = apply(Z.obs.unstacked, FUN = max, MARGIN = 1)
  Z.obs.min = apply(Z.obs.unstacked, FUN = min, MARGIN = 1)
  odds = exp(gamma * (Z.obs.max - Z.obs.min))
  return(odds)
}

# Create a dataframe containing Z, D, R with N=2I rows. The first I rows correspond to one unit within each matched pair and the second I rows are the corresponding matched units
create_data = function(Z.ind, D, R) {
  df = data.frame(Z.ind, D, R)
  colnames(df) = c("Z.ind", "D", "R")
  return(df)
}

## The lower limit of the (1-alpha)-level one-sided lower CI for (LB - K0)
transformed.tau.bar.lower.limit = function(data, K0, alpha, R, Gamma) {
  I = nrow(data) / 2
  tau_hat =  R[which(data$Z==1)]- K0*data$D[which(data$Z==1)] - (R[which(data$Z==0)]- K0*data$D[which(data$Z==0)])
  tau_hat_Gamma = (Gamma + 1)*((Gamma + 1)*tau_hat - (Gamma - 1)*abs(tau_hat))/(4*Gamma)
  tau_bar_Gamma = mean(tau_hat_Gamma)
  S_square_Gamma = sum((tau_hat_Gamma - tau_bar_Gamma)^2)/(I*(I-1)) 
  transformed.tau.bar.lower.limit = tau_bar_Gamma - qnorm(1-alpha)*sqrt(S_square_Gamma)
  return(transformed.tau.bar.lower.limit)
} 

## The lower limit of the (1-alpha/2)-level one-sided lower CI for LB
LB.lower.limit = function(data, K0, alpha = 0.05/2, R, Gamma) {
  LB.lower.limit = transformed.tau.bar.lower.limit(data = data, K0 = K0, alpha = alpha, R = R, Gamma = Gamma) + K0 
  return(LB.lower.limit)
} 

## The upper limit of the (1-alpha/2)-level one-sided upper CI for UB
UB.upper.limit = function(data, K1, alpha = 0.05/2, R, Gamma) {
  data$Z =  1-data$Z
  UB.upper.limit = -transformed.tau.bar.lower.limit(data = data, K0 = K1, alpha = alpha, R = R, Gamma = Gamma) + K1
  return(UB.upper.limit)
} 

# Write a function that computes the width of the confidence interval for partial identification bounds
compute_width_CI_PI = function(gamma, Z.obs.matched, pi, K, data) {
  K0 = K[1]
  K1 = K[2]
  R = data$R
  alpha = 0.05 / 2
  
  Z.obs.unstacked = gen_obs_IV_dose_unstacked(Z.obs.matched, pi)
  Z.ind = create_binary_Z(Z.obs.unstacked)
  Gamma = compute_IV_odds(gamma, Z.obs.unstacked)
  
  # CI for LB
  CI.LB = LB.lower.limit(data, K0, alpha, R, Gamma)
  
  # CI for UB
  CI.UB = UB.upper.limit(data, K1, alpha, R, Gamma)
  
  # Compute width of the CI for the partial identification bound
  width = CI.UB - CI.LB
  
  return(width) 
}

################################################################################
# PERFORM SIMULATION #
################################################################################
# Perform one iteration of the simulation study
one_iter = function(n) {
  # Generate data
  x = gen_covariates(n)
  Z = gen_Z_obs(n)
  r_d = gen_r(n)
  
  # Perform matching
  matches = nbp_match(x, Z, penalty = 10, caliper.1 = 6, caliper.2 = 9, 
                      caliper.3 = 12, caliper.4 = 15)
  
  matched.idx.caliper.0 = extract_match(matches, caliper.idx = 1)
  matched.idx.caliper.1 = extract_match(matches, caliper.idx = 2)
  matched.idx.caliper.2 = extract_match(matches, caliper.idx = 3)
  matched.idx.caliper.3 = extract_match(matches, caliper.idx = 4)
  matched.idx.caliper.4 = extract_match(matches, caliper.idx = 5)
  
  Z.obs.matched.caliper.0 = extract_IV_doses(matched.idx.caliper.0, Z)
  Z.obs.matched.caliper.1 = extract_IV_doses(matched.idx.caliper.1, Z)
  Z.obs.matched.caliper.2 = extract_IV_doses(matched.idx.caliper.2, Z)
  Z.obs.matched.caliper.3 = extract_IV_doses(matched.idx.caliper.3, Z)
  Z.obs.matched.caliper.4 = extract_IV_doses(matched.idx.caliper.4, Z)
  
  # Obtain compliance statuses
  S.0 = gen_compliance(Z.obs.matched.caliper.0) 
  S.1 = gen_compliance(Z.obs.matched.caliper.1) 
  S.2 = gen_compliance(Z.obs.matched.caliper.2) 
  S.3 = gen_compliance(Z.obs.matched.caliper.3) 
  S.4 = gen_compliance(Z.obs.matched.caliper.4) 
  
  # Compute compliance rates 
  iota.C0 = calc_compl_rate(S.0)
  iota.C1 = calc_compl_rate(S.1)
  iota.C2 = calc_compl_rate(S.2)
  iota.C3 = calc_compl_rate(S.3)
  iota.C4 = calc_compl_rate(S.4)
  
  # Compute average IV dose difference
  ivdiff.0 = compute_IV_dose_avg_diff(Z.obs.matched.caliper.0)
  ivdiff.1 = compute_IV_dose_avg_diff(Z.obs.matched.caliper.1)
  ivdiff.2 = compute_IV_dose_avg_diff(Z.obs.matched.caliper.2)
  ivdiff.3 = compute_IV_dose_avg_diff(Z.obs.matched.caliper.3)
  ivdiff.4 = compute_IV_dose_avg_diff(Z.obs.matched.caliper.4)
  
  # Generate d
  d.0 = gen_poten_trt(S.0)
  d.1 = gen_poten_trt(S.1)
  d.2 = gen_poten_trt(S.2)
  d.3 = gen_poten_trt(S.3)
  d.4 = gen_poten_trt(S.4)
  
  # Generate r
  r.0 = gen_poten_out_enc(S.0, r_d)
  r.1 = gen_poten_out_enc(S.1, r_d)
  r.2 = gen_poten_out_enc(S.2, r_d)
  r.3 = gen_poten_out_enc(S.3, r_d)
  r.4 = gen_poten_out_enc(S.4, r_d)
  
  # Generate pi
  pi.0.001 = gen_pi(Z.obs.matched.caliper.0, gamma = 0.001)
  pi.0.005 = gen_pi(Z.obs.matched.caliper.0, gamma = 0.005)
  pi.0.01 = gen_pi(Z.obs.matched.caliper.0, gamma = 0.01)
  pi.0.02 = gen_pi(Z.obs.matched.caliper.0, gamma = 0.02)

  pi.1.001 = gen_pi(Z.obs.matched.caliper.1, gamma = 0.001)
  pi.1.005 = gen_pi(Z.obs.matched.caliper.1, gamma = 0.005)
  pi.1.01 = gen_pi(Z.obs.matched.caliper.1, gamma = 0.01)
  pi.1.02 = gen_pi(Z.obs.matched.caliper.1, gamma = 0.02)

  pi.2.001 = gen_pi(Z.obs.matched.caliper.2, gamma = 0.001)
  pi.2.005 = gen_pi(Z.obs.matched.caliper.2, gamma = 0.005)
  pi.2.01 = gen_pi(Z.obs.matched.caliper.2, gamma = 0.01)
  pi.2.02 = gen_pi(Z.obs.matched.caliper.2, gamma = 0.02) 
  
  pi.3.001 = gen_pi(Z.obs.matched.caliper.3, gamma = 0.001)
  pi.3.005 = gen_pi(Z.obs.matched.caliper.3, gamma = 0.005)
  pi.3.01 = gen_pi(Z.obs.matched.caliper.3, gamma = 0.01)
  pi.3.02 = gen_pi(Z.obs.matched.caliper.3, gamma = 0.02)

  pi.4.001 = gen_pi(Z.obs.matched.caliper.4, gamma = 0.001)
  pi.4.005 = gen_pi(Z.obs.matched.caliper.4, gamma = 0.005)
  pi.4.01 = gen_pi(Z.obs.matched.caliper.4, gamma = 0.01)
  pi.4.02 = gen_pi(Z.obs.matched.caliper.4, gamma = 0.02)
  
  # Generate Z
  Z.obs.unstacked.m0.001 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.0, pi.0.001)
  Z.obs.unstacked.m0.005 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.0, pi.0.005)
  Z.obs.unstacked.m0.01 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.0, pi.0.01)
  Z.obs.unstacked.m0.02 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.0, pi.0.02)
  
  Z.ind.m0.001 = create_binary_Z(Z.obs.unstacked.m0.001)
  Z.ind.m0.005 = create_binary_Z(Z.obs.unstacked.m0.005)
  Z.ind.m0.01 = create_binary_Z(Z.obs.unstacked.m0.01)
  Z.ind.m0.02 = create_binary_Z(Z.obs.unstacked.m0.02)
  
  Z.obs.stacked.m0.001 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m0.001)
  Z.obs.stacked.m0.005 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m0.005)
  Z.obs.stacked.m0.01 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m0.01)
  Z.obs.stacked.m0.02 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m0.02)
  
  Z.obs.unstacked.m1.001 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.1, pi.1.001)
  Z.obs.unstacked.m1.005 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.1, pi.1.005)
  Z.obs.unstacked.m1.01 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.1, pi.1.01)
  Z.obs.unstacked.m1.02 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.1, pi.1.02)

  Z.ind.m1.001 = create_binary_Z(Z.obs.unstacked.m1.001)
  Z.ind.m1.005 = create_binary_Z(Z.obs.unstacked.m1.005)
  Z.ind.m1.01 = create_binary_Z(Z.obs.unstacked.m1.01)
  Z.ind.m1.02 = create_binary_Z(Z.obs.unstacked.m1.02)

  Z.obs.stacked.m1.001 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m1.001)
  Z.obs.stacked.m1.005 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m1.005)
  Z.obs.stacked.m1.01 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m1.01)
  Z.obs.stacked.m1.02 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m1.02)

  Z.obs.unstacked.m2.001 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.2, pi.2.001)
  Z.obs.unstacked.m2.005 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.2, pi.2.005)
  Z.obs.unstacked.m2.01 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.2, pi.2.01)
  Z.obs.unstacked.m2.02 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.2, pi.2.02)
  
  Z.ind.m2.001 = create_binary_Z(Z.obs.unstacked.m2.001)
  Z.ind.m2.005 = create_binary_Z(Z.obs.unstacked.m2.005)
  Z.ind.m2.01 = create_binary_Z(Z.obs.unstacked.m2.01)
  Z.ind.m2.02 = create_binary_Z(Z.obs.unstacked.m2.02)

  Z.obs.stacked.m2.001 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m2.001)
  Z.obs.stacked.m2.005 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m2.005)
  Z.obs.stacked.m2.01 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m2.01)
  Z.obs.stacked.m2.02 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m2.02)
  
  Z.obs.unstacked.m3.001 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.3, pi.3.001)
  Z.obs.unstacked.m3.005 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.3, pi.3.005)
  Z.obs.unstacked.m3.01 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.3, pi.3.01)
  Z.obs.unstacked.m3.02 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.3, pi.3.02)

  Z.ind.m3.001 = create_binary_Z(Z.obs.unstacked.m3.001)
  Z.ind.m3.005 = create_binary_Z(Z.obs.unstacked.m3.005)
  Z.ind.m3.01 = create_binary_Z(Z.obs.unstacked.m3.01)
  Z.ind.m3.02 = create_binary_Z(Z.obs.unstacked.m3.02)
  
  Z.obs.stacked.m3.001 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m3.001)
  Z.obs.stacked.m3.005 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m3.005)
  Z.obs.stacked.m3.01 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m3.01)
  Z.obs.stacked.m3.02 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m3.02)
  
  Z.obs.unstacked.m4.001 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.4, pi.4.001)
  Z.obs.unstacked.m4.005 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.4, pi.4.005)
  Z.obs.unstacked.m4.01 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.4, pi.4.01)
  Z.obs.unstacked.m4.02 = gen_obs_IV_dose_unstacked(Z.obs.matched.caliper.4, pi.4.02)
  
  Z.ind.m4.001 = create_binary_Z(Z.obs.unstacked.m4.001)
  Z.ind.m4.005 = create_binary_Z(Z.obs.unstacked.m4.005)
  Z.ind.m4.01 = create_binary_Z(Z.obs.unstacked.m4.01)
  Z.ind.m4.02 = create_binary_Z(Z.obs.unstacked.m4.02)
  
  Z.obs.stacked.m4.001 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m4.001)
  Z.obs.stacked.m4.005 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m4.005)
  Z.obs.stacked.m4.01 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m4.01)
  Z.obs.stacked.m4.02 = gen_obs_IV_dose_stacked(Z.obs.unstacked.m4.02)
  
  # Generate D
  D.m0.001 = gen_obs_trt(Z.ind.m0.001, d.0)
  D.m0.005 = gen_obs_trt(Z.ind.m0.005, d.0)
  D.m0.01 = gen_obs_trt(Z.ind.m0.01, d.0)
  D.m0.02 = gen_obs_trt(Z.ind.m0.02, d.0)
  
  D.m1.001 = gen_obs_trt(Z.ind.m1.001, d.1)
  D.m1.005 = gen_obs_trt(Z.ind.m1.005, d.1)
  D.m1.01 = gen_obs_trt(Z.ind.m1.01, d.1)
  D.m1.02 = gen_obs_trt(Z.ind.m1.02, d.1)
  
  D.m2.001 = gen_obs_trt(Z.ind.m2.001, d.2)
  D.m2.005 = gen_obs_trt(Z.ind.m2.005, d.2)
  D.m2.01 = gen_obs_trt(Z.ind.m2.01, d.2)
  D.m2.02 = gen_obs_trt(Z.ind.m2.02, d.2)
  
  D.m3.001 = gen_obs_trt(Z.ind.m3.001, d.3)
  D.m3.005 = gen_obs_trt(Z.ind.m3.005, d.3)
  D.m3.01 = gen_obs_trt(Z.ind.m3.01, d.3)
  D.m3.02 = gen_obs_trt(Z.ind.m3.02, d.3)
  
  D.m4.001 = gen_obs_trt(Z.ind.m4.001, d.4)
  D.m4.005 = gen_obs_trt(Z.ind.m4.005, d.4)
  D.m4.01 = gen_obs_trt(Z.ind.m4.01, d.4)
  D.m4.02 = gen_obs_trt(Z.ind.m4.02, d.4)
  
  # Generate R
  R.m0.001 = gen_obs_outcome(Z.ind.m0.001, r.0)
  R.m0.005 = gen_obs_outcome(Z.ind.m0.005, r.0)
  R.m0.01 = gen_obs_outcome(Z.ind.m0.01, r.0)
  R.m0.02 = gen_obs_outcome(Z.ind.m0.02, r.0)
  
  R.m1.001 = gen_obs_outcome(Z.ind.m1.001, r.1)
  R.m1.005 = gen_obs_outcome(Z.ind.m1.005, r.1)
  R.m1.01 = gen_obs_outcome(Z.ind.m1.01, r.1)
  R.m1.02 = gen_obs_outcome(Z.ind.m1.02, r.1)

  R.m2.001 = gen_obs_outcome(Z.ind.m2.001, r.2)
  R.m2.005 = gen_obs_outcome(Z.ind.m2.005, r.2)
  R.m2.01 = gen_obs_outcome(Z.ind.m2.01, r.2)
  R.m2.02 = gen_obs_outcome(Z.ind.m2.02, r.2)
  
  R.m3.001 = gen_obs_outcome(Z.ind.m3.001, r.3)
  R.m3.005 = gen_obs_outcome(Z.ind.m3.005, r.3)
  R.m3.01 = gen_obs_outcome(Z.ind.m3.01, r.3)
  R.m3.02 = gen_obs_outcome(Z.ind.m3.02, r.3)

  R.m4.001 = gen_obs_outcome(Z.ind.m4.001, r.4)
  R.m4.005 = gen_obs_outcome(Z.ind.m4.005, r.4)
  R.m4.01 = gen_obs_outcome(Z.ind.m4.01, r.4)
  R.m4.02 = gen_obs_outcome(Z.ind.m4.02, r.4)
  
  # Create data
  df.m0.001 = create_data(Z.ind.m0.001, D.m0.001, R.m0.001)
  df.m0.005 = create_data(Z.ind.m0.005, D.m0.005, R.m0.005)
  df.m0.01 = create_data(Z.ind.m0.01, D.m0.01, R.m0.01)
  df.m0.02 = create_data(Z.ind.m0.02, D.m0.02, R.m0.02)
  
  df.m1.001 = create_data(Z.ind.m1.001, D.m1.001, R.m1.001)
  df.m1.005 = create_data(Z.ind.m1.005, D.m1.005, R.m1.005)
  df.m1.01 = create_data(Z.ind.m1.01, D.m1.01, R.m1.01)
  df.m1.02 = create_data(Z.ind.m1.02, D.m1.02, R.m1.02)

  df.m2.001 = create_data(Z.ind.m2.001, D.m2.001, R.m2.001)
  df.m2.005 = create_data(Z.ind.m2.005, D.m2.005, R.m2.005)
  df.m2.01 = create_data(Z.ind.m2.01, D.m2.01, R.m2.01)
  df.m2.02 = create_data(Z.ind.m2.02, D.m2.02, R.m2.02)

  df.m3.001 = create_data(Z.ind.m3.001, D.m3.001, R.m3.001)
  df.m3.005 = create_data(Z.ind.m3.005, D.m3.005, R.m3.005)
  df.m3.01 = create_data(Z.ind.m3.01, D.m3.01, R.m3.01)
  df.m3.02 = create_data(Z.ind.m3.02, D.m3.02, R.m3.02)

  df.m4.001 = create_data(Z.ind.m4.001, D.m4.001, R.m4.001)
  df.m4.005 = create_data(Z.ind.m4.005, D.m4.005, R.m4.005)
  df.m4.01 = create_data(Z.ind.m4.01, D.m4.01, R.m4.01)
  df.m4.02 = create_data(Z.ind.m4.02, D.m4.02, R.m4.02)
  
  K = c(4,6)
  
  ###############################
  # MATCH 0: NO IV DOSE CALIPER #
  ###############################
  ### gamma = 0.001
  CI.width.m0.001 = compute_width_CI_PI(gamma = 0.001, Z.obs.matched.caliper.0, pi.0.001, K, df.m0.001)
  Gamma.m0.001 = compute_IV_odds(gamma = 0.001, Z.obs.unstacked.m0.001)
  Gamma.m0.001.max = max(Gamma.m0.001)
  Gamma.m0.001.med = median(Gamma.m0.001)

  ### gamma = 0.005
  CI.width.m0.005 = compute_width_CI_PI(gamma = 0.005, Z.obs.matched.caliper.0, pi.0.005, K, df.m0.005)
  Gamma.m0.005 = compute_IV_odds(gamma = 0.005, Z.obs.unstacked.m0.005)
  Gamma.m0.005.max = max(Gamma.m0.005)
  Gamma.m0.005.med = median(Gamma.m0.005)

  ### gamma = 0.01
  CI.width.m0.01 = compute_width_CI_PI(gamma = 0.01, Z.obs.matched.caliper.0, pi.0.01, K, df.m0.01)
  Gamma.m0.01 = compute_IV_odds(gamma = 0.01, Z.obs.unstacked.m0.01)
  Gamma.m0.01.max = max(Gamma.m0.01)
  Gamma.m0.01.med = median(Gamma.m0.01)

  ### gamma = 0.02
  CI.width.m0.02 = compute_width_CI_PI(gamma = 0.02, Z.obs.matched.caliper.0, pi.0.02, K, df.m0.02)
  Gamma.m0.02 = compute_IV_odds(gamma = 0.02, Z.obs.unstacked.m0.02)
  Gamma.m0.02.max = max(Gamma.m0.02)
  Gamma.m0.02.med = median(Gamma.m0.02)
  
  ##################################
  # MATCH 1: SMALL IV DOSE CALIPER #
  ##################################
  ### gamma = 0.001
  CI.width.m1.001 = compute_width_CI_PI(gamma = 0.001, Z.obs.matched.caliper.1, pi.1.001, K, df.m1.001)
  Gamma.m1.001 = compute_IV_odds(gamma = 0.001, Z.obs.unstacked.m1.001)
  Gamma.m1.001.max = max(Gamma.m1.001)
  Gamma.m1.001.med = median(Gamma.m1.001)
  
  ### gamma = 0.005
  CI.width.m1.005 = compute_width_CI_PI(gamma = 0.005, Z.obs.matched.caliper.1, pi.1.005, K, df.m1.005)
  Gamma.m1.005 = compute_IV_odds(gamma = 0.005, Z.obs.unstacked.m1.005)
  Gamma.m1.005.max = max(Gamma.m1.005)
  Gamma.m1.005.med = median(Gamma.m1.005)
  
  ### gamma = 0.01
  CI.width.m1.01 = compute_width_CI_PI(gamma = 0.01, Z.obs.matched.caliper.1, pi.1.01, K, df.m1.01)
  Gamma.m1.01 = compute_IV_odds(gamma = 0.01, Z.obs.unstacked.m1.01)
  Gamma.m1.01.max = max(Gamma.m1.01)
  Gamma.m1.01.med = median(Gamma.m1.01)
  
  ### gamma = 0.02
  CI.width.m1.02 = compute_width_CI_PI(gamma = 0.02, Z.obs.matched.caliper.1, pi.1.02, K, df.m1.02)
  Gamma.m1.02 = compute_IV_odds(gamma = 0.02, Z.obs.unstacked.m1.02)
  Gamma.m1.02.max = max(Gamma.m1.02)
  Gamma.m1.02.med = median(Gamma.m1.02)
  
  ##################################
  # MATCH 2: LARGE IV DOSE CALIPER #
  ##################################
  ### gamma = 0.001
  CI.width.m2.001 = compute_width_CI_PI(gamma = 0.001, Z.obs.matched.caliper.2, pi.2.001, K, df.m2.001)
  Gamma.m2.001 = compute_IV_odds(gamma = 0.001, Z.obs.unstacked.m2.001)
  Gamma.m2.001.max = max(Gamma.m2.001)
  Gamma.m2.001.med = median(Gamma.m2.001)

  ### gamma = 0.005
  CI.width.m2.005 = compute_width_CI_PI(gamma = 0.005, Z.obs.matched.caliper.2, pi.2.005, K, df.m2.005)
  Gamma.m2.005 = compute_IV_odds(gamma = 0.005, Z.obs.unstacked.m2.005)
  Gamma.m2.005.max = max(Gamma.m2.005)
  Gamma.m2.005.med = median(Gamma.m2.005)

  ### gamma = 0.01
  CI.width.m2.01 = compute_width_CI_PI(gamma = 0.01, Z.obs.matched.caliper.2, pi.2.01, K, df.m2.01)
  Gamma.m2.01 = compute_IV_odds(gamma = 0.01, Z.obs.unstacked.m2.01)
  Gamma.m2.01.max = max(Gamma.m2.01)
  Gamma.m2.01.med = median(Gamma.m2.01)

  ### gamma = 0.02
  CI.width.m2.02 = compute_width_CI_PI(gamma = 0.02, Z.obs.matched.caliper.2, pi.2.02, K, df.m2.02)
  Gamma.m2.02 = compute_IV_odds(gamma = 0.02, Z.obs.unstacked.m2.02)
  Gamma.m2.02.max = max(Gamma.m2.02)
  Gamma.m2.02.med = median(Gamma.m2.02)

  ####################################
  # MATCH 3: X-LARGE IV DOSE CALIPER #
  ####################################
  ### gamma = 0.001
  CI.width.m3.001 = compute_width_CI_PI(gamma = 0.001, Z.obs.matched.caliper.3, pi.3.001, K, df.m3.001)
  Gamma.m3.001 = compute_IV_odds(gamma = 0.001, Z.obs.unstacked.m3.001)
  Gamma.m3.001.max = max(Gamma.m3.001)
  Gamma.m3.001.med = median(Gamma.m3.001)
  
  ### gamma = 0.005
  CI.width.m3.005 = compute_width_CI_PI(gamma = 0.005, Z.obs.matched.caliper.3, pi.3.005, K, df.m3.005)
  Gamma.m3.005 = compute_IV_odds(gamma = 0.005, Z.obs.unstacked.m3.005)
  Gamma.m3.005.max = max(Gamma.m3.005)
  Gamma.m3.005.med = median(Gamma.m3.005)
  
  ### gamma = 0.01
  CI.width.m3.01 = compute_width_CI_PI(gamma = 0.01, Z.obs.matched.caliper.3, pi.3.01, K, df.m3.01)
  Gamma.m3.01 = compute_IV_odds(gamma = 0.01, Z.obs.unstacked.m3.01)
  Gamma.m3.01.max = max(Gamma.m1.01)
  Gamma.m3.01.med = median(Gamma.m1.01)
  
  ### gamma = 0.02
  CI.width.m3.02 = compute_width_CI_PI(gamma = 0.02, Z.obs.matched.caliper.3, pi.3.02, K, df.m3.02)
  Gamma.m3.02 = compute_IV_odds(gamma = 0.02, Z.obs.unstacked.m3.02)
  Gamma.m3.02.max = max(Gamma.m3.02)
  Gamma.m3.02.med = median(Gamma.m3.02)
  
  ####################################
  # MATCH 4: XX-LARGE IV DOSE CALIPER #
  ####################################
  ### gamma = 0.001
  CI.width.m4.001 = compute_width_CI_PI(gamma = 0.001, Z.obs.matched.caliper.4, pi.4.001, K, df.m4.001)
  Gamma.m4.001 = compute_IV_odds(gamma = 0.001, Z.obs.unstacked.m4.001)
  Gamma.m4.001.max = max(Gamma.m4.001)
  Gamma.m4.001.med = median(Gamma.m4.001)
  
  ### gamma = 0.005
  CI.width.m4.005 = compute_width_CI_PI(gamma = 0.005, Z.obs.matched.caliper.4, pi.4.005, K, df.m4.005)
  Gamma.m4.005 = compute_IV_odds(gamma = 0.005, Z.obs.unstacked.m4.005)
  Gamma.m4.005.max = max(Gamma.m4.005)
  Gamma.m4.005.med = median(Gamma.m4.005)
  
  ### gamma = 0.01
  CI.width.m4.01 = compute_width_CI_PI(gamma = 0.01, Z.obs.matched.caliper.4, pi.4.01, K, df.m4.01)
  Gamma.m4.01 = compute_IV_odds(gamma = 0.01, Z.obs.unstacked.m4.01)
  Gamma.m4.01.max = max(Gamma.m4.01)
  Gamma.m4.01.med = median(Gamma.m4.01)
  
  ### gamma = 0.02
  CI.width.m4.02 = compute_width_CI_PI(gamma = 0.02, Z.obs.matched.caliper.4, pi.4.02, K, df.m4.02)
  Gamma.m4.02 = compute_IV_odds(gamma = 0.02, Z.obs.unstacked.m4.02)
  Gamma.m4.02.max = max(Gamma.m4.02)
  Gamma.m4.02.med = median(Gamma.m4.02)
  
  # Record outputs
  m0.output = data.frame(iota = c(iota.C0, iota.C0, iota.C0, iota.C0),
                         ci.width = c(CI.width.m0.001, CI.width.m0.005, CI.width.m0.01, CI.width.m0.02),
                         ivdiff = c(ivdiff.0, ivdiff.0, ivdiff.0, ivdiff.0),
                         Gamma.max = c(Gamma.m0.001.max, Gamma.m0.005.max, Gamma.m0.01.max, Gamma.m0.02.max),
                         Gamma.median = c(Gamma.m0.001.med, Gamma.m0.005.med, Gamma.m0.01.med, Gamma.m0.02.med))

  m1.output = data.frame(iota = c(iota.C1, iota.C1, iota.C1, iota.C1),
                         ci.width = c(CI.width.m1.001, CI.width.m1.005, CI.width.m1.01, CI.width.m1.02),
                         ivdiff = c(ivdiff.1, ivdiff.1, ivdiff.1, ivdiff.1),
                         Gamma.max = c(Gamma.m1.001.max, Gamma.m1.005.max, Gamma.m1.01.max, Gamma.m1.02.max),
                         Gamma.median = c(Gamma.m1.001.med, Gamma.m1.005.med, Gamma.m1.01.med, Gamma.m1.02.med))
  
  m2.output = data.frame(iota = c(iota.C2, iota.C2, iota.C2, iota.C2),
                         ci.width = c(CI.width.m2.001, CI.width.m2.005, CI.width.m2.01, CI.width.m2.02),
                         ivdiff = c(ivdiff.2, ivdiff.2, ivdiff.2, ivdiff.2),
                         Gamma.max = c(Gamma.m2.001.max, Gamma.m2.005.max, Gamma.m2.01.max, Gamma.m2.02.max),
                         Gamma.median = c(Gamma.m2.001.med, Gamma.m2.005.med, Gamma.m2.01.med, Gamma.m2.02.med))

  m3.output = data.frame(iota = c(iota.C3, iota.C3, iota.C3, iota.C3),
                         ci.width = c(CI.width.m3.001, CI.width.m3.005, CI.width.m3.01, CI.width.m3.02),
                         ivdiff = c(ivdiff.3, ivdiff.3, ivdiff.3, ivdiff.3),
                         Gamma.max = c(Gamma.m3.001.max, Gamma.m3.005.max, Gamma.m3.01.max, Gamma.m3.02.max),
                         Gamma.median = c(Gamma.m3.001.med, Gamma.m3.005.med, Gamma.m3.01.med, Gamma.m3.02.med))
  
  m4.output = data.frame(iota = c(iota.C4, iota.C4, iota.C4, iota.C4),
                         ci.width = c(CI.width.m4.001, CI.width.m4.005, CI.width.m4.01, CI.width.m4.02),
                         ivdiff = c(ivdiff.4, ivdiff.4, ivdiff.4, ivdiff.4),
                         Gamma.max = c(Gamma.m4.001.max, Gamma.m4.005.max, Gamma.m4.01.max, Gamma.m4.02.max),
                         Gamma.median = c(Gamma.m4.001.med, Gamma.m4.005.med, Gamma.m4.01.med, Gamma.m4.02.med))
  
  return(list(m0.output, m1.output, m2.output, m3.output, m4.output))
}

one_iter(1000)

# Perform simulation
sim.results = replicate(n = 1000, one_iter(1000), simplify = T) 

# Extract results
m0.iota = 0
m0.ivdiff = 0
m0.001.width = 0
m0.005.width = 0
m0.01.width = 0
m0.02.width = 0
m0.001.Gamma.max = 0
m0.001.Gamma.median = 0
m0.005.Gamma.max = 0
m0.005.Gamma.median = 0
m0.01.Gamma.max = 0
m0.01.Gamma.median = 0
m0.02.Gamma.max = 0
m0.02.Gamma.median = 0

m1.iota = 0
m1.ivdiff = 0
m1.001.width = 0
m1.005.width = 0
m1.01.width = 0
m1.02.width = 0
m1.001.Gamma.max = 0
m1.001.Gamma.median = 0
m1.005.Gamma.max = 0
m1.005.Gamma.median = 0
m1.01.Gamma.max = 0
m1.01.Gamma.median = 0
m1.02.Gamma.max = 0
m1.02.Gamma.median = 0

m2.iota = 0
m2.ivdiff = 0
m2.001.width = 0
m2.005.width = 0
m2.01.width = 0
m2.02.width = 0
m2.001.Gamma.max = 0
m2.001.Gamma.median = 0
m2.005.Gamma.max = 0
m2.005.Gamma.median = 0
m2.01.Gamma.max = 0
m2.01.Gamma.median = 0
m2.02.Gamma.max = 0
m2.02.Gamma.median = 0

m3.iota = 0
m3.ivdiff = 0
m3.001.width = 0
m3.005.width = 0
m3.01.width = 0
m3.02.width = 0
m3.001.Gamma.max = 0
m3.001.Gamma.median = 0
m3.005.Gamma.max = 0
m3.005.Gamma.median = 0
m3.01.Gamma.max = 0
m3.01.Gamma.median = 0
m3.02.Gamma.max = 0
m3.02.Gamma.median = 0

m4.iota = 0
m4.ivdiff = 0
m4.001.width = 0
m4.005.width = 0
m4.01.width = 0
m4.02.width = 0
m4.001.Gamma.max = 0
m4.001.Gamma.median = 0
m4.005.Gamma.max = 0
m4.005.Gamma.median = 0
m4.01.Gamma.max = 0
m4.01.Gamma.median = 0
m4.02.Gamma.max = 0
m4.02.Gamma.median = 0

for (i in 1:1000) {
  m0.iota = m0.iota + sim.results[ ,i][[1]][1,1]
  m1.iota = m1.iota + sim.results[ ,i][[2]][1,1]
  m2.iota = m2.iota + sim.results[ ,i][[3]][1,1]
  m3.iota = m3.iota + sim.results[ ,i][[4]][1,1]
  m4.iota = m4.iota + sim.results[ ,i][[5]][1,1]
  
  m0.001.width = m0.001.width + sim.results[ ,i][[1]][1,2]
  m0.005.width = m0.005.width + sim.results[ ,i][[1]][2,2]
  m0.01.width = m0.01.width + sim.results[ ,i][[1]][3,2]
  m0.02.width = m0.02.width + sim.results[ ,i][[1]][4,2]
  
  m0.001.Gamma.max = m0.001.Gamma.max + sim.results[ ,i][[1]][1,4]
  m0.001.Gamma.median = m0.001.Gamma.median + sim.results[ ,i][[1]][1,5]

  m0.005.Gamma.max = m0.005.Gamma.max + sim.results[ ,i][[1]][2,4]
  m0.005.Gamma.median = m0.005.Gamma.median + sim.results[ ,i][[1]][2,5]

  m0.01.Gamma.max = m0.01.Gamma.max + sim.results[ ,i][[1]][3,4]
  m0.01.Gamma.median = m0.01.Gamma.median + sim.results[ ,i][[1]][3,5]

  m0.02.Gamma.max = m0.02.Gamma.max + sim.results[ ,i][[1]][4,4]
  m0.02.Gamma.median = m0.02.Gamma.median + sim.results[ ,i][[1]][4,5]

  m1.001.width = m1.001.width + sim.results[ ,i][[2]][1,2]
  m1.005.width = m1.005.width + sim.results[ ,i][[2]][2,2]
  m1.01.width = m1.01.width + sim.results[ ,i][[2]][3,2]
  m1.02.width = m1.02.width + sim.results[ ,i][[2]][4,2]
  
  m1.001.Gamma.max = m1.001.Gamma.max + sim.results[ ,i][[2]][1,4]
  m1.001.Gamma.median = m1.001.Gamma.median + sim.results[ ,i][[2]][1,5]
  m1.005.Gamma.max = m1.005.Gamma.max + sim.results[ ,i][[2]][2,4]
  m1.005.Gamma.median = m1.005.Gamma.median + sim.results[ ,i][[2]][2,5]
  m1.01.Gamma.max = m1.01.Gamma.max + sim.results[ ,i][[2]][3,4]
  m1.01.Gamma.median = m1.01.Gamma.median + sim.results[ ,i][[2]][3,5]
  m1.02.Gamma.max = m1.02.Gamma.max + sim.results[ ,i][[2]][4,4]
  m1.02.Gamma.median = m1.02.Gamma.median + sim.results[ ,i][[2]][4,5]
  
  m2.001.width = m2.001.width + sim.results[ ,i][[3]][1,2]
  m2.005.width = m2.005.width + sim.results[ ,i][[3]][2,2]
  m2.01.width = m2.01.width + sim.results[ ,i][[3]][3,2]
  m2.02.width = m2.02.width + sim.results[ ,i][[3]][4,2]
  
  m2.001.Gamma.max = m2.001.Gamma.max + sim.results[ ,i][[3]][1,4]
  m2.001.Gamma.median = m2.001.Gamma.median + sim.results[ ,i][[3]][1,5]

  m2.005.Gamma.max = m2.005.Gamma.max + sim.results[ ,i][[3]][2,4]
  m2.005.Gamma.median = m2.005.Gamma.median + sim.results[ ,i][[3]][2,5]

  m2.01.Gamma.max = m2.01.Gamma.max + sim.results[ ,i][[3]][3,4]
  m2.01.Gamma.median = m2.01.Gamma.median + sim.results[ ,i][[3]][3,5]

  m2.02.Gamma.max = m2.02.Gamma.max + sim.results[ ,i][[3]][4,4]
  m2.02.Gamma.median = m2.02.Gamma.median + sim.results[ ,i][[3]][4,5]

  m3.001.width = m3.001.width + sim.results[ ,i][[4]][1,2]
  m3.005.width = m3.005.width + sim.results[ ,i][[4]][2,2]
  m3.01.width = m3.01.width + sim.results[ ,i][[4]][3,2]
  m3.02.width = m3.02.width + sim.results[ ,i][[4]][4,2]
  
  m3.001.Gamma.max = m3.001.Gamma.max + sim.results[ ,i][[4]][1,4]
  m3.001.Gamma.median = m3.001.Gamma.median + sim.results[ ,i][[4]][1,5]
  m3.005.Gamma.max = m3.005.Gamma.max + sim.results[ ,i][[4]][2,4]
  m3.005.Gamma.median = m3.005.Gamma.median + sim.results[ ,i][[4]][2,5]
  m3.01.Gamma.max = m3.01.Gamma.max + sim.results[ ,i][[4]][3,4]
  m3.01.Gamma.median = m3.01.Gamma.median + sim.results[ ,i][[4]][3,5]
  m3.02.Gamma.max = m3.02.Gamma.max + sim.results[ ,i][[4]][4,4]
  m3.02.Gamma.median = m3.02.Gamma.median + sim.results[ ,i][[4]][4,5]
  
  m4.001.width = m4.001.width + sim.results[ ,i][[5]][1,2]
  m4.005.width = m4.005.width + sim.results[ ,i][[5]][2,2]
  m4.01.width = m4.01.width + sim.results[ ,i][[5]][3,2]
  m4.02.width = m4.02.width + sim.results[ ,i][[5]][4,2]
  
  m4.001.Gamma.max = m4.001.Gamma.max + sim.results[ ,i][[5]][1,4]
  m4.001.Gamma.median = m4.001.Gamma.median + sim.results[ ,i][[5]][1,5]
  m4.005.Gamma.max = m4.005.Gamma.max + sim.results[ ,i][[5]][2,4]
  m4.005.Gamma.median = m4.005.Gamma.median + sim.results[ ,i][[5]][2,5]
  m4.01.Gamma.max = m4.01.Gamma.max + sim.results[ ,i][[5]][3,4]
  m4.01.Gamma.median = m4.01.Gamma.median + sim.results[ ,i][[5]][3,5]
  m4.02.Gamma.max = m4.02.Gamma.max + sim.results[ ,i][[5]][4,4]
  m4.02.Gamma.median = m4.02.Gamma.median + sim.results[ ,i][[5]][4,5]
  
  m0.ivdiff = m0.ivdiff + sim.results[ ,i][[1]][1,3]
  m1.ivdiff = m1.ivdiff + sim.results[ ,i][[2]][1,3]
  m2.ivdiff = m2.ivdiff + sim.results[ ,i][[3]][1,3]
  m3.ivdiff = m3.ivdiff + sim.results[ ,i][[4]][1,3]
  m4.ivdiff = m4.ivdiff + sim.results[ ,i][[5]][1,3]
  
}

compliance.rates = c(m0.iota / 1000, m1.iota / 1000, m2.iota / 1000, m3.iota / 1000, m4.iota / 1000)
m0.widths = c(m0.001.width / 1000, m0.005.width / 1000, m0.01.width / 1000, m0.02.width / 1000)
m1.widths = c(m1.001.width / 1000, m1.005.width / 1000, m1.01.width / 1000, m1.02.width / 1000)
m2.widths = c(m2.001.width / 1000, m2.005.width / 1000, m2.01.width / 1000, m2.02.width / 1000)
m3.widths = c(m3.001.width / 1000, m3.005.width / 1000, m3.01.width / 1000, m3.02.width / 1000)
m4.widths = c(m4.001.width / 1000, m4.005.width / 1000, m4.01.width / 1000, m4.02.width / 1000)

ivdifferences = c(m0.ivdiff / 1000, m1.ivdiff / 1000, m2.ivdiff / 1000, m3.ivdiff / 1000, m4.ivdiff / 1000)

Gamma.max.m0 = c(m0.001.Gamma.max / 1000, m0.005.Gamma.max/ 1000, m0.01.Gamma.max/ 1000, m0.02.Gamma.max/ 1000)
Gamma.median.m0 = c(m0.001.Gamma.median/ 1000, m0.005.Gamma.median/ 1000, m0.01.Gamma.median/ 1000, m0.02.Gamma.median/ 1000)

Gamma.max.m1 = c(m1.001.Gamma.max/ 1000, m1.005.Gamma.max/ 1000, m1.01.Gamma.max/ 1000, m1.02.Gamma.max/ 1000)
Gamma.median.m1 = c(m1.001.Gamma.median/ 1000, m1.005.Gamma.median/ 1000, m1.01.Gamma.median/ 1000, m1.02.Gamma.median/ 1000)

Gamma.max.m2 = c(m2.001.Gamma.max/ 1000, m2.005.Gamma.max/ 1000, m2.01.Gamma.max/ 1000, m2.02.Gamma.max/ 1000)
Gamma.median.m2 = c(m2.001.Gamma.median/ 1000, m2.005.Gamma.median/ 1000, m2.01.Gamma.median/ 1000, m2.02.Gamma.median/ 1000)

Gamma.max.m3 = c(m3.001.Gamma.max/ 1000, m3.005.Gamma.max/ 1000, m3.01.Gamma.max/ 1000, m3.02.Gamma.max/ 1000)
Gamma.median.m3 = c(m3.001.Gamma.median/ 1000, m3.005.Gamma.median/ 1000, m3.01.Gamma.median/ 1000, m3.02.Gamma.median/ 1000)

Gamma.max.m4 = c(m4.001.Gamma.max/ 1000, m4.005.Gamma.max/ 1000, m4.01.Gamma.max/ 1000, m4.02.Gamma.max/ 1000)
Gamma.median.m4 = c(m4.001.Gamma.median/ 1000, m4.005.Gamma.median/ 1000, m4.01.Gamma.median/ 1000, m4.02.Gamma.median/ 1000)

data.frame(compliance.rates, ivdifferences)
data.frame(m0.widths, m1.widths, m2.widths, m3.widths, m4.widths,
           Gamma.max.m0, Gamma.median.m0,
           Gamma.max.m1, Gamma.median.m1,
           Gamma.max.m2, Gamma.median.m2,
           Gamma.max.m3, Gamma.median.m3,
           Gamma.max.m4, Gamma.median.m4)












