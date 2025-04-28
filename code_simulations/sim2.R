# IMPLEMENT SIMULATION STUDY 2
################################################################################
set.seed(373)
################################################################################
# Import libraries
library(Rlab)
library(dplyr)
library(tidyr)
library(xtable)
library(locfit)
################################################################################
# EXTRACT IV DOSES FROM REAL DATA #
################################################################################
# Read data
m1.age35.long = read.csv("match1_long_age35_2003-4.csv")
m1.age35.short = read.csv("match1_short_age35_2003-4.csv")
# Extract IV doses
Z.max0 = m1.age35.long$dif.travel.time
Z.min0 = m1.age35.short$dif.travel.time
Z0 = c(Z.max0, Z.min0)
################################################################################
# GENERATE POTENTIAL OUTCOMES DATA #
################################################################################
# Write a function that generates IV doses for matched pairs
gen_IV_dose = function(I, Z0) {
  # Sort and divide Z into two non-overlapping sets
  Z0 = sort(Z0, decreasing = T)
  n = length(Z0) / 2
  Z.max0 = Z0[1:n]
  Z.min0 = Z0[(n+1):length(Z0)]
  # Randomly sample minimum continuous IV dose
  Z.min = sample(Z.min0, size = I, replace = T)
  # Randomly sample maximum continuous IV dose
  Z.max = sample(Z.max0, size = I, replace = T)
  
  return(cbind(Z.min, Z.max))
}

# Write a function that generates compliance status for each unit 
gen_compliance = function(Z, beta) {
  # Find the difference between maximum and minimum IV doses
  Z.min = Z[, 1]
  Z.max = Z[, 2]
  Z.diff = Z.max - Z.min
  # Replicate the IV dose for each unit in matched pairs
  Z.diff2 = rep(Z.diff, each = 2)
  # Create a linear predictor that's a function of the difference in IV dose
  lp = beta * Z.diff2
  # Create another linear predictor that's a function of Z.max
  I = length(Z.diff)
  alpha = runif(I, min = 0, max = 0.4)
  lp2 = -alpha * Z.max
  lp2 = rep(lp2, each = 2)
  # Generate probabilities of compliance status
  denom = (1 + exp(lp) + exp(lp2)) 
  prob.1 = 1 / denom 
  prob.2 = exp(lp) / denom 
  prob.3 = exp(lp2) / denom
  prob = cbind(prob.1, prob.2, prob.3) # within-matched-pair units have same probabilities
  
  # Randomly sample compliance status from the said probabilities
  compliance = apply(prob, MARGIN = 1, function(x) sample(x = 1:3, size = 1, prob = x))
  return(data.frame(compliance))
}

# Write a function that computes compliance rate
calc_compl_rate = function(S) {
  n = nrow(S)
  n.compliers = length(S[S$compliance == 1, ])
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

# Generate IV dose assignment probability for each matched pair
gen_pi = function(Z, gamma) { # If gamma == 0, randomization i.e., Prop. 1
  I = nrow(Z)
  Z.min = Z[, 1]
  Z.max = Z[, 2]
  Z.diff = Z.max - Z.min
  exp.term = exp(gamma * Z.diff)
  LB = 0.5 # Mechanism 1
  UB = exp.term / (exp.term + 1)
  # LB = UB # Mechanism 2
  pi = runif(n = I, min = LB, max = UB)
  
  return(data.frame(pi))
}

# Generate potential outcome based on receipt of treatment or control for one unit 
gen_poten_out_agn_indi = function(Sij, cts = T) {
  epsilon = runif(n = 1, min = -1, max = 1) 
  if (cts == T) { # Continuous potential outcomes
    r_d0 = rnorm(n = 1, mean = 0, sd = 1)
    r_d1 = r_d0 + epsilon
  } else { # Binary potential outcomes
    p = expit(0.5)
    r_d0 = rbern(n = 1, p)
    p.r_d1 = expit(0.5 + epsilon)
    r_d1 = rbern(n = 1, p.r_d1)
  }
  return(c(r_d1, r_d0))
}

# Generate potential outcome based on receipt of treatment or control
gen_poten_out_agn = function(S, cts) {
  df = t(apply(S, MARGIN=1, FUN = gen_poten_out_agn_indi, cts = cts))
  colnames(df) = c("r_d1", "r_d0")
  df = data.frame(df)
  df$pairID = rep(1:(nrow(df)/2), each = 2)
  df = df %>%
    group_by(pairID) %>%
    arrange(desc(r_d0), .by_group = T) %>%
    ungroup(pairID) %>%
    select(-pairID)
  df = as.matrix(df)
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
  df = mapply(FUN = gen_poten_out_enc_indi, S$compliance, r_d[ ,1], r_d[ ,2])
  df = t(df)
  colnames(df) = c("r_T", "r_C")
  return(df)
}

# Create a dataframe of potential outcome variables
gen_df_pot = function(Z, pi, S, d, r, r_d) {
  # Replicate Z and pi for each matched unit within a pair
  Z.min = Z[, 1]
  Z.max = Z[, 2]
  Z.min = rep(Z.min, each = 2)
  Z.max = rep(Z.max, each = 2)
  Z = data.frame(cbind(Z.min, Z.max))
  pi = rep(pi$pi, each = 2)
  pi = data.frame(pi)
  
  return(cbind(Z, pi, S, d, r, r_d))
}

################################################################################
# GENERATE OBSERVED OUTCOMES DATA #
################################################################################
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

# Write a function that creates an indicator of which matched unit within a pair has larger IV dose
create_Z_ind = function(Z.obs.unstacked, pi) {
  pi1 = pi[[1]]
  Z.i1.ind = as.numeric(Z.obs.unstacked[,1] >= Z.obs.unstacked[,2]) 
  Z.i2.ind = 1 - Z.i1.ind
  Z.ind = data.frame(cbind(Z.i1.ind, Z.i2.ind))
  return(Z.ind)
}

# Write a function that creates a stacked indicator of which matched unit within a pair has larger IV dose
create_Z_ind_stacked = function(Z.ind) {
  Z.ind.stacked = data.frame(data = c(t(Z.ind)))
  return(Z.ind.stacked)
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

################################################################################
# COMPUTE ESTIMANDS BASED ON POTENTIAL OUTCOMES DATA #
################################################################################
# Write a function that computes the SATE
compute_SATE = function(r_d) {
  n = nrow(r_d)
  r_d1 = r_d[ ,1]
  r_d0 = r_d[ ,2]
  kappa = sum(r_d1 - r_d0) / n
  return(kappa)
}

# Write a function that computes the lower and upper bounds
compute_bounds = function(r, d, K) {
  r_T = r[ ,1]
  r_C = r[ ,2]
  d_T = d[ ,1]
  d_C = d[ ,2]
  n = nrow(d)
  K0 = K[1]
  K1 = K[2]
  
  LB = (sum(r_T - r_C) / n) - K0 * (sum(d_T - d_C) / n) + K0
  UB = (sum(r_T - r_C) / n) - K1 * (sum(d_T - d_C) / n) + K1
  return(c(LB, UB))
}

################################################################################
# COMPUTE ESTIMATES BASED ON OBSERVED DATA: PROPOSITION 1 #
################################################################################
# Write a function that computes the test statistic
compute_test_stat = function(K, R, D, Z.ind.stacked, l) {
  I = nrow(R) / 2
  K0 = K[1]
  K1 = K[2]
  test.stat.lb = (sum(Z.ind.stacked * (R - K0*D) - (1 - Z.ind.stacked) * (R - K0*D)) / I) - (l - K0)
  return(test.stat.lb)
}

# Write a function that computes within-matched-pair ATE
compute_within_pair_ATE = function(Z.ind, K, R, D) {
  n = nrow(R)
  I = n / 2
  K0 = K[1]
  # Extract the first unit in each matched pair i.e., the odd rows
  odd.rows = seq_len(n) %% 2 
  R.first = R[odd.rows == 1, ] 
  D.first = D[odd.rows == 1, ]
  Z.first = Z.ind[ ,1]
  # Extract the second unit in each matched pair i.e., the even rows
  R.second = R[odd.rows == 0, ]
  D.second = D[odd.rows == 0, ]
  Z.second = Z.ind[ ,2]
  
  tau.i.hat = (2*Z.first - 1) * (R.first - K0 * D.first) + (2*Z.second - 1) * (R.second - K0 * D.second)
  return(tau.i.hat)
}

# Write a function that computes the sample variance
compute_sample_var = function(tau.i.hat) {
  I = length(tau.i.hat)
  sample.var = stats::var(tau.i.hat) / I
  return(sample.var)
}

# Write a function that checks whether the confidence interval covers LB (prop.1)
check_ci_cover = function(K, R, D, Z.ind, Z.ind.stacked, l, alpha) {
  test.stat.lb = compute_test_stat(K, R, D, Z.ind.stacked, l)
  tau.i.hat = compute_within_pair_ATE(Z.ind, K, R, D)
  sample.var = compute_sample_var(tau.i.hat)
  z = qnorm(1 - alpha / 2)
  ci.lb = test.stat.lb + l - z * sqrt(sample.var)
  ci.ub = test.stat.lb + l + z * sqrt(sample.var)
  return(as.numeric(ci.lb <= l && ci.ub >= l))
}

################################################################################
# COMPUTE ESTIMATES BASED ON OBSERVED DATA: PROPOSITION 3 #
################################################################################
# Write a function that computes IV dose-dependent odds
compute_IV_odds = function(gamma, Z.obs.unstacked) {
  Z.obs.max = apply(Z.obs.unstacked, FUN = max, MARGIN = 1)
  Z.obs.min = apply(Z.obs.unstacked, FUN = min, MARGIN = 1)
  odds = exp(gamma * (Z.obs.max - Z.obs.min))
  return(odds)
}

# Write a function that computes scaled within-matched-pair ATE
compute_scaled_within_pair_ATE = function(tau.i.hat, Gamma) {
  constant = (Gamma + 1) / (4 * Gamma)
  output = constant * ((Gamma + 1) * tau.i.hat - (Gamma - 1) * abs(tau.i.hat))
  return(output)
}

# Write a function that tests proposition 3
test_prop3 = function(gamma, Z.obs.unstacked, tau.i.hat, alpha, l, K) {
  Gamma = compute_IV_odds(gamma, Z.obs.unstacked)
  tau.i.hat.biased = compute_scaled_within_pair_ATE(tau.i.hat, Gamma)
  tau.bar.biased = mean(tau.i.hat.biased)
  sample.var.biased = compute_sample_var(tau.i.hat.biased)
  z = qnorm(1 - alpha)
  K0 = K[1]
  return(as.numeric(tau.bar.biased - (l - K0) < z * sqrt(sample.var.biased)))
}

################################################################################
# PERFORM SIMULATION #
################################################################################
# Write a function that performs one iteration of the simulation
one_iter = function(I, Z0, beta, gamma, alpha, cts) {
  # Generate potential outcomes data
  Z = gen_IV_dose(I, Z0)
  S = gen_compliance(Z, beta)
  iota.C = calc_compl_rate(S)
  d = gen_poten_trt(S)
  pi = gen_pi(Z, gamma)
  r_d = gen_poten_out_agn(S, cts = cts)
  r = gen_poten_out_enc(S, r_d)
  
  # Generate observed data
  pi1 = pi[[1]]
  Z.obs.unstacked = gen_obs_IV_dose_unstacked(Z, pi1)
  Z.ind = create_Z_ind(Z.obs.unstacked, pi)
  Z.ind.stacked = create_Z_ind_stacked(Z.ind)
  D = gen_obs_trt(Z.ind.stacked, d)
  R = gen_obs_outcome(D, r)
  
  # Compute estimands based on potential outcomes data
  if (cts == T) {
    K = c(-1,1) 
  } else {
    K = c(-1,1) 
  }
  l = compute_bounds(r, d, K)[1]
  # Compute estimates based on observed data (propositions 1 & 3)
  tau.i.hat = compute_within_pair_ATE(Z.ind, K, R, D)
  
  return(c(check_ci_cover(K, R, D, Z.ind, Z.ind.stacked, l, alpha),
           test_prop3(gamma, Z.obs.unstacked, tau.i.hat, alpha, l, K), 
           iota.C))
}

# Perform simulation for a combination of different parameter values
I = c(100, 500, 1000, 2000)
beta = c(-0.2)
gamma = c(0, 0.025, 0.05)
params = expand.grid(I = I, beta = beta, gamma = gamma)

# Write a function that performs simulation study for a combination of parameter values
perform_sim = function(n, Z0, alpha, params, cts) {
  I = params[[1]]
  beta = params[[2]]
  gamma = params[[3]]
  sim.results = replicate(n = n, one_iter(I, Z0, beta, gamma, alpha, cts), simplify = T)
  results = apply(sim.results, MARGIN = 1, mean)
  print(paste("I: ", I, "; beta: ", beta, "; gamma: ", gamma))
  return(results)
}

# Continuous case
sim.results.fin = apply(params, MARGIN = 1, perform_sim, n = 1000, Z0 = Z0, alpha = 0.05, cts = T)
sim.results.fin2 = t(sim.results.fin)
colnames(sim.results.fin2) = c("coverage1", "coverage3", "compliance")
results = cbind(params, sim.results.fin2)
print(xtable(results, digits = 3), include.rownames = F)

# Binary case
sim.results.fin3 = apply(params, MARGIN = 1, perform_sim, n = 1000, Z0 = Z0, alpha = 0.05, cts = F)
sim.results.fin4 = t(sim.results.fin3)
colnames(sim.results.fin4) = c("coverage1", "coverage3", "compliance")
results2 = cbind(params, sim.results.fin4)
print(xtable(results2, digits = 3), include.rownames = F)


