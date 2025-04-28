# IMPLEMENT SIMULATION STUDY 1
################################################################################
# Import libraries
library(nbpMatching)
library(dplyr)
library(smd)
################################################################################
set.seed(333)
################################################################################
# Write a function that generates covariates
gen_covariates = function(n) {
  x1 = rnorm(n, mean = 0, sd = 1)
  x2 = rnorm(n, mean = 2, sd = sqrt(5))
  x3 = runif(n, min = 1, max = 3)
  x4 = runif(n, min = -2, max = 0)
  x5 = rbinom(n, size = 1, p = 0.5)
  x = data.frame(x1, x2, x3, x4, x5)
  return(x)
}

# Write a function that generates thresholds
gen_threshold = function(n) {
  threshold = runif(n, min = 20, max = 30)
  return(threshold)
}

# Write a function that generates observed IV dose
gen_Z_obs = function(n, x, perfect.IV = T) {
  if (perfect.IV == T) { # Factor 1: IV is randomized
    Z.obs = runif(n, min = 5, max = 50)
  } else { # Factor 2: IV is correlated with some covariates
    x2 = x[ ,2]
    x3 = x[ ,3]
    noise = runif(n, min = 0, max = 2)
    Z.obs = 4*x2 + 6*x3 + noise
  }
  return(Z.obs)
}

# Write a function that generates observed D for an individual
gen_D_obs_ind = function(Z.obs.ij, Tij) {
  if (Z.obs.ij <= Tij) {
    return(0)
  } else {
    return(1)
  }
}

# Write a function that generates observed D
gen_D_obs = function(Z.obs, threshold) {
  mapply(FUN = gen_D_obs_ind, Z.obs, threshold)
}

# Create a dataframe of potential outcome variables
gen_df_pot = function(x, thresholds, Z.obs, Z.obs.corr) {
  return(cbind(x, thresholds, Z.obs, Z.obs.corr))
}

# Write a function that adds a dose caliper to a distance matrix
add_dose_caliper = function(Z.obs, penalty, caliper, distmat) {
  distmat.Z.obs = abs(outer(Z.obs, Z.obs, "-"))
  penalty.dose.caliper = penalty * (distmat.Z.obs <= caliper)
  distmat = distmat + penalty.dose.caliper
  return(distmat)
}

# Write a function that performs nbp matching
nbp_match = function(x, Z.obs, penalty, caliper.sm, caliper.big) {
  # Construct a distance matrix
  dist.mat.list = gendistance(x)
  dist.mat = dist.mat.list$dist
  # Match 2: Add a small IV dose caliper to the distance matrix
  dist.mat.caliper.sm = add_dose_caliper(Z.obs, penalty, caliper.sm, dist.mat)
  # Match 3: Add a large IV dose caliper to the distance matrix
  dist.mat.caliper.big = add_dose_caliper(Z.obs, penalty, caliper.big, dist.mat)
  
  # Match 1: Do matching without IV dose caliper
  dist.mat.no.caliper = distancematrix(dist.mat.list)
  matching.no.caliper = nonbimatch(dist.mat.no.caliper)$halves
  # Match 2: Do matching with the small IV dose caliper
  dist.mat.caliper.sm = distancematrix(dist.mat.caliper.sm)
  matching.caliper.sm = nonbimatch(dist.mat.caliper.sm)$halves
  # Match 3: Do matching with the large IV dose caliper
  dist.mat.caliper.big = distancematrix(dist.mat.caliper.big)
  matching.caliper.big = nonbimatch(dist.mat.caliper.big)$halves
  
  return(list("matching.no.caliper" = matching.no.caliper, 
              "matching.caliper.sm" = matching.caliper.sm,
              "matching.caliper.big" = matching.caliper.big))
}

# Write a function that extracts the indices of matched pairs
extract_match = function(matches, caliper.idx) { # Caliper.idx takes values in {1,2,3}
  matched.pairs = matches[[caliper.idx]] # Caliper.idx == 1 means no IV dose caliper
  first.idx = matched.pairs$Group1.Row # Caliper.idx == 3 means largest IV dose caliper
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

# Find the compliance status of each matched unit
find_compliance_indi = function(D.obs.min, D.obs.max) {
  Sij = case_when(D.obs.min == 0 & D.obs.max == 1 ~ 1, # complier
                  D.obs.min == 1 & D.obs.max == 1 ~ 2, # always-taker
                  D.obs.min == 0 & D.obs.max == 0 ~ 3) # never-taker
  return(Sij)
}

# Find the compliance status of all the matched units
find_compliance = function(Z.obs.matched, thresholds, matched.idx) {
  first.idx = matched.idx[, 1]
  second.idx = matched.idx[, 2]
  # Extract data for the first unit in each matched pair
  df.first = cbind(first.idx, Z.obs.matched, "thresholds" = thresholds[first.idx])
  # Generate D.obs for the first unit for the smaller IV dose
  Z.obs.min.first = df.first[,2]
  thresholds.first = df.first[,4]
  D.obs.first.min = gen_D_obs(Z.obs.min.first, thresholds.first)
  # Generate D.obs for the first unit for the larger IV dose
  Z.obs.max.first = df.first[,3]
  D.obs.first.max = gen_D_obs(Z.obs.max.first, thresholds.first)
  D.obs.first = cbind("first.idx" = df.first[,1], D.obs.first.min, D.obs.first.max)
  
  # Extract data for the second unit in each matched pair
  df.second = cbind(second.idx, Z.obs.matched, "thresholds" = thresholds[second.idx])
  # Generate D.obs for the second unit for the smaller IV dose
  Z.obs.min.second = df.second[,2]
  thresholds.second = df.second[,4]
  D.obs.second.min = gen_D_obs(Z.obs.min.second, thresholds.second)
  # Generate D.obs for the second unit for the larger IV dose
  Z.obs.max.second = df.second[,3]
  D.obs.second.max = gen_D_obs(Z.obs.max.second, thresholds.second)
  D.obs.second = cbind("second.idx" = df.second[,1], D.obs.second.min, D.obs.second.max)
  
  df = rbind(D.obs.first, D.obs.second)
  colnames(df) = c("unit.idx", "D.obs.min.IV", "D.obs.max.IV")
  # Find the compliance status of each matched unit
  compliance = mapply(FUN = find_compliance_indi, df[,2], df[,3])
  
  return(cbind("idx" = df[,1], compliance))
}

# Write a function that computes compliance rate
calc_compl_rate = function(S) {
  n = nrow(S)
  S.comp = S[,2]
  n.compliers = length(S.comp[S.comp == 1])
  iota.C = n.compliers / n
  return(iota.C)
}

# Generate potential outcome r_d1 based on compliance status for an individual
gen_rij_compliance = function(Sij) {
  if (Sij == 1) { # Unit ij is a complier
    r_d1 = runif(n = 1, min = 2, max = 5)
  } else if (Sij == 2) { # Unit ij is an always-taker
    r_d1 = runif(n = 1, min = 4, max = 6)
  } else { # Unit ij is a never-taker
    r_d1 = runif(n = 1, min = 1, max = 3)
  }
  return(r_d1)
}

# Write a function that generates potential outcomes for all the units
gen_r = function(n, complier = F, S = NA) {
  r_d0 = rnorm(n, mean = 0, sd = 1)
  if (complier == F) { # Potential outcomes are independent of compliance status
    psi = runif(n, min = 4, max = 6)
    r_d1 = r_d0 + psi
  } else { # Potential outcomes depend on compliance status
    r_d1 = sapply(S[ ,2], gen_rij_compliance)
  }
  return(cbind(r_d0, r_d1))
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
  # Order compliance rates by increasing index
  S = S[order(S[ ,1], decreasing = F), ]
  df = mapply(FUN = gen_poten_out_enc_indi, S[ ,2], r_d[ ,2], r_d[ ,1])
  df = t(df)
  colnames(df) = c("r_T", "r_C")
  return(df)
}

# Write a function that generates potential treatments for each unit
gen_poten_trt_indi = function(S_ij) {
  if (S_ij == 1) { # complier
    d_T = 1
    d_C = 0
  } else if (S_ij == 2) { # always-taker
    d_T = 1
    d_C = 1
  } else { # never-taker
    d_T = 0
    d_C = 0
  }
  return(c(d_T, d_C))
}

# Write a function that generates potential treatments for all the units
gen_poten_trt = function(S) {
  # Order compliance rate by increasing index
  S = S[order(S[ ,1], decreasing = F), ]
  df = t(sapply(S[ ,2], FUN = gen_poten_trt_indi))
  colnames(df) = c("d_T", "d_C")
  return(df)
}

# Create a dataframe of potential outcome variables for matched pairs
gen_df_pot_matched = function(df_pot, matched.idx, S, r) {
  first.idx = matched.idx[,1]
  second.idx = matched.idx[,2]
  first.units = df_pot[first.idx, ]
  second.units = df_pot[second.idx, ]
  compliance = S[,2]
  
  df = rbind(first.units, second.units)
  df = cbind(df, compliance, r)
  
  return(df) # first half is matched to second half
}

# Compute lower and upper bounds
compute_bounds = function(r_enc, d, K) {
  r_T = r_enc[ ,1]
  r_C = r_enc[ ,2]
  d_T = d[ ,1]
  d_C = d[ ,2]
  n = nrow(d)
  K0 = K[1]
  K1 = K[2]
  
  LB = (sum(r_T - r_C) / n) - K0 * (sum(d_T - d_C) / n) + K0
  UB = (sum(r_T - r_C) / n) - K1 * (sum(d_T - d_C) / n) + K1
  return(c(LB, UB))
}

# Find the difference UB - LB
diff_in_bounds = function(bounds) {
  LB = bounds[1]
  UB = bounds[2]
  return(UB - LB)
}

# Based on matching, create remaining data and compute bounds
dgp_bounds = function(matched.idx, Z.obs, thresholds, r_d, K) {
  # Extract IV doses of matched units
  Z.obs.matched = extract_IV_doses(matched.idx, Z.obs)
  # Find compliance status of matched units
  S = find_compliance(Z.obs.matched, thresholds, matched.idx)
  # Find potential outcomes based on encouragement of treatment or control
  r_enc = gen_poten_out_enc(S, r_d)
  d = gen_poten_trt(S)
  
  # Compute bounds
  bounds = compute_bounds(r_enc, d, K)
  # Compute difference between upper and lower bounds
  diff.in.bounds = diff_in_bounds(bounds)
  
  return(diff.in.bounds)
  #return(list(Z.obs.matched = Z.obs.matched, S = S, r_enc = r_enc, d = d))
}

# Initialize empty dataframes
### M1 ###
m1.f1.df = data.frame()
m1.f2.df = data.frame()
m1.f3.df = data.frame()
m1.f4.df = data.frame()

### M2 ###
m2.f1.df = data.frame()
m2.f2.df = data.frame()
m2.f3.df = data.frame()
m2.f4.df = data.frame()

### M3 ###
m3.f1.df = data.frame()
m3.f2.df = data.frame()
m3.f3.df = data.frame()
m3.f4.df = data.frame()

# Write a function that performs one iteration of the simulation study
one_iter = function(n, penalty) {
  # Generate data
  x = gen_covariates(n)
  thresholds = gen_threshold(n)
  # Factor 1
  Z.obs = gen_Z_obs(n, x, perfect.IV = T) # 1(a): perfectly randomized IV
  Z.obs.corr = gen_Z_obs(n, x, perfect.IV = F) # 1(b): IV correlated with covariates
  # Factor 2
  r.no.comp = gen_r(n, complier = F) # 2(a): r_d independent of compliance status
  
  # Generate preliminary potential outcomes data
  df_pot = gen_df_pot(x, thresholds, Z.obs, Z.obs.corr)
  
  # Perform matching
  ## Perfect IV 
  matches = nbp_match(x, Z.obs, penalty, caliper.sm = 7, caliper.big = 15)
  matched.idx.no.caliper = extract_match(matches, 1) # no IV dose caliper
  matched.idx.caliper.sm = extract_match(matches, 2) # small IV dose caliper
  matched.idx.caliper.big = extract_match(matches, 3) # large IV dose caliper
  
  ## Imperfect IV 
  matches.corr = nbp_match(x, Z.obs.corr, penalty, caliper.sm = 7, caliper.big = 15)
  matched.idx.no.caliper.corr = extract_match(matches.corr, 1) # no IV dose caliper
  matched.idx.caliper.sm.corr = extract_match(matches.corr, 2) # small IV dose caliper
  matched.idx.caliper.big.corr = extract_match(matches.corr, 3) # large IV dose caliper
  
  # For each matching design, compute UB - LB
  ############################
  ## M1: no IV dose caliper ##
  ############################
  ### Perfect IV & r_d independent of S
  m1.f1 = dgp_bounds(matched.idx.no.caliper, Z.obs, thresholds, r.no.comp, c(4,6))
  Z.obs.matched.no.caliper = extract_IV_doses(matched.idx.no.caliper, Z.obs) 
  S.no.caliper = find_compliance(Z.obs.matched.no.caliper, thresholds, matched.idx.no.caliper)
  m1.f1.data = gen_df_pot_matched(df_pot, matched.idx.no.caliper, S.no.caliper, r.no.comp)
  m1.f1.df <<- rbind(m1.f1.df, m1.f1.data)
  
  ### Perfect IV & r_d dependent of S
  r.comp = gen_r(n, complier = T, S = S.no.caliper) # 2(b)
  m1.f2 = dgp_bounds(matched.idx.no.caliper, Z.obs, thresholds, r.comp, c(1,6))
  m1.f2.data = gen_df_pot_matched(df_pot, matched.idx.no.caliper, S.no.caliper, r.comp)
  m1.f2.df <<- rbind(m1.f2.df, m1.f2.data)
  
  ### Imperfect IV & r_d independent of S
  m1.f3 = dgp_bounds(matched.idx.no.caliper.corr, Z.obs.corr, thresholds, r.no.comp, c(4,6))
  Z.obs.matched.no.caliper.corr = extract_IV_doses(matched.idx.no.caliper.corr, Z.obs.corr) 
  S.no.caliper.corr = find_compliance(Z.obs.matched.no.caliper.corr, thresholds, matched.idx.no.caliper.corr)
  m1.f3.data = gen_df_pot_matched(df_pot, matched.idx.no.caliper.corr, S.no.caliper.corr, r.no.comp)
  m1.f3.df <<- rbind(m1.f3.df, m1.f3.data)
  
  ### Imperfect IV & r_d dependent of S
  r.comp.corr = gen_r(n, complier = T, S = S.no.caliper.corr)
  m1.f4 = dgp_bounds(matched.idx.no.caliper.corr, Z.obs.corr, thresholds, r.comp.corr, c(1,6))
  m1.f4.data = gen_df_pot_matched(df_pot, matched.idx.no.caliper.corr, S.no.caliper.corr, r.comp.corr)
  m1.f4.df <<- rbind(m1.f4.df, m1.f4.data)
  
  m1 = c(m1.f1, m1.f2, m1.f3, m1.f4)
  
  ###############################
  ## M2: small IV dose caliper ##
  ###############################
  ### Perfect IV & r_d independent of S
  m2.f1 = dgp_bounds(matched.idx.caliper.sm, Z.obs, thresholds, r.no.comp, c(4,6))
  Z.obs.matched.caliper.sm = extract_IV_doses(matched.idx.caliper.sm, Z.obs) 
  S.caliper.sm = find_compliance(Z.obs.matched.caliper.sm, thresholds, matched.idx.caliper.sm)
  m2.f1.data = gen_df_pot_matched(df_pot, matched.idx.caliper.sm, S.caliper.sm, r.no.comp)
  m2.f1.df <<- rbind(m2.f1.df, m2.f1.data)
  
  ### Perfect IV & r_d dependent of S
  r.comp = gen_r(n, complier = T, S = S.caliper.sm)
  m2.f2 = dgp_bounds(matched.idx.caliper.sm, Z.obs, thresholds, r.comp, c(1,6))
  m2.f2.data = gen_df_pot_matched(df_pot, matched.idx.caliper.sm, S.caliper.sm, r.comp)
  m2.f2.df <<- rbind(m2.f2.df, m2.f2.data)
  
  ### Imperfect IV & r_d independent of S
  m2.f3 = dgp_bounds(matched.idx.caliper.sm.corr, Z.obs.corr, thresholds, r.no.comp, c(4,6))
  Z.obs.matched.caliper.sm.corr = extract_IV_doses(matched.idx.caliper.sm.corr, Z.obs.corr) 
  S.caliper.sm.corr = find_compliance(Z.obs.matched.caliper.sm.corr, thresholds, matched.idx.caliper.sm.corr)
  m2.f3.data = gen_df_pot_matched(df_pot, matched.idx.caliper.sm.corr, S.caliper.sm.corr, r.comp)
  m2.f3.df <<- rbind(m2.f3.df, m2.f3.data)
  
  ### Imperfect IV & r_d dependent of S
  r.comp.corr = gen_r(n, complier = T, S = S.caliper.sm.corr)
  m2.f4 = dgp_bounds(matched.idx.caliper.sm.corr, Z.obs.corr, thresholds, r.comp.corr, c(1,6))
  m2.f4.data = gen_df_pot_matched(df_pot, matched.idx.caliper.sm.corr, S.caliper.sm, r.comp.corr)
  m2.f4.df <<- rbind(m2.f4.df, m2.f4.data)
  
  m2 = c(m2.f1, m2.f2, m2.f3, m2.f4)
  
  
  ###############################
  ## M3: large IV dose caliper ##
  ###############################
  ### Perfect IV & r_d independent of S
  m3.f1 = dgp_bounds(matched.idx.caliper.big, Z.obs, thresholds, r.no.comp, c(4,6))
  Z.obs.matched.caliper.big = extract_IV_doses(matched.idx.caliper.big, Z.obs) 
  S.caliper.big = find_compliance(Z.obs.matched.caliper.big, thresholds, matched.idx.caliper.big)
  m3.f1.data = gen_df_pot_matched(df_pot, matched.idx.caliper.big, S.caliper.big, r.no.comp)
  m3.f1.df <<- rbind(m3.f1.df, m3.f1.data)
  
  ### Perfect IV & r_d dependent of S
  
  r.comp = gen_r(n, complier = T, S = S.caliper.big)
  m3.f2 = dgp_bounds(matched.idx.caliper.big, Z.obs, thresholds, r.comp, c(1,6))
  m3.f2.data = gen_df_pot_matched(df_pot, matched.idx.caliper.big, S.caliper.big, r.comp)
  m3.f2.df <<- rbind(m3.f2.df, m3.f2.data)
  
  ### Imperfect IV & r_d independent of S
  m3.f3 = dgp_bounds(matched.idx.caliper.big.corr, Z.obs.corr, thresholds, r.no.comp, c(4,6))
  Z.obs.matched.caliper.big.corr = extract_IV_doses(matched.idx.caliper.big.corr, Z.obs.corr) 
  S.caliper.big.corr = find_compliance(Z.obs.matched.caliper.big.corr, thresholds, matched.idx.caliper.big.corr)
  m3.f3.data = gen_df_pot_matched(df_pot, matched.idx.caliper.big.corr, S.caliper.big.corr, r.no.comp)
  m3.f3.df <<- rbind(m3.f3.df, m3.f3.data)
  
  ### Imperfect IV & r_d dependent of S
  
  r.comp.corr = gen_r(n, complier = T, S = S.caliper.big.corr)
  m3.f4 = dgp_bounds(matched.idx.caliper.big.corr, Z.obs.corr, thresholds, r.comp.corr, c(1,6))
  m3.f4.data = gen_df_pot_matched(df_pot, matched.idx.caliper.big.corr, S.caliper.big.corr, r.comp.corr)
  m3.f4.df <<- rbind(m3.f4.df, m3.f4.data)
  
  m3 = c(m3.f1, m3.f2, m3.f3, m3.f4)
  
  return(cbind(m1, m2, m3))
}

################################################################################
# Perform simulation
sim.results = replicate(n = 1000, one_iter(n=1000, penalty = 10), simplify = F)
sim.table = Reduce("+", sim.results) / length(sim.results)

################################################################################
# COMPUTE the SMD OF MATCHED PAIRS
################################################################################
# Write a function that computes the average SMD of matched pairs
compute_smd = function(df) {
  # Subset data to only include covariates
  df = df[, 2:6]
  # Create ID for matched pairs
  df$id = rep(1:2000, each = 500)
  # Compute average absolute SMD of each covariate across 1,000 simulations
  df.smd = df %>% 
    group_by(grp = as.integer(gl(n(), 1000, n()))) %>%
    summarise_at(
      .vars = vars(dplyr::matches("^x")),
      .funs = list(smd = ~ abs(smd(., g = id)$estimate))) %>% # absolute value of SMD
    summarise_all(mean)
  return(df.smd)
}

h = compute_smd(m1.f1.df)
summary(h[,2])

compute_smd(m1.f2.df)
compute_smd(m1.f3.df)
compute_smd(m1.f4.df)
compute_smd(m2.f1.df)
compute_smd(m2.f2.df)
compute_smd(m2.f3.df)
compute_smd(m2.f4.df)
compute_smd(m3.f1.df)
compute_smd(m3.f2.df)
compute_smd(m3.f3.df)
compute_smd(m3.f4.df)

################################################################################
# COMPUTE SAMPLE AVERAGE TREATMENT EFFECT AMONG COMPLIERS 
################################################################################
# Write a function that computes SATE among compliers
compute_complier_SATE = function(df) {
  SATE = df %>%
    # Divide data according to 1,000 simulations
    group_by(grp = as.integer(gl(n(), 1000, n()))) %>%
    # Only include compliers
    filter(compliance == 1) %>% 
    mutate(r.diff = (r_d1 - r_d0)) %>%
    # Compute complier SATE for each simulation
    summarize(SATE = mean(r.diff)) %>%
    # Compute the average complier SATE across 1,000 simulations
    summarize(SATE.final = mean(SATE))
  return(SATE)
}

compute_complier_SATE(m1.f1.df)
compute_complier_SATE(m1.f2.df)
compute_complier_SATE(m1.f3.df)
compute_complier_SATE(m1.f4.df)
compute_complier_SATE(m2.f1.df)
compute_complier_SATE(m2.f2.df)
compute_complier_SATE(m2.f3.df)
compute_complier_SATE(m2.f4.df)
compute_complier_SATE(m3.f1.df)
compute_complier_SATE(m3.f2.df)
compute_complier_SATE(m3.f3.df)
compute_complier_SATE(m3.f4.df)
