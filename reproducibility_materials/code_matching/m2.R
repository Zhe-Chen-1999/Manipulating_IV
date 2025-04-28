# IMPLEMENT MATCH 2
################################################################################
# Load packages
library(dplyr)
library(nbpMatching)
library(ggplot2)
library(xtable)
################################################################################
# Source helper functions
source('./util_matching.R')
################################################################################
# Load data
data.names = paste("preemies", sep="", 1995:2005, ".csv")
# Set seed
set.seed(321)
################################################################################
# Define helper function
mt = function(df, sink_prop = 0.15, i) {
  df$highrisk = case_when(df$gestage <= 30 & df$ageM >= 35 ~ 1, # highest risk
                          df$gestage <= 30 & df$ageM < 35 ~ 2,
                          df$gestage > 30 & df$ageM >= 35 ~ 3,
                          df$gestage > 30 & df$ageM < 35 ~ 4) # lowest risk
  
  # Collapse "other.race" and "race.mis" variables into one
  df$race.OTHER = ifelse(df$other.race == 1 | df$race.mis == 1, 1, 0)
  
  num.units = dim(df)[1]
  num.sinks = comp_sink_size(num.units, sink_prop)
  template.idx = sample(nrow(df), size = num.sinks, replace = F)
  print(paste("N.sinks for year", i+1994, ":", num.sinks))
  print(paste("N.expected matched pairs for year", i+1994, ":", (num.units-num.sinks)/2))
  
  # Define template data
  template.full = df[sort(template.idx), ]
  rownames(template.full) = NULL # Re-index rows
  
  # Select covariates to match on
  template = template.full %>% dplyr::select(c(bthwght, parity, fee.for.service, below_poverty,
                                               ageM, medu, singlebirth, highrisk))
  template.rest = template.full %>% dplyr::select(-c(bthwght, parity, fee.for.service, below_poverty,
                                                     ageM, medu, singlebirth, highrisk))
  
  # Define observational data
  obs.full = df
  
  # Select covariates to match on
  obs = obs.full %>% dplyr::select(c(bthwght, parity, fee.for.service, below_poverty,
                                     ageM, medu, singlebirth, highrisk))
  obs.rest = obs.full %>% dplyr::select(-c(bthwght, parity, fee.for.service, below_poverty,
                                           ageM, medu, singlebirth, highrisk))
  num.obs = dim(obs)[1]
  IV = obs.full$dif.travel.time
  
  # Construct X
  X = rbind(template, obs)
  X$s = NA
  
  # Create a dummy variable indicating membership in template
  X[1:num.sinks, "s"] = 1
  X[(num.sinks+1): nrow(X), "s"] = 0
  
  # Estimate probabilities of participation
  prob.part = est_prob_part(X)
  
  # Create an array of est. prob. of part. for template data
  s.hat.temp = prob.part[X$s==1]
  
  # Construct a distance matrix btw obs. units & templates
  distmat = construct_dist_mat(obs, template) # Rank-based Maha distance matrix
  
  # Implement generalizability caliper v2
  distmat = gen_caliper2(prob.part, X, s.hat.temp, gamma.mult = 1, D=500,
                         num.obs, distmat)
  
  # Implement dose caliper
  distmat = dose_caliper(IV, C=250, tau=sd(IV), distmat)
  
  # Prevent sinks from being matched to other sinks
  distmat = no_sink_match(distmat, num.obs, max(distmat)*3)
  
  # Do the matching
  distmat = distancematrix(distmat)
  matching = nonbimatch(distmat)$halves
  
  ### Subset matches we want to keep
  pair_kept = which(matching$Group1.Row <= num.obs & matching$Group2.Row <= num.obs)
  matching_keep = matching[matching$Group1.Row <= num.obs & matching$Group2.Row <= num.obs, ]
  
  ### Subset obs.--sinks matches
  pair_discard = which(matching$Group1.Row <= num.obs & matching$Group2.Row > num.obs)
  matching_discard = matching[matching$Group1.Row <= num.obs & matching$Group2.Row > num.obs, ]
  discarded.full = obs.full[matching_discard$Group1.Row, ]
  
  # Subset moms in matched pairs
  obs.pair1 = obs.full[matching_keep$Group1.Row, ]
  obs.pair2 = obs.full[matching_keep$Group2.Row, ]
  
  # Extract IVs within matched pairs
  obs.pair1$IV = ifelse(obs.pair1$dif.travel.time <= obs.pair2$dif.travel.time, "near", "far")
  obs.pair2$IV = ifelse(obs.pair1$dif.travel.time > obs.pair2$dif.travel.time, "near", "far")
  
  # Append the matched observational data
  obs.pairfull = rbind(obs.pair1, obs.pair2)
  match_shorttime = obs.pairfull[obs.pairfull$IV == "near", ]
  match_longtime = obs.pairfull[obs.pairfull$IV == "far", ]
  dif.travel.time.near = obs.pairfull[obs.pairfull$IV == "near", ]$dif.travel.time
  dif.travel.time.far = obs.pairfull[obs.pairfull$IV == "far", ]$dif.travel.time
  
  print(paste("N. matched pairs for year", i+1994, ":", nrow(obs.pair1)))
  
  # Compute abs. avg. difference in excess travel time
  df.tt = sum(abs(obs.pair1$dif.travel.time - obs.pair2$dif.travel.time)) / nrow(obs.pair1)
  return(list(obs.full = obs.full,
              template.full = template.full,
              match_shorttime = match_shorttime,
              match_longtime = match_longtime,
              dif.travel.time.near = dif.travel.time.near,
              dif.travel.time.far = dif.travel.time.far,
              discarded.full = discarded.full,
              match.distances = matching_keep$Distance,
              discarded.distances = matching_discard$Distance))
}

# Define helper function
mt2 = function(df, sink_prop = 0.15, i) {
  df$highrisk = case_when(df$gestage <= 30 & df$ageM >= 35 ~ 1, # highest risk
                          df$gestage <= 30 & df$ageM < 35 ~ 2,
                          df$gestage > 30 & df$ageM >= 35 ~ 3,
                          df$gestage > 30 & df$ageM < 35 ~ 4) # lowest risk
  
  # Collapse "other.race" and "race.mis" variables into one
  df$race.OTHER = ifelse(df$other.race == 1 | df$race.mis == 1, 1, 0)
  
  num.units = dim(df)[1]
  num.sinks = comp_sink_size(num.units, sink_prop)
  template.idx = sample(nrow(df), size = num.sinks, replace = F)
  print(paste("N.sinks for year", i+1994, ":", num.sinks))
  print(paste("N.expected matched pairs for year", i+1994, ":", (num.units-num.sinks)/2))
  
  # Define template data
  template.full = df[sort(template.idx), ]
  rownames(template.full) = NULL # Re-index rows
  
  # Select covariates to match on
  template = template.full %>% dplyr::select(c(bthwght, parity, fee.for.service, below_poverty,
                                               ageM, medu, singlebirth, highrisk, race.OTHER))
  template.rest = template.full %>% dplyr::select(-c(bthwght, parity, fee.for.service, below_poverty,
                                                     ageM, medu, singlebirth, highrisk, race.OTHER))
  
  # Define observational data
  obs.full = df
  
  # Select covariates to match on
  obs = obs.full %>% dplyr::select(c(bthwght, parity, fee.for.service, below_poverty,
                                     ageM, medu, singlebirth, highrisk, race.OTHER))
  obs.rest = obs.full %>% dplyr::select(-c(bthwght, parity, fee.for.service, below_poverty,
                                           ageM, medu, singlebirth, highrisk, race.OTHER))
  num.obs = dim(obs)[1]
  IV = obs.full$dif.travel.time
  
  # Construct X
  X = rbind(template, obs)
  X$s = NA
  
  # Create a dummy variable indicating membership in template
  X[1:num.sinks, "s"] = 1
  X[(num.sinks+1): nrow(X), "s"] = 0
  
  # Estimate probabilities of participation
  prob.part = est_prob_part(X)
  
  # Create an array of est. prob. of part. for template data
  s.hat.temp = prob.part[X$s==1]
  
  # Construct a distance matrix btw obs. units & templates
  distmat = construct_dist_mat(obs, template) # Rank-based Maha distance matrix
  
  # Implement generalizability caliper v2
  distmat = gen_caliper2(prob.part, X, s.hat.temp, gamma.mult = 1, D=500,
                         num.obs, distmat)
  
  # Implement dose caliper
  distmat = dose_caliper(IV, C=250, tau=sd(IV), distmat)
  
  # Prevent sinks from being matched to other sinks
  distmat = no_sink_match(distmat, num.obs, max(distmat)*3)
  
  # Do the matching
  distmat = distancematrix(distmat)
  matching = nonbimatch(distmat)$halves
  
  ### Subset matches we want to keep
  pair_kept = which(matching$Group1.Row <= num.obs & matching$Group2.Row <= num.obs)
  matching_keep = matching[matching$Group1.Row <= num.obs & matching$Group2.Row <= num.obs, ]
  
  ### Subset obs.--sinks matches
  pair_discard = which(matching$Group1.Row <= num.obs & matching$Group2.Row > num.obs)
  matching_discard = matching[matching$Group1.Row <= num.obs & matching$Group2.Row > num.obs, ]
  discarded.full = obs.full[matching_discard$Group1.Row, ]
  
  # Subset moms in matched pairs
  obs.pair1 = obs.full[matching_keep$Group1.Row, ]
  obs.pair2 = obs.full[matching_keep$Group2.Row, ]
  
  # Extract IVs within matched pairs
  obs.pair1$IV = ifelse(obs.pair1$dif.travel.time <= obs.pair2$dif.travel.time, "near", "far")
  obs.pair2$IV = ifelse(obs.pair1$dif.travel.time > obs.pair2$dif.travel.time, "near", "far")
  
  # Append the matched observational data
  obs.pairfull = rbind(obs.pair1, obs.pair2)
  match_shorttime = obs.pairfull[obs.pairfull$IV == "near", ]
  match_longtime = obs.pairfull[obs.pairfull$IV == "far", ]
  dif.travel.time.near = obs.pairfull[obs.pairfull$IV == "near", ]$dif.travel.time
  dif.travel.time.far = obs.pairfull[obs.pairfull$IV == "far", ]$dif.travel.time
  
  print(paste("N. matched pairs for year", i+1994, ":", nrow(obs.pair1)))
  
  # Compute abs. avg. difference in excess travel time
  df.tt = sum(abs(obs.pair1$dif.travel.time - obs.pair2$dif.travel.time)) / nrow(obs.pair1)
  return(list(obs.full = obs.full,
              template.full = template.full,
              match_shorttime = match_shorttime,
              match_longtime = match_longtime,
              dif.travel.time.near = dif.travel.time.near,
              dif.travel.time.far = dif.travel.time.far,
              discarded.full = discarded.full,
              match.distances = matching_keep$Distance,
              discarded.distances = matching_discard$Distance))
}

################################################################################
obs.pool = data.frame()
template.pool = data.frame()
match_shorttime.pool = data.frame()
match_longtime.pool = data.frame()
discard.pool = data.frame()
match.distances.pool = c() 
discarded.distances.pool = c()

# Perform matching
for (i in 1:10) {
  df = process_obs(data.names, i)
  
  #########################################
  ###### MATCHING FOR BLACK SUBGROUP ###### 
  #########################################
  # Perform matching for the black subgroup
  df_black = df %>%
    filter(black == 1)
  
  res_black = mt(df_black, sink_prop = 0.2, i)
  
  # Save matched data for the black subgroup
  write.csv(res_black$match_longtime, paste("m2_black_long", 1994 + i, ".csv", sep = "_"))
  write.csv(res_black$match_shorttime, paste("m2_black_short", 1994 + i, ".csv", sep = "_"))
  write.csv(res_black$template.full, paste("m2_black_template", 1994 + i, ".csv", sep = "_"))
  
  # Append matched data to the dataframe
  obs.pool = rbind(obs.pool, res_black$obs.full)
  template.pool = rbind(template.pool, res_black$template.full)
  match_shorttime.pool = rbind(match_shorttime.pool, res_black$match_shorttime)
  match_longtime.pool = rbind(match_longtime.pool, res_black$match_longtime)
  discard.pool = rbind(discard.pool, res_black$discarded.full)
  match.distances.pool = append(match.distances.pool, res_black$match.distances)
  discarded.distances.pool = append(discarded.distances.pool, res_black$discarded.distances)
  
  #############################################
  ###### MATCHING FOR NON-BLACK SUBGROUP ###### 
  #############################################
  # Subset non-black subgroup
  df_nonblack = df %>%
    filter(black == 0)
  
  df_nonblack$highrisk = case_when(df_nonblack$gestage <= 30 & df_nonblack$ageM >= 35 ~ 1, # highest risk
                                   df_nonblack$gestage <= 30 & df_nonblack$ageM < 35 ~ 2,
                                   df_nonblack$gestage > 30 & df_nonblack$ageM >= 35 ~ 3,
                                   df_nonblack$gestage > 30 & df_nonblack$ageM < 35 ~ 4)
  
  # Stratify non-black subgroup 1
  df_nonblack1 = df_nonblack %>%
    filter(highrisk == 1)
  
  # Perform matching for non-black subgroup 1
  res_nonblack1 = mt2(df_nonblack1, sink_prop = 0.2, i)
  
  # Save matched data for non-black subgroup 1
  write.csv(res_nonblack1$match_longtime, paste("m2_nonblack1_long", 1994 + i, ".csv", sep = "_"))
  write.csv(res_nonblack1$match_shorttime, paste("m2_nonblack1_short", 1994 + i, ".csv", sep = "_"))
  write.csv(res_nonblack1$template.full, paste("m2_nonblack1_template", 1994 + i, ".csv", sep = "_"))
  
  # Append matched data to the dataframe
  obs.pool = rbind(obs.pool, res_nonblack1$obs.full)
  template.pool = rbind(template.pool, res_nonblack1$template.full)
  match_shorttime.pool = rbind(match_shorttime.pool, res_nonblack1$match_shorttime)
  match_longtime.pool = rbind(match_longtime.pool, res_nonblack1$match_longtime)
  discard.pool = rbind(discard.pool, res_nonblack1$discarded.full)
  match.distances.pool = append(match.distances.pool, res_nonblack1$match.distances)
  discarded.distances.pool = append(discarded.distances.pool, res_nonblack1$discarded.distances)
  
  # Stratify non-black subgroup 2
  df_nonblack2 = df_nonblack %>%
    filter(highrisk == 2)
  
  # Perform matching for non-black subgroup 2
  res_nonblack2 = mt2(df_nonblack2, sink_prop = 0.2, i)
  
  # Save matched data for non-black subgroup 2
  write.csv(res_nonblack2$match_longtime, paste("m2_nonblack2_long", 1994 + i, ".csv", sep = "_"))
  write.csv(res_nonblack2$match_shorttime, paste("m2_nonblack2_short", 1994 + i, ".csv", sep = "_"))
  write.csv(res_nonblack2$template, paste("m2_nonblack2_template", 1994 + i, ".csv", sep = "_"))
  
  # Append matched data to the dataframe
  obs.pool = rbind(obs.pool, res_nonblack2$obs.full)
  template.pool = rbind(template.pool, res_nonblack2$template.full)
  match_shorttime.pool = rbind(match_shorttime.pool, res_nonblack2$match_shorttime)
  match_longtime.pool = rbind(match_longtime.pool, res_nonblack2$match_longtime)
  discard.pool = rbind(discard.pool, res_nonblack2$discarded.full)
  match.distances.pool = append(match.distances.pool, res_nonblack2$match.distances)
  discarded.distances.pool = append(discarded.distances.pool, res_nonblack2$discarded.distances)
  
  # Stratify non-black subgroups 3 & 4
  df_nonblack3 = df_nonblack %>%
    filter(highrisk == 3 | highrisk == 4)
  
  # Perform matching for non-black subgroups 3 & 4
  res_nonblack3 = mt2(df_nonblack3, sink_prop = 0.2, i)
  
  # Save matched data for non-black subgroup 2
  write.csv(res_nonblack3$match_longtime, paste("m2_nonblack34_long", 1994 + i, ".csv", sep = "_"))
  write.csv(res_nonblack3$match_shorttime, paste("m2_nonblack34_short", 1994 + i, ".csv", sep = "_"))
  write.csv(res_nonblack3$template, paste("m2_nonblack34_template", 1994 + i, ".csv", sep = "_"))
  
  # Append matched data to the dataframe
  obs.pool = rbind(obs.pool, res_nonblack3$obs.full)
  template.pool = rbind(template.pool, res_nonblack3$template.full)
  match_shorttime.pool = rbind(match_shorttime.pool, res_nonblack3$match_shorttime)
  match_longtime.pool = rbind(match_longtime.pool, res_nonblack3$match_longtime)
  discard.pool = rbind(discard.pool, res_nonblack3$discarded.full)
  match.distances.pool = append(match.distances.pool, res_nonblack3$match.distances)
  discarded.distances.pool = append(discarded.distances.pool, res_nonblack3$discarded.distances)
}

# Plot rank-based Mahalanobis distances of two types of matched pairs
boxplot(match.distances.pool, discarded.distances.pool,
        ylab = "Rank-Based Mahalanobis Distance", 
        names = c("Obs-Obs", "Obs-Template"))

# Check covariate balance
mean_template = apply(template.pool, 2, mean, na.rm = T)

match_longtime.pool$death_total = match_longtime.pool$death + match_longtime.pool$death_fetal
mean_far = apply(match_longtime.pool[, which(names(match_longtime.pool) != "IV")], 2, mean, na.rm = T)

match_shorttime.pool$death_total = match_shorttime.pool$death + match_shorttime.pool$death_fetal
mean_near = apply(match_shorttime.pool[, which(names(match_shorttime.pool) != "IV")], 2, mean, na.rm = T)

sd_far = apply(match_longtime.pool[, which(names(match_longtime.pool) != "IV")], 2, stats::sd, na.rm = T)
sd_near = apply(match_shorttime.pool[, which(names(match_shorttime.pool) != "IV")], 2, stats::sd, na.rm = T)
sd_pooled = sqrt((sd_far^2 + sd_near^2) / 2)
asd = abs((mean_near - mean_far) / sd_pooled)

balance_table = round(data.frame(mean_near, mean_far, asd), 4)
xtable(balance_table, digits = 4)

mean_discard = apply(discard.pool, 2, mean, na.rm = T)
balance_table_template = round(data.frame(mean_template, mean_discard), 4)
xtable(balance_table_template, digits = 4)

