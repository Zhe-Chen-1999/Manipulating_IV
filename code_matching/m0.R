# IMPLEMENT MATCH 0
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
############################# M0: RETAIN ALL MOMS ##############################
################################################################################
# Store data
match_shorttime.pool = data.frame()
match_longtime.pool = data.frame()
dif.travel.time.near.pool = c()
dif.travel.time.far.pool = c()

for (i in 1:10) {
  # Read and process data for each year
  df = process_obs(data.names, i)
  
  # Define a categorical variable based on gestage & ageM
  df$highrisk = case_when(df$gestage <= 30 & df$ageM >= 35 ~ 1, # highest risk
                          df$gestage <= 30 & df$ageM < 35 ~ 2,
                          df$gestage > 30 & df$ageM >= 35 ~ 3,
                          df$gestage > 30 & df$ageM < 35 ~ 4) # lowest risk
  IV = df$dif.travel.time
  
  # Select covariates to match on
  obs = df %>% dplyr::select(c(bthwght, parity, fee.for.service, below_poverty, 
                               ageM, medu, singlebirth, black, highrisk))
  obs.rest = df %>% dplyr::select(-c(bthwght, parity, fee.for.service, below_poverty, 
                                     ageM, medu, singlebirth, black, highrisk))
  num.obs = dim(obs)[1]
  
  # Construct X
  X = obs
  
  # Construct a distance matrix between observational units
  distmat = smahal(X) # Rank-based Maha distance in the order of [0, 27]
  distmat = round(distmat * 10) # Now the distmat in the order of [0, 270]
  
  # Do the matching
  distmat = distancematrix(distmat)
  matching = nonbimatch(distmat)
  matches.halves = matching$halves
  
  # Subset moms in matched pairs
  obs.pair1 = df[matches.halves$Group1.Row, ] 
  obs.pair2 = df[matches.halves$Group2.Row, ]
  
  # Extract IVs within matched pairs
  obs.pair1$IV = ifelse(obs.pair1$dif.travel.time <= obs.pair2$dif.travel.time, "near", "far")
  obs.pair2$IV = ifelse(obs.pair1$dif.travel.time > obs.pair2$dif.travel.time, "near", "far")
  
  # Append the matched observational data
  obs.pairfull = rbind(obs.pair1, obs.pair2)
  match_shorttime = obs.pairfull[obs.pairfull$IV == "near", ]
  match_longtime = obs.pairfull[obs.pairfull$IV == "far", ]
  
  dif.travel.time.near = obs.pairfull[obs.pairfull$IV == "near", ]$dif.travel.time
  dif.travel.time.far = obs.pairfull[obs.pairfull$IV == "far", ]$dif.travel.time
  
  # Compute abs. avg. difference in excess travel time within matched pairs
  diff.tt = sum(abs(obs.pair1$dif.travel.time - obs.pair2$dif.travel.time)) / nrow(obs.pair1)
  
  # Store data
  match_shorttime.pool = rbind(match_shorttime.pool, match_shorttime)
  match_longtime.pool = rbind(match_longtime.pool, match_longtime)
  dif.travel.time.near.pool = c(dif.travel.time.near.pool, dif.travel.time.near)
  dif.travel.time.far.pool = c(dif.travel.time.far.pool, dif.travel.time.far)
  
  print(paste("N. matched pairs for year ", i+1994, ":", nrow(obs.pair1)))
  print(paste("Avg. abs. diff in IV for year ", i+1994, ":", diff.tt))
}

# Save matched data
write.csv(match_shorttime.pool, "m0_short_1995-2004.csv", row.names=F)
write.csv(match_longtime.pool, "m0_long_1995-2004.csv", row.names=F)

