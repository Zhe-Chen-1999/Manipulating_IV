# IMPLEMENT MATCH 1
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
# Read and process data for each year
df = process_obs(data.names, 1)

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

# Randomly split the data into two samples
num.obs1 = floor(dim(obs)[1] / 2)
obs.1.idx = sample(nrow(obs), num.obs1)
obs.1 = obs[obs.1.idx, ]
obs.rest.1 = obs.rest[obs.1.idx, ]

obs.2 = obs[-obs.1.idx, ]
obs.rest.2 = obs.rest[-obs.1.idx, ]
num.obs2 = nrow(obs.2)

IV.1 = IV[obs.1.idx]
IV.2 = IV[-obs.1.idx]

# Define number of sinks to be half of the number of observations in each of the two samples
num.sinks1 = floor(num.obs1 / 2)
num.sinks2 = floor(num.obs2 / 2)

# Construct a matrix for sinks
X.sink1 = matrix(0, nrow = num.sinks1, ncol = dim(obs.1)[2])
colnames(X.sink1) = colnames(obs.1)

X.sink2 = matrix(0, nrow = num.sinks2, ncol = dim(obs.2)[2])
colnames(X.sink2) = colnames(obs.2)

# Construct X
X1 = rbind(obs.1, X.sink1)

X2 = rbind(obs.2, X.sink2)

# Construct a distance matrix between observational units
distmat.1 = smahal(X1) # Rank-based Maha distance in the order of [0, 27]
distmat.1 = round(distmat.1 * 10) # Now the distmat in the order of [0, 270]

distmat.2 = smahal(X2) # Rank-based Maha distance in the order of [0, 27]
distmat.2 = round(distmat.2 * 10) # Now the distmat in the order of [0, 270]

# Implement dose caliper
distmat.1 = dose_caliper(IV.1, C=250, tau=sd(IV.1), distmat.1)

distmat.2 = dose_caliper(IV.2, C=250, tau=sd(IV.2), distmat.2)

# Make each sink - baby pair at zero discrepancy
distmat.1[(num.obs1+1):(num.obs1+num.sinks1), 1:num.obs1] = 0
distmat.1[1:num.obs1, (num.obs1+1):(num.obs1+num.sinks1)] = 0

distmat.2[(num.obs2+1):(num.obs2+num.sinks2), 1:num.obs2] = 0
distmat.2[1:num.obs2, (num.obs2+1):(num.obs2+num.sinks2)] = 0

# Make each sink - sink pair at infinite discrepancy
lambda1 = 10 * max(distmat.1)
distmat.1[(num.obs1+1):(num.obs1+num.sinks1), (num.obs1+1):(num.obs1+num.sinks1)] = lambda1

lambda2 = 10 * max(distmat.2)
distmat.2[(num.obs2+1):(num.obs2+num.sinks2), (num.obs2+1):(num.obs2+num.sinks2)] = lambda2

# Do the matching
distmat.1 = distancematrix(distmat.1)
matching.1 = nonbimatch(distmat.1)

distmat.2 = distancematrix(distmat.2)
matching.2 = nonbimatch(distmat.2)

# Among the matches, retain babies NOT matched to sinks
matches.halves1 = matching.1$halves
matches.halves1 = matches.halves1[matches.halves1$Group2.Row <= num.obs1, ]

matches.halves2 = matching.2$halves
matches.halves2 = matches.halves2[matches.halves2$Group2.Row <= num.obs2, ]

# Subset moms in matched pairs
obs.pair1.s1 = obs.1[matches.halves1$Group1.Row, ] 
obs.pair2.s1 = obs.1[matches.halves1$Group2.Row, ]
obs.pair1.s1.rest = obs.rest.1[matches.halves1$Group1.Row, ]
obs.pair2.s1.rest = obs.rest.1[matches.halves1$Group2.Row, ]

obs.pair1.s1.full = cbind(obs.pair1.s1, obs.pair1.s1.rest)
obs.pair2.s1.full = cbind(obs.pair2.s1, obs.pair2.s1.rest)

obs.pair1.s2 = obs.2[matches.halves2$Group1.Row, ] 
obs.pair2.s2 = obs.2[matches.halves2$Group2.Row, ]
obs.pair1.s2.rest = obs.rest.2[matches.halves2$Group1.Row, ]
obs.pair2.s2.rest = obs.rest.2[matches.halves2$Group2.Row, ]

obs.pair1.s2.full = cbind(obs.pair1.s2, obs.pair1.s2.rest)
obs.pair2.s2.full = cbind(obs.pair2.s2, obs.pair2.s2.rest)

# Extract IVs within matched pairs
obs.pair1.s1.full$IV = ifelse(obs.pair1.s1.full$dif.travel.time <= obs.pair2.s1.full$dif.travel.time, "near", "far")
obs.pair2.s1.full$IV = ifelse(obs.pair1.s1.full$dif.travel.time > obs.pair2.s1.full$dif.travel.time, "near", "far")

obs.pair1.s2.full$IV = ifelse(obs.pair1.s2.full$dif.travel.time <= obs.pair2.s2.full$dif.travel.time, "near", "far")
obs.pair2.s2.full$IV = ifelse(obs.pair1.s2.full$dif.travel.time > obs.pair2.s2.full$dif.travel.time, "near", "far")

# Append the matched observational data
obs.pairfull.s1 = rbind(obs.pair1.s1.full, obs.pair2.s1.full)
obs.pairfull.s2 = rbind(obs.pair1.s2.full, obs.pair2.s2.full)

match_shorttime = rbind(obs.pairfull.s1[obs.pairfull.s1$IV == "near", ],
                        obs.pairfull.s2[obs.pairfull.s2$IV == "near", ])

match_longtime = rbind(obs.pairfull.s1[obs.pairfull.s1$IV == "far", ],
                       obs.pairfull.s2[obs.pairfull.s2$IV == "far", ])

write.csv(match_shorttime, "m1_short_1995.csv", row.names=F)
write.csv(match_longtime, "m1_long_1995.csv", row.names=F)
