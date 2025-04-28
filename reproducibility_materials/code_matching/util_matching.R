# HELPER FUNCTIONS FOR MATCHING
################################################################################
dat <- function(preemiedata){ 
  preemiedata=read.csv(preemiedata)
  
  bthwght=preemiedata$bthwght;
  bthwghtreg=lm(bthwght~preemiedata$gestage_weeks);
  bthwght[is.na(bthwght)]=coef(bthwghtreg)[1]+coef(bthwghtreg)[2]*preemiedata$gestage_weeks[is.na(bthwght)];
  gestage=preemiedata$gestage_weeks;
  gestagereg=lm(gestage~preemiedata$bthwght);
  gestage[is.na(gestage)]=coef(gestagereg)[1]+coef(gestagereg)[2]*preemiedata$bthwght[is.na(gestage)];
  gestdiabetes=preemiedata$Gestational_DiabetesM  
  precare = preemiedata$precare
  precare[is.na(precare)]=mean(precare,na.rm=TRUE)
  singlebirth=as.numeric(preemiedata$multiples==1); 
  parity=preemiedata$parity;
  parity[is.na(parity)]=mean(parity,na.rm=TRUE);
  
  #Mother's covariates
  ageM=preemiedata$ageM; 
  medu=preemiedata$medu;
  medu[is.na(medu)]=mean(medu, na.rm = TRUE)
  race=preemiedata$raceM;
  
  # combine fee.for.service and other.insur: Federal/state/other
  insurance=preemiedata$insuranceM
  fee.for.service <- as.numeric(insurance==1) 
  fee.for.service[is.na(fee.for.service)]=0
  HMO <- as.numeric(insurance==2)
  HMO[is.na(HMO)]=0
  Federal.state.other.insur <- as.numeric(insurance==3|insurance==4|insurance==5|is.na(insurance))
  Federal.state.other.insur[is.na(Federal.state.other.insur)]=0
  
  white <- as.numeric(race==1)
  white[is.na(white)]=0 
  black <- as.numeric(race==2)
  black[is.na(black)]=0
  other.race <- as.numeric(race!=1 & race!=2)
  other.race[is.na(other.race)]=0
  race.mis <- as.numeric(is.na(race))
  
  #Mother's neighborhood
  below_poverty=preemiedata$below_poverty
  below_poverty[is.na(below_poverty)]=mean(below_poverty, na.rm = TRUE)
  
  # IV
  dif.travel.time = preemiedata$diff_travel_2500
  
  #create high level NICU indicator
  high_level_NICU <- preemiedata$vol_level_big_2500
  losI = preemiedata$losI
  death = preemiedata$death
  death_fetal = preemiedata$death_fetal
  
  # all covariates + IV
  X.reduced=cbind(bthwght,gestage,gestdiabetes,precare,singlebirth,parity,
                  ageM,medu,white,black,other.race,race.mis,below_poverty,
                  fee.for.service,HMO,Federal.state.other.insur,
                  dif.travel.time,
                  high_level_NICU, losI, death, death_fetal);
  return(X.reduced)
}

# Process obs. data to include only relevant variables
process_obs = function(data.names, num) {
  df = dat(data.names[num])
  df = as.data.frame(df)
  return(df)
}
 
# Compute robust Mahalanobis distance 
smahal=function(X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  k<-dim(X)[2]
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,n,n)
  library(MASS)
  icov<-ginv(cv)
  for (i in 1:n) out[i,]<-mahalanobis(X,X[i,],icov,inverted=T)
  out
}

# Compute number of sinks
comp_sink_size = function(num.obs, p) {
  num.sinks = num.obs * p
  num.sinks = 2*ceiling((num.obs + num.sinks) / 2) - num.obs
  return(num.sinks)
}

# Construct rank-based Mahalanobis distance matrix of obs. units & template
construct_dist_mat = function(obs, template) {
  X = rbind(obs, template)
  distmat = smahal(X) # Rank-based Maha distance in the order of [0, 27]
  distmat = round(distmat * 10) # Now the distmat in the order of [0, 270]
  return(distmat)
}

# Implement dose caliper
dose_caliper = function(IV, C, tau, fullmat) {
  # Construct a matrix of abs. diff. in cts. IV doses for obs. units
  IV.distmat.obs = abs(outer(IV, IV, "-"))
  # Add penalty if doses aren't far apart enough by pre-specified 'tau'
  penalty.dose_caliper = C * (IV.distmat.obs <= tau)
  n.obs = dim(penalty.dose_caliper)[1]
  # Subset distance matrix for obs. units
  distmat.obs_obs = fullmat[1:n.obs, 1:n.obs]
  distmat.obs_obs = distmat.obs_obs + penalty.dose_caliper
  # Update distances for obs. units
  fullmat[1:n.obs, 1:n.obs] = distmat.obs_obs
  return (fullmat)
}

# Construct X
construct_X = function(d05, idx, obs) {
  # Check which indices in template haven't been sampled
  num.temp.uniq = nrow(d05)
  missing.idx = sort(setdiff(1:num.temp.uniq, unique(idx)))
  # Create a template dataframe whose rows are unique & have been sampled
  if (!is.na(missing.idx[1])) { # Check if there're indices that haven't been sampled
    d05 = d05[-c(missing.idx), ] # Remove those indices from template
  }
  X = rbind(d05, obs)
  X$s = NA
  # Create a dummy variable indicating membership in template
  num.obs = nrow(obs)
  num.sinks.orig = nrow(X) - num.obs
  X[1:num.sinks.orig, "s"] = 1
  X[(num.sinks.orig+1): nrow(X), "s"] = 0
  
  return(list(X=X, missing.idx=missing.idx))
}

# Estimate probabilities of participation i.e., generalizability scores
est_prob_part = function(X) {
  prob.part = glm(s ~ bthwght + fee.for.service + below_poverty + medu,
                  data = X, family = binomial)$fitted.values
  return(prob.part)
}

# Create an array of est. prob. of part. for template data
prob_part_temp = function(prob.part, X, missing.idx, num.temp.uniq, idx) {
  s.temp = prob.part[X$s==1]
  # If there's a missing index, append zero prob. to missing index
  # Total length of est. prob. of part. for template should = num.temp.uniq
  while(!is.na(missing.idx[1])) {
    missing.first = missing.idx[1]
    s.temp = c(s.temp[1:missing.first - 1], 0, s.temp[missing.first:length(s.temp)])
    missing.idx = missing.idx[-1]
  }
  # Construct a vector of est. prob. of part. for template data
  s.hat.temp = c()
  for (i in 1:num.temp.uniq) {
    n = length(which(idx == i))
    s.hat.temp = c(s.hat.temp, rep(s.temp[i], n))
  }
  return(s.hat.temp)
}

# Generalizability caliper 2
gen_caliper2 = function(prob.part, X, s.hat.temp, gamma.mult, D, num.obs, distmat) {
  # Construct a matrix of |S_i - S_j|: each row is an obs. unit, each column 
  # a template unit
  prob.part.distmat = abs(outer(prob.part[X$s==0], s.hat.temp, "-"))
  # Determine caliper width
  gamma = gamma.mult * sd(prob.part.distmat)
  # Add a penalty if units have similar probabilities of participation
  total = dim(prob.part.distmat)[1] + dim(prob.part.distmat)[2]
  penalty.gen_caliper = D * (prob.part.distmat <= gamma)
  distmat.obs_temp = distmat[1:num.obs, (num.obs+1):total]
  C.tilde = max(distmat.obs_temp) + 1
  distmat.obs_temp = C.tilde - distmat.obs_temp + penalty.gen_caliper
  
  # Update full matrix with new distances
  distmat[1:num.obs, (num.obs+1):total] = distmat.obs_temp
  distmat[(num.obs+1):total, 1:num.obs] = t(distmat.obs_temp)
  
  return(distmat)
}

# Prevent sinks from being matched to other sinks, including themselves
no_sink_match = function(fullmat, num.obs, penalty) {
  total = dim(fullmat)[1]
  distmat.temp_temp = fullmat[(num.obs+1):total, (num.obs+1):total]
  distmat.temp_temp = penalty
  
  # Update full matrix
  fullmat[(num.obs+1):total, (num.obs+1):total] = distmat.temp_temp
  return(fullmat)
}



