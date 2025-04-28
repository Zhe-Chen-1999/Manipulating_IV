# IMPLEMENT CLASSIFICATION PERMUTATION TEST (CPT)
################################################################################
# Define helper functions
## Compute test statistic: in-sample classification accuracy rate
compute_accuracy = function(t, t.hat) {
  sum1 = sum(t * t.hat)
  
  t2 = ifelse(t==0, 1, 0)
  t.hat2 = ifelse(t.hat == 0, 1, 0)
  sum2 = sum(t2 * t.hat2)
  
  output = (sum1 + sum2) / length(t.hat)
  return(output)
}

## Plot null distribution
plot_null_dist = function(rmse.obs, rmse.vec, plottitle) {
  xmax = max(c(rmse.obs,rmse.vec)) + min(rmse.vec)/20
  xmin = min(c(rmse.obs,rmse.vec)) - max(rmse.vec)/20
  p = ggplot(data.frame(rmse = rmse.vec),aes(x=rmse)) + 
    geom_density(alpha=0.75, fill="grey") + xlim(xmin,xmax) + 
    labs(title = plottitle, y = "", x = "\n Test Statistic Under the Null") +
    geom_vline(xintercept = rmse.obs, col = "black", lwd=1) + theme_bw()+ 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(panel.border = element_blank(), 
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.text.y = element_blank(),           
          axis.title = element_text(size = 12))
  return(p)
}

## Compute p-value
compute_pval = function(rmse.vec, rmse.obs, perm.N) {
  length(which(rmse.vec >= rmse.obs)) / perm.N
}

# Compute bounding p-value
pvalue_bound_ss_cpt = function(t, I, Gamma, score_i1, score_i2){
  max_score = pmax(score_i1, score_i2)
  min_score = pmin(score_i1, score_i2)
  T_mean = sum(Gamma/(1+Gamma)*max_score + 1/(1+Gamma)*min_score)
  T_var = sum(Gamma/(1+Gamma)^2 * (max_score - min_score)^2)
  
  T_null = (t - T_mean)/sqrt(T_var)
  return(1 - pnorm(T_null))
}

# Compute RSV for SS-CPT
RV_Calculate_SS_CPT = function(dt_match, d, score_type = 0){
  dt_match = dt_match[order(dt_match$pair, dt_match$IV.binary),]
  
  n = dim(dt_match)[1]
  if(round(n/2) %% 2 == 0){
    k = round(n/2)
  }else{
    k = round(n/2) - 1
  }
  
  df1 = dt_match[(1:k), (1:(d+2))]
  df2 = dt_match[((k+1):n), (1:(d+2))]
  
  model_logistic_1 = glm(IV.binary ~ ., data = df1, family = "binomial")
  
  propen2 = predict(model_logistic_1, newdata = df2, 
                    type = "response") 
  
  model_logistic_2 = glm(IV.binary ~ ., data = df2, family = "binomial")
  
  propen1 = predict(model_logistic_2, newdata = df1, 
                    type = "response") 
  
  if (score_type == 0){
    score2 = cbind(dt_match[(k+1):n,], propen2 = propen2) %>% group_by(pair) %>%
      mutate(pred = ifelse(propen2 == max(propen2), 1, 0)) %>%
      ungroup() 
    
    
    score2 = c(score2)$pred
    
    
    score1 = cbind(dt_match[1:k,], propen1 = propen1) %>% group_by(pair) %>%
      mutate(pred = ifelse(propen1 == max(propen1), 1, 0)) %>%
      ungroup() 
    
    
    score1 = c(score1)$pred
    
    
  } else {
    score2 = propen2
    score1 = propen1
  }
  
  score2_i1 = score2[seq(1, length(score2)-1, 2)]
  score2_i2 = score2[seq(2, length(score2), 2)]
  test_statistic_2 = sum(score2 * df2$IV.binary)
  
  score1_i1 = score1[seq(1, length(score1)-1, 2)]
  score1_i2 = score1[seq(2, length(score1), 2)]
  test_statistic_1 = sum(score1 * df1$IV.binary)
  
  
  if (min(pvalue_bound_ss_cpt(t = test_statistic_2, I = length(score2)/2, Gamma = 1,
                              score_i1 = score2_i1, score_i2 = score2_i2),
          pvalue_bound_ss_cpt(t = test_statistic_1, I = length(score1)/2, Gamma = 1,
                              score_i1 = score1_i1, score_i2 = score1_i2)) >= 0.025){
    return(1)
  }
  else {
    RV=1+0.01
    while (min(pvalue_bound_ss_cpt(t = test_statistic_2, I = length(score2)/2, Gamma = RV,
                                   score_i1 = score2_i1, score_i2 = score2_i2),
               pvalue_bound_ss_cpt(t = test_statistic_1, I = length(score1)/2, Gamma = RV,
                                   score_i1 = score1_i1, score_i2 = score1_i2)) < 0.025) {
      RV=RV+0.01
    }
    return(RV)
  }
}

# Load data
m1.long = data.frame()
m1.short = data.frame()

for (i in 1:10) {
  # Read far data
  df.long = read.csv(paste("m1_black_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  df.long = read.csv(paste("m1_nonblack1_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  df.long = read.csv(paste("m1_nonblack2_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  df.long = read.csv(paste("m1_nonblack34_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  
  # Read near data
  df.short = read.csv(paste("m1_black_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
  df.short = read.csv(paste("m1_nonblack1_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
  df.short = read.csv(paste("m1_nonblack2_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
  df.short = read.csv(paste("m1_nonblack34_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
}

# Choose covariates
m1.long = m1.long %>% dplyr::select(-c(X, gestage, gestdiabetes, precare, white, 
                                       other.race, race.mis, HMO, 
                                       Federal.state.other.insur, IV))

m1.short = m1.short %>% dplyr::select(-c(X, gestage, gestdiabetes, precare, white, 
                                         other.race, race.mis, HMO, 
                                         Federal.state.other.insur, IV))

# Reorder columns s.t. the covariates are together
m1.long = m1.long[, c(1:8, 14:15, 9:13)]
m1.short = m1.short[, c(1:8, 14:15, 9:13)]

# Create a pair ID variable
m1.long$pair = 1:nrow(m1.long)
m1.short$pair = 1:nrow(m1.short)

# Create merged data for Z
m1 = rbind(m1.short, m1.long)
m1 = m1[order(m1$pair), ]

# Dichotomize IV: within each matched pair, the larger IV = 1 & smaller IV = 0.
m1 = m1 %>% 
  group_by(pair) %>%
  mutate(IV.binary = ifelse(dif.travel.time == max(dif.travel.time), 1, 0))

# Check if there are any pairs with the same IV values
m1 %>% 
  group_by(pair) %>%
  mutate(n.dichotomizedIV = n_distinct(IV.binary) != 2) %>%
  ungroup %>%
  count(n.dichotomizedIV) 

# Separate the matched pairs into 2 groups by whether or not there are ties in the IVs
m1.ties = m1 %>% 
  group_by(pair) %>%
  filter(n_distinct(IV.binary) != 2) %>%
  ungroup()

m1.no.ties = m1 %>% 
  group_by(pair) %>%
  filter(n_distinct(IV.binary) == 2) %>%
  ungroup()

# Break ties for any matched pairs with the same IV values  
m1.ties = m1.ties %>% 
  group_by(pair) %>%
  mutate(IV.binary = sample(0:1, n(), replace = F)) %>% # break ties
  ungroup() 

# Append the two sets of matched sub-samples
m1 = rbind(m1.no.ties, m1.ties)

# Fit a simple logistic regression model 
lr = glm(IV.binary ~ bthwght + parity + fee.for.service + below_poverty + 
           ageM + medu + singlebirth + black + highrisk + race.OTHER, data = m1, family = "binomial")

# Create a column of predicted treatment assignments: within each matched pair, 
# the mother with a higher predicted probability gets a 1.
m1$lr.pred = predict(lr, newdata = dplyr::select(m1, -dif.travel.time, -pair, -IV.binary), 
                     type = "response")

m1 = m1 %>% 
  group_by(pair) %>%
  mutate(lr.pred = ifelse(lr.pred == max(lr.pred), 1, 0))

# Check if there are any pairs with the same predicted IV assignments 
m1 %>% 
  group_by(pair) %>%
  mutate(n.sameprob = n_distinct(lr.pred) != 2) %>%
  ungroup %>%
  count(n.sameprob) 

# Create a sequence of pairs in which to perform random permutations
n.pairs = 1:nrow(m1)
pair.list = list()
for (i in 1:(nrow(m1) / 2)) {
  pair.list[[i]] = n.pairs[(1+((i-1)*2)) : (2+((i-1)*2))]
}
block = how(blocks = rep(seq_along(pair.list), sapply(pair.list, length)))

# Compute test statistic: in-sample classification accuracy rate
test.stat2.m1 = compute_accuracy(m1$IV.binary, m1$lr.pred)

# Define number of permutations
perm.N = 500
# Define vector of RMSEs for null distribution
accuracy.vec.m1 = matrix(NA, ncol = 1, nrow = perm.N)

# Create a null distribution of in-sample classification accuracy rates
for (s in c(1:perm.N)) {
  if(s %% 50 == 0){cat("Iteration: ", s, "\n")}
  # Randomly permute binary IV within matched pairs
  permuted.idx = shuffle(length(block$blocks), block) 
  # Extract randomly permuted IV doses
  m1$permuted.IV = m1[permuted.idx, "dif.travel.time"]
  # Dichotomize IV: within each matched pair, the larger IV = 1 & smaller IV = 0.
  m1 = m1 %>% 
    group_by(pair) %>%
    mutate(permuted.IV = ifelse(permuted.IV == max(permuted.IV), 1, 0))
  # Fit a simple logistic regression model 
  lr = glm(permuted.IV ~ bthwght + parity + fee.for.service + below_poverty +
             ageM + medu + singlebirth + black + highrisk + race.OTHER, data = m1, family = "binomial")
  # Create a column of predicted treatment assignments: within each matched pair, 
  # the mother with a higher predicted probability gets a 1.
  m1$lr.pred.perm = predict(lr, newdata = m1[, 1:10], type = "response")
  m1 = m1 %>% 
    group_by(pair) %>%
    mutate(lr.pred.perm = ifelse(lr.pred.perm == max(lr.pred.perm), 1, 0))
  # Compute test statistic: in-sample classification accuracy rate
  accuracy.vec.m1[s] = compute_accuracy(m1$permuted.IV, m1$lr.pred.perm)
}

# Plot distribution
nulldist.m1 = plot_null_dist(test.stat2.m1, accuracy.vec.m1, "M1 \n")

# Compute p-value
pval.m1 = 1 - compute_pval(test.stat2.m1, accuracy.vec.m1, perm.N)

# Subset relevant columns for SS-CPT computation
m1.sscpt = m1[, c(1:10, 17, 16)]

# Compute RSV for SS-CPT accuracy
RV_Calculate_SS_CPT(m1.sscpt, d=9, score_type = 0) 
