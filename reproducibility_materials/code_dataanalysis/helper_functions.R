#####################################################################################
############## Helper functions for biased randomization inference in Theorem 1 #####
#####################################################################################
## Parameters:
# data: a dataframe containing Z, D, R with N=2I rows. The first I rows are matched treated units and the second I rows are the corresponding matched control units
# K0: the assumed lower bound for treatment effects across all units
# alpha: one-sided confidence level
# R: outcome vector of length 2I
# Gamma: Capital Gamma_i vector of length I

## The lower limit of the (1-alpha)-level one-sided lower CI for (LB - K0)
transformed.tau.bar.lower.limit<-function(data = dat_all, K0 = 0, alpha = 0.05, R = dat_all$death_total, Gamma = NULL){
  tau_hat =  R[which(data$Z==1)]- K0*data$D[which(data$Z==1)] - (R[which(data$Z==0)]- K0*data$D[which(data$Z==0)])
  tau_hat_Gamma = (Gamma + 1)*((Gamma + 1)*tau_hat - (Gamma - 1)*abs(tau_hat))/(4*Gamma)
  tau_bar_Gamma = mean(tau_hat_Gamma)
  S_square_Gamma = sum((tau_hat_Gamma - tau_bar_Gamma)^2)/(I*(I-1)) 
  transformed.tau.bar.lower.limit = tau_bar_Gamma - qnorm(1-alpha)*sqrt(S_square_Gamma)
  return(transformed.tau.bar.lower.limit)
} 

## The lower limit of the (1-alpha/2)-level one-sided lower CI for LB
LB.lower.limit<-function(data = dat_all, K0 = 0, alpha = 0.05/2, R = dat_all$death_total, Gamma = NULL){
  LB.lower.limit = transformed.tau.bar.lower.limit(data = data, K0 = K0, alpha = alpha, R = R, Gamma = Gamma) + K0 
  return(LB.lower.limit)
} 

## The upper limit of the (1-alpha/2)-level one-sided upper CI for UB
UB.upper.limit<-function(data = dat_all, K1 = 1, alpha = 0.05/2, R = dat_all$death_total, Gamma = NULL){
  data$Z =  1-data$Z
  UB.upper.limit = -transformed.tau.bar.lower.limit(data = data, K0 = K1, alpha = alpha, R = R, Gamma = Gamma) + K1
  return(UB.upper.limit)
} 

## Combine inference results for SATE and effect ratio in one function
biased.RI <- function(CPT.max.Gamma = 1.17, gamma = NULL,
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.01,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, R = dat_all$death_total, alpha = 0.05){
  
  if(is.null(gamma)){
    gamma = log(CPT.max.Gamma)/max(data$dif.travel.time[which(data$Z==1)] - data$dif.travel.time[which(data$Z==0)])
  }
  Gamma = exp(gamma*(data$dif.travel.time[which(data$Z==1)] - data$dif.travel.time[which(data$Z==0)]))
  
  I = dim(data)[1]/2
  
  ### Biased ramdomization inference for SATE ### 
  UB.est = c()
  UB.upper.limit = c()
  LB.est = c()
  LB.lower.limit = c()
  
  if(!is.null(K0)){
    ## 97.5% lower confidence limit of LB with K0 = 0
    LB.lower.limit = LB.lower.limit(data = data, K0 = K0, alpha = alpha/2, R = R, Gamma = Gamma)
  }
  
  for (K in seq(K.lb, K.ub, by = step) ) {
    
    if(is.null(K0)){
      
      ## 97.5% lower confidence limit of LB with -K.ub <= K0 = -K <= -K.lb
      LB.lower.limit = c(LB.lower.limit, LB.lower.limit(data = data, K0 = -K, alpha = alpha/2, R = R, Gamma = Gamma))
    }
    
    ## 97.5% upper confidence limit of UB with K.lb <= K1 = K <= K.ub
    UB.upper.limit = c(UB.upper.limit, UB.upper.limit(data = data, K1 = K, alpha = alpha/2, R = R, Gamma = Gamma))
    
  }
  
  ## 95% CI for SATE (kappa)
  ci.kappa = data.frame(K = seq(K.lb, K.ub, by = step), 
                        LB.lower.limit = LB.lower.limit,
                        UB.upper.limit = UB.upper.limit)
  
  ### Inference for effect ratio ### 
  
  ## point estimate under ramdomization
  lambda.estimate = sum(R[which(data$Z==1)] - R[which(data$Z==0)]) / sum(data$D[which(data$Z==1)] - data$D[which(data$Z==0)])
  
  ##  CI for effect ratio (lambda)
  effect.ratio.ci = function(lambda, data, R, Gamma, alpha){
    tau_hat = R[which(data$Z==1)] - lambda*data$D[which(data$Z==1)] - (R[which(data$Z==0)] - lambda*data$D[which(data$Z==0)])
    tau_hat_Gamma = (Gamma + 1)*((Gamma + 1)*tau_hat - (Gamma - 1)*abs(tau_hat))/(4*Gamma)
    tau_bar_Gamma = mean(tau_hat_Gamma)
    I = length(tau_hat)
    S_square_Gamma = sum((tau_hat_Gamma - tau_bar_Gamma)^2)/(I*(I-1)) 
    return( tau_bar_Gamma - qnorm(1-alpha)*sqrt(S_square_Gamma) ) 
  }
  
  ## Upper One-Sided Confidence Interval 
  upper.ci.lambda = c()
  for (lambda in seq(lambda.test.lb, lambda.test.ub, 0.001)) {
    if(effect.ratio.ci(lambda, data, R, Gamma, alpha/2)<0){
      upper.ci.lambda = c(upper.ci.lambda, lambda)
    }
  }
  
  ## Lower One-Sided Confidence Interval
  lower.ci.lambda = c()
  for (lambda in seq(lambda.test.lb, lambda.test.ub, 0.001)) {
    if(effect.ratio.ci(lambda, data, -R, Gamma, alpha/2)<0){
      lower.ci.lambda = c(lower.ci.lambda, -lambda)
    }
  }
  
  ### Return all results
  return(list(kappa.ci = ci.kappa, 
              effect.ratio.upper.ci =  upper.ci.lambda, 
              effect.ratio.lower.ci = lower.ci.lambda,
              effect.ratio.est = lambda.estimate
  ))
}