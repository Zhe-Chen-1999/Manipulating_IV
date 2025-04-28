################################################################################
##  This file contains the M1 primary analysis performed under two assumptions:
##  1. A randomization assumption 
##  2. A biased randomization assumption
################################################################################

source(helper_functions.R)
library(ggplot2)
library(tidyr)
library(ggpubr)

####################### Data Processing ################################

## read merged data for all years and all subgroups
dat_all = read.csv("dat_all.csv")

## IV
sum(is.na(dat_all$IV))
dat_all$Z = ifelse(dat_all$IV == "far", 1, 0)

## Treatment Received
dat_all$D = 1-dat_all$high_level_NICU

## Binary outcome: death_total = death (death after delivery) + death_fetal (death in the womb)
dat_all$death_total = dat_all$death_fetal + dat_all$death

## Composite outcome U = LOS * ind\{not die\} + C * ind\{die\}, where C is the top 1% among all survivors' LOS
C = quantile(dat_all$losI[which(dat_all$death_total!=1)], 0.99) 
dat_all$composite_outcome = ifelse(dat_all$death_total==0, dat_all$losI, C)

I = dim(dat_all)[1]/2 

## Estimated compliance rate
com_rate = sum(dat_all$D[which(dat_all$Z==1)] - dat_all$D[which(dat_all$Z==0)])/I 

## Table comparing group far from a high-level NICU (treatment) and near a high-level NICU (control) 
mean_treated = apply(dat_all[which(dat_all$Z==1), -which(names(dat_all) %in% c("IV","Filename","Filename.1"))], 2, mean)
mean_control = apply(dat_all[which(dat_all$Z==0), -which(names(dat_all) %in% c("IV","Filename","Filename.1"))], 2, mean)
mean_diff = mean_treated - mean_control

sd_treated = apply(dat_all[which(dat_all$Z==1), -which(names(dat_all) %in% c("IV","Filename","Filename.1"))], 2, stats::sd) 
sd_control = apply(dat_all[which(dat_all$Z==0), -which(names(dat_all) %in% c("IV","Filename","Filename.1"))], 2, stats::sd) 
pooled_sd = sqrt((sd_treated^2 + sd_control^2)/2)
std = mean_diff/pooled_sd

##  mean difference for survivors' losI
losI.dat = na.omit(dat_all)
losI_mean_treated = apply(losI.dat[which(losI.dat$Z==1), -which(names(losI.dat) %in% c("IV","Filename","Filename.1"))], 2, mean)
losI_mean_control = apply(losI.dat[which(losI.dat$Z==0), -which(names(losI.dat) %in% c("IV","Filename","Filename.1"))], 2, mean)

losI_mean_diff = losI_mean_treated - losI_mean_control
sd_treated = apply(losI.dat[which(losI.dat$Z==1), -which(names(losI.dat) %in% c("IV","Filename","Filename.1"))], 2, stats::sd) 
sd_control = apply(losI.dat[which(losI.dat$Z==0), -which(names(losI.dat) %in% c("IV","Filename","Filename.1"))], 2, stats::sd) 
pooled_sd = sqrt((sd_treated^2 + sd_control^2)/2)
losI_std = losI_mean_diff/pooled_sd


####################################################################
################### Randomization Inference ########################
####################################################################

## Mortality Outcome ##
M1_RI_death_res = biased.RI(gamma = 0, 
                            K0 = 0, 
                            K.lb = 0, K.ub = 0.03, step = 0.0005,
                            lambda.test.lb = -1, lambda.test.ub = 1,
                            data = dat_all, alpha = 0.05, R = dat_all$death_total)

saveRDS(M1_RI_death_res, "M1_RI_death_res.rds")

# Inference results for effect ratio 
M1_RI_death_res$effect.ratio.est # point estimate
c(range(M1_RI_death_res$effect.ratio.upper.ci)[1], range(M1_RI_death_res$effect.ratio.lower.ci)[2]) # 95% CI

## Composite Outcome ##
M1_RI_composite_res = biased.RI(gamma = 0, 
                                K0 = NULL, 
                                K.lb = 0, K.ub = 2, 
                                lambda.test.lb = -10, lambda.test.ub = 10,
                                data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

saveRDS(M1_RI_composite_res, "M1_RI_composite_res.rds")

# Inference results for effect ratio 
M1_RI_composite_res$effect.ratio.est # point estimate:-0.3763209
c(range(M1_RI_composite_res$effect.ratio.upper.ci)[1], range(M1_RI_composite_res$effect.ratio.lower.ci)[2]) # 95% CI

####################################################################
################### Biased Randomization Inference #################
####################################################################

## Mortality Outcome ##
death_res = biased.RI(CPT.max.Gamma = 1.17, 
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.0005,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, alpha = 0.05, R = dat_all$death_total)

saveRDS(death_res, "M1_BI_death_res.rds")

# Inference results for effect ratio 
death_res$effect.ratio.est # point estimate: 0.01278858
c(range(death_res$effect.ratio.upper.ci)[1], range(death_res$effect.ratio.lower.ci)[2]) # 95% CI


## Composite Outcome ##
composite_res = biased.RI(CPT.max.Gamma = 1.17, 
                          K0 = NULL, 
                          K.lb = 0, K.ub = 2, 
                          lambda.test.lb = -10, lambda.test.ub = 10,
                          data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

saveRDS(composite_res, "M1_BI_composite_res.rds")

# Inference results for effect ratio (lambda)
composite_res$effect.ratio.est 
c(range(composite_res$effect.ratio.upper.ci)[1], range(composite_res$effect.ratio.lower.ci)[2]) 

###########################################################################
### Confidence Interval (CI) Plot for SATE under Both Inference Schemes ###
###########################################################################

RI_plot_df = M1_RI_death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
              names_to = 'Bound') %>%
  mutate(Scheme = 'Randomization')

BI_plot_df = death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') %>%
  mutate(Scheme = 'Biased Randomization')

plot_df = rbind(RI_plot_df, BI_plot_df)

plot_df = plot_df %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

death_plot = ggplot(plot_df, aes(x = K, y = value, color = Bound, 
                    linetype = Scheme)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d(guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_linetype_manual('', values = c('Randomization' = 'dotted',
                                       'Biased Randomization' = 'solid'),
                        guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(x = expression(paste(K[1])),
       y = expression(paste("Mortality rate"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')


RI_plot_df2 = M1_RI_composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') %>%
  mutate(Scheme = 'Randomization')

BI_plot_df2 = composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') %>%
  mutate(Scheme = 'Biased Randomization')

plot_df2 = rbind(RI_plot_df2, BI_plot_df2)

plot_df2 = plot_df2 %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

composite_plot = ggplot(plot_df2, aes(x = K, y = value, color = Bound, 
                        linetype = Scheme)) + 
  geom_line(size = 2) +
  scale_color_viridis_d(guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_linetype_manual('', values = c('Randomization' = 'dotted',
                                       'Biased Randomization' = 'solid'),
                        guide = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "K",
       y = expression(paste("Composite outcome"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

death.ci.plot = ggarrange(death_plot, composite_plot, common.legend = TRUE, ncol=2)

ggsave(filename="M1_CI_plots.pdf", height = 7, width = 14, plot=death.ci.plot)

############################################################################
########### Confidence Interval (CI) Plots for Subgroup Analysis ###########
############################################################################
### BLACK SUBGROUP ###
## merge data
m1.long = data.frame()
m1.short = data.frame()

for (i in 1:10) {
  # Read far data
  df.long = read.csv(paste("m1_black_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  
  # Read near data
  df.short = read.csv(paste("m1_black_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
}

## read merged data
dat_all = rbind(m1.long, m1.short)

## IV
dat_all$Z = ifelse(dat_all$IV == "far", 1, 0)

## Treatment Received
dat_all$D = 1-dat_all$high_level_NICU

## Binary outcome: death_total = death (death after delivery) + death_fetal (death in the womb)
dat_all$death_total = dat_all$death_fetal + dat_all$death

## Composite outcome U = LOS * ind\{not die\} + C * ind\{die\}, where C is the top 1% among all survivors' LOS
C = quantile(dat_all$losI[which(dat_all$death_total!=1)], 0.99) 
dat_all$composite_outcome = ifelse(dat_all$death_total==0, dat_all$losI, C)

death_res = biased.RI(CPT.max.Gamma = 1.17, 
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.001,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, alpha = 0.05, R = dat_all$death_total)

## CI Plot for mortality rate
plot_df = death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df = plot_df %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.black.mortality = ggplot(plot_df, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(-0.01, 0.03)) +
  labs(x = expression(paste(K[1])),
       y = expression(paste("Mortality rate"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

composite_res = biased.RI(CPT.max.Gamma = 1.17, 
                          K0 = NULL, 
                          K.lb = 0, K.ub = 2, 
                          lambda.test.lb = -10, lambda.test.ub = 10,
                          data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

## CI Plot for composite outcome
plot_df2 = composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df2 = plot_df2 %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.black.comp = ggplot(plot_df2, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  labs(x = "K",
       y = expression(paste("Composite outcome"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

### HIGH-RISK NON-BLACK SUBGROUP ###
## merge data for all years and all subgroups
m1.long = data.frame()
m1.short = data.frame()

for (i in 1:10) {
  # Read far data
  df.long = read.csv(paste("m1_nonblack1_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  
  # Read near data
  df.short = read.csv(paste("m1_nonblack1_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
}

## read merged data
dat_all = rbind(m1.long, m1.short)

## IV
dat_all$Z = ifelse(dat_all$IV == "far", 1, 0)

## Treatment Received
dat_all$D = 1-dat_all$high_level_NICU

## Binary outcome: death_total = death (death after delivery) + death_fetal (death in the womb)
dat_all$death_total = dat_all$death_fetal + dat_all$death

## Composite outcome U = LOS * ind\{not die\} + C * ind\{die\}, where C is the top 1% among all survivors' LOS
C = quantile(dat_all$losI[which(dat_all$death_total!=1)], 0.99) 
dat_all$composite_outcome = ifelse(dat_all$death_total==0, dat_all$losI, C)

death_res = biased.RI(CPT.max.Gamma = 1.17, 
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.001,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, alpha = 0.05, R = dat_all$death_total)

## CI Plot for mortality rate
plot_df = death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df = plot_df %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.nonblack1.mortality = ggplot(plot_df, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(x = expression(paste(K[1])),
       y = expression(paste("Mortality rate"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

composite_res = biased.RI(CPT.max.Gamma = 1.17, 
                          K0 = NULL, 
                          K.lb = 0, K.ub = 2, 
                          lambda.test.lb = -10, lambda.test.ub = 10,
                          data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

## CI Plot for composite outcome
plot_df2 = composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df2 = plot_df2 %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.nonblack1.comp = ggplot(plot_df2, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  labs(x = "K",
       y = expression(paste("Composite outcome"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

### MEDIUM RISK NON-BLACK SUBGROUP ###
## merge data for all years and all subgroups
m1.long = data.frame()
m1.short = data.frame()

for (i in 1:10) {
  # Read far data
  df.long = read.csv(paste("m1_nonblack2_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  
  # Read near data
  df.short = read.csv(paste("m1_nonblack2_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
}

## read merged data
dat_all = rbind(m1.long, m1.short)

## IV
dat_all$Z = ifelse(dat_all$IV == "far", 1, 0)

## Treatment Received
dat_all$D = 1-dat_all$high_level_NICU

## Binary outcome: death_total = death (death after delivery) + death_fetal (death in the womb)
dat_all$death_total = dat_all$death_fetal + dat_all$death

## Composite outcome U = LOS * ind\{not die\} + C * ind\{die\}, where C is the top 1% among all survivors' LOS
C = quantile(dat_all$losI[which(dat_all$death_total!=1)], 0.99)
dat_all$composite_outcome = ifelse(dat_all$death_total==0, dat_all$losI, C)

death_res = biased.RI(CPT.max.Gamma = 1.17, 
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.001,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, alpha = 0.05, R = dat_all$death_total)

## CI Plot for mortality rate
plot_df = death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df = plot_df %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.nonblack2.mortality = ggplot(plot_df, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(x = expression(paste(K[1])),
       y = expression(paste("Mortality rate"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

composite_res = biased.RI(CPT.max.Gamma = 1.17, 
                          K0 = NULL, 
                          K.lb = 0, K.ub = 2, 
                          lambda.test.lb = -10, lambda.test.ub = 10,
                          data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

## CI Plot for composite outcome
plot_df2 = composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df2 = plot_df2 %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.nonblack2.comp = ggplot(plot_df2, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  labs(x = "K",
       y = expression(paste("Composite outcome"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

### LOW RISK NON-BLACK SUBGROUP
## merge data for all years and all subgroups
m1.long = data.frame()
m1.short = data.frame()

for (i in 1:10) {
  # Read far data
  df.long = read.csv(paste("m1_nonblack34_long_", 1994 + i, "_032724.csv", sep = ""))
  m1.long = rbind(m1.long, df.long)
  
  # Read near data
  df.short = read.csv(paste("m1_nonblack34_short_", 1994 + i, "_032724.csv", sep = ""))
  m1.short = rbind(m1.short, df.short)
}

## read merged data
dat_all = rbind(m1.long, m1.short)

## IV
dat_all$Z = ifelse(dat_all$IV == "far", 1, 0)

## Treatment Received
dat_all$D = 1-dat_all$high_level_NICU

## Binary outcome: death_total = death (death after delivery) + death_fetal (death in the womb)
dat_all$death_total = dat_all$death_fetal + dat_all$death

## Composite outcome U = LOS * ind\{not die\} + C * ind\{die\}, where C is the top 1% among all survivors' LOS
C = quantile(dat_all$losI[which(dat_all$death_total!=1)], 0.99)
dat_all$composite_outcome = ifelse(dat_all$death_total==0, dat_all$losI, C)

death_res = biased.RI(CPT.max.Gamma = 1.17, 
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.001,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, alpha = 0.05, R = dat_all$death_total)

## CI Plot for mortality rate
plot_df = death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df = plot_df %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.nonblack34.mortality = ggplot(plot_df, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(x = expression(paste(K[1])),
       y = expression(paste("Mortality rate"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

composite_res = biased.RI(CPT.max.Gamma = 1.17, 
                          K0 = NULL, 
                          K.lb = 0, K.ub = 2, 
                          lambda.test.lb = -10, lambda.test.ub = 10,
                          data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

## CI Plot for composite outcome
plot_df2 = composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df2 = plot_df2 %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p.nonblack34.comp = ggplot(plot_df2, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  labs(x = "K",
       y = expression(paste("Composite outcome"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

# Merge the plots
mortality.plots = ggarrange(M0_p.black.mortality, M0_p.nonblack1.mortality, 
                            M0_p.nonblack2.mortality, M0_p.nonblack34.mortality, 
                            common.legend = TRUE, nrow = 2, ncol=2)

ggsave(filename="CI_plots.mortality.pdf", height = 14, width = 14, plot=mortality.plots)

comp.plots = ggarrange(M0_p.black.comp, M0_p.nonblack1.comp, 
                       M0_p.nonblack2.comp, M0_p.nonblack34.comp, 
                       common.legend = TRUE, nrow = 2, ncol=2)

ggsave(filename="CI_plots.comp.pdf", height = 14, width = 14, plot=comp.plots)

###########################################################################
########################## Plot of Gamma Distribution #####################
###########################################################################
data = dat_all
gamma = log(1.17)/max(data$dif.travel.time[which(data$Z==1)] - data$dif.travel.time[which(data$Z==0)])
Gamma = exp(gamma*(data$dif.travel.time[which(data$Z==1)] - data$dif.travel.time[which(data$Z==0)]))
Gamma_df = data.frame(Gamma = Gamma)
Gamma.dist.plot = ggplot(Gamma_df, aes(x=Gamma)) + 
  geom_density(size = 2)+
  xlab( expression(paste(Gamma)))+
  theme_bw(base_size = 25)
ggsave(filename="Gamma_dist_plot.pdf", height = 5, width = 6, plot=Gamma.dist.plot)
