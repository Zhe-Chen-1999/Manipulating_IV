################################################################################
##  This file contains the M0 primary analysis under a randomization assumption 
################################################################################
library(ggplot2)
library(tidyr)
library(dplyr)

## read merged data for all years and all subgroups
dat_all = read.csv("data_all.csv") 

dat_all = dat_all %>%
  filter(!is.na(IV))

# IV dose
dat_all$Z = ifelse(dat_all$IV == "far", 1, 0)

## Treatment received
dat_all$D = 1-dat_all$high_level_NICU

## Binary outcome: death_total = death (death after delivery) + death_fetal (death in the womb)
dat_all$death_total = dat_all$death_fetal + dat_all$death

## Composite outcome U = LOS * ind\{not die\} + C * ind\{die\}, where C is the top 1% among all survivors' LOS
C = quantile(dat_all$losI[which(dat_all$death_total!=1)], 0.99) 
dat_all$composite_outcome = ifelse(dat_all$death_total==0, dat_all$losI, C)

I = dim(dat_all)[1]/2 

## Estimated compliance rate
com_rate = sum(dat_all$D[which(dat_all$Z==1)] - dat_all$D[which(dat_all$Z==0)])/I 

##############################################################################################
################## Randomization inference for the mortality outcome #########################
##############################################################################################

death_res = biased.RI(gamma = 0,
                      K0 = 0, 
                      K.lb = 0, K.ub = 0.03, step = 0.001,
                      lambda.test.lb = -1, lambda.test.ub = 1,
                      data = dat_all, alpha = 0.05, R = dat_all$death_total)

saveRDS(death_res, "M0_RI_death_res.rds")

## Inference results for effect ratio 
death_res$effect.ratio.est # point estimate
c(range(death_res$effect.ratio.upper.ci)[1], range(death_res$effect.ratio.lower.ci)[2]) # 95% CI

##############################################################################################
################ Randomization Inference for the composite outcome ###########################
##############################################################################################

composite_res = biased.RI(gamma = 0, 
                          K0 = NULL, 
                          K.lb = 0, K.ub = 2, 
                          lambda.test.lb = -10, lambda.test.ub = 10,
                          data = dat_all, alpha = 0.05, R = dat_all$composite_outcome)

saveRDS(composite_res, "M0_RI_composite_res.rds")

## Inference results for effect ratio
composite_res$effect.ratio.est # point estimate
c(range(composite_res$effect.ratio.upper.ci)[1], range(composite_res$effect.ratio.lower.ci)[2]) # 95% CI

##############################################################################################
################ SATE CI plot for M0 under a randomization assumption (Figure S4) ############
##############################################################################################

plot_df = death_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df = plot_df %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p = ggplot(plot_df, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  labs(x = expression(paste(K[1])),
       y = expression(paste("Mortality rate"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

plot_df2 = composite_res$kappa.ci %>% 
  pivot_longer(cols = c('LB.lower.limit', 'UB.upper.limit'),
               names_to = 'Bound') 

plot_df2 = plot_df2 %>%
  mutate(Bound = ifelse(Bound == 'UB.upper.limit', 'Upper confidence limit',
                        'Lower confidence limit'))

M0_p2 = ggplot(plot_df2, aes(x = K, y = value, color = Bound)) + 
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  labs(x = "K",
       y = expression(paste("Composite outcome"))) +
  theme_bw(base_size = 25) +
  theme(legend.position = 'top')

death.ci.plot = ggarrange(M0_p, M0_p2, common.legend = TRUE, ncol=2)
ggsave(filename="M0_CI_plots.pdf", height = 7, width = 14, plot=death.ci.plot)
