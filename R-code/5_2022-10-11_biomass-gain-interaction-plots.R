#
# Effect size plots - Bayesian CIs 
#
#Model including interaction

library(patchwork)
library(DataGLMRepeat)
library(tidyverse)
library(mgcv)
library(viridis)
library(visreg)
library(maps)
library(parallel)
library(patchwork)

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*%  matrix(rnorm(m * n), m, n))
}

load("rls-dontshare/biomass-gain-models.rda")
load("rls-dontshare/rls-tour_step_three_data-R1.RDA")

# ------------ 
# Predictions from additive model 
# ------------ 

#1000 draws of the parameters from covariance matrix 
br <- rmvn(1000, coef(tour_mod), tour_mod$Vp)

#
# Predictions for tourism
#
xfit2 <-  visreg(tour_mod, xvar = "TI", partial = FALSE, plot = FALSE)
newdat2 <- xfit2$fit
newdat2$NEOLI_Total <- 1
lp2 <- predict(tour_mod, newdata = newdat2, type = "lpmatrix")
#set RE to zero 
lp2[,28] <- 0
B20mult2 <- matrix(NA, nrow = nrow(newdat2),
                  ncol = 1000)
B20diff <- rep(NA, 1000)
#Simulate posterior
for (i in 1:1000){
  val1 <- lp2 %*% br[i,]  # start of decade
  B20mult2[, i] <- 10^val1
  B20diff[i] <- 10^(val1[nrow(newdat2)] - val1[1])
}


CItrend2 <- data.frame(t(apply(B20mult2, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat2) %>% filter((TI %% 1) == 0)
sum(B20mult2[1,]>B20mult2[nrow(newdat2),])/1000
quantile(B20diff,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))


g3 <- ggplot(CItrend2) +
  aes(x = TI, y = X50.) +
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) +
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2, 
                 position = position_dodge(0.25)) +
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1, 
                 position = position_dodge(0.25)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_hline(yintercept = 1) +
  ylab(NULL) +
  xlab("Tourism index") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))+
  scale_y_continuous(breaks = seq(0,10, 2), limits = c(0, 10.5))
g3

# ------------ 
# Predictions from interactive model 
# ------------ 

#1000 draws of the parameters from covariance matrix 
br <- rmvn(1000, coef(tour_mod), tour_mod$Vp)

#visreg - easy way to get prediction  dataframe
# with all values held at most common/median values
xfit <-  visreg(tour_mod, xvar = "NEOLI_Total", partial = FALSE, 
                plot = FALSE, by = "has_fee")
newdat1 <- xfit$fit %>% filter((NEOLI_Total %% 1) == 0)
lp1 <- predict(tour_mod, newdata = newdat1, type = "lpmatrix")
#set RE to zero 
lp1[,grepl("MPA", colnames(lp1))] <- 0
#find matching rows so can see B20 is conditionally higher with a fee
ifee0 <- match(subset(newdat1, has_fee == 0)$NEOLI_Total, 
               subset(newdat1, has_fee == 1)$NEOLI_Total)
ifee1 <-  which(newdat1$has_fee == 1)
#Simulate posterior
B20mult <-  matrix(NA, nrow = nrow(newdat1),
                   ncol = 1000)
res <- res2 <- matrix(NA, nrow = nrow(newdat1)/2,
                      ncol = 1000)
res3 <- res4 <- rep(NA, nrow(newdat1))

for (i in 1:1000){
  val1 <- lp1 %*% br[i,]  
  B20mult[, i] <- 10^val1
  #Fees
  res[,i] <- val1[ifee0] < val1[ifee1] #Is it higher when there is a fee?
  res2[,i] <- 10^(val1[ifee1] - val1[ifee0]) #how much higher? 
  #NEOLI
  res3[i] <- val1[1] < val1[5] #prob its higher 
  res4[i] <- 10^(val1[5])/10^(val1[1]) #how much higher? 
}

CItrend <- data.frame(t(apply(B20mult, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat1)
CItrend$MPAfee <- factor(CItrend$has_fee, levels = c(0,1), labels = c("No fee", "Entry fee"))
fee_probCIs <- data.frame(prob_fee = apply(res, 1, sum)/1000) %>%
  cbind(newdat1)
fee_diff_magCIs <- data.frame(t(apply(res2, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat1)
neoli_probCIs <- data.frame(prob_neoli = sum(res3)/1000) %>%
  cbind(newdat1)
fee_diff_dat <- inner_join(fee_probCIs, fee_diff_magCIs) %>%
  inner_join(neoli_probCIs)

#effect multiples for NEOLI
quantile(res4, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

#Biomass benefit at NEOLI 1 and 5, both with fees
quantile(B20mult[6,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(B20mult[10,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

#Biomass benefit at NEOLI 1 and 5, both without fees
quantile(B20mult[1,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(B20mult[5,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# Probability NEOLI 5 > NEOLI 1
sum(B20mult[5,] > B20mult[1,])/1000

g1 <- ggplot(CItrend) +
  aes(x = NEOLI_Total, y = X50., colour = MPAfee) +
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2, 
                 position = position_dodge(0.25)) +
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1, 
                 position = position_dodge(0.25)) +
  geom_point(size = 3, position = position_dodge(0.25)) +
  geom_hline(yintercept = 1)+ 
  theme_classic() + labs(colour="", fill="") +
  ylab("Biomass gain \n (multiple of no MPA baseline)") + 
  xlab("NEOLI") + scale_y_continuous(breaks = seq(0,40, 5), 
                                     limits = c(0, 20)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), 
        legend.text = element_text(size=14),
        legend.position = c(.80, .90))  +
  scale_colour_manual(values = c("black", "grey50"))
g1


g2 <- ggplot(filter(CItrend, NEOLI_Total == 1)) + 
  aes(x = MPAfee, y = X50.) +
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2, 
                 position = position_dodge(0.25)) +
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1, 
                 position = position_dodge(0.25)) +
  geom_point(size = 3, position = position_dodge(0.25)) +
  geom_hline(yintercept = 1)+ 
  theme_classic() + labs(colour="", fill="") +
  ylab(NULL) + 
  xlab("Fee type") + scale_y_continuous(breaks = seq(0,10, 2), 
                                        limits = c(0, 10)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), 
        legend.position = c(.90, .90)) 
g2


## Effect multiples for fee at each NEOLI level
sum(res[1,])/1000
sum(res[2,])/1000
sum(res[3,])/1000
sum(res[4,])/1000
sum(res[5,])/1000

#Comparison of with and without fees for each NEOLI level
 thquantile(res2[1,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(res2[2,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(res2[3,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(res2[5,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))


# ------------ 
# Figures 
# ------------ 
pw <- g1 |  g3

pw2 <- pw  + plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 2, widths = c(1.5,1.5)) & 
  theme(plot.tag = element_text(size = 16)) 
pw2

ggsave(file = "figures/2022-10-11_RLS-tourism-effects.png",
       width = 8, height = 4, dpi = 300)

#
# Sample sizes 
#

MPAdat %>%
  ungroup() %>% select(MPA, has_fee) %>%
  distinct() %>%
  group_by(has_fee) %>%
  summarise(n())
length(unique(MPAdat$MPA))  

#
# Probabilities 
#

