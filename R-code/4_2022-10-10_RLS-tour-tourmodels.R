#Run biomass gain models

library(tidyverse)
library(mgcv)
library(viridis)
library(visreg)
library(maps)
library(parallel)
library(patchwork)


# Location of data and model folder for analysis (add your own)
data_folder <- "rls-dontshare"

load(file = paste0(data_folder,  "/rls-tour_step_three_data-R1.RDA" ))

# Make sure all models have the same data
colSums(is.na(MPAdat))
MPAdat <- MPAdat %>% filter(!is.na(TI)) 


# Run the gam models ----------------------------------------------------------

#Make sure bio_resids is numeric, not an array, the way it was
# created from a list was causing errors, as it was an array

## Compare model SOS
tour_gam1 <- gam(bio_resids ~ s(NEOLI_Total,  k = 5, by = has_fee) + 
                    s(TI, k = 5, by = has_fee) +
                    has_fee +
                    s(MPA, bs = "re"), data = MPAdat, method = "ML")

tour_gam1.0 <- update(tour_gam1, . ~ . + s(Lat, Lon, bs = "sos", k = 50)) # This is enough
tour_gam1.1 <- update(tour_gam1, . ~ . + s(Lat, Lon, bs = "sos", k = 100))
tour_gam1.1k <- update(tour_gam1, . ~ . + s(Lat, Lon, bs = "sos", k = 200))

AIC(tour_gam1, tour_gam1.0, tour_gam1.1, tour_gam1.1k)#, #tour_gam1.1k, tour_gam1.1m, tour_gam1.1l)
#all same same!

## Stepwise selection to find the best model

tour_gam1.2 <- gam(bio_resids ~ s(NEOLI_Total,  k = 5, by = has_fee) + 
                     s(TI, k = 5) +
                     has_fee +
                     s(MPA, bs = "re"), data = MPAdat, method = "ML")
AIC(tour_gam1, tour_gam1.2) #Use linear TI
tour_gam1.3 <- update(tour_gam1.2, . ~ . - TI ) 
AIC(tour_gam1.2, tour_gam1.3) # Keep linear TI


tour_gam1.2 <- update(tour_gam1, . ~ . - s(TI, k = 5)+ TI) 
AIC(tour_gam1, tour_gam1.2) #Use linear TI
tour_gam1.3 <- update(tour_gam1.2, . ~ . - TI ) 
AIC(tour_gam1.2, tour_gam1.3) # Keep linear TI

tour_gam1.4 <- gam(bio_resids ~ 
                     NEOLI_Total*has_fee + 
                     TI +
                     s(MPA, bs = "re"), data = MPAdat, method = "ML")
AIC(tour_gam1.3, tour_gam1.4) #drop non-linear interaction

#additive model
tour_gam1.5 <- gam(bio_resids ~ 
                     NEOLI_Total + 
                     has_fee + 
                     TI +
                     s(MPA, bs = "re"), data = MPAdat, method = "ML")
AIC(tour_gam1.5, tour_gam1.4)
#Interaction slightly preferred

tour_mod <- tour_gam1.4
gam.check(tour_mod)
summary(tour_mod)

visreg(tour_mod, rug = FALSE)
plot(tour_mod)


model_decription <- c("All terms" = "tour_gam1",
                      "Linear tourism" = "tour_gam1.2",
                      "No tourism" = "tour_gam1.3",
                      "Linear NEOLI, fee interaction" = "tour_gam1.4",
                      "All terms additive" = "tour_gam1.5"
)

model_selection_tourmodel <- NULL

# Extract the information that we want out of the biomodels
for (i in seq_along(model_decription)){
  
  temp_model <- get(model_decription[[i]])
  # Add results
  model_selection_tourmodel <- rbind(model_selection_tourmodel,
                                    data.frame(Model = names(model_decription)[[i]],
                                               df = summary(temp_model)[["edf"]] %>% sum() %>% round(2),
                                               ML = - summary(temp_model)[["sp.criterion"]] %>% as.vector() %>% round(2),
                                               AIC = AIC(temp_model) %>% round(2),
                                               Adjusted_R_square = round(summary(temp_model)$r.sq, 4))
  )
}

write_csv(model_selection_tourmodel,
          file = "figures/model_selection_tourmodel-R1.csv")
          # path = "~/Dropbox/Dropbox_active_projects/RLS_gams_and_data/model_selection_tourmodel.csv")


#
# Compare models ----------------------------------------------------------
#

visreg(tour_mod, xvar = "NEOLI_Total", by = "has_fee")
visreg(tour_gam1.5, xvar = "NEOLI_Total", by = "has_fee")

visreg(tour_gam1.5, xvar = "has_fee")
visreg(tour_mod, xvar = "has_fee")

# Check for spatial autocorrelation ---------------------------------------


# Look as the residuals in space
p1 <- ggplot(MPAdat) + aes(x = Lon, y = Lat) + theme_bw() +
  borders("world",fill="black",colour="black") +
  geom_jitter(mapping = aes(col = sign(resid(tour_mod)),
                           #size = abs(resid(tour_mod))
                           ),
             alpha = 0.5 , width = 5,
             height = 5) +
  scale_colour_gradient(low = "purple",
                        high = "green", space = "Lab",
                        na.value = "red")

p1

## Look at Moran's I
site_dist <- with(MPAdat, {
  
  Lon <- as.numeric(Lon)
  Lat <- as.numeric(Lat)
  
  geosphere::distm(cbind(Lon, Lat),
                   cbind(Lon, Lat),
                   fun = geosphere::distHaversine)/1000
})

#This Requires CB's semivariance function - doesnt seem like much spatial AC
sp_ac_df1 <- DataGLMRepeat::semivariance(xdists = site_dist, yresp = resid(tour_mod), 150)
#sp_ac_df2 <- semivariance(xdists = site_dist, yresp = resid(tour_mod_comp), 150)

p3 <- ggplot(sp_ac_df1, aes(distances, moransI)) +
  geom_point()  +
  geom_line() +
  theme_bw() +
  ggtitle("Tour_mod")

p3


# Residual plots ----------------------------------------------------------

# QQplot
png(filename = "figures/2022-10-10_tourism_qqplot.png", width = 12,
    height = 11, units = "cm", res = 300)
par(mar = c(4, 4, 1,1))
qqnorm(residuals(tour_mod, type = "deviance"), main = NULL)
qqline(residuals(tour_mod, type = "deviance"))
dev.off()

# Residual vs fitted plot
png(filename = "figures/2022-10-10_tourism_residual-vs-fitted.png", width = 12,
    height = 11, units = "cm", res = 300)
par(mar = c(4, 4, 1,1))
plot(y = residuals(tour_mod, type = "deviance"), x = fitted(tour_mod), 
     xlab = "Linear predictor", ylab = "Deviance residuals")
dev.off()

#
# Save results for plotting later 
#

save(tour_mod, tour_gam1.5, file = "rls-dontshare/biomass-gain-models.rda")

#
# Effect size plots - Bayesian CIs ----------------------------------------
#

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*%  matrix(rnorm(m * n), m, n))
}
#1000 draws of the parameters from covariance matrix 
br <- rmvn(1000, coef(tour_mod), tour_mod$Vp)

#visreg - easy way to get prediction  dataframe
# with all values held at most common/median values
xfit <-  visreg(tour_mod, xvar = "NEOLI_Total", partial = FALSE, plot = FALSE, by = "has_fee")
newdat1 <- xfit$fit %>% filter((NEOLI_Total %% 1) == 0)
lp1 <- predict(tour_mod, newdata = newdat1, type = "lpmatrix")
B20mult <- matrix(NA, nrow = nrow(newdat1),
                  ncol = 1000)

#Simulate posterior
for (i in 1:1000){
  val1 <- lp1 %*% br[i,]  # start of decade
  B20mult[, i] <- 10^val1
}

CItrend <- data.frame(t(apply(B20mult, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat1)
CItrend$MPAfee <- factor(CItrend$has_fee, levels = c(0,1), labels = c("No fee", "Entry fee"))

g1 <- ggplot(filter(CItrend, MPAfee == "No fee")) + 
  aes(x = NEOLI_Total, y = X50.) +
  #geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.5) + 
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) + 
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2, 
                 position = position_dodge(0.25)) +
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1, 
                 position = position_dodge(0.25)) +
  geom_point(size = 3, position = position_dodge(0.25)) +
  geom_hline(yintercept = 1)+ 
  theme_classic() + labs(colour="", fill="") +
  ylab("Multiple of no MPA baseline") + 
  xlab("NEOLI") + scale_y_continuous(breaks = seq(0,40, 5), limits = c(0, 22)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), 
        legend.position = c(.90, .90)) 
g1

g2 <- ggplot(filter(CItrend, NEOLI_Total == 3)) + 
  aes(x = MPAfee, y = X50.) +
  #geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.5) + 
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5) + 
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2, 
                 position = position_dodge(0.25)) +
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1, 
                 position = position_dodge(0.25)) +
  geom_point(size = 3, position = position_dodge(0.25)) +
  geom_hline(yintercept = 1)+ 
  theme_classic() + labs(colour="", fill="") +
  ylab(NULL) + 
  xlab("Fee type") + scale_y_continuous(breaks = seq(0,10, 2), limits = c(0, 10)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), 
        legend.position = c(.90, .90)) 
g2


# g2 <- ggplot(filter(CItrend, NEOLI_Total == 4)) + 
#   aes(x = MPAfee, y = X50., color = NULL) +
#   geom_linerange(aes(ymin = X25., ymax = X75.), size = 2) + 
#   geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1) + 
#   geom_point(size = 4) +
#   theme_classic() + 
#   geom_hline(yintercept =1)+ 
#   ylab("Multiple of no MPA baseline") + 
#   xlab("")+ 
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=16))+ 
#   ylim(0, 130)
# g2

#visreg - easy way to get prediction  dataframe
# with all values held at most common/median values
xfit2 <-  visreg(tour_mod, xvar = "TI", partial = FALSE, plot = FALSE)
newdat2 <- xfit2$fit
lp2 <- predict(tour_mod, newdata = newdat2, type = "lpmatrix")
B20mult2 <- matrix(NA, nrow = nrow(newdat2),
                  ncol = 1000)

#Simulate posterior
for (i in 1:1000){
  val1 <- lp2 %*% br[i,]  # start of decade
  B20mult2[, i] <- 10^val1
}

CItrend2 <- data.frame(t(apply(B20mult2, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat2) %>% filter((TI %% 1) == 0)

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

library(patchwork)

pw <- g1 | g2 | g3

pw2 <- pw  + plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 3, widths = c(1.5,1,1.5)) & 
  theme(plot.tag = element_text(size = 16)) 
pw2

ggsave(file = "figures/2021-04-23_RLS-tourism-effects.png",
       width = 12, height = 6, dpi = 300)

#
# Sample sizes 
#

MPAdat %>%
  ungroup() %>% select(MPA, has_fee) %>%
  distinct() %>%
  group_by(has_fee) %>%
  summarize(n())
length(unique(MPAdat$MPA))  

#
# Explore MPAs with Fees 
#
