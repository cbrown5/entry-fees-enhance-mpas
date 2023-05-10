#Plot tourism analysis

library(visreg)
library(ggplot2)
library(mgcv)
library(DataGLMRepeat)
library(tmap)
library(sf)

rm(list = ls())
load("rls-dontshare/rls-tour_step_one_data-R1.RDA")

#
# Maps
#
#Fix one MPA that was coded incorrectly
MPAdat$has_fee[MPAdat$MPA == "Tulamben"] <- "1"

sMPAdat <- st_as_sf(MPAdat, coords = c("Lon", "Lat"))
# sMPAdat$sbio <- as.factor(sign(sMPAdat$bio_resids))
data("World")

#keep refining map

tm1 <- tm_shape(World) + 
  tm_fill(col = "grey80") +
  # tm_shape(World) +
  # tm_borders(col = "white") + 
  tm_shape(sMPAdat) + 
  tm_symbols(col = "has_fee", 
             shape = "has_fee",
             size = 0.1,
             alpha = 0.3,
             legend.shape.show = FALSE,
             legend.col.show = FALSE,
             palette = c("grey30", "red"),
             border.lwd = 0.2) +
  tm_add_legend(type = c("symbol"),
                labels = c("No fee", "Fee"),
                col = c("black", "red"),
                shape = c(21, 22))

tmap_save(tm1, "figures/fee-map.png",
          width = 4, height = 2.5,
          dpi = 600)

# Check out this to see if has_fee is clustered in space 
# which explains why when we add the SOS has_fee isn't important anymore
feemap <- ggplot(MPAdat) + aes(x = Lon, y = Lat) + theme_bw() + 
  borders("world",fill="black",colour="black") +
  geom_point(mapping = aes(col = has_fee,
                           shape = as.factor(sign(bio_resids))), 
             size = 4,
             alpha = 0.4)

ggsave("figures/feemap.png", feemap)

### Check fo spatial autocorrelation

# Look as the residuals in space
ggplot(mpadat2) + aes(x = Lon, y = Lat) + theme_bw() +
  borders("world",fill="black",colour="black") +
  geom_point(mapping = aes(col = sign(resid(tour_best$gam)),
                           size = abs(resid(tour_best$gam))),
             alpha = 0.2) +
  scale_colour_gradient(low = "purple",
                        high = "green", space = "Lab",
                        na.value = "red")

site_dist <- with(mpadat2, {
  
  Lon <- as.numeric(Lon)
  Lat <- as.numeric(Lat)
  
  geosphere::distm(cbind(Lon, Lat),
                   cbind(Lon, Lat),
                   fun = geosphere::distHaversine)/1000
})

#This Requires CB's semivariance function - doesnt seem like much spatial AC
sp_ac_df2 <- DataGLMRepeat::semivariance(xdists = site_dist, yresp = resid(tour_best$gam), 150)

ggplot(sp_ac_df2, aes(distances, moransI)) +
  geom_point()  +
  geom_line() +
  theme_bw() +
  ylab("moransI correlation coef")

# Effect size plots - Bayesian CIs ----------------------------------------

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*%  matrix(rnorm(m * n), m, n))
}
#1000 draws of the parameters from covariance matrix 
br <- rmvn(1000, coef(tour_best$gam), tour_best$gam$Vp)

#visreg - easy way to get prediction  dataframe
# with all values held at most common/median values
xfit <-  visreg(tour_best$gam, xvar = "NEOLI_Total", 
                partial = FALSE, plot = FALSE, by = "has_fee",
                data = mpadat2)

newdat1 <- xfit$fit
lp1 <- predict(tour_best$gam, newdata = newdat1, type = "lpmatrix")
B20mult <- matrix(NA, nrow = nrow(newdat1),
                  ncol = 1000)

#Simulate posterior
for (i in 1:1000){
  val1 <- lp1 %*% br[i,]  # start of decade
  B20mult[, i] <- 10^val1
}

CItrend <- data.frame(t(apply(B20mult, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  cbind(newdat1) %>%
  dplyr::filter(NEOLI_Total %in% c(1, 2,3,4,5))


g1 <- ggplot(filter(CItrend, has_fee == "0")) + 
  aes(x = NEOLI_Total, y = X50., color = NULL) +
  geom_linerange(aes(ymin = X25., ymax = X75.), alpha = 0.5, size = 5) + 
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.5, size = 1) + 
  geom_point(size = 2) +
  geom_hline(yintercept = 1)+ 
  theme_classic() + 
  ylab("Multiple of no MPA baseline") + 
  xlab("NEOLI") + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) + 
  ylim(0, 30)
g1

CItrend$MPAfee <- factor(CItrend$has_fee, labels = c("No fee", "Entry fee"))
g2 <- ggplot(filter(CItrend, NEOLI_Total == 4)) + 
  aes(x = MPAfee, y = X50., color = NULL) +
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 2) + 
  geom_linerange(aes(ymin = X2.5., ymax = X97.5.), size = 1) + 
  geom_point(size = 4) +
  theme_classic() + 
  geom_hline(yintercept =1)+ 
  ylab("Multiple of no MPA baseline") + 
  xlab("")+ 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))+ 
  ylim(0, 21)
g2

library(patchwork)

pw <- g1 | g2

pw2 <- pw  + plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 2, widths = c(2,1)) & 
  theme(plot.tag = element_text(size = 16)) 
pw2

ggsave(pw2, file = "figures/2020-10-29_RLS-tourism-effects.png",
       width = 8, height = 6)

#
# Sample sizes 
#

mpadat2 %>%
  ungroup() %>% select(MPA, has_fee) %>%
  distinct() %>%
  group_by(has_fee) %>%
  summarize(n())
length(unique(mpadat2$MPA))  

summary(tour_best$gam)

