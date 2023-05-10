#Runs step two, just for the best model
# see 2021-09-23_RLS-gam-biomodels_step2.R 
# for model selection 

library(DataGLMRepeat)
library(tidyverse)
library(mgcv)
library(viridis)
library(visreg)
library(maps)
library(parallel)
library(patchwork)

# Location of data and model folder for analysis (add your own)
data_folder <- "../rls-dontshare" #case_when(Sys.info()["user"] == "" ~ "",
               #          TRUE ~ "")



load(file = paste0(data_folder,  "/rls-tour_step_one_data-R1.RDA" ))


par(mfrow = c(2,2))

with(noMPAdat,{
  
  hist(x = B20)
  hist(x = sqrt(B20))
  hist(x = B20^(1/4))
  hist(x = log10_B20)
  
})

with(noMPAdat,{
  
  hist(x = Kernel2015)
  hist(x = sqrt(Kernel2015))
  hist(x = Kernel2015^(1/4))
  hist(x = log(Kernel2015 + min(Kernel2015)/2))
  
})


# Construct GAMs biodata ----------------------------------------------


knotts_base <- 15
m2.7 <- gam(log10_B20 ~ 
            Kernel2015_4throot +
            salinty +
            Phosphate + 
            Nitrate + 
            s(SSTrange, k = knotts_base) + 
            s(SSTMean, k = knotts_base) + 
            s(Calcite, k = knotts_base) + 
            s(Location, bs = "re"), method = "ML",
          data = noMPAdat)

### Check fo spatial autocorrelation

# Look as the residuals in space
p1 <- ggplot(noMPAdat) + aes(x = Lon, y = Lat) + theme_bw() +
  borders("world",fill="black",colour="black") +
  geom_point(mapping = aes(col = sign(resid(m2.7)),
                           size = abs(resid(m2.7))),
             alpha = 0.2) +
  scale_colour_gradient(low = "purple",
                        high = "darkgreen", space = "Lab",
                        na.value = "red")
p1

site_dist <- with(noMPAdat, {
  
  Lon <- as.numeric(Lon)
  Lat <- as.numeric(Lat)
  
  geosphere::distm(cbind(Lon, Lat), 
                   cbind(Lon, Lat), 
                   fun = geosphere::distHaversine)/1000
})

#This Requires CB's semivariance function - doesnt seem like much spatial AC
sp_ac_df <- semivariance(xdists = site_dist, yresp = resid(m2.7), 100)

p3 <- ggplot(sp_ac_df, aes(distances, moransI)) +
  geom_point()  +
  geom_line() +
  theme_bw() +
  ylab("moransI correlation coef")
p3


# Model summary
bio_mod <- m2.7

gam.check(bio_mod)
summary(bio_mod)

par(mfrow = c(2,2))
plot(bio_mod)

# save(bio_mod, noMPAdat, MPAdat,
  # file = paste0(data_folder,  "/rls-tour_step_one_data-R1.RDA" ))


# Predict MPA data --------------------------------------------------------


# Predict MPA log10_b20 based on the biovariable model
predicted_log10_B20 <- predict(bio_mod, newdata = MPAdat ,
                               exclude= c(#'s(Location, bs = "re")',
                                          's(Lat, Lon, bs = "sos", k = 100)'
                               ), 
                               se = TRUE)

bio_resids <- MPAdat$log10_B20 - predicted_log10_B20[["fit"]]

# Look at the distribution of the bio-residuals
par(mfrow = c(1,1))
hist(bio_resids) # looks normally distributed 

MPAdat$bio_resids <- as.numeric(bio_resids)

save("MPAdat", "noMPAdat", "bio_mod" , file = paste0(data_folder,  "/rls-tour_step_two_data-R1.RDA" ))
