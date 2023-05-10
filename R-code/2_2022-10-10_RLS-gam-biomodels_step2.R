#Fit models for no MPA data

library(tidyverse)
library(mgcv)
library(viridis)
library(visreg)
library(maps)
library(parallel)
library(patchwork)

# Location of data and model folder for analysis (add your own)
data_folder <- "rls-dontshare"

load(file = paste0(data_folder,  "/rls-tour_step_one_data-R1.RDA" ))




# 4th root transformation seems like a good idea 
#because that downweights the few extreme outliers

with(noMPAdat,{
  
  hist(x = Kernel2015)
  hist(x = sqrt(Kernel2015))
  hist(x = Kernel2015^(1/4))
  hist(x = log(Kernel2015 + min(Kernel2015)/2))
  
})


# Construct GAMs biodata ----------------------------------------------

# Not including dissolved oxygen, parmean or min and max SST due to their confounded
# relationships

knotts_base <- 15
  
## Model to compare correlation structure
f1 <- function (){
  m2 <- gam(log10_B20 ~ s(Kernel2015_4throot, k = knotts_base) + 
                 s(salinty, k = knotts_base) +
                 s(chlomean, k = knotts_base) + 
                 s(Phosphate, k = knotts_base) + s(Silicate, k = knotts_base) + 
                 s(SSTrange, k = knotts_base) + s(SSTMean, k = knotts_base) + 
                 s(Nitrate, k = knotts_base) + 
                 s(Calcite, k = knotts_base) + 
                 s(Location, bs = "re"), method = "ML",
               data = noMPAdat)
}

f2 <- function (){ 
  m2.0 <- gam(log10_B20 ~ s(Kernel2015_4throot, k = knotts_base) + 
                 s(salinty, k = knotts_base) +
                 s(chlomean, k = knotts_base) + 
                 s(Phosphate, k = knotts_base) + s(Silicate, k = knotts_base) + 
                 s(SSTrange, k = knotts_base) + s(SSTMean, k = knotts_base) + 
                 s(Nitrate, k = knotts_base) + 
                 s(Calcite, k = knotts_base) + 
                 s(Location, bs = "re") +
                 s(Lat, Lon, bs = "sos", k = 50), method = "ML",
               data = noMPAdat)
}

f3 <- function (){ # more knotts
  m2.1 <- gam(log10_B20 ~ s(Kernel2015_4throot, k = knotts_base) + 
                 s(salinty, k = knotts_base) +
                 s(chlomean, k = knotts_base) + 
                 s(Phosphate, k = knotts_base) + s(Silicate, k = knotts_base) + 
                 s(SSTrange, k = knotts_base) + s(SSTMean, k = knotts_base) + 
                 s(Nitrate, k = knotts_base) + 
                 s(Calcite, k = knotts_base) + 
                 s(Location, bs = "re") +
                 s(Lat, Lon, bs = "sos", k = 100), method = "ML",
               data = noMPAdat)
}

# How many cores are available
Ncores <- detectCores()
Ncores

# Formatting to evaluate models in Parallel
ExpressionVect <- c(substitute(f1()),
                    substitute(f2()),
                    substitute(f3()),
                    substitute(f4()),
                    substitute(f5()))

## Run Models in Parallel
system.time({model.list <- mclapply(ExpressionVect, eval, mc.cores = (Ncores - 2))})

# Manually unlist models
m2 <- model.list[[1]]
m2.0 <- model.list[[2]]
m2.1 <- model.list[[3]]

# Check which correlation structure is preferred
AIC(m2, m2.0, m2.1)

# Results appear to be very similar
summary(m2)
summary(m2.0)
summary(m2.1)

gam.check(m2)

summary(m2) 

# # Stepwise selection
m2.2 <- update(m2, . ~ . - s(chlomean, k = knotts_base))
m2.3 <- update(m2.2, . ~ . - s(Silicate, k = knotts_base))
m2.4 <- update(m2.3, . ~ . - s(Phosphate, k = knotts_base))
m2.5 <- update(m2.3, . ~ . - s(Nitrate, k = knotts_base))
m2.6 <- update(m2.3, . ~ . - s(Calcite, k = knotts_base)) # Worsens AIC so stop

AIC(m2, m2.2, m2.3, m2.4, m2.5, m2.6)

# Convert 1df terms to parametric terms
m2.7 <- update(m2.3, . ~ . - s(Kernel2015_4throot, k = knotts_base) +  
                 Kernel2015_4throot -
                 s(salinty, k = knotts_base) + salinty - 
                 s(Phosphate, k = knotts_base) +  Phosphate -
                 s(Nitrate, k = knotts_base) +  Nitrate )

AIC(m2.3, m2.7)

# Increase the number of knotts
m2.8 <- update(m2.7, . ~ . - s(SSTrange, k = knotts_base) + s(SSTrange, k = 25) - 
                 s(SSTMean, k = knotts_base) + s(SSTMean, k = 25) -
                 s(Calcite, k = knotts_base) + s(Calcite, k = 25) )

gam.check(m2.7)
gam.check(m2.8)

AIC(m2.7, m2.8) # seems like we have enough already

# Add SOS
m2.9 <- update(m2.7, . ~ . + s(Lat, Lon, bs = "sos", k = 100))

AIC(m2.7, m2.9) # same same

model_decription <- c("base model" = "m2",
                      "without Chlorophyll" = "m2.2",
                      "without Chlorophyll and Silicate" = "m2.3",
                      "without Chlorophyll, Silicate and Phosphate" = "m2.4",
                      "without Chlorophyll, Silicate and Nitrate" = "m2.5",
                      "without Chlorophyll, Silicate and Calcite" = "m2.6",
                      "with SOS and without Chlorophyll and Silicate" = "m2.9")

model_selection_biomodel <- NULL

# Extract the information that we want out of the biomodels
for (i in seq_along(model_decription)){
  
  temp_model <- get(model_decription[[i]])
  # Add results
  model_selection_biomodel <- rbind(model_selection_biomodel,
    data.frame(Model = names(model_decription)[[i]],
             df = summary(temp_model)[["edf"]] %>% sum(),
             ML = - summary(temp_model)[["sp.criterion"]] %>% as.vector() %>% round(2),
             AIC = AIC(temp_model),
             Adjusted_R_square = round(summary(temp_model)$r.sq, 4))
  )
}

write_csv(model_selection_biomodel,
          file = "figures/model_selection_biomodel-R1.csv")

### Check fo spatial autocorrelation
# Look as the residuals in space
p1 <- ggplot(noMPAdat) + aes(x = Lon, y = Lat) + theme_bw() +
  borders("world",fill="black",colour="black") +
  geom_point(mapping = aes(col = sign(resid(m2.7))),
                           #size = abs(resid(m2.7))),
             alpha = 0.5, size = 0.5) +
  scale_colour_gradient(low = "purple",
                        high = "darkgreen", space = "Lab",
                        na.value = "red") +
  guides(col=guide_legend(title="Residual"))

p2 <- ggplot(noMPAdat) + aes(x = Lon, y = Lat) + theme_bw() +
  borders("world",fill="black",colour="black") +
  geom_point(mapping = aes(col = sign(resid(m2.9)),
                           size = abs(resid(m2.9))),
             alpha = 0.2) +
  scale_colour_gradient(low = "purple",
                        high = "darkgreen", space = "Lab",
                        na.value = "red")

# Practically the same
p1/p2

#For supplemntal
ggsave(p1, filename = "figures/resids-stage1-gamm.png",
       width = 8, height = 4)

site_dist <- with(noMPAdat, {
  
  Lon <- as.numeric(Lon)
  Lat <- as.numeric(Lat)
  
  geosphere::distm(cbind(Lon, Lat), 
                   cbind(Lon, Lat), 
                   fun = geosphere::distHaversine)/1000
})

#This Requires CB's semivariance function - doesnt seem like much spatial AC
sp_ac_df <- DataGLMRepeat::semivariance(xdists = site_dist, yresp = resid(m2.7), 100)
sp_ac_df2 <- DataGLMRepeat::semivariance(xdists = site_dist, yresp = resid(m2.9), 100)

p3 <- ggplot(sp_ac_df, aes(distances, moransI)) +
  geom_point()  +
  geom_line() +
  theme_bw() +
  ylab("moransI correlation coef")

p4 <- ggplot(sp_ac_df2, aes(distances, moransI)) +
  geom_point()  +
  geom_line() +
  theme_bw() +
  ylab("moransI correlation coef")

p3/p4
ggsave(p3, filename = "figures/Morans I.png",
       width = 8, height = 4)

# Model summary
bio_mod <- m2.7

gam.check(bio_mod)
summary(bio_mod)

par(mfrow = c(3,2))
plot(bio_mod)

save(bio_mod, noMPAdat, 
  file = paste0(data_folder,  "/rls-tour_step_two_data-R1.RDA" ))


