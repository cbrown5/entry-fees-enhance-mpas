#Predict biomass gain 
rm(list = ls())

data_folder <- "rls-dontshare"
load(file = paste0(data_folder,  "/rls-tour_step_one_data-R1.RDA" ))
load(file = paste0(data_folder,  "/rls-tour_step_two_data-R1.RDA" ))

# Predict MPA log10_b20 based on the biovariable model
predicted_log10_B20 <- predict(bio_mod, newdata = MPAdat ,
                               exclude= c('s(Location, bs = "re")',
                                          's(Lat, Lon, bs = "sos", k = 100)'
                               ), 
                               se = TRUE)

bio_resids <- MPAdat$log10_B20 - predicted_log10_B20[["fit"]]

# Look at the distribution of the bio-residuals
hist(bio_resids) # looks normally distributed

MPAdat$bio_resids <- as.numeric(bio_resids)

save("MPAdat", "noMPAdat", "bio_mod" , file = paste0(data_folder,  "/rls-tour_step_three_data-R1.RDA" ))
