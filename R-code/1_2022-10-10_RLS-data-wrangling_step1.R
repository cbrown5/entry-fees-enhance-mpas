#Data wrangling to prepare databases 

library(tidyverse)

# Location of rls-data folder (add your own)
data_fp <- ""
# Location of data and model folder for analysis (add your own)
data_folder <- "rls-dontshare"

# Get data into R
#Covariates 
biodat <- read.csv(paste0(data_fp,"rls-dontshare/RLS_with_BIO.csv"), header = TRUE)

#Tourism status
tdat <- read.csv(paste0(data_fp,"rls-dontshare/MPAs-for-tourism-review-GE_NAD_RSS_GE.csv"))

#Fee status 
feedat <- read.csv(paste0(data_fp,"rls-dontshare/MPAs-for-tourism-review_CAB-CB.csv"))
#B20 metric
dat <- read.csv(paste0(data_fp,"rls-dontshare/Metrics_global.csv"))
#Site level covariates (inc NEOLI)
sdat <- read_csv(paste0(data_fp,"rls-dontshare/Copy of Sites-master_consolidated 4June2020GS_GE.csv"))
humdat <- read.csv(paste0(data_fp,"rls-dontshare/RLSPoPOutput.txt"), header = TRUE)

# fee data
feedat2 <- select(feedat, LocationName, MPA, Tourism.fee.amount) %>% distinct() %>%
  mutate(fee = (Tourism.fee.amount > 0))

# clean tourist intensity index data
tdat2 <-  select(tdat, LocationName, MPA, Tourism.intensity.index.consensus, Zone.name) %>%
  distinct()

# Join site data to fee data and tourism intensity - 
# Joined by Zone.name too which removed duplicates (a large number) 
sdat2 <- select(sdat, SiteCode, SiteCode_rev, LocationName, 
                MPA, `NEOLI Total`:Age, Zone.name = `Zone name`) %>% 
  distinct() %>% #removes one duplicate
  left_join(tdat2, by = c("LocationName", "MPA", "Zone.name")) %>% 
  left_join(feedat2, by = c("LocationName", "MPA")) %>%
  mutate(fee = ifelse(is.na(fee), 0, fee))
nrow(sdat)
nrow(sdat2)

colSums(is.na(sdat2))

#
# MPA NEOLI summary
#

MPAsum <- sdat2 %>% filter(!is.na(`NEOLI Total`)) %>%
  select(MPA, `NEOLI Total`:Age) %>%
  distinct()
nMPA <- length(unique(MPAsum$MPA))
table(MPAsum$`NEOLI Total`)/nMPA
table(MPAsum$`No-take`)/nMPA
table(MPAsum$Effectiveness_NEOLI)/nMPA
table(MPAsum$area)/nMPA
table(MPAsum$Isolation_NEOLI)/nMPA
table(MPAsum$Age)/nMPA

table(with(MPAsum, NEOLI_Total - as.numeric(as.character(notake)) - 
             as.numeric(as.character(Effectiveness_NEOLI)) -
             as.numeric(as.character(area)) - 
             as.numeric(as.character(Isolation_NEOLI))))
#replace missing sitecodes
# dat has sitecode, whereas
# sdat has sitecode_rev
# in dat, sitecodes are same as sitecode_rev in sdat

missing_site_codes <- which(!(sdat2$SiteCode %in% dat$SiteCode))
sdat3 <- sdat2
sdat3$SiteCode[missing_site_codes] <- sdat2$SiteCode_rev[missing_site_codes]


# Join tourism/site/fee data with biomass data
fdat <- left_join(dat, sdat3, by = "SiteCode")

# Combine bio variable data and human data
biohumdat <- humdat %>% select(ID, Country, SiteLatitu, SiteLongit, 
                               KernelD202, KernelD201, Kernel2020, Kernel2015) %>% 
  left_join(select(biodat, - Name),
            by = c("SiteLatitu" = "SiteLatitude", "SiteLongit" = "SiteLongitude", 
                   "Country", "ID")) %>%  
  select(-ID) %>% distinct() %>% rename(Lat = SiteLatitu, Lon = SiteLongit) # Some double ups in long and lat based on ID

# Join bio/human data with tourism/site/fee data and edit vars
alldata <- left_join(fdat, select(biohumdat, -Country), by = c("Lon", "Lat")) %>% 
  mutate(TI = Tourism.intensity.index.consensus,
         notake = ifelse(is.na(`No-take`), 0, `No-take`)) %>% 
  select(SurveyID, SiteCode, Lat, Lon, SurveyDate, Country, Location, Year, 
         B20, LocationName, MPA, NEOLI_Total = `NEOLI Total`, Effectiveness_NEOLI, 
         Isolation_NEOLI, area,  Zone_name = Zone.name, TI, notake,
         fee_amount = Tourism.fee.amount, has_fee = fee, KernelD202, KernelD201, 
         Kernel2020, Kernel2015, salinty, chlomean, dissox, SSTMin, Phosphate, 
         Silicate, ParMean, SSTrange, SSTMax, SSTMean, Nitrate, Calcite)

nrow(alldata)

# Get non MPA data with all biovars
noMPAdat <- alldata %>% filter(is.na(MPA), !is.na(salinty))

# Get MPA data with NEOLI scores and all biovars
MPAdat <- alldata %>% filter(!is.na(MPA)) %>% # Get MPA data
  filter(!is.na(NEOLI_Total), !is.na(salinty)) # Remove those samples which we are missing data for
colSums(is.na(MPAdat))

## Function to filter the closest years to 2015
rm_other_years <- function(data){
  
  # Calculate the differences
  data$diff2015 <- abs(data$Year - 2015)
  
  # Find which differences we want to keep for each site
  which_difs <- data %>% group_by(Lon, Lat) %>% 
    summarise(SiteCode = SiteCode[which.min(diff2015)], # To remove mislabeled SiteCode (duplicates) 
              min_dif = min(diff2015)) %>% 
    mutate(keep = TRUE) %>% ungroup() %>% select(-Lon, -Lat)
  
  # Join to data and filter by the differences we want to keep
  output <- full_join(data, which_difs, by = c("SiteCode", "diff2015" = "min_dif")) %>% 
    filter(!is.na(keep)) %>% select(-keep, -diff2015)
  
  return(output)
}

# Filter to only the closest years to 2015
noMPAdat <- rm_other_years(noMPAdat)
MPAdat <- rm_other_years(MPAdat)

# Average samples from the same site
# This groups by the columns we want to  keep so we can average the years when there are two
noMPAdat <- noMPAdat %>% 
  group_by_at(setdiff(names(noMPAdat), c("SurveyID", "Year", "SurveyDate", "B20"))) %>%
  summarise(B20 = mean(B20), Years = paste(unique(Year), collapse = " "))
MPAdat <- MPAdat %>% 
  group_by_at(setdiff(names(MPAdat), c("SurveyID", "Year", "SurveyDate", "B20"))) %>%
  summarise(B20 = mean(B20), Years = paste(unique(Year), collapse = " "))



# Formatting data for the modelling ---------------------------------------


# Transform the response
plusmin <- min(rbind(noMPAdat, MPAdat)$B20[rbind(noMPAdat, MPAdat)$B20>0])/2
noMPAdat <- mutate(.data = noMPAdat, log10_B20 = log10(B20 + plusmin))
MPAdat <- mutate(.data = MPAdat, log10_B20 = log10(B20 + plusmin))

# Formatting Non MPA data for the modelling
noMPAdat <- within(noMPAdat,{
  SiteCode <- factor(SiteCode)
  Years <- factor(Years) 
  LocationName <- factor(LocationName)
  Location <- factor(Location)
  Country <- factor(Country)
  Kernel2015_4throot<- Kernel2015^(1/4)
  log10_B20_more <-  log10(B20 + plusmin)
  
}) %>% filter(Kernel2015 != -9999)


# Formatting MPA data for the modelling
MPAdat <- within(MPAdat,{
  
  SiteCode <- factor(SiteCode)
  Years <- factor(Years) 
  LocationName <- factor(LocationName)
  Location <- factor(Location)
  Country <- factor(Country)
  Kernel2015_4throot<- Kernel2015^(1/4)
  log10_B20_more <-  log10(B20 + plusmin)
  fee_amount[fee_amount == -1] <- NA # replace missing values
  notake<- factor(notake)
  has_fee<- factor(has_fee)
  MPA <- factor(MPA)
  
}) %>% filter(Kernel2015 != -9999)


# Save for the next step
save("MPAdat", "noMPAdat", file = paste0(data_folder,  "/rls-tour_step_one_data-R1.RDA" ))





