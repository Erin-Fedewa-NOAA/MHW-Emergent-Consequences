# Calculate Bering Sea (EBS and NBS) average sea ice concentration

#Author: Erin Fedewa

#Notes:
#In 2025, NSIDC sea ice extent index was replaced with average sea ice concentration from ERA5 reanalysis
#To do in 2026: Automate ERA5 data pull in script: https://cds.climate.copernicus.eu/how-to-api 

# load ----
library(tidyverse)
library(tidync)
library(lubridate)
library(magrittr)
library(akgfmaps)

###########################################################
#To download and process ice cover data from ERA 5 monthly averaged data
# 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview
# 2) Click on "Download data" tab
# 3) Click on the product type you’d like. We have been using “Monthly averaged reanalysis”
# 4) For ice cover, click on the “Other” drop down, and then check the “Sea-ice cover” box. 
# 5) Select years you’d like the data to cover. Initial pull has been divided into two time stanzas, as
#a larger request can result in an error with processing 
# 6) Select months you’d like the data to cover (Jan-Apr in this case) 
# 7) Select geographical area for the data. For the data below, the region has been set to:
#North 65°, West -178°, South 56°, and East -165°. 
# 8) Select NetCDF(experimental) as the data format and unarchived as download format, and submit form to query/download data. 

################################################################
### Process Ice Data
early_dat <- tidync("./Data/ERA5_ice_1975_1999.nc") %>% 
  hyper_tibble() %>% 
  separate(valid_time, into = c("Year", "Month", "Time"), sep = "-") %>%
  select(-Time) %>%
  mutate(Year = as.numeric(Year),
         Month = as.numeric(Month))

recent_dat <- tidync("./Data/ERA5_ice_2000_2025.nc") %>% 
  hyper_tibble() %>% 
  separate(valid_time, into = c("Year", "Month", "Time"), sep = "-") %>%
  select(-Time) %>%
  mutate(Year = as.numeric(Year),
         Month = as.numeric(Month))

rbind(early_dat, recent_dat) %>%
  group_by(Year) %>%
  summarize(ice_avg = mean(siconc,  na.rm=T)*100) -> ice_extent

#Plot timeseries
ice_extent %>%
  ggplot(aes(Year, ice_avg)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(ice_avg, na.rm=TRUE)), linetype = 5) +
  geom_hline(yintercept = 15, color= "red") + #ice-free threshold
  theme_bw()
#Less than 15% ice cover in Bering Sea in spring 2018-2019- below the threshold
#for ice covered area

#write csv 
write.csv(ice_extent, "./output/seaice_output.csv", row.names=F)
