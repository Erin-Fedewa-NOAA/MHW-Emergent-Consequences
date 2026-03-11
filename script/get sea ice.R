# Calculate EBS average sea ice concentration, masked for survey grid
  #Source: ERA5 ice cover data, monthly averaged reanalysis
      #https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview

#Original script author: E. Ryznar
  #modifictions: E Fedewa


# load ----
library(tidyverse)
library(tidync)
library(lubridate)
library(magrittr)
library(akgfmaps)
library(ecmwfr)
library(purrr)
library(akgfmaps)

#-----------------------------------------------
# Pull ice data using Climate Data Store API
#-----------------------------------------------

# specify years (run requests sequentially, otherwise file size is too big)
ice.years <- 1980:1988
ice.years <-1989:2000
ice.years <- 2001:2013
ice.years <- 2014:2025

# specify login credentials for the climate data store
user_id = "erin.fedewa@noaa.gov"
  api_key = "3e25eba0-9dda-48d8-9d53-df47aeae8620"
  
#set key
  wf_set_key(user = user_id,
             key = api_key) 

#specify request for current year
request <- list(
  "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
  "product_type" = "monthly_averaged_reanalysis",
  "variable" = c("sea_ice_cover"),    
  "year" = ice.years,                     
  "month" = sprintf("%02d", 1:12),
  "day" = sprintf("%02d", 1:31),
  "time" = sprintf("%02d:00", 0:23),
  "area" = c(65, -178, 50, -154),      "format" = "netcdf",                  
  "target" = paste0("ERA5_ice_", min(ice.years), "-", max(ice.years), ".nc") # target file name
)

# run request (you may need to manually click accept license on website --> follow link in error message if it appears)
wf_request(
  user     = user_id,
  request  = request,
  transfer = TRUE,
  path     = paste0("./data/"), # where do you want the data to be saved?
  verbose = TRUE
)

#-----------------------------
# Process Ice Data
#-----------------------------

#pull all ice files
ice_files <- list.files("./data", pattern = "ERA5_ice")

ice_all <- map_dfr(ice_files, \(file) {
  
  tidync(file.path("./data", file)) %>%
    hyper_tibble() %>%
    separate(valid_time, into = c("year", "month", "time"), sep = "-") %>%
    select(-time) %>%
    mutate(across(c(year, month, latitude, longitude), as.numeric)) %>%
    filter(month %in% 1:4)
})

#calculate monthly averages for entire area specified in data pull
ice_means <- ice_all %>%
  group_by(year, month) %>%
  summarise(ice_conc = mean(siconc), .groups = "drop")

#plot
ice_means %>%
  ggplot(aes(year, ice_conc)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = .15, color= "red") + #ice-free threshold
  facet_wrap(~month)

#calculate spatial (grid-cell) averages for Jan-April period
ice_spatial <- ice_all %>%
  group_by(year, latitude, longitude) %>%
  summarise(ice_conc = mean(siconc), .groups = "drop")

#-------------------------------------------------------------
# Subset spatial extent of ice data to EBS survey grid area
#-------------------------------------------------------------

#Get EBS grid for masking
region_layers <- akgfmaps::get_base_layers("sebs")
survey_area <- region_layers$survey.area 

# Convert ice data to spatial points and keep only those within survey area
ice <- ice_spatial %>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = "epsg:4326") %>%
  st_transform(., st_crs(survey_area)) %>%
  st_intersection(., survey_area) 

#Calculate yearly mean ice concentration within the survey area
ebs_ice <- ice %>%
  st_drop_geometry() %>% 
  drop_na(ice_conc) %>%
  group_by(year) %>%
  summarise(ice_conc = mean(ice_conc), .groups = "drop")

#plot
ebs_ice %>%
  ggplot(aes(year, ice_conc)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = .15, color= "red")  #ice-free threshold
  
#write output
write.csv(ebs_ice, "./output/seaice_output.csv", row.names=F)
