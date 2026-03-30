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

# specify login credentials for the climate data store
user_id = "erin.fedewa@noaa.gov"
  api_key = "3e25eba0-9dda-48d8-9d53-df47aeae8620"
  
#set key
  wf_set_key(user = user_id,
             key = api_key) 
  
#specify years (run below requests sequentially, otherwise file size is too big)
  ice.years <- 1980:1988
  ice.years <-1989:2000
  ice.years <- 2001:2013
  ice.years <- 2014:2025

#specify request for current year
request <- list(
  "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
  "product_type" = "monthly_averaged_reanalysis",
  "variable" = c("sea_ice_cover"),    
  "year" = ice.years,                     
  "month" = sprintf("%02d", 1:12),
  "day" = sprintf("%02d", 1:31),
  "time" = sprintf("%02d:00", 0:23),
  "area" = c(66, -179, 54, -158),      "format" = "netcdf",                  
  "target" = paste0("ERA5_ice_", min(ice.years), "-", max(ice.years), ".nc") # target file name
  )

# run request (you may need to manually click accept license on website --> follow link in error message if it appears)
wf_request(
  user     = user_id,
  request  = request,
  transfer = TRUE,
  path     = paste0("./data/"), # where do you want the data to be saved?
  verbose = TRUE)

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

#-----------------------------
# Plot Ice Data
#-----------------------------

#Load survey layers
region_layers <- akgfmaps::get_base_layers("ebs")
survey_area <- region_layers$survey.area
akland <- region_layers$akland

#compute annual grid averages across Jan-Apr period
siconc_yearly <- ice_all %>%
  group_by(longitude, latitude, year) %>%
  summarise(siconc = mean(siconc, na.rm = TRUE),
            .groups = "drop")

#plot faceted by year 
ggplot() +
  # sea ice raster
  geom_raster(
    data = siconc_yearly,
    aes(longitude, latitude, fill = siconc)) +
  # land
  geom_sf(data = akland,
          fill = "grey85",
          color = NA) +
  # survey boundary
  geom_sf(data = survey_area,
          fill = NA,
          color = "white",
          linewidth = 0.4) +
  scale_fill_viridis_c(
    name = "Sea Ice\nConcentration",
    limits = c(0,1),
    option = "C") +
  coord_sf(
    xlim = region_layers$plot.boundary$x,
    ylim = region_layers$plot.boundary$y,
    expand = FALSE) +
  facet_wrap(~year, ncol = 8) +
  theme_bw(base_size = 11) +
  labs(title = "Annual Mean Sea Ice Concentration in the Bering Sea",
    x = "Longitude",
    y = "Latitude") +
  theme(plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm"))

#API data request uses EBS/NBS survey grid max/min lat and longs, but because 
  #this box ends up being quite a bit bigger than the survey grid, we'll 
  #spatially subset the sea ice data to the EBS/NBS grid 

#-------------------------------------------------------------
# Subset spatial extent of ice data to survey grid area
#-------------------------------------------------------------

#Get survey grids for masking
ebs_survey_area <- akgfmaps::get_base_layers("sebs")$survey.area
bs_survey_area  <- akgfmaps::get_base_layers("ebs")$survey.area

#function for calculating seasonal spatial averages 
  #(b/c we'll evaluate over several time periods)
calc_spatial_mean <- function(data, months = NULL) {
  data %>%
    { if (!is.null(months)) filter(., month %in% months) else . } %>%
    group_by(year, latitude, longitude) %>%
    summarise(ice_conc = mean(siconc), .groups = "drop")
}

#compute time periods of interest 
ice_spatial <- list(
  Jan_Feb = calc_spatial_mean(ice_all, 1:2),
  Mar_Apr = calc_spatial_mean(ice_all, 3:4),
  Jan_Apr = calc_spatial_mean(ice_all))

#function for spatial masking + calculating yearly mean ice 
mask_and_average <- function(data, survey_area, name) {
    data %>%
    st_as_sf(coords = c("longitude","latitude"), crs = "epsg:4326") %>%
    st_transform(st_crs(survey_area)) %>%
    st_intersection(survey_area) %>%
    st_drop_geometry() %>%
    drop_na(ice_conc) %>%
    group_by(year) %>%
    summarise(value = mean(ice_conc), .groups="drop") %>%
    rename(!!name := value)
}

#EBS-only yearly averages
Jan_Feb_ice_EBS <- mask_and_average(ice_spatial$Jan_Feb, ebs_survey_area, "Jan_Feb_ice_EBS")
Mar_Apr_ice_EBS <- mask_and_average(ice_spatial$Mar_Apr, ebs_survey_area, "Mar_Apr_ice_EBS")
Jan_Apr_ice_EBS <- mask_and_average(ice_spatial$Jan_Apr, ebs_survey_area, "Jan_Apr_ice_EBS")

#EBS+NBS yearly averages
Jan_Feb_ice_EBS_NBS <- mask_and_average(ice_spatial$Jan_Feb, bs_survey_area, "Jan_Feb_ice_EBS_NBS")
Mar_Apr_ice_EBS_NBS <- mask_and_average(ice_spatial$Mar_Apr, bs_survey_area, "Mar_Apr_ice_EBS_NBS")
Jan_Apr_ice_EBS_NBS <- mask_and_average(ice_spatial$Jan_Apr, bs_survey_area, "Jan_Apr_ice_EBS_NBS")

#-------------------------------------------
# Combine all time series and output
#-------------------------------------------
ice_timeseries <- Mar_Apr_ice_EBS %>%
  full_join(Jan_Feb_ice_EBS) %>%
  full_join(Jan_Apr_ice_EBS) %>%
  full_join(Mar_Apr_ice_EBS_NBS) %>%
  full_join(Jan_Feb_ice_EBS_NBS) %>%
  full_join(Jan_Apr_ice_EBS_NBS)

#plot
#jan-april ice
ice_timeseries %>%
  pivot_longer(-year, names_to = "metric", values_to = "ice_conc") %>%
  filter(metric %in% c("Jan_Apr_ice_EBS", "Jan_Apr_ice_EBS_NBS")) %>%
  ggplot(aes(year, ice_conc, color=metric)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = .15, color= "red") -> a   #ice-free threshold

#jan-feb ice
ice_timeseries %>%
  pivot_longer(-year, names_to = "metric", values_to = "ice_conc") %>%
  filter(metric %in% c("Jan_Feb_ice_EBS", "Jan_Feb_ice_EBS_NBS")) %>%
  ggplot(aes(year, ice_conc, color=metric)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = .15, color= "red") -> b   #ice-free threshold

#march-april ice
ice_timeseries %>%
  pivot_longer(-year, names_to = "metric", values_to = "ice_conc") %>%
  filter(metric %in% c("Mar_Apr_ice_EBS", "Mar_Apr_ice_EBS_NBS")) %>%
  ggplot(aes(year, ice_conc, color=metric)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = .15, color= "red") -> c   #ice-free threshold

#combined plot 
a/b/c
  
#write output
write.csv(ice_timeseries, "./output/seaice_output.csv", row.names=F)












