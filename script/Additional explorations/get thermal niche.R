#Calculate hybrid temperatures of occupancy (CPUE weighted) 

# Author: Erin Fedewa

# load ----
library(tidyverse)

# Set years ----
current_year <- 2025
years <- 1988:current_year

## Pull specimen data from crabpack ----
hybrid <- get_specimen_data(species = "HYBRID",
                            region = "EBS",
                            channel = "KOD")

haul <- hybrid$haul

#Calculate station-level CPUE ----
hybrid_cpue <- calc_cpue(crab_data = hybrid,
                         species = "HYBRID",
                         years = years)

##########################################
#Temperature of Occupancy Calculations ----
#Note that this is ignoring NA's in bottom temp data!!!
hybrid_cpue %>%
  left_join(., haul %>% select(YEAR, GEAR_TEMPERATURE, STATION_ID)) %>%
  group_by(YEAR) %>% 
  summarise(temp_occ = weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T)) %>%
  print(n=50) -> temp_occ

#plot
temp_occ %>%
  ggplot(aes(x = YEAR, y = temp_occ))+
  geom_point(size=3)+
  geom_line() +
  theme_bw()



