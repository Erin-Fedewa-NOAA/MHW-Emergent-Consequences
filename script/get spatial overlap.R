#Calculate spatial overlap of snow lg male/Tanner mature female and Tanner 
  #lg male/snow mature female as a proxy for interspecies mating opportunities

#Original Author: E. Ryznar
  #modifications by E. Fedewa

## Load packages
library(crabpack)
library(tidyverse)
library(patchwork)

#------------------------------
# #data wrangling
#------------------------------

## Pull specimen data from crabpack ----
snow <- get_specimen_data(species = "SNOW",
                          region = "EBS",
                          channel = "KOD")

tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")

# Set years ----
current_year <- 2025
years <- 1988:current_year

#define corner stations so we can drop
stations <- read.csv("./data/station_lookup.csv")

corners <- stations %>% 
  filter(REGION == "EBS", STATION_TYPE == "MTCA_CORNER") %>%
  pull(STATION_ID)

#----------------------------------------------------
# #spatial overlap functions from Carroll et al 2019
#------------------------------------------------------

## Bhattacharyya's coefficient for two species
## input is station-level density estimates for two species
## measures whether two species use space independently
bhatta_coeffn <- function(prey, pred) {
  #calculate station-level cpue/total cpue within a given year for each species
  p_prey <- prey/sum(prey, na.rm = T) 
  p_pred <- pred/sum(pred, na.rm = T)
  #multiply station-level species probabilities by year
  sum(sqrt(p_prey*p_pred), na.rm = T)
}

## Bhattacharyya's coefficient for four species
## input is station-level density estimates for two species
## measures whether two species use space independently
bhatta4_coeffn <- function(sp1, sp2, sp3, sp4) {
  #calculate station-level cpue/total cpue within a given year for each species
  p_sp1 <- sp1/sum(sp1, na.rm = T) 
  p_sp2 <- sp1/sum(sp2, na.rm = T)
  p_sp3 <- sp1/sum(sp3, na.rm = T)
  p_sp4 <- sp1/sum(sp4, na.rm = T)
  #multiply station-level species probabilities by year
  sum(sqrt(p_sp1*p_sp2*p_sp3*p_sp4), na.rm = T)
}

#----------------------------------------------------
# #Pull CPUE data from crabpack
#------------------------------------------------------

#Station-level mature female snow crab CPUE
snow_fem_cpue <- crabpack::calc_cpue(crab_data = snow,
                                     region = "EBS",
                                     species = "SNOW",
                                     years = years,
                                     crab_category = "mature_female")

#Station-level mature female tanner crab CPUE
tanner_fem_cpue <- crabpack::calc_cpue(crab_data = tanner,
                                       region = "EBS",
                                       species = "TANNER",
                                       years = years,
                                       crab_category = "mature_female")

#Station-level large male snow crab CPUE
snow_male_cpue <- crabpack::calc_cpue(crab_data = snow,
                                      region = "EBS",
                                      species = "SNOW",
                                      years = years,
                                      crab_category = "large_male")

#Station-level large male snow crab CPUE
tanner_male_cpue <- crabpack::calc_cpue(crab_data = tanner,
                                        region = "EBS",
                                        species = "TANNER",
                                        years = years,
                                        crab_category = "large_male") 

#----------------------------------------------------
# #Bhattacharyya's Overlap Calculations
#------------------------------------------------------

#Calculate bhatt overlap between snow mature females and tanner large males
bhatt_snowfem_tanmale <- snow_fem_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, cpue) %>%
  rename(snow_female_cpue = cpue) %>%
  full_join(tanner_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(tanner_male_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) %>% #exclude corner stations
  group_by(year) %>%
  summarise(bhatta_snowfem_tanmale = bhatta_coeffn(snow_female_cpue, tanner_male_cpue)) 

#plot  
bhatt_snowfem_tanmale %>%
  ggplot(aes(year, bhatta_snowfem_tanmale)) +
  geom_point() +
  geom_line() -> snowfem_tanmale_plot

#Calculate bhatt overlap between Tanner mature females and snow large males
bhatt_tannerfem_snowmale <- tanner_fem_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, cpue) %>%
  rename(tanner_female_cpue = cpue) %>%
  full_join(snow_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(snow_male_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) %>% #exclude corner stations
  group_by(year) %>%
  summarise(bhatta_tannerfem_snowmale = bhatta_coeffn(tanner_female_cpue, snow_male_cpue)) 

#plot  
bhatt_tannerfem_snowmale %>%
  ggplot(aes(year, bhatta_tannerfem_snowmale)) +
  geom_point() +
  geom_line()

#Calculate global bhatt overlap: Tanner mature females, snow large males
  #snow mature females and tanner large males
bhatt_global <- tanner_fem_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, cpue) %>%
  rename(tanner_female_cpue = cpue) %>%
  full_join(snow_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(snow_male_cpue = cpue)) %>%
  full_join(snow_fem_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(snow_female_cpue = cpue)) %>%
  full_join(tanner_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(tanner_male_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) %>% #exclude corner stations
  group_by(year) %>%
  summarise(bhatta_global = bhatta4_coeffn(tanner_female_cpue, snow_male_cpue,
                                          snow_female_cpue, tanner_male_cpue)) 

#plot  
bhatt_global %>%
  ggplot(aes(year, bhatta_global)) +
  geom_point() +
  geom_line()

#----------------------------------------
# #Join datasets and produce output
#----------------------------------------

#join datasets
bhatt_global %>%
  full_join(bhatt_tannerfem_snowmale) %>%
  full_join(bhatt_snowfem_tanmale) -> bhatt_dat

#plot all bhatt overlap indices
bhatt_dat %>%
  pivot_longer(2:4, names_to = "method", values_to = "overlap") %>%
  ggplot(aes(year, overlap, color=method)) +
  geom_point() +
  geom_line() +
  theme_bw()

#Write output 
missing <- data.frame(year = 2020)

bhatt_dat %>%
  bind_rows(missing) %>%
  arrange(year) %>%
  write.csv(file="./output/overlap_output.csv", row.names = F)

#----------------------------------------------------
# #Quick look at a few other overlap metrics:
#------------------------------------------------------

#Functions from Carroll et al:

## Schoener's D
## density or probability of occurrence data
## measures how equally predator and prey share available resources
schoeners_overlapfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  1 - 0.5 * (sum(abs(p_prey-p_pred), na.rm = T))
}

## local index of collocation
## estimates correlation of predator and prey densities
loc_collocfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(p_prey*p_pred, na.rm = T)/(sqrt(sum(p_prey^2, na.rm = T)*sum(p_pred^2, na.rm = T)))
}
#################################

#Calculate Schoener's D overlap between snow mature females and tanner large males
schoeD_snowfem_tanmale <- snow_fem_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, cpue) %>%
  rename(snow_female_cpue = cpue) %>%
  full_join(tanner_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(tanner_male_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) %>% #exclude corner stations
  group_by(year) %>%
  summarise(schoeD = schoeners_overlapfn(snow_female_cpue, tanner_male_cpue))   

#compare to bhatt
schoeD_snowfem_tanmale %>%
  full_join(bhatt_snowfem_tanmale) %>%
  pivot_longer(2:3, names_to = "overlap_metric", values_to = "overlap") %>%
  ggplot(aes(year, overlap, color=overlap_metric)) +
  geom_point() +
  geom_line()
#trends very similar! 

##Calculate local index of collocation for snow mature females and tanner large males
colloc_snowfem_tanmale <- snow_fem_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, cpue) %>%
  rename(snow_female_cpue = cpue) %>%
  full_join(tanner_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue) %>%
              rename(tanner_male_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) %>% #exclude corner stations
  group_by(year) %>%
  summarise(colloc = loc_collocfn(snow_female_cpue, tanner_male_cpue)) 

#compare to bhatt
colloc_snowfem_tanmale %>%
  full_join(bhatt_snowfem_tanmale) %>%
  pivot_longer(2:3, names_to = "overlap_metric", values_to = "overlap") %>%
  ggplot(aes(year, overlap, color=overlap_metric)) +
  geom_point() +
  geom_line()

#Simple presence/absence of mature female snow crab/lg male tanner crab 
  #this should be similar to bhatt, but doesn't consider density
snow_fem_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, cpue, latitude, longitude) %>%
  rename(snow_female_cpue = cpue) %>%
  full_join(tanner_male_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, cpue, latitude, longitude) %>%
              rename(tanner_male_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) %>%
  group_by(year, station_id, latitude, longitude) %>%
  filter(snow_female_cpue > 0 & tanner_male_cpue > 0) -> pres_snowfem_tannermale

#plot spatially 
pres_snowfem_tannermale %>%
  ggplot(aes(latitude, longitude)) +
  geom_point() +
  facet_wrap(~year)
  
#And plot number of stations alongside bhatt metric 
pres_snowfem_tannermale %>%
  group_by(year) %>%
  summarize(count = n_distinct(station_id)) %>%
  ggplot(aes(year, count)) +
  geom_point() +
  geom_line() -> presab_snowfem_tanmale

snowfem_tanmale_plot / presab_snowfem_tanmale








