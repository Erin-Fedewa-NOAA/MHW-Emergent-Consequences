#calculate upper percentile of mature female size distribution

#Author: EJF

## Load packages
library(crabpack)
library(tidyverse)
library(reldist)

########################################

# Set years ----
current_year <- 2025
years <- 1988:current_year

## Pull specimen data from crabpack ----
snow <- get_specimen_data(species = "SNOW",
                          region = "EBS",
                          channel = "KOD")

##########################################

#First- 90th percentile of mature female population, unweighted
perc_90 <- snow$specimen %>%
  filter((SEX == 2 & CLUTCH_SIZE > 0 & SHELL_CONDITION <= 2),
         YEAR %in% years) %>%
  group_by(YEAR) %>%
  summarise(upper_90th = quantile(SIZE, 0.90, na.rm = TRUE)) 

#now weighted quantile
perc_90_weight <- snow$specimen %>%
  filter((SEX == 2 & CLUTCH_SIZE > 0 & SHELL_CONDITION <= 2),
          YEAR %in% years) %>%
  group_by(YEAR, SIZE, STATION_ID, AREA_SWEPT) %>%
  mutate(cpue = sum(SAMPLING_FACTOR) / AREA_SWEPT) %>%
  group_by(YEAR) %>% 
  summarize(upper_90th_weighted = wtd.quantile(SIZE, weight = cpue)) 

#weighted mean size of mature females (i.e. avg SAM)
mean_size <- snow$specimen %>%
  filter((SEX == 2 & CLUTCH_SIZE > 0 & SHELL_CONDITION <= 2),
         YEAR %in% years) %>%
  group_by(YEAR, SIZE, STATION_ID, AREA_SWEPT) %>%
  mutate(cpue = sum(SAMPLING_FACTOR) / AREA_SWEPT) %>%
  ungroup() %>%
  group_by(YEAR) %>% 
  summarize(female_size = weighted.mean(SIZE, w = cpue, na.rm = TRUE))

#relative proportion of mature females larger than timeseries mean 

#average size of mature females 
snow$specimen %>%
  filter((SEX == 2 & CLUTCH_SIZE > 0 & SHELL_CONDITION <= 2),
         YEAR %in% years) %>%
  group_by(YEAR, SIZE, STATION_ID, AREA_SWEPT) %>%
  mutate(cpue = sum(SAMPLING_FACTOR) / AREA_SWEPT) %>%
  ungroup() %>%
  summarize(female_size = weighted.mean(SIZE, weight = cpue))

avg_size <- 55.7

prop_large <- calc_bioabund(crab_data=snow,
              species="SNOW",
              year = years,
              crab_category = "mature_female",
              bin_1mm = TRUE) %>%
     group_by(YEAR) %>%
     summarise(total_fem_abun = sum(ABUNDANCE),
                large_fem_abun = sum(ABUNDANCE[SIZE_1MM >= avg_size], na.rm = TRUE),
                proportion_large = (large_fem_abun/total_fem_abun)) %>%
  select(YEAR, proportion_large)

#combine datasets
perc_90 %>%
  full_join(perc_90_weight) %>%
  full_join(mean_size) %>%
  full_join(prop_large) %>%
  pivot_longer(cols= 2:5, names_to = "method", values_to = "estimate") -> all_dat

#plot different methods 
all_dat %>%
  ggplot(aes(YEAR, estimate)) +
  geom_point() +
  geom_line() +
  facet_wrap(~method, scales="free_y")

#lets go with proportional size distribution for now, it seems to best reflect 
  #size structure of the mature female population
write_csv(prop_large, "./output/female_size.csv")

#Follow up: might be worth considering a metric for immatures b/c we'd be able to 
  #pick up on the 2024-2025 signal earlier - but that gets tricky b/c you'd need 
  #to control for cohort effects, vrs. primiparous mature females as a presumed
  #single cohort 

#Also to follow up- why does this not produce the same data as the avg size method above??
  #Rounding down to 1mm bin should produce consistent underestimation
snow_abund <- calc_bioabund(crab_data = snow,
                            species = "SNOW",
                            region = "EBS",
                            years = c(1988:2025),
                            shell_condition = c("soft_molting","new_hardshell"),
                            crab_category = "mature_female",
                            bin_1mm = TRUE) %>%
  group_by(YEAR) %>%
  summarise(fem_abund = weighted.mean(SIZE_1MM, w = ABUNDANCE))
