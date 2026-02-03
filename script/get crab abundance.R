#Produce crab abundance survey output for 
#1) hybrids 
#2) snow and tanner population
#2) immature/mature female snow crab
#4) large male snow crab (95mm +)
#note: Emily's sdmtmb maturity workflow is not yet built into crabpack, so we're 
  #not generating male abundance estimates by maturity here 

#Author: EJF

## Load packages
library(crabpack)
library(tidyverse)

########################################

# Set years ----
current_year <- 2025
years <- 1988:current_year

## Pull specimen data from crabpack ----
hybrid <- get_specimen_data(species = "HYBRID",
                          region = "EBS",
                          channel = "KOD")

snow <- get_specimen_data(species = "SNOW",
                            region = "EBS",
                            channel = "KOD")

tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")

#################################

#calculate hybrid abundance timeseries 
hyb_abun_pop <- calc_bioabund(crab_data = hybrid,
                              species = "HYBRID",
                              region = "EBS",
                              years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population",
         species = "hybrid")

#Plot
hyb_abun_pop %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#now let's subset hybrids for sizes that we know we're able to detect 
  #first females > 65mm
hyb_sub_female <- calc_bioabund(crab_data = hybrid,
                              species = "HYBRID",
                              region = "EBS",
                              sex = "female", 
                              size_min = 65,
                              years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  rename(abundance_fem = abundance, abundance_female_ci=abundance_ci)

#and then males > 79mm
calc_bioabund(crab_data = hybrid,
                          species = "HYBRID",
                          region = "EBS",
                          sex = "male", 
                          size_min = 79,
                          years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  rename(abundance_male = abundance, abundance_male_ci=abundance_ci) %>%
#join to females and sum for hybrid abundance
full_join(hyb_sub_female) %>%
  group_by(year) %>%
  mutate(abundance = sum(abundance_male, abundance_fem)) %>%
  select(year, abundance) %>%
  mutate(category = "population_subset",
         species = "hybrid") %>%
#join to full hybrid population time series 
full_join(hyb_abun_pop) -> hyb_abun


#Plot both timeseries
hyb_abun %>%
  ggplot(aes(x = year, y = abundance, color=category)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate snow crab abundance timeseries 
snow_abun_pop <- calc_bioabund(crab_data = snow,
                          species = "SNOW",
                          region = "EBS",
                          years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population",
         species = "snow")

#Plot
snow_abun_pop %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate tanner crab abundance timeseries 
tanner_abun_pop <- calc_bioabund(crab_data = tanner,
                           species = "TANNER",
                           region = "EBS",
                           spatial_level = "region",
                           years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population",
         species = "tanner")

#Plot
tanner_abun_pop %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate snow crab mature/immature abundance timeseries 
snow_abun_female <- calc_bioabund(crab_data = snow,
                           species = "SNOW",
                           region = "EBS",
                           sex = "female",
                           crab_category = c("mature_female", "immature_female"),
                           years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI, CATEGORY) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(species = "snow")

#Plot
snow_abun_female %>%
  filter(category != "NA") %>%
  ggplot(aes(x = year, y = abundance, color=category)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate snow crab large male abundance timeseries 
snow_abun_lgmale <- calc_bioabund(crab_data = snow,
                                  species = "SNOW",
                                  region = "EBS",
                                  sex = "male",
                                  size_min = 95,
                                  years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(species = "snow",
         category = "large_male")

#Plot
snow_abun_lgmale  %>%
  filter(category != "NA") %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#join all abundance datasets
hyb_abun %>%
  full_join(snow_abun_pop) %>%
  full_join(tanner_abun_pop) %>%
  full_join(snow_abun_female) %>%
  full_join(snow_abun_lgmale) %>%
  select(-abundance_ci) -> abun_dat

#Write output 
write_csv(abun_dat, "./output/crab_abundance.csv")
