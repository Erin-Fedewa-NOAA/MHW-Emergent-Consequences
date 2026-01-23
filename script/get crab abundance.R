#Produce crab abundance survey output for 
#1) hybrids (not yet split out by sex/maturity categories)
#2) immature/mature female snow crab
#3) immature/mature female tanner crab
#4) snow and tanner population
#note: Emily's sdmtmb maturity workflow is not yet built into crabpack, so we're 
  #not generating male abundance estimates here 

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
hyb_abun <- calc_bioabund(crab_data = hybrid,
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
hyb_abun %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate snow crab abundance timeseries 
snow_abun <- calc_bioabund(crab_data = snow,
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
snow_abun %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate tanner crab abundance timeseries 
tanner_abun <- calc_bioabund(crab_data = tanner,
                           species = "TANNER",
                           region = "EBS",
                           district = "ALL",
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
tanner_abun %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#stop here and ask shannon about ebs wide Tanner estimates and then do
#female estimates


#calculate snow crab mature/immature abundance timeseries 
snow_fem_abun <- calc_bioabund(crab_data = snow,
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
snow_fem_abun %>%
  filter(category != "NA") %>%
  ggplot(aes(x = year, y = abundance, color=category)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()


#join all abundance datasets




#Write output 
write_csv(hyb_abun, "./output/hybrid_abundance.csv")
