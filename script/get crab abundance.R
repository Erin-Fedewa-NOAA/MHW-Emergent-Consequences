#Produce crab abundance survey output for 
#1) hybrid population
#2) snow and tanner population
#3) immature snow/tanner
#4) mature female and large male tanner/snow crab (large male cutoffs
  #represent size at functional maturity b/c we're interested in portion
  #of the population participating in the mate)

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

#-----------------------------
# Hybrid abundance timeseries
#-----------------------------

#calculate hybrid population abundance timeseries 
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

#now let's subset hybrids for a representative psuedocohort
hyb_sub <- calc_bioabund(crab_data = hybrid,
                              species = "HYBRID",
                              region = "EBS",
                              size_min = 50,
                              size_max = 65,
                              years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "50_65mm",
         species = "hybrid") 

#and now mature females/large males only
hybrid_lg_mat <- calc_bioabund(crab_data = hybrid,
                             species = "HYBRID",
                             region = "EBS",
                             crab_category = c("mature_female",
                                               "large_male"),
                             years = years) %>%
  select(YEAR, ABUNDANCE, CATEGORY) %>%
  pivot_wider(names_from = CATEGORY, values_from = ABUNDANCE) %>%
  group_by(YEAR) %>%
  mutate(ABUNDANCE = sum(large_male + mature_female)) %>%
  select(YEAR, ABUNDANCE) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "maturefem_lgmale",
         species = "hybrid") %>%
  #join to full hybrid population/psuedo-cohort time series 
  full_join(hyb_abun_pop) %>%
  full_join(hyb_sub) -> hyb_abun

#Plot timeseries
hyb_abun %>%
  ggplot(aes(x = year, y = abundance, color=category)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#Lastly, to inform a hybrid population abundance metric that best represents 
  #the peak seen in 2024-2025, let's look at the size interval that contains 
  #the central 90% of the abundance for 2025

size_range <- calc_bioabund(crab_data = hybrid,
                      species = "HYBRID",
                      region = "EBS",
                      years = 2025,
                      bin_1mm = T) %>%
  select(YEAR, ABUNDANCE, SIZE_1MM) %>%
  arrange(SIZE_1MM, .by_group = TRUE) %>%
  mutate(total_abundance = sum(ABUNDANCE),
        cum_abundance = cumsum(ABUNDANCE),
        cum_prop = cum_abundance / total_abundance) %>%
  summarise(lower_90 = SIZE_1MM[which(cum_prop >= 0.05)[1]],
             upper_90 = SIZE_1MM[which(cum_prop >= 0.95)[1]], .groups = "drop")

size_range #48 - 103mm

#let's plot this now
calc_bioabund(crab_data = hybrid,
              species = "HYBRID",
              region = "EBS",
              years = 2025,
              bin_1mm = T) %>%
  select(YEAR, ABUNDANCE, SIZE_1MM) %>%
  ggplot(aes(x = SIZE_1MM, y = ABUNDANCE)) +
  geom_col(width = 1, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 48, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = 103, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Carapace width (mm)", y = "Hybrid Abundance") +
  theme_minimal()

#And then calculate an abundance timeseries using these cutoffs
hybrid_90perc <- calc_bioabund(crab_data = hybrid,
                               species = "HYBRID",
                               region = "EBS",
                               size_min = 48,
                               size_max = 103,
                               years = years) %>%
  select(YEAR, ABUNDANCE) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population_subset",
         species = "hybrid")

#and join to full hybrid population/psuedo-cohort time series 
hybrid_90perc %>%
full_join(hyb_abun) -> hyb_abun_dat

#-----------------------------
# Snow Crab abundance timeseries
#----------------------------- 

#calculate snow crab population abundance timeseries 
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

#------------------------------------------------
# Immature snow and tanner crab abundance timeseries
#-------------------------------------------------------

#immature female snow/tanner abundance timeseries 
snow_imm_female <- calc_bioabund(crab_data = snow,
                           species = "SNOW",
                           region = "EBS",
                           size_min = 30,
                           size_max = 60,
                           crab_category = "immature_female",
                           years = years) %>%
  select(YEAR, ABUNDANCE, CATEGORY) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE_IMMFEM = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  select(-category, -abundance) %>%
  mutate(species = "snow")

tanner_imm_female <- calc_bioabund(crab_data = tanner,
                                 species = "TANNER",
                                 spatial_level = "region",
                                 size_min = 30,
                                 size_max = 60,
                                 crab_category = "immature_female",
                                 years = years) %>%
  select(YEAR, ABUNDANCE, CATEGORY) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE_IMMFEM = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  select(-category, -abundance) %>%
  mutate(species = "tanner")


#add female abundance to small male abundance for both species 
snow_immature <- calc_bioabund(crab_data = snow,
                                 species = "SNOW",
                                 region = "EBS",
                                 sex = "male",
                                 size_min = 30,
                                 size_max = 60,
                                 years = years) %>%
  select(YEAR, ABUNDANCE) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE_SMMALE = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  select(-abundance) %>%
  mutate(species = "snow") %>%
  full_join(snow_imm_female) %>%
  group_by(year) %>%
  mutate(abundance = abundance_smmale + abundance_immfem,
         category = "immature") %>%
  select(-abundance_smmale, -abundance_immfem)

#Plot
snow_immature %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

tanner_immature <- calc_bioabund(crab_data = tanner,
                               species = "TANNER",
                               spatial_level = "region",
                               sex = "male",
                               size_min = 30,
                               size_max = 60,
                               years = years) %>%
  select(YEAR, ABUNDANCE) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE_SMMALE = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  select(-abundance) %>%
  mutate(species = "tanner") %>%
  full_join(tanner_imm_female) %>%
  group_by(year) %>%
  mutate(abundance = abundance_smmale + abundance_immfem,
         category = "immature") %>%
  select(-abundance_smmale, -abundance_immfem)

#Plot
tanner_immature %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#-----------------------------------------------
# Mature female/large male abundance timeseries
#--------------------------------------------------

#calculate mature female/lg male snow crab abundance 
snow_lg_mat <- calc_bioabund(crab_data = snow,
                                species = "SNOW",
                                region = "EBS",
                                crab_category = c("mature_female",
                                                  "large_male"),
                                years = years) %>%
  select(YEAR, ABUNDANCE, CATEGORY) %>%
  pivot_wider(names_from = CATEGORY, values_from = ABUNDANCE) %>%
  group_by(YEAR) %>%
  mutate(ABUNDANCE = sum(large_male + mature_female)) %>%
  select(YEAR, ABUNDANCE) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "maturefem_lgmale",
         species = "snow")

#Plot
snow_lg_mat  %>%
  filter(category != "NA") %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#calculate mature female/lg male tanner crab abundance 
tanner_lg_mat <- calc_bioabund(crab_data = tanner,
                             species = "TANNER",
                             region = "EBS",
                             spatial_level = "region",
                             crab_category = c("mature_female",
                                               "large_male"),
                             years = years) %>%
  select(YEAR, ABUNDANCE, CATEGORY) %>%
  pivot_wider(names_from = CATEGORY, values_from = ABUNDANCE) %>%
  group_by(YEAR) %>%
  mutate(ABUNDANCE = sum(large_male + mature_female)) %>%
  select(YEAR, ABUNDANCE) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "maturefem_lgmale",
         species = "tanner")

#Plot
tanner_lg_mat  %>%
  filter(category != "NA") %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#-----------------------------
# Join datasets and output
#-----------------------------

hyb_abun_dat %>%
  full_join(snow_abun_pop) %>%
  full_join(tanner_abun_pop) %>%
  full_join(snow_lg_mat) %>%
  full_join(tanner_lg_mat) %>%
  full_join(snow_immature) %>%
  full_join(tanner_immature) %>%
  select(-abundance_ci) -> abun_dat

#Write output 
write_csv(abun_dat, "./output/crab_abundance.csv")
