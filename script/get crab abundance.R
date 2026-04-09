#Produce crab abundance survey time series from crabpack for 
#1) hybrids
#2) snow crab
#3) tanner crab


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

#calculate hybrid population abundance timeseries with minimum size cutoff 
hybrid_50mmplus <- calc_bioabund(crab_data = hybrid,
                               species = "HYBRID",
                               region = "EBS",
                               size_min = 50,
                               years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population_50mm_plus",
         species = "hybrid")

#Plot
hybrid_50mmplus %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#Mature females/large males 
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
         species = "hybrid") 
 
#Plot timeseries
hybrid_lg_mat %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#Pre-recruits, i.e. 65-80mm males, newshell
hybrid_prerecruit <- calc_bioabund(crab_data = hybrid,
                               species = "HYBRID",
                               region = "EBS",
                               sex = "male", 
                               size_min = 65, 
                               size_max = 80, 
                               shell_condition = "new_hardshell",
                               years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "prerecruit",
         species = "hybrid") 

#Plot timeseries
hybrid_prerecruit %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#join all hybrid time series 
hyb_abun_pop %>%
  full_join(hybrid_50mmplus) %>%
  full_join(hybrid_lg_mat) %>%
  full_join(hybrid_prerecruit) -> hyb_abun_dat

#and plot on a single plot
hyb_abun_dat %>%
  ggplot(aes(year, abundance, color=category)) +
  geom_point() +
  geom_line() +
  theme_bw() 

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

#calculate snow crab population abundance timeseries with minimum size cutoff 
snow_50mmplus <- calc_bioabund(crab_data = snow,
                                 species = "SNOW",
                                 region = "EBS",
                                 size_min = 50,
                                 years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population_50mm_plus",
         species = "snow")

#Plot
snow_50mmplus %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#Pre-recruits, i.e. 65-80mm males, newshell
snow_prerecruit <- calc_bioabund(crab_data = snow,
                                   species = "SNOW",
                                   region = "EBS",
                                   sex = "male", 
                                   size_min = 65, 
                                   size_max = 80, 
                                   shell_condition = "new_hardshell",
                                   years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "prerecruit",
         species = "snow") 

#Plot timeseries
snow_prerecruit %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

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

#join all snow crab time series 
snow_abun_pop %>%
  full_join(snow_50mmplus) %>%
  full_join(snow_lg_mat) %>%
  full_join(snow_prerecruit) -> snow_abun_dat

#and plot on a single plot
snow_abun_dat %>%
  ggplot(aes(year, abundance, color=category)) +
  geom_point() +
  geom_line() +
  theme_bw()

#-----------------------------
# Tanner Crab abundance timeseries
#----------------------------- 

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

#calculate tanner crab population abundance timeseries with minimum size cutoff 
tanner_50mmplus <- calc_bioabund(crab_data = tanner,
                               species = "TANNER",
                               region = "EBS",
                               spatial_level = "region",
                               size_min = 50,
                               years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "population_50mm_plus",
         species = "tanner")

#Plot
tanner_50mmplus %>%
  ggplot(aes(x = year, y = abundance)) +
  geom_point() +
  geom_line()+
  labs(y = "Number of crab (millions)", x = "") +
  theme_bw()

#Pre-recruits, i.e. 65-80mm males, newshell
tanner_prerecruit <- calc_bioabund(crab_data = tanner,
                                 species = "TANNER",
                                 region = "EBS",
                                 spatial_level = "region",
                                 sex = "male", 
                                 size_min = 65, 
                                 size_max = 80, 
                                 shell_condition = "new_hardshell",
                                 years = years) %>%
  select(YEAR, ABUNDANCE, ABUNDANCE_CI) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  mutate(ABUNDANCE = as.numeric(ABUNDANCE/1e6),
         ABUNDANCE_CI = as.numeric(ABUNDANCE_CI/1e6)) %>%
  rename_with(tolower) %>%
  mutate(category = "prerecruit",
         species = "tanner") 

#Plot timeseries
tanner_prerecruit %>%
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

#join all tanner crab time series 
tanner_abun_pop %>%
  full_join(tanner_50mmplus) %>%
  full_join(tanner_lg_mat) %>%
  full_join(tanner_prerecruit) -> tanner_abun_dat

#and plot on a single plot
tanner_abun_dat %>%
  ggplot(aes(year, abundance, color=category)) +
  geom_point() +
  geom_line() +
  theme_bw()

#-------------------------------------------------------------
# Proportion of legal males that are hybrids/tanners/snows 
#------------------------------------------------------------- 

# Calculate total abundance of 101mm snow/tanner/hybrid males
snow_legal <- calc_bioabund(crab_data = snow,
                     species ="SNOW",
                     region = "EBS",
                     sex = "male",
                     size_min = 101,
                     years = years) %>%
              select(SPECIES, YEAR, ABUNDANCE, ABUNDANCE_CI)

hybrid_legal <- calc_bioabund(crab_data = hybrid,
                            species ="HYBRID",
                            region = "EBS",
                            sex = "male",
                            size_min = 101,
                            years = years) %>%
                select(SPECIES, YEAR, ABUNDANCE, ABUNDANCE_CI)

tanner_legal <- calc_bioabund(crab_data = tanner,
                            species ="TANNER",
                            region = "EBS",
                            spatial_level = "region",
                            sex = "male",
                            size_min = 101,
                            years = years) %>%
                select(SPECIES, YEAR, ABUNDANCE, ABUNDANCE_CI)

#join and calculate ratio
ratio <- snow_legal %>% 
  bind_rows(hybrid_legal) %>%
  bind_rows(tanner_legal) %>%
  pivot_wider(names_from = SPECIES, values_from = c(ABUNDANCE, ABUNDANCE_CI)) %>%
  mutate(total = ABUNDANCE_SNOW + ABUNDANCE_HYBRID + ABUNDANCE_TANNER,
         p_snow   = ABUNDANCE_SNOW / total,
         p_hybrid = ABUNDANCE_HYBRID / total,
         p_tanner = ABUNDANCE_TANNER / total) %>%
  rename_with(tolower) %>%
  select(year, p_snow, p_hybrid, p_tanner) %>%
  pivot_longer(cols = starts_with("p_"),
    names_to = "stock", values_to = "prop") %>%
  mutate(stock = recode(stock, "p_snow" = "Snow",
                          "p_tanner" = "Tanner", "p_hybrid" = "Hybrid"))

#plot
ratio %>%
  ggplot(aes(year, prop, color = stock, fill = stock)) +
  geom_point() +
  geom_line() 

#-----------------------------
# Join datasets and output
#-----------------------------

hyb_abun_dat %>%
  full_join(snow_abun_dat) %>%
  full_join(tanner_abun_dat) -> abun_dat

#Write output for abundance 
write_csv(abun_dat, "./output/crab_abundance.csv")

#write output for proportion legal
write_csv(ratio, "./output/proportion_legal.csv")
