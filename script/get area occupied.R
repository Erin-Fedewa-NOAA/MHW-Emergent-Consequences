# Calculate "D95" for each size sex group in EBS:
#area of stations that make up 95% of the cumulative cpue

#Author: Erin Fedewa

#Follow ups: 
#move toward spatiotemporal modeling approach to account for changing footprint
#current approach is only using stations sampled every year

# load ----
library(tidyverse)

## Read in setup
source("./Scripts/get_crab_data.R")

##########################################
#Assign maturity to specimen data; calculate CPUE
cpue <- snow$specimen %>% 
  left_join(., mat_size) %>%
  mutate(CATEGORY = case_when((SEX == 1 & SIZE >= MAT_SIZE) ~ "mature_male",
                              (SEX == 1 & SIZE < MAT_SIZE) ~ "immature_male",
                              (SEX == 2 & CLUTCH_SIZE >= 1) ~ "mature_female",
                              (SEX == 2 & CLUTCH_SIZE == 0) ~ "immature_female",
                              TRUE ~ NA)) %>%
  filter(YEAR %in% years,
         !is.na(CATEGORY)) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE, AREA_SWEPT, CATEGORY) %>%
  summarise(COUNT = round(sum(SAMPLING_FACTOR))) %>%
  pivot_wider(names_from = CATEGORY, values_from = COUNT) %>%
  mutate(population = sum(immature_male, mature_male, immature_female, mature_female, na.rm = T)) %>%
  pivot_longer(c(6:10), names_to = "CATEGORY", values_to = "COUNT") %>%
  filter(CATEGORY != "NA") %>%
  mutate(COUNT = replace_na(COUNT, 0),
         CPUE = COUNT / AREA_SWEPT) %>%
  ungroup() 

# compute D95 by each size and sex category ----
# i.e. the number of stations contributing to 95% of cumulative cpue

# function to compute D95
f_d95_est <- function(x){
  x %>%
    arrange(-CPUE) %>% #sort by cpue (large:small)
    mutate(prop_cpue = CPUE/sum(CPUE),  #calculate the proportion of total cpue for each station
           cum_cpue = cumsum(prop_cpue)) %>%  
    filter(cum_cpue <= 0.95) %>% #T if in d95, F if not
    count() %>%
    mutate(d95 = (n + 1) * 401) %>% #add 1 station to n to push over 95%, multiply by 401 nm
    pull(d95)
}

# do the estimation
cpue %>%
  filter(!(STATION_ID %in% corners)) %>% #exclude corner stations
  nest(data = c(-YEAR, -CATEGORY)) %>%
  mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
  unnest(cols = c(data)) %>%
  group_by(YEAR, CATEGORY) %>%
  summarise(mean_cpue = mean(CPUE), # add a column for mean cpue of each group in each year
            d95 = mean(d95)) -> d95 # take 'mean' just to get one value (they are all the same)

#plot by size/sex
d95 %>%
  select(YEAR, CATEGORY, d95) %>%
  ggplot(aes(x = YEAR, y = d95, group= CATEGORY, color = CATEGORY))+
  geom_point(size=3)+
  geom_line() +
  theme_bw() +
  facet_wrap(~CATEGORY)

#plot just mature males, our indicator
d95 %>%
  select(YEAR, CATEGORY, d95) %>%
  filter(CATEGORY == "mature_male") %>%
  ggplot(aes(x = YEAR, y = d95))+
  geom_point(size=3)+
  geom_line() +
  geom_hline(aes(yintercept = mean(d95, na.rm=TRUE)), linetype = 5) +
  theme_bw()

#d95 vs. abund plot
d95 %>%
  filter(!CATEGORY == "population") %>%
  ggplot(aes(x = mean_cpue, y = d95, group = CATEGORY, color = CATEGORY)) +
  geom_point() +
  # geom_line() +
  geom_smooth(method = 'lm') +
  labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~CATEGORY, scales = "free")
#interesting male vrs female relationship!

#d95 vs. bottom temperature plot
haul %>%
  filter(!HAUL_TYPE == 17) %>%
  distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
  group_by(YEAR) %>%
  summarise(summer_bt = mean(GEAR_TEMPERATURE, na.rm = T)) %>%
  right_join(d95, by="YEAR") %>%
  ggplot(aes(x = summer_bt, y = d95, group = CATEGORY, color = CATEGORY)) +
  geom_point() +
  # geom_line() +
  geom_smooth(method = 'lm') +
  labs(x = "Bottom Temperature (C)", y = expression("Area Occupied ("~nmi^2~")")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~CATEGORY, scales = "free")


#Write output for D95 indicator
missing <- data.frame(YEAR = 2020)

d95 %>%
  select(-mean_cpue) %>%
  pivot_wider(names_from = "CATEGORY", values_from = "d95") %>%
  rename(immature_female_d95=immature_female,
         immature_male_d95=immature_male,
         mature_female_d95=mature_female,
         mature_male_d95=mature_male,
         population_d95=population) %>%
  bind_rows(missing) %>%
  arrange(YEAR) %>%
  write.csv(file="./Output/D95_output.csv", row.names = F)

########################################################################
#compare area occupied to core range


#Calculate station-level snow crab CPUE ----
cpue <- crabpack::calc_cpue(crab_data = snow,
                            species = "SNOW",
                            years = years)

#Define snow crab core area for spatially subsetting indicators ----
#i.e. stations containing 50th percentile of snow crab CPUE
cpue50_core <- cpue %>%
  filter(!STATION_ID %in% corners) %>% #remove corner stations
  group_by(STATION_ID) %>%
  summarise(MEAN_CPUE = mean(CPUE)) %>%
  filter(MEAN_CPUE > quantile(MEAN_CPUE, 0.50)) %>%
  left_join(., stations %>% select(STATION_ID, LATITUDE, LONGITUDE))

#Quick plot for visual check 
ggplot(ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf() +
  geom_point(data = cpue50_core, aes(x = LONGITUDE, y = LATITUDE), size = 2, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-180, -159), ylim = c(54, 64), expand = FALSE) +
  theme_bw() #looks about right! 

#Write csv for stations in 50th percentile of avg CPUE  
write.csv(cpue50_core, file="./Output/snow_core_area.csv", row.names = FALSE)

#Pull stations for filtering in indicator script
snow_core <- cpue50_core %>%
  pull(STATION_ID)