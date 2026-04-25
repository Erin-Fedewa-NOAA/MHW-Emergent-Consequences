# Calculate "D95" for Tanner, hybrid and snow crab:
#area of stations that make up 95% of the cumulative cpue

#Note that these are population-wide metrics only- large males is probably most
  #meaningful 

#Author: Erin Fedewa

# load ----
library(tidyverse)
library(crabpack)

# Set years ----
current_year <- 2025
years <- 1988:current_year

## Pull specimen data from crabpack ----
snow <- get_specimen_data(species = "SNOW",
                          region = "EBS",
                          channel = "KOD")

tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")

hybrid <- get_specimen_data(species = "HYBRID",
                            region = "EBS",
                            channel = "KOD")

#Calculate station-level CPUE ----
snow_cpue <- calc_cpue(crab_data = snow,
                            species = "SNOW",
                            years = years)

tanner_cpue <- calc_cpue(crab_data = tanner,
                       species = "TANNER",
                       years = years,
                       district = "ALL")

hybrid_cpue <- calc_cpue(crab_data = hybrid,
                species = "HYBRID",
                years = years)

#define corner stations
stations <- read.csv("Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/lookup_tables/station_lookup.csv")

corners <- stations %>% 
  filter(STATION_TYPE == "MTCA_CORNER") %>%
  pull(STATION_ID)


##########################################
# compute D95 for snow crab ----
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
snow_cpue %>%
  filter(!(STATION_ID %in% corners)) %>% #exclude corner stations
  nest(data = -YEAR) %>%
  mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
  unnest(cols = c(data)) %>%
  group_by(YEAR) %>%
  summarise(mean_cpue = mean(CPUE), # add a column for mean cpue of each group in each year
            d95 = mean(d95)) -> d95_snow # take 'mean' just to get one value (they are all the same)

#plot 
d95_snow %>%
  ggplot(aes(x = YEAR, y = d95))+
  geom_point(size=3)+
  geom_line() +
  theme_bw() 

#d95 vs. abund plot
d95_snow %>%
  ggplot(aes(x = mean_cpue, y = d95)) +
  geom_point() +
  # geom_line() +
  geom_smooth(method = 'lm') +
  labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
  theme_bw() +
  theme(legend.title = element_blank()) 

##########################################
# now compute D95 for tanner crab ----

tanner_cpue %>%
  filter(!(STATION_ID %in% corners)) %>% #exclude corner stations
  nest(data = -YEAR) %>%
  mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
  unnest(cols = c(data)) %>%
  group_by(YEAR) %>%
  summarise(mean_cpue = mean(CPUE), # add a column for mean cpue of each group in each year
            d95 = mean(d95)) -> d95_tanner # take 'mean' just to get one value (they are all the same)

#plot 
d95_tanner %>%
  ggplot(aes(x = YEAR, y = d95))+
  geom_point(size=3)+
  geom_line() +
  theme_bw() 

#d95 vs. abund plot
d95_tanner %>%
  ggplot(aes(x = mean_cpue, y = d95)) +
  geom_point() +
  # geom_line() +
  geom_smooth(method = 'lm') +
  labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
  theme_bw() +
  theme(legend.title = element_blank()) 

#combine datasets
d95_tanner %>%
  select(YEAR, d95) %>%
  rename(tanner_area=d95) %>%
  full_join(d95_snow %>%
              select(YEAR, d95) %>%
              rename(snow_area=d95)) -> area_occupied

#plot
area_occupied %>%
  ggplot(aes(tanner_area, snow_area)) +
  geom_point() +
  geom_smooth(method = 'lm')

##########################################
# and compute D95 for hybrids ----

hybrid_cpue %>%
  filter(!(STATION_ID %in% corners)) %>% #exclude corner stations
  nest(data = -YEAR) %>%
  mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
  unnest(cols = c(data)) %>%
  group_by(YEAR) %>%
  summarise(mean_cpue = mean(CPUE), # add a column for mean cpue of each group in each year
            d95 = mean(d95)) -> d95_hybrid # take 'mean' just to get one value (they are all the same)

#plot 
d95_hybrid %>%
  ggplot(aes(x = YEAR, y = d95))+
  geom_point(size=3)+
  geom_line() +
  theme_bw() 

#d95 vs. abund plot
d95_hybrid %>%
  ggplot(aes(x = mean_cpue, y = d95)) +
  geom_point() +
  # geom_line() +
  geom_smooth(method = 'lm') +
  labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
  theme_bw() +
  theme(legend.title = element_blank()) 

############################################

#combine datasets
d95_tanner %>%
  select(YEAR, d95) %>%
  rename(tanner_area=d95) %>%
  full_join(d95_hybrid %>%
  select(YEAR, d95) %>%
  rename(hybrid_area=d95)) %>%
  full_join(d95_snow %>%
              select(YEAR, d95) %>%
              rename(snow_area=d95)) -> area_occupied

#Write output 
missing <- data.frame(YEAR = 2020)

area_occupied %>%
  bind_rows(missing) %>%
  arrange(YEAR) %>%
  write.csv(file="./output/area_occupied_output.csv", row.names = F)

