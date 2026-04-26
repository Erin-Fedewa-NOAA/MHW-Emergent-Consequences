#Hybrid zone exploration

#Author: EJF

## Load packages
library(crabpack)
library(tidyverse)
library(patchwork)
library(ggrepel)

# Set years ----
current_year <- 2025
years <- 1988:current_year

#define corner stations so we can drop
stations <- read.csv("./data/station_lookup.csv")

corners <- stations %>% 
  filter(REGION == "EBS", STATION_TYPE == "MTCA_CORNER") %>%
  pull(STATION_ID)

#------------------------------
# #data wrangling
#------------------------------

## Pull specimen data from crabpack 
snow <- get_specimen_data(species = "SNOW",
                          region = "EBS",
                          channel = "KOD")

tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")

hybrid <- get_specimen_data(species = "HYBRID",
                            region = "EBS",
                            channel = "KOD")

## Pull station-level cpue data from crabpack 
snow_cpue <- crabpack::calc_cpue(crab_data = snow,
                                     region = "EBS",
                                     species = "SNOW",
                                     years = years)

tanner_cpue <- crabpack::calc_cpue(crab_data = tanner,
                                       region = "EBS",
                                       species = "TANNER",
                                       years = years) 

hybrid_cpue <- crabpack::calc_cpue(crab_data = hybrid,
                                   region = "EBS",
                                   species = "HYBRID",
                                   years = years) 

#join
crab_cpue <- snow_cpue %>%
  rename_with(tolower) %>%
  select(year, station_id, latitude, longitude, cpue) %>%
  rename(snow_cpue = cpue) %>%
  full_join(tanner_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, latitude, longitude, cpue) %>%
              rename(tanner_cpue = cpue)) %>%
  full_join(hybrid_cpue %>%
              rename_with(tolower) %>%
              select(year, station_id, latitude, longitude, cpue) %>%
              rename(hybrid_cpue = cpue)) %>%
  filter(!(station_id %in% corners)) 

#------------------------------
# #Spatial Plots
#------------------------------

#hybrid fraction
crab_cpue2 <- crab_cpue %>%
  mutate(total = snow_cpue + tanner_cpue + hybrid_cpue,
         hybrid_frac = ifelse(total > 0, hybrid_cpue / total, NA))

ggplot(crab_cpue2, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = hybrid_frac, size = hybrid_cpue)) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  facet_wrap(~year) +
  coord_fixed() +
  theme_minimal()

#hybrid presence 
crab_cpue2 <- crab_cpue2 %>%
  mutate(hybrid_present = hybrid_cpue > 0)

ggplot(crab_cpue2, aes(longitude, latitude)) +
  geom_point(aes(color = hybrid_present)) +
  facet_wrap(~year)

#centroid of hybrid occurrence 
centroids <- crab_cpue %>%
  filter(hybrid_cpue > 0) %>%   # only where hybrids occur
  group_by(year) %>%
  summarize(
    centroid_lat = weighted.mean(latitude, hybrid_cpue, na.rm = TRUE),
    centroid_lon = weighted.mean(longitude, hybrid_cpue, na.rm = TRUE),
    total_hybrid = sum(hybrid_cpue, na.rm = TRUE))

ggplot(centroids, aes(x = centroid_lon, y = centroid_lat)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = year), size = 3) +
  coord_fixed() +
  theme_minimal()

#------------------------------
# #Hybrid zone metrics
#------------------------------

#calculate hybrid fraction
crab2 <- crab_cpue %>%
  mutate(
    total_cpue = snow_cpue + tanner_cpue + hybrid_cpue,
    hybrid_frac = ifelse(total_cpue > 0, hybrid_cpue / total_cpue, NA)) %>%
  filter(total_cpue > 0)   # remove empty stations

#define spatial axis - defining as long here so easier to intepret, but can 
  #use a PCA with lat/long
crab2 <- crab2 %>%
  mutate(x = longitude)

#fit clines for each year
fit_cline <- function(df) {
  glm(cbind(hybrid_cpue, snow_cpue + tanner_cpue) ~ x,
      family = binomial,
      data = df)
}

models <- crab2 %>%
  group_by(year) %>%
  group_split() %>%
  setNames(unique(crab2$year)) %>%
  map(fit_cline)

#extract hybrid zone metrics 
  #center = where hybrid fraction = 0.5
  #width = distance between 10% and 90%

extract_metrics <- function(model) {
  b0 <- coef(model)[1]
  b1 <- coef(model)[2]
  
  # logistic inverse
  logit <- function(p) log(p / (1 - p))
  
  center <- -b0 / b1
  x10 <- (logit(0.1) - b0) / b1
  x90 <- (logit(0.9) - b0) / b1
  
  width <- abs(x90 - x10)
  
  data.frame(center = center, width = width)
}

#------------------------------
# #Hybrid zone visualization
#------------------------------

# Plot clines for specific years
plot_cline <- function(year_val) {
  df <- crab2 %>% filter(year == year_val)
  mod <- models[[as.character(year_val)]]
  
  newx <- data.frame(x = seq(min(df$x), max(df$x), length.out = 100))
  preds <- predict(mod, newdata = newx, type = "response")
  
  ggplot(df, aes(x, hybrid_frac)) +
    geom_point(alpha = 0.4) +
    geom_line(data = data.frame(x = newx$x, y = preds),
              aes(x, y), color = "red") +
    ggtitle(paste("Year:", year_val))
}

plot_cline(2024)
#We really don't see a typical s-shaped cline in any years- suggests a localized 
  #region of hybrid occurrence rather than a well-defined clinal transition

metrics <- map_df(models, extract_metrics, .id = "year") %>%
  mutate(year = as.numeric(year))

#hybrid zone centroid mvmt
ggplot(metrics, aes(year, center)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(y = "Hybrid zone center")

#hybrid zone width
ggplot(metrics, aes(year, width)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(y = "Hybrid zone width")
