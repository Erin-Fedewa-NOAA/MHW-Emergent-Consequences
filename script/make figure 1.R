#Goals ----
#Create Figure 1 for manuscript 

#Author: EJF

## Load packages
library(RColorBrewer)
library(tidyverse)
library(crabpack)
library(patchwork)
library(akgfmaps)
library(sf)
library(sp)
library(adehabitatHR)

## Load Alaska Region layers (FROM AKGFMAPS R package) 
ebs_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")
ebs_survey_areas <- ebs_layers$survey.area

#--------------------------------#
#Range map ----
#--------------------------------#
#Pull specimen data for each species
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
snow_cpue <- crabpack::calc_cpue(crab_data = snow,
                            species = "SNOW",
                            years = 1988:2025)

tanner_cpue <- crabpack::calc_cpue(crab_data = tanner,
                                 species = "TANNER",
                                 years = 1988:2025)

hybrid_cpue <- crabpack::calc_cpue(crab_data = hybrid,
                                 species = "HYBRID",
                                 years = 1988:2025)

#combine and convert to spatial points
pts <- snow_cpue %>%
  bind_rows(tanner_cpue) %>%
  bind_rows(hybrid_cpue) %>%
  filter(CPUE > 0) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
  st_transform(3338)

# function to calculate 50% KDE polygon for each species
  #using the 
calc_kde <- function(x){
  sp_obj <- as(x, "Spatial")
    #estimate kernel density surface using adehabitatHR package
    kde <- kernelUD(sp_obj, h = "href")
    #extract 50% contour
  getverticeshr(kde, percent = 50) %>%
    st_as_sf()
}

# calculate polygons species-by-species
core50 <- pts %>%
  group_split(SPECIES) %>%
  map_dfr( ~ calc_kde(.x) %>%
             mutate(SPECIES = unique(.x$SPECIES)))

#plot
map <- core50 %>%
  mutate(SPECIES = factor(SPECIES, levels = c("SNOW", "TANNER", "HYBRID"))) %>%
ggplot() +
  geom_sf(aes(fill = SPECIES), alpha = 0.3, color = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey85",
          color = "grey50", linewidth = 0.2) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  coord_sf(crs = 3338) +
  scale_fill_manual(values = c(SNOW = "#0072B2",
                                TANNER = "#D55E00",
                                HYBRID = "#009E73"),
    breaks = c("SNOW", "TANNER", "HYBRID"),
    labels = c("Snow crab", "Tanner crab",
      "Snow × Tanner hybrids")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
      plot.margin = margin(0, 0, 0, 0))
       #legend.title = element_blank(),
       #legend.position = c(0.07, 0.01),
       #legend.justification = c(0, 0),
       #legend.text = element_text(size = 6)) 

#--------------------------------#
#Abundance plots ----
#--------------------------------#

#data
abundance <- read.csv("./output/crab_abundance.csv") %>%
                filter(category == "population") %>%
                mutate(species = recode(species, "hybrid" = "Snow × Tanner hybrids",
                           "tanner" = "Tanner crab", "snow" = "Snow crab"))
  
#faceted plot 
abun_plot <- ggplot(abundance, aes(year, abundance, color = species)) +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci,
                  fill = species, group = species),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8, lineend = "round") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
  facet_wrap(~factor(species, levels = c("Snow crab", "Tanner crab", "Snow × Tanner hybrids")),
    scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Tanner crab" = "#D55E00",
                                "Snow crab"   = "#0072B2",
                                "Snow × Tanner hybrids" = "#4DA02C")) +
  scale_fill_manual(values = c("Tanner crab" = "#D55E00",
                               "Snow crab"   = "#0072B2",
                               "Snow × Tanner hybrids" = "#4DA02C")) +
  labs(y = "Abundance (millions)", x = "") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(panel.spacing.y = unit(0.15, "lines"),
    strip.background = element_blank(),
    strip.clip = "off",
    strip.text = element_text(size = 11, margin = margin(b = -2)),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 5))

#snow crab plot
snow_plot <- abundance %>%
  filter(species == "Snow crab") %>%
  ggplot(aes(year, abundance)) +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci),
                  fill = "#0072B2", alpha = 0.12) +
  geom_line(linewidth = 0.9, lineend = "round", color = "#0072B2") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
 labs(y = "Abundance (millions)", x = "") +
  theme_classic(base_size = 10) +
  annotate("text", x = 2008, y = Inf,
              label = "Snow crab", hjust = 0.5, vjust = 1.3,
              size = 3.5, fontface = "bold", color="grey50") +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.4),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "plain"),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey96", linewidth = 0.3))

#tanner crab plot
tanner_plot <- abundance %>%
  filter(species == "Tanner crab") %>%
  ggplot(aes(year, abundance)) +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci),
              fill = "#D55E00", alpha = 0.12) +
  geom_line(linewidth = 0.9, lineend = "round", color = "#D55E00") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
  labs(y = "Abundance (millions)", x = "") +
  theme_classic(base_size = 10) +
  annotate("text", x = 2008, y = Inf,
           label = "Tanner crab", hjust = 0.5, vjust = 1,
           size = 3.5, fontface = "bold", color="grey50") +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.4),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "plain"),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey96", linewidth = 0.3))

#hybrid plot
hybrid_plot <- abundance %>%
  filter(species == "Snow × Tanner hybrids") %>%
  ggplot(aes(year, abundance)) +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci),
              fill = "#4DA02C", alpha = 0.12) +
  geom_line(linewidth = 0.9, lineend = "round", color = "#4DA02C") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
  labs(y = "Abundance (millions)", x = "") +
  theme_classic(base_size = 10) +
  annotate("text", x = 2008, y = Inf,
           label = "Snow × Tanner hybrids", hjust = 0.5, vjust = 1.3,
           size = 3.5, fontface = "bold", color="grey50") +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.4),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "plain"),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey96", linewidth = 0.3))


#--------------------------------------#
#Sea Ice plot ----
#--------------------------------------#

#data
ice <- read.csv("./output/seaice_output.csv") %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
  filter(year >= 1988) %>%
  mutate(ice = Mar_Apr_ice_EBS_NBS *100)

#plot
ice_plot <- ggplot(ice, aes(year, ice)) +
  geom_line(linewidth = 0.8, lineend = "round", color="#0072B2") +
  labs(y = "Sea ice extent (%)", x = "") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
            alpha=0.1, fill= "grey60") +
  annotate("text", x = 2007, y = 76,
           label = "Bering Sea ice extent", hjust = 0.5, vjust = 0.3,
           size = 3.5, fontface = "bold", color="grey50") +
  theme_classic(base_size = 10) +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.4),
        panel.grid.major = element_line(color = "grey96", linewidth = 0.3),
        panel.grid.minor = element_blank())
    
#--------------------------------------#
#Proportion 4 inch males plot ----
#--------------------------------------#

#data
ratio <- read.csv("./output/proportion_hybrid.csv")

#plot
ratio_plot <- ratio %>%
  ggplot(aes(year, prop_hybrid)) +
  geom_line(linewidth = 0.8, lineend = "round", color = "#4DA02C") +
  theme_classic(base_size = 10) +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
  geom_ribbon(aes(ymin = prop_lower, ymax = prop_upper), alpha = 0.3, fill = "#4FA32D") +
  annotate("text", x = 2007, y = 76,
           label = "Commercial-sized hybrid fraction", hjust = 0.5, vjust = 0.3,
           size = 3.5, fontface = "bold", color="grey50") +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.4),
        panel.grid.major = element_line(color = "grey96", linewidth = 0.3),
        panel.grid.minor = element_blank()) +
  labs(y="Hybrid proportion (%)", x="") 

#--------------------------------#
#Combine and save figures ----
#--------------------------------#

fig1 <- (free(map) + ice_plot + snow_plot + tanner_plot + hybrid_plot + ratio_plot) +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("./figures/fig1.jpeg", fig1,
  width = 7.2, height = 8.5, units = "in", dpi = 600)

#--------------------------------------------------------------#
#Population increase calculations (Results paragraph 1)  ----
#--------------------------------------------------------------#

#tanner crab
tanner <- read.csv("./output/crab_abundance.csv") %>%
  filter(category == "population" & species == "tanner")

historic_baseline <- tanner %>%
  filter(year >= 1988, year <= 2023) %>%
  summarise(mean_baseline = mean(abundance, na.rm = TRUE)) %>%
  pull(mean_baseline)

recent_mean <- tanner %>%
  filter(year %in% c(2024, 2025)) %>%
  summarise(mean_recent = mean(abundance, na.rm = TRUE)) %>%
  pull(mean_recent)

fold_change <- recent_mean / historic_baseline

#hybrids
hybrid <- read.csv("./output/crab_abundance.csv") %>%
  filter(category == "population" & species == "hybrid")

historic_baseline <- hybrid %>%
  filter(year >= 1988, year <= 2024) %>%
  summarise(mean_baseline = mean(abundance, na.rm = TRUE)) %>%
  pull(mean_baseline)

prior_max <- hybrid %>%
  filter(year <= 2024) %>%
  slice_max(abundance, n = 1) %>%
  pull(abundance)

recent_mean <- hybrid %>%
  filter(year == 2025) %>%
  summarise(mean_recent = mean(abundance, na.rm = TRUE)) %>%
  pull(mean_recent)

recent_mean / historic_baseline
recent_mean / prior_mean
pct_of_prior_max <- (100 * recent_mean / prior_max) - 100

#snow crab % decline during collapse
read.csv("./output/crab_abundance.csv") %>%
  filter(category == "population" & species == "snow") %>%
  filter(year %in% c(2018, 2021)) %>%
  summarise(
    a2019 = abundance[year == 2018],
    a2021 = abundance[year == 2021],
    pct_decline = (1 - a2021 / a2019) * 100) %>%
  pull(pct_decline)

#proportion of Chionoecetes population comprised of hybrids pre-MHW
abundance %>%
  filter(year < 2021) %>%
  group_by(year) %>%
  mutate(proportion = abundance / sum(abundance)) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE))

#proportion of Chionoecetes population comprised of hybrids in 2025
abundance %>%
  filter(year == 2025) %>%
  group_by(year) %>%
  mutate(proportion = abundance / sum(abundance)) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE))



