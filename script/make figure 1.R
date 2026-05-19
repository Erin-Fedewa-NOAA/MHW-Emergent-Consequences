#Goals ----
#Create Figure 1 for manuscript 

#Author: EJF

## Load packages
library(RColorBrewer)
library(tidyverse)
library(crabpack)
library(patchwork)

#--------------------------------#
#Abundance plots ----
#--------------------------------#

#data
abundance <- read.csv("./output/crab_abundance.csv") %>%
                filter(category == "population") %>%
                mutate(species = recode(species, "hybrid" = "Snow x Tanner hybrids",
                           "tanner" = "Tanner Crab", "snow" = "Snow Crab"))
  
#faceted plot 
abun_plot <- ggplot(abundance, aes(year, abundance, color = species)) +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci,
                  fill = species, group = species),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8, lineend = "round") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
  facet_wrap(~factor(species, levels = c("Snow Crab", "Tanner Crab", "Snow x Tanner hybrids")),
    scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Tanner Crab" = "#D55E00",
                                "Snow Crab"   = "#0072B2",
                                "Snow x Tanner hybrids" = "#009E73")) +
  scale_fill_manual(values = c("Tanner Crab" = "#D55E00",
                               "Snow Crab"   = "#0072B2",
                               "Snow x Tanner hybrids" = "#009E73")) +
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
  labs(y = "Spring Sea Ice Extent (%)", x = "") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
            alpha=0.1, fill= "grey60") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank())
    
#--------------------------------------#
#Proportion 4 inch males plot ----
#--------------------------------------#

#data
ratio <- read.csv("./output/proportion_hybrid.csv")

#plot
ratio_plot <- ratio %>%
  ggplot(aes(year, prop_hybrid)) +
  geom_line(linewidth = 0.8, lineend = "round", color = "#009E73") +
  theme_minimal(base_size = 11) +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
           alpha=0.1, fill= "grey60") +
  geom_ribbon(aes(ymin = prop_lower, ymax = prop_upper), alpha = 0.3, fill = "#009E73") +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  labs(y="Hybrid proportion (%)", x="") 

#--------------------------------#
#Combine and save figures ----
#--------------------------------#

fig1 <- ((ice_plot / ratio_plot) | abun_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11),
    plot.margin = margin(5, 5, 5, 5))

ggsave("./figures/fig1.jpeg", fig1,
  width = 17.8, height = 14, units = "cm", dpi = 600)

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
  filter(year >= 1988, year <= 2023) %>%
  summarise(mean_baseline = mean(abundance, na.rm = TRUE)) %>%
  pull(mean_baseline)

recent_mean <- hybrid %>%
  filter(year %in% c(2024, 2025)) %>%
  summarise(mean_recent = mean(abundance, na.rm = TRUE)) %>%
  pull(mean_recent)

fold_change <- recent_mean / historic_baseline

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
  group_by(year) %>%
  mutate(proportion = abundance / sum(abundance)) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE))



