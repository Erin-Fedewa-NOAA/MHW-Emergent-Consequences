#Goals ----
#Create Figure 1 for manuscript 

#Author: EJF

## Load packages
library(RColorBrewer)
library(tidyverse)
library(crabpack)

#--------------------------------#
#Abundance plots ----
#--------------------------------#

#data
abundance <- read.csv("./output/crab_abundance.csv") %>%
                filter(category == "population") %>%
                mutate(species = recode(species, "hybrid" = "Snow x Tanner hybrids",
                           "tanner" = "Tanner Crab", "snow" = "Snow Crab"))
  
#plot with shared y axis
abundance %>%
  ggplot(aes(year, abundance, color = species)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci,
                  fill = species), alpha = 0.2, color = NA, na.rm = TRUE) +
  #scale_color_manual(values = c("#F09E70", "#A3C5ED", "#B3E5A0")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(y="Abundance (millions)", x="") +
  #theme(legend.position=c(.2,.13)) +
  theme(legend.background = element_rect(color="transparent", fill="transparent")) +
  theme(axis.title.y = element_text(size=10)) +
  theme(legend.title=element_blank()) 

#faceted plot 
ggplot(abundance, aes(year, abundance, color = species)) +
  geom_ribbon(aes(ymin = abundance - abundance_ci,
                  ymax = abundance + abundance_ci,
                  fill = species, group = species),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~factor(species, levels = c("Snow Crab", "Tanner Crab", "Snow x Tanner hybrids")),
    scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Tanner Crab" = "#D55E00",
                                "Snow Crab"   = "#0072B2",
                                "Snow x Tanner hybrids" = "#009E73")) +
  scale_fill_manual(values = c("Tanner Crab" = "#D55E00",
                               "Snow Crab"   = "#0072B2",
                               "Snow x Tanner hybrids" = "#009E73")) +
  labs(y = "Abundance (millions)", x = "") +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing = unit(1.2, "lines"),
    strip.text = element_text(face = "plain", color = "grey40", size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.border = element_blank())

#--------------------------------------#
#Sea Ice plot ----
#--------------------------------------#

#data
ice <- read.csv("./output/seaice_output.csv") %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
  filter(year >= 1988)

#plot
ggplot(ice, aes(year, Mar_Apr_ice_EBS_NBS)) +
  geom_line(linewidth = 0.7, color="#0072B2") +
  labs(y = "Spring Sea Ice Extent", x = "") +
  annotate("rect", xmin= 2017.5, xmax=2019.5 ,ymin=-Inf , ymax=Inf, 
            alpha=0.1, fill= "#FF9B9B") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank())
    
#--------------------------------------#
#Proportion 4 inch males plot ----
#--------------------------------------#

#data
ratio <- read.csv("./output/proportion_legal.csv")

#plot
ratio %>%
  ggplot(aes(year, prop, color = stock, fill = stock)) +
  geom_line() +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  labs(y="Proportion of legal males", x="") 

#--------------------------------#
#Lag effect figures ----
#--------------------------------#

#Data
tanner_lags <- read.csv("./output/tanner_lags.csv")
hybrid_lags <- read.csv("./output/hybrid_lags.csv")

# Tanner Plot
ggplot(tanner_lags, aes(x = Causal_Pathway, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
            width = 0.65, color = "black",linewidth = 0.3) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error),
                      position = position_dodge(width = 0.65),
                      width = 0.15, linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black",linewidth = 0.3) +
  scale_fill_manual(values = c("#F7FBFF","#DEEBF7","#C6DBEF",
                               "#9ECAE1","#6BAED6","#3182BD","#08519C")) +
  labs(x = "", y = "Causal estimate (± SE)", fill = "Lag") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9))

#hybrid plot
ggplot(hybrid_lags, aes(x = Causal_Pathway, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.65, color = "black",linewidth = 0.3) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error),
                position = position_dodge(width = 0.65),
                width = 0.15, linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black",linewidth = 0.3) +
  scale_fill_manual(values = c("#F7FBFF","#DEEBF7","#C6DBEF",
                               "#9ECAE1","#6BAED6","#3182BD","#08519C")) +
  labs(x = "", y = "Causal estimate (± SE)", fill = "Lag") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9))

#option 2 faceted plot, tanner
tanner_lags_plot <- tanner_lags %>%
  mutate(Causal_Pathway = factor(Causal_Pathway,
                                 levels = unique(Causal_Pathway)))

ggplot(tanner_lags_plot, aes(x = lag, y = Estimate, fill = lag)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Estimate - Std_Error,
                    ymax = Estimate + Std_Error),
                width = 0.15, size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_fill_manual(values = lag_colors, name = "Lag") +
  facet_wrap(~Causal_Pathway, nrow = 1) +  # single row = cleaner comparison
  labs(x = "", y = "Effect size (± SE)") +
  theme_bw(base_size = 12) +
  theme(legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(face = "bold", size = 14, hjust = 0))

#option 2 faceted plot, hybrids
hybrid_lags_plot <- hybrid_lags %>%
  mutate(Causal_Pathway = factor(Causal_Pathway,
                                 levels = unique(Causal_Pathway)))

ggplot(hybrid_lags_plot, aes(x = lag, y = Estimate, fill = lag)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Estimate - Std_Error,
                    ymax = Estimate + Std_Error),
                width = 0.15, size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_fill_manual(values = lag_colors, name = "Lag") +
  facet_wrap(~Causal_Pathway, nrow = 1) +  # single row = cleaner comparison
  labs(x = "", y = "Effect size (± SE)") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(face = "bold", size = 14, hjust = 0))

#--------------------------------#
#Combine and save figures ----
#--------------------------------#

#follow up once figures are finalized!
ggsave("figure1.tiff", width = 3.5, height = 3, dpi = 600, compression = "lzw")

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



