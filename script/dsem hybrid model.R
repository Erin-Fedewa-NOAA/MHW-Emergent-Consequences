#Goals ----
#Fit a dsem model to examine causal linkages between the
#MHW and hybrid abundance increase

#Specifically, we'll test the following causal pathways:
#sea ice -> snow crab population abundance 
  #Mechanism testing: sea ice loss/mass mortality events 
  #Using the 1 yr lag tested in "dsem tanner model.R"
#sea ice -> hybrid population abundance 
  #Mechanism testing: sea ice effects on juvenile survival/recruitment
  #Hypothesized lags to test: 4-6 years
#snow crab population abundance -> hybrid population abundance
  #Mechanism: snow crab collapse results in competitive release/niche & habitat expansion
  #Hypothesized lags to test: 1-3 years 
#sea ice -> snow/tanner spatial overlap
  #Mechanism: temperature effects on spatial distributions
  #Hypothesized lags to test: 0-1 years
#snow/tanner spatial overlap -> hybrid abundance
  #Mechanism: increased interspecies mating opportunities
  #Hypothesized lags to test: 5-7 years

#Author: EJF

#load 
library(tidyverse)
library(dsem)
library(ggplot2)
library(dplyr)
library(dagitty)
library(ggdag)
library(knitr)
library(corrplot)
library(patchwork)
library(TMB)
library(knitr)
library(kableExtra)
library(phylopath)

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")
overlap <- read.csv("./output/overlap_output.csv")

#----------------------------------
# Define covariates and standardize
#----------------------------------
start_year = 1988

#Note: because we're testing a max lag of 7 years, this workflow was originally tested 
#with a longer sea ice timeseries (1982+). We report sensitivities to this model specification
#in the paper, noting that the shorter sea ice timeseries was selected in order to 
#optimize latent state estimation of crab abundance timeseries 

dat_hybrid <- sea_ice %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
  filter(year >= start_year) %>% 
  rename(sea_ice = Mar_Apr_ice_EBS_NBS) %>%
  full_join(overlap %>%
              select(-bhatta_global)) %>%
  full_join(crab_abund %>%
              filter(category == "population_50mm_plus" & species != "tanner") %>%
              select(year, abundance, species) %>%
              pivot_wider(names_from = "species", values_from = "abundance") %>%
              rename(snow_abundance=snow, hybrid_abundance=hybrid)) 
#plot 
dat_hybrid %>%
  pivot_longer(-year, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#Distributions of covariates
plot_histo <- function(data) {
  plots <- data %>%
    imap(~ ggplot(data, aes(x = .data[[.y]])) +
           geom_histogram(aes(y = after_stat(density)),
                          bins = 30,
                          fill = "skyblue", color = "black") +
           geom_density(color = "red") +
           labs(x = .y, y = "Density") +
           theme_minimal())
  patchwork::wrap_plots(plots)
}

plot_histo(dat_hybrid %>% select(-year))

#Normalize and standardize all variables 
vars <- c("sea_ice", "bhatta_snowfem_tanmale",
          "snow_abundance",
          "hybrid_abundance", "log_hybrid_abundance", "bhatta_tannerfem_snowmale")

#And center covariates
hybrid_data <- dat_hybrid %>%
  mutate(log_hybrid_abundance = log(hybrid_abundance)) %>%
  mutate(across(all_of(c(vars)), ~ as.numeric(scale(.)))) %>%
  select(-hybrid_abundance)

#Assess collinearity b/w variables
hybrid_data %>% 
  select(-year) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(method="number")

#plot all standardized variables
hybrid_data %>%
  pivot_longer(cols = -year, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#-----------------------------
#Prep data for dsem models
#-----------------------------

data <- hybrid_data %>%
  select(-year) %>%
  ts()

family <- rep("normal", ncol(data))










