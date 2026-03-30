#Test a hybrid model to examine causal linkages between the
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

#Note: because we're testing a max lag of 6 years, this workflow was originally tested 
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
           geom_histogram(bins = 30, fill = "skyblue", color = "black") +
           labs(x = .y, y = "Frequency") +
           theme_minimal())
  wrap_plots(plots)
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

#--------------------------------
# Test for DAG-data consistency
#--------------------------------



#-------------------------------------------
#Lag testing: sea ice -> tanner causal pathway
#-------------------------------------------

#Define SEMs for each sea ice -> tanner lag
sem_lag4 <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> tanner_abundance, 4, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 1, snowtotanner"

sem_lag5 <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> tanner_abundance, 5, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 1, snowtotanner"

sem_lag6 <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> tanner_abundance, 6, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 1, snowtotanner"

# -----------------------------
# Fit Lag 4
#-----------------------------
#build model without running it (needed so we can modify TMB inputs)
fit_build4 <- dsem(
  sem = sem_lag4,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

#Extract parameters & map from full model
pars4 <- fit_build4$tmb_inputs$parameters
map4  <- fit_build4$tmb_inputs$map

#Set reasonable starting values
n_vars <- ncol(data) #used to index parameters
pars4$delta0_j <- rep(0.01, n_vars) #small positive numbers improve numerical stability during optimization
pars4$lnsigma_j <- rep(log(0.5), n_vars) #observation error fixed at 0.5
map4$lnsigma_j <- factor(rep(NA, n_vars)) #tells TMB to fix observation error
n_beta4 <- length(pars4$beta_z) #total # of coefficients in the SEM (beta_z)
pars4$beta_z[(n_vars+1):n_beta4] <- 0.05  # small lag starting values

fit_lag4 <- dsem(
  sem = sem_lag4,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars4,
    map = map4,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Fit Lag 5
#-----------------------------
fit_build5 <- dsem(
  sem = sem_lag5,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars5 <- fit_build5$tmb_inputs$parameters
map5  <- fit_build5$tmb_inputs$map

#set starting parameters
pars5$delta0_j <- rep(0.01, n_vars)
pars5$lnsigma_j <- rep(log(0.5), n_vars)
map5$lnsigma_j <- factor(rep(NA, n_vars))
n_beta5 <- length(pars5$beta_z)
pars5$beta_z[(n_vars+1):n_beta5] <- 0.05

fit_lag5 <- dsem(
  sem = sem_lag5,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars5,
    map = map5,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Fit Lag 6
#-----------------------------
fit_build6 <- dsem(
  sem = sem_lag6,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars6 <- fit_build6$tmb_inputs$parameters
map6  <- fit_build6$tmb_inputs$map

pars6$delta0_j <- rep(0.01, n_vars)
pars6$lnsigma_j <- rep(log(0.5), n_vars)
map6$lnsigma_j <- factor(rep(NA, n_vars))
n_beta6 <- length(pars6$beta_z)
pars6$beta_z[(n_vars+1):n_beta6] <- 0.05

fit_lag6 <- dsem(
  sem = sem_lag6,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars6,
    map = map6,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Compare models
#-----------------------------
aic_values_icetotanner <- c(
  lag4 = AIC(fit_lag4),
  lag5 = AIC(fit_lag5),
  lag6 = AIC(fit_lag6))

loglik_values_icetotanner <- c(
  lag4 = logLik(fit_lag4),
  lag5 = logLik(fit_lag5),
  lag6 = logLik(fit_lag6))

aic_values_icetotanner
loglik_values_icetotanner

#----------------------------------
# Inspect causal pathways and plot
#----------------------------------
summary(fit_lag4)
summary(fit_lag5)
summary(fit_lag6)

#extract icetotanner estimates from each model
sm_lag4 <- as.data.frame(summary(fit_lag4))
sm_lag5 <- as.data.frame(summary(fit_lag5))
sm_lag6 <- as.data.frame(summary(fit_lag6))

est_lag4 <- sm_lag4 %>% filter(name == "icetotanner") %>% mutate(lag = "Lag 4")
est_lag5 <- sm_lag5 %>% filter(name == "icetotanner") %>% mutate(lag = "Lag 5")
est_lag6 <- sm_lag6 %>% filter(name == "icetotanner") %>% mutate(lag = "Lag 6")

# Combine into one data frame
est_all_icetotanner <- bind_rows(est_lag4, est_lag5, est_lag6)

# Plot effect sizes with error bars
ggplot(est_all_icetotanner, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 4" = "#1f78b4", 
                               "Lag 5" = "#33a02c", 
                               "Lag 6" = "#e31a1c")) +
  labs(title = "Sea Ice → Tanner Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#Lag 5 best supported by AIC

#############################################################
#-------------------------------------------
#Lag testing: sea ice -> snow causal pathway
#-------------------------------------------

#Define SEMs for each sea ice -> snow crab lag
sem_lag1 <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> tanner_abundance, 5, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 1, snowtotanner"

sem_lag2 <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

   # Causal pathways
  sea_ice -> tanner_abundance, 5, icetotanner
  sea_ice -> snow_abundance, 2, icetosnow
  snow_abundance -> tanner_abundance, 1, snowtotanner"

# -----------------------------
# Fit Lag 1
#-----------------------------
fit_build1 <- dsem(
  sem = sem_lag1,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars1 <- fit_build1$tmb_inputs$parameters
map1  <- fit_build1$tmb_inputs$map

# Safe starting values
n_vars <- ncol(data)
pars1$delta0_j <- rep(0.01, n_vars) #delta0
pars1$lnsigma_j <- rep(log(0.5), n_vars) #observation error
map1$lnsigma_j <- factor(rep(NA, n_vars))
n_beta1 <- length(pars1$beta_z)
pars1$beta_z[(n_vars+1):n_beta1] <- 0.05  # small lag starting values

fit_lag1 <- dsem(
  sem = sem_lag1,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars1,
    map = map1,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Fit Lag 2
#-----------------------------
fit_build2 <- dsem(
  sem = sem_lag2,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars2 <- fit_build2$tmb_inputs$parameters
map2  <- fit_build2$tmb_inputs$map

pars2$delta0_j <- rep(0.01, n_vars)
pars2$lnsigma_j <- rep(log(0.5), n_vars)
map2$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2 <- length(pars2$beta_z)
pars2$beta_z[(n_vars+1):n_beta2] <- 0.05

fit_lag2 <- dsem(
  sem = sem_lag2,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars2,
    map = map2,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Compare models
#-----------------------------
aic_values_icetosnow <- c(
  lag1 = AIC(fit_lag1),
  lag2 = AIC(fit_lag2))

loglik_values_icetosnow <- c(
  lag1 = logLik(fit_lag1),
  lag2 = logLik(fit_lag2))

aic_values_icetosnow
loglik_values_icetosnow

#----------------------------------
# Inspect causal pathways and plot
#----------------------------------
summary(fit_lag1)
summary(fit_lag2)

#extract icetosnow estimates from each model
sm_lag1 <- as.data.frame(summary(fit_lag1))
sm_lag2 <- as.data.frame(summary(fit_lag2))

est_lag1 <- sm_lag1 %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 1")
est_lag2 <- sm_lag2 %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 2")

# Combine into one data frame
est_all_icetosnow <- bind_rows(est_lag1, est_lag2)

# Plot effect sizes with error bars
ggplot(est_all_icetosnow, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 1" = "#1f78b4", 
                               "Lag 2" = "#33a02c")) +
  labs(title = "Sea Ice → Snow Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#lag 1 best supported by AIC

#############################################################
#-------------------------------------------
#Lag testing: snow -> tanner causal pathway
#-------------------------------------------

#Define SEMs for each snow to tanner (stt) lag being tested
sem_lag1_stt <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> tanner_abundance, 5, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 1, snowtotanner"

sem_lag2_stt <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

   # Causal pathways
  sea_ice -> tanner_abundance, 5, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 2, snowtotanner"

sem_lag3_stt <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  tanner_abundance -> tanner_abundance, 1, ar_tanner

   # Causal pathways
  sea_ice -> tanner_abundance, 5, icetotanner
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> tanner_abundance, 3, snowtotanner"

# -----------------------------
# Fit Lag 1
#-----------------------------
fit_build1_stt <- dsem(
  sem = sem_lag1_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars1_stt <- fit_build1_stt$tmb_inputs$parameters
map1_stt  <- fit_build1_stt$tmb_inputs$map

# Safe starting values
n_vars <- ncol(data)
pars1_stt$delta0_j <- rep(0.01, n_vars) #delta0
pars1_stt$lnsigma_j <- rep(log(0.5), n_vars) #observation error
map1_stt$lnsigma_j <- factor(rep(NA, n_vars))
n_beta1_stt <- length(pars1_stt$beta_z)
pars1_stt$beta_z[(n_vars+1):n_beta1_stt] <- 0.05  # small lag starting values

fit_lag1_stt <- dsem(
  sem = sem_lag1_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars1_stt,
    map = map1_stt,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Fit Lag 2
#-----------------------------
fit_build2_stt <- dsem(
  sem = sem_lag2_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars2_stt <- fit_build2_stt$tmb_inputs$parameters
map2_stt  <- fit_build2_stt$tmb_inputs$map

pars2_stt$delta0_j <- rep(0.01, n_vars)
pars2_stt$lnsigma_j <- rep(log(0.5), n_vars)
map2_stt$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2 <- length(pars2_stt$beta_z)
pars2_stt$beta_z[(n_vars+1):n_beta2] <- 0.05

fit_lag2_stt <- dsem(
  sem = sem_lag2_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars2_stt,
    map = map2_stt,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Fit Lag 3
#-----------------------------
fit_build3_stt <- dsem(
  sem = sem_lag3_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars3_stt <- fit_build3_stt$tmb_inputs$parameters
map3_stt  <- fit_build3_stt$tmb_inputs$map

pars3_stt$delta0_j <- rep(0.01, n_vars)
pars3_stt$lnsigma_j <- rep(log(0.5), n_vars)
map3_stt$lnsigma_j <- factor(rep(NA, n_vars))
n_beta3 <- length(pars3_stt$beta_z)
pars3_stt$beta_z[(n_vars+1):n_beta3] <- 0.05

fit_lag3_stt <- dsem(
  sem = sem_lag3_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars3_stt,
    map = map3_stt,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------
# Compare models
#-----------------------------
aic_values_snowtotanner <- c(
  lag1 = AIC(fit_lag1_stt),
  lag2 = AIC(fit_lag2_stt),
  lag3 = AIC(fit_lag3_stt))

loglik_values_snowtotanner <- c(
  lag1 = logLik(fit_lag1_stt),
  lag2 = logLik(fit_lag2_stt),
  lag3 = logLik(fit_lag3_stt))

aic_values_snowtotanner
loglik_values_snowtotanner

#----------------------------------
# Inspect causal pathways and plot
#----------------------------------
summary(fit_lag1_stt)
summary(fit_lag2_stt)
summary(fit_lag3_stt)

#extract snow to tanner estimates from each model
sm_lag1_stt <- as.data.frame(summary(fit_lag1_stt))
sm_lag2_stt <- as.data.frame(summary(fit_lag2_stt))
sm_lag3_stt <- as.data.frame(summary(fit_lag3_stt))

est_lag1_stt <- sm_lag1_stt %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 1")
est_lag2_stt <- sm_lag2_stt %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 2")
est_lag3_stt <- sm_lag3_stt %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 3")

# Combine into one data frame
est_all_snowtotanner <- bind_rows(est_lag1_stt, est_lag2_stt, est_lag3_stt)

# Plot effect sizes with error bars
ggplot(est_all_snowtotanner, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 1" = "#1f78b4", 
                               "Lag 2" = "#33a02c",
                               "Lag 3" = "#e31a1c")) +
  labs(title = "Snow Abundance → Tanner Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#lag 3 best supported by AIC

#----------------------------------------------------------------
# Produce AIC table and combined plots for all 3 causal pathways 
#----------------------------------------------------------------

#create dataframe with all AIC scores
aic_long <- bind_rows(
  tibble(Pathway = "Sea Ice → Tanner", Lag = c(4,5,6), AIC = aic_values_icetotanner),
  tibble(Pathway = "Sea Ice → Snow",   Lag = c(1,2),   AIC = aic_values_icetosnow),
  tibble(Pathway = "Snow → Tanner",    Lag = c(1,2,3), AIC = aic_values_snowtotanner))

# Compute delta AIC
aic_table <- aic_long %>%
  group_by(Pathway) %>%
  mutate(Delta_AIC = AIC - min(AIC)) %>%
  arrange(Pathway, AIC) %>%
  ungroup() 

#Print table
aic_table %>%
  kable("html", digits = 2,
        col.names = c("Causal Pathway", "Best Lag", "AIC", "ΔAIC"),
        caption = "Best-supported lag for each causal pathway") %>%
  kable_styling(full_width = FALSE, position = "center",
                bootstrap_options = c("striped", "hover")) %>%
  row_spec(2, extra_css = "border-bottom: 2px solid black;") %>%
  row_spec(5, extra_css = "border-bottom: 2px solid black;")

# Combine effect sizes
est_all <- bind_rows(
  est_all_icetotanner %>% mutate(Causal_Pathway = "Sea Ice → Tanner"),
  est_all_icetosnow %>% mutate(Causal_Pathway = "Sea Ice → Snow"),
  est_all_snowtotanner %>% mutate(Causal_Pathway = "Snow → Tanner"))

# Color palette (colorblind-friendly)
lag_colors <- c(
  "Lag 1" = "#E69F00",
  "Lag 2" = "#56B4E9",
  "Lag 3" = "#009E73",
  "Lag 4" = "#F0E442",
  "Lag 5" = "#0072B2",
  "Lag 6" = "#D55E00")

# Plot
ggplot(est_all, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = lag_colors) +
  facet_wrap(~Causal_Pathway, scales = "free_y") +
  labs(title = "Tanner Crab Candidate Lags",
       x = "",
       y = "Effect Size (± Std Error)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10))

##########################################################################

#----------------------------------------------------------------
# Fit final Tanner crab model 
#----------------------------------------------------------------

#define SEM based on best lag structure from above
sem_final <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  snow_abundance -> snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_tanner
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

  # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> snow_abundance, 1, icetosnow
  snow_abundance -> log_hybrid_abundance, 3, snowtohybrid
  sea_ice -> bhatta_snowfem_tanmale, 1, icetooverlap
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 4, overlaptohybrid"

#two-stage modeling fitting procedure: 
##initial first model run without delta0 (to improve starting values)
fit0 <- dsem(sem = sem_final, tsdata = data,
             family = family, estimate_delta0 = FALSE,
             control = dsem_control(quiet = FALSE, getsd = FALSE))

# extract starting parameters
parameters <- fit0$obj$env$parList()

# add starting values for delta0
parameters$delta0_j <- rep(0, ncol(data))

#build model without running it
#(needed so we can modify TMB inputs)
fit_build_final <- dsem(sem = sem_final, tsdata = data,
                        family = family,
                        estimate_delta0 = TRUE,
                        control = dsem_control(
                          run_model = FALSE,
                          parameters = parameters))

pars_final <- fit_build_final$tmb_inputs$parameters
map_final  <- fit_build_final$tmb_inputs$map

# fix observation SD = 0.1
pars_final$lnsigma_j <- rep(log(0.1), ncol(data))

# prevent estimation of observation SD
map_final$lnsigma_j <- factor(rep(NA, ncol(data)))

#run final model fit with Delta0 and fixed observation error
#this final model estimates AR1 coefficients, lagged effects, 
#process variances and delta0 while keeping observation SD fixed

fit_dsem <- dsem(sem=sem_final, tsdata=data, 
                 family=family,
                 estimate_delta0=TRUE,
                 control=dsem_control(parameters = pars_final,
                                      map = map_final,
                                      quiet = TRUE,
                                      getsd = TRUE))

summary(fit_dsem)

#----------------------------------------
# Final Tanner crab model diagnostics
#----------------------------------------
#Convergence diagnostics:

#Hessian/SE - should be no hessian warnings or NA SE estimates
fit_dsem$sdrep

#check gradient
max(abs(fit_dsem$opt$gradient))

#Identifiablity/overfitting
fit_dsem$sdrep$pdHess # TRUE = good here

#Residual diagnostics:

resid <- residuals(fit_dsem)

#plot fitted vrs residuals
plot(fitted(fit_dsem), resid)
abline(h = 0, col = "red")

#autocorrelation of residuals
acf(resid[, "tanner_abundance"])
acf(resid[, "snow_abundance"])
acf(resid[, "sea_ice"])
acf(residuals(fit_dsem))


sim <- simulate(fit_dsem)
plot(data[, "tanner_abundance"], type = "l")
lines(sim[, "tanner_abundance"], col = "blue")

# sample-based quantile residuals
samples = loo_residuals(fit_dsem, what="samples", track_progress=FALSE)
which_use = which(!is.na(tanner_data))
fitResp = loo_residuals( fit_dsem, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]

# Build and display DHARMa object
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(tanner_data)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)
#looks good

# Calculate total effects
effect = total_effect(fit_dsem, n_lags = 6)
#specifying n+1 longest lag so sparse matrix doesn't have NAs

# Plot total effect
ggplot( effect) + 
  geom_bar(aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )
#if we test an iid model structure above (ie switch 1's in AR1 terms to 0) there shouldn't
#be n-1 lags, only the specified lag
#note that for response variables with direct and indirect causal links, we see total
#effects for each specified lag/causal path

#relative importance of variables as predictors of female size structure
partition_variance(fit_dsem,
                   which_response = "hybrid_abundance",
                   n_times = 10 )
#hmm maybe this is in development version of dsem only? 

#-----------------------------------------------------------
#Plot final Tanner crab model fit
#-----------------------------------------------------------

#build function to plot final model fit
#function revised from C. Monnahan plot_fit function
plot_fit <- function(Y, fit, start_year = NULL){
  
  # Extract latent states (model estimates of underlying system state)
  ParHat <- fit$obj$env$parList()
  pred <- ParHat$x_tj
  
  # Extract standard errors around latent states
  sd_list <- as.list(fit$sdrep, what = "Std.")
  SD <- sd_list$x_tj
  
  #this is needed to prevent NA's from being offset due to missing data
  if (is.null(start_year)) {
    if (!is.null(rownames(Y))) {
      years <- as.numeric(rownames(Y))
    } else {
      years <- seq_len(nrow(Y))
    }
  } else {
    years <- start_year + seq_len(nrow(Y)) - 1
  }
  
  #build a dataset for each causal pathway
  out <- lapply(seq_len(ncol(Y)), function(i) {
    tmp <- data.frame(
      year = years,
      variable = colnames(Y)[i],
      obs = as.numeric(Y[, i]),
      pred = as.numeric(pred[, i]),
      sd = as.numeric(SD[, i]))
    
    tmp %>%
      mutate(lower = pred - ifelse(is.na(sd), 0, sd),
             upper = pred + ifelse(is.na(sd), 0, sd))
    
    #combine all variables
  }) %>% bind_rows()
  
  out$variable <- factor(out$variable, levels = colnames(Y))
  
  #plot output
  ggplot(out, aes(x = year)) +
    #uncertainty around latent state
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "blue", alpha = 0.2) +
    #latent state (model estimate)
    geom_line(aes(y = pred), color = "blue", linewidth = 1) +
    #observed data
    geom_point(aes(y = obs), color = "red", size = 2, na.rm = TRUE) +
    geom_line(aes(y = obs), color = "red", linetype = "dashed") +
    facet_wrap(~variable, scales = "free_y", ncol = 1) +
    labs(x = "", y = "Value",
         title = "Hybrid model DSEM fit: observed vs estimated latent state") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold"))
}

#plot
plot_fit(data, fit_dsem, start_year = start_year)

#and now plot fitted DAG
plot(as_fitted_DAG(fit_dsem, what = "path_coefficient"))
plot(as_fitted_DAG(fit_dsem_auto, what = "p_value"), type = "width")
#need to follow up on this- dsem output doesn't retain path names that 
#as_fitted_DAG() function needs since we manually specified map and par

#-----------------------------------------------------------
# Compute % variance explained for each causal pathway
#-----------------------------------------------------------




#-----------------------------------------------------------
# Test simulations
#-----------------------------------------------------------






#-----------------------------------------------------------
# Evaluate sensitivity to changes in fixed observation error
#-----------------------------------------------------------

#test plausible range of fixed values
obs_sd_values <- c(0.1, 0.2, 0.3, 0.4)

results_list <- list()

#loop through values
for (i in seq_along(obs_sd_values)) {
  
  sd_val <- obs_sd_values[i]
  
  # update observation error
  pars_tmp <- pars_final
  map_tmp  <- map_final
  
  pars_tmp$lnsigma_j <- rep(log(sd_val), ncol(data))
  map_tmp$lnsigma_j  <- factor(rep(NA, ncol(data)))
  
  #fit model with each observation error value
  fit_tmp <- dsem(
    sem = sem_final,
    tsdata = data,
    family = family,
    estimate_delta0 = TRUE,
    control = dsem_control(
      parameters = pars_tmp,
      map = map_tmp,
      quiet = TRUE,
      getsd = TRUE))
  
  sm <- as.data.frame(summary(fit_tmp))
  
  # extract key pathways
  effects <- sm %>%
    filter(name %in% c("icetotanner", "icetosnow", "snowtotanner")) %>%
    mutate(obs_sd = sd_val)
  
  results_list[[i]] <- effects
}

sensitivity_results <- bind_rows(results_list)

#plot sensitivity
ggplot(sensitivity_results, aes(x = obs_sd, y = Estimate, color = name)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - Std_Error,
                    ymax = Estimate + Std_Error),
                width = 0.02) +
  labs(x = "Fixed Observation SD",
       y = "Effect Size",
       color = "Path") +
  theme_minimal()

#The sea ice -> snow crab linkage seems most sensitive to changes in fixed 
#observation error- a positive effect exists at all tested values, but 
#the effect is weaker and non-significant when accounting for higher uncertainty (>0.2)
#Other two pathways are robust to changes in observation SD




































































