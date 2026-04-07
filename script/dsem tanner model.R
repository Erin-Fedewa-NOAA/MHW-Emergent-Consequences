#Goals -----
#Fit a dsem model to examine causal linkages between the
#MHW and tanner abundance increase

#Specifically, we'll test the following causal pathways:
#sea ice -> snow crab population abundance 
  #Mechanism testing: sea ice loss/mass mortality events 
  #Hypothesized lags to test: 1-2 years
#sea ice -> tanner crab population abundance 
  #Mechanism testing: sea ice effects on survival/productivity
  #Hypothesized lags to test: 1-2 years
#snow crab population abundance -> tanner crab population abundance
  #Mechanism: snow crab collapse results in competitive release/niche & habitat expansion
  #Hypothesized lags to test: 1-3 years 

#Approach: because each dsem model fit produces new parameters (dynamic mapping), if
  #you don't reset starting values each iteration, optimizers reuse previous parameter 
  #values, which causes NA or singular matrix errors when looping through lags and 
  #refitting models. It's long & clunky, but instead of a loop we'll 1) evaluate lags for each 
  #causal link first by iteratively refitting models and comparing with AIC, 2) do the same 
  #for each causal link, and 3) use best lag structure to fit a final Tanner crab model  


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
library(broom)

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")

#-----------------------------------------#
# Define covariates and standardize ---- 
#-----------------------------------------#
start_year = 1988

#Note: because we're testing a max lag of 6 years, this workflow was originally tested 
  #with a longer sea ice timeseries (1982+). We report sensitivities to this model,
  #specification, noting that the shorter sea ice timeseries was selected in order to 
  #optimize latent state estimation of crab abundance timeseries 

dat_tanner <- sea_ice %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
  filter(year >= start_year) %>% 
  rename(sea_ice = Mar_Apr_ice_EBS_NBS) %>%
  full_join(crab_abund %>%
              filter(category == "population_50mm_plus" & species != "hybrid") %>%
              select(year, abundance, species) %>%
              pivot_wider(names_from = "species", values_from = "abundance") %>%
              rename(snow_abundance=snow, tanner_abundance=tanner)) 
#plot 
dat_tanner %>%
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

plot_histo(dat_tanner %>% select(-year))

#Scale all variables and transform abundance variables 
vars <- c("sea_ice",
          "snow_abundance", "log_snow_abundance",
          "tanner_abundance", "log_tanner_abundance")

tanner_data <- dat_tanner %>%
  mutate(log_tanner_abundance = log(tanner_abundance),
         log_snow_abundance = log(snow_abundance)) %>%
  mutate(across(all_of(c(vars)), ~ as.numeric(scale(.)))) %>%
  select(-tanner_abundance, -snow_abundance)

#check distributions now
plot_histo(tanner_data %>% select(-year))

#Assess collinearity b/w variables
tanner_data %>% 
  select(-year) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(method="number")

#plot all standardized/transformed variables
tanner_data %>%
  pivot_longer(cols = -year, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#-------------------------------#
#Prep data for dsem models ----
#-------------------------------#

data <- tanner_data %>%
  select(-year) %>%
  ts()

family <- rep("normal", ncol(data))

#---------------------------------------#
# Test for DAG-data consistency ----
#---------------------------------------#

#download specified DAG from dagitty.net - conditional independencies are 
#identified based on structure of DAG drawn in dagitty
tanner_dag <- dagitty('dag {
"Sea Ice" [exposure,pos="-1.800,0.078"]
"Snow crab abundance" [pos="-0.487,-0.799"]
"Tanner crab abundance" [outcome,pos="-0.454,0.832"]
"Sea Ice" -> "Snow crab abundance"
"Sea Ice" -> "Tanner crab abundance"
"Snow crab abundance" -> "Tanner crab abundance"
}
')

#plot DAG
ggdag(tanner_dag, layout = "nicely") +
  theme_dag()

plot(tanner_dag) 

ggdag_status(tanner_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(tanner_dag)
#2 open causal pathways 

#and plot paths
ggdag_paths(tanner_dag, text = FALSE, use_labels = "name") +
  theme_dag()

#find adjustment sets for response variable 
adjustmentSets(tanner_dag, exposure="Sea Ice", outcome="Tanner crab abundance")
#This tells us that no covariate adjustment is necessary to identify the causal effect
  #ie there are no open backdoor paths 

#and visualize adjustment sets, if there are any
ggdag_adjustment_set(tanner_dag, shadow = TRUE) +
  theme_dag()

#find conditional independencies- i.e. two variables that are implied to 
# be independent and not correlated shouldn't be connected by a node
impliedConditionalIndependencies(tanner_dag)
#our DAG doesn't imply any conditional independence relationships, and no
  #d-separation claims exist 

#----------------------------------------------------#
#Lag testing: sea ice -> tanner causal pathway ----
#----------------------------------------------------#

#Define SEMs for each sea ice -> tanner lag
sem_lag1_tanner <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 1, snowtotanner"

sem_lag2_tanner <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> log_tanner_abundance, 2, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 1, snowtotanner"

# ----------------------------#
# Fit Lag 1 ----
#-----------------------------#
#build model without running it (needed so we can modify TMB inputs)
fit_build1_tanner <- dsem(
  sem = sem_lag1_tanner,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

#Hessian warnings are expected from fit_build b/c parameters are at
  #default/initial value

#Extract parameters & map from full model
pars1_tanner <- fit_build1_tanner$tmb_inputs$parameters
map1_tanner  <- fit_build1_tanner$tmb_inputs$map

#Set reasonable starting values
n_vars <- ncol(data) #used to index parameters
pars1_tanner$delta0_j <- rep(0.01, n_vars) #small positive numbers improve numerical stability during optimization
pars1_tanner$lnsigma_j <- rep(log(0.5), n_vars) #observation error fixed at 0.5
map1_tanner$lnsigma_j <- factor(rep(NA, n_vars)) #tells TMB to fix observation error
n_beta1_tanner <- length(pars1_tanner$beta_z) #total # of coefficients in the SEM (beta_z)
pars1_tanner$beta_z[(n_vars+1):n_beta1_tanner] <- 0.05  # small lag starting values

fit_lag1_tanner <- dsem(
  sem = sem_lag1_tanner,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars1_tanner,
    map = map1_tanner,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 2 ----
#-----------------------------#
fit_build2_tanner <- dsem(
  sem = sem_lag2_tanner,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars2_tanner <- fit_build2_tanner$tmb_inputs$parameters
map2_tanner  <- fit_build2_tanner$tmb_inputs$map

#set starting parameters
pars2_tanner$delta0_j <- rep(0.01, n_vars)
pars2_tanner$lnsigma_j <- rep(log(0.5), n_vars)
map2_tanner$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2_tanner <- length(pars2_tanner$beta_z)
pars2_tanner$beta_z[(n_vars+1):n_beta2_tanner] <- 0.05

fit_lag2_tanner <- dsem(
  sem = sem_lag2_tanner,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars2_tanner,
    map = map2_tanner,
    quiet = TRUE,
    getsd = TRUE))

#---------------------------------#
# Compare ice/tanner models ----
#---------------------------------#
aic_values_icetotanner <- c(
  lag1 = AIC(fit_lag1_tanner),
  lag2 = AIC(fit_lag2_tanner))

loglik_values_icetotanner <- c(
  lag1 = logLik(fit_lag1_tanner),
  lag2 = logLik(fit_lag2_tanner))

aic_values_icetotanner
loglik_values_icetotanner

summary(fit_lag1_tanner)
summary(fit_lag2_tanner)

#extract icetotanner estimates from each model
sm_lag1_tanner <- as.data.frame(summary(fit_lag1_tanner))
sm_lag2_tanner <- as.data.frame(summary(fit_lag2_tanner))

est_lag1_tanner <- sm_lag1_tanner %>% filter(name == "icetotanner") %>% mutate(lag = "Lag 1")
est_lag2_tanner <- sm_lag2_tanner %>% filter(name == "icetotanner") %>% mutate(lag = "Lag 2")

# Combine into one data frame
est_all_icetotanner <- bind_rows(est_lag1_tanner, est_lag2_tanner)

# Plot effect sizes with error bars
ggplot(est_all_icetotanner, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 1" = "#1f78b4", 
                               "Lag 2" = "#33a02c")) +
  labs(title = "Sea Ice → Tanner Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#Can't differentiate between models using AIC, and effect is small and insignificant 
  #using either lag. Let's stick with lag 1 for consistency with snow crab mechanism

#------------------------------------------------#
#Lag testing: ice -> snow causal pathway ----
#------------------------------------------------#

#Define SEMs for each sea ice -> snow crab lag
sem_lag1_snow <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 1, snowtotanner"

sem_lag2_snow <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

   # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 2, icetosnow
  log_snow_abundance -> log_tanner_abundance, 1, snowtotanner"

#-----------------------------#
# Fit Lag 1 ----
#-----------------------------#
fit_build1_snow <- dsem(
  sem = sem_lag1_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars1_snow <- fit_build1_snow$tmb_inputs$parameters
map1_snow  <- fit_build1_snow$tmb_inputs$map

# Safe starting values
n_vars <- ncol(data)
pars1_snow$delta0_j <- rep(0.01, n_vars) #delta0
pars1_snow$lnsigma_j <- rep(log(0.5), n_vars) #observation error
map1_snow$lnsigma_j <- factor(rep(NA, n_vars))
n_beta1_snow <- length(pars1_snow$beta_z)
pars1_snow$beta_z[(n_vars+1):n_beta1_snow] <- 0.05  # small lag starting values

fit_lag1_snow <- dsem(
  sem = sem_lag1_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars1_snow,
    map = map1_snow,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 2 ----
#-----------------------------#
fit_build2_snow <- dsem(
  sem = sem_lag2_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars2_snow <- fit_build2_snow$tmb_inputs$parameters
map2_snow  <- fit_build2_snow$tmb_inputs$map

pars2_snow$delta0_j <- rep(0.01, n_vars)
pars2_snow$lnsigma_j <- rep(log(0.5), n_vars)
map2_snow$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2_snow <- length(pars2_snow$beta_z)
pars2_snow$beta_z[(n_vars+1):n_beta2_snow] <- 0.05

fit_lag2_snow <- dsem(
  sem = sem_lag2_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(
    parameters = pars2_snow,
    map = map2_snow,
    quiet = TRUE,
    getsd = TRUE))

#---------------------------------#
# Compare ice/snow models ----
#---------------------------------#
aic_values_icetosnow <- c(
  lag1 = AIC(fit_lag1_snow),
  lag2 = AIC(fit_lag2_snow))

loglik_values_icetosnow <- c(
  lag1 = logLik(fit_lag1_snow),
  lag2 = logLik(fit_lag2_snow))

aic_values_icetosnow
loglik_values_icetosnow

summary(fit_lag1_snow)
summary(fit_lag2_snow)

#extract icetosnow estimates from each model
sm_lag1_snow <- as.data.frame(summary(fit_lag1_snow))
sm_lag2_snow <- as.data.frame(summary(fit_lag2_snow))

est_lag1_snow <- sm_lag1_snow %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 1")
est_lag2_snow <- sm_lag2_snow %>% filter(name == "icetosnow") %>% mutate(lag = "Lag 2")

# Combine into one data frame
est_all_icetosnow <- bind_rows(est_lag1_snow, est_lag2_snow)

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

#Again, lags are nearly identical, suggesting that lag 1 is effectively
  #capturing cumulative/lagged effects at longer lag 2

#--------------------------------------------------#
#Lag testing: snow -> tanner causal pathway ----
#--------------------------------------------------#

#Define SEMs for each snow to tanner (stt) lag being tested
sem_lag1_stt <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

  # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 1, snowtotanner"

sem_lag2_stt <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

   # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 2, snowtotanner"

sem_lag3_stt <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner

   # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 3, snowtotanner"

# -----------------------------#
# Fit Lag 1 ----
#-----------------------------#
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
pars1_stt$delta0_j <- rep(0, n_vars) #delta0
pars1_stt$lnsigma_j <- rep(log(0.4), n_vars) #slightly smaller than above to improve stability
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

#-----------------------------#
# Fit Lag 2 ----
#-----------------------------#
fit_build2_stt <- dsem(
  sem = sem_lag2_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars2_stt <- fit_build2_stt$tmb_inputs$parameters
map2_stt  <- fit_build2_stt$tmb_inputs$map

pars2_stt$delta0_j <- rep(0, n_vars)
pars2_stt$lnsigma_j <- rep(log(0.4), n_vars)
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

#-----------------------------#
# Fit Lag 3 ----
#-----------------------------#
fit_build3_stt <- dsem(
  sem = sem_lag3_stt,
  tsdata = data,
  family = family,
  estimate_delta0 = TRUE,
  control = dsem_control(run_model = FALSE))

pars3_stt <- fit_build3_stt$tmb_inputs$parameters
map3_stt  <- fit_build3_stt$tmb_inputs$map

pars3_stt$delta0_j <- rep(0, n_vars)
pars3_stt$lnsigma_j <- rep(log(0.4), n_vars)
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

#---------------------------------#
# Compare snow/tanner models ----
#---------------------------------#
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

summary(fit_lag1_stt)
summary(fit_lag2_stt)
summary(fit_lag3_stt)

#extract snow to tanner estimates from each model
sm_lag1_stt <- as.data.frame(summary(fit_lag1_stt))
sm_lag2_stt <- as.data.frame(summary(fit_lag2_stt))
sm_lag3_stt <- as.data.frame(summary(fit_lag3_stt))

est_lag1_stt <- sm_lag1_stt %>% filter(name == "snowtotanner") %>% mutate(lag = "Lag 1")
est_lag2_stt <- sm_lag2_stt %>% filter(name == "snowtotanner") %>% mutate(lag = "Lag 2")
est_lag3_stt <- sm_lag3_stt %>% filter(name == "snowtotanner") %>% mutate(lag = "Lag 3")

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

#---------------------------------------------------#
# Model comparison for causal pathways ----
#---------------------------------------------------#

#create dataframe with all AIC scores
aic_long <- bind_rows(
  tibble(Pathway = "Sea Ice → Tanner", Lag = c(1,2), AIC = aic_values_icetotanner),
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
  row_spec(4, extra_css = "border-bottom: 2px solid black;")

# Combine effect sizes
est_all <- bind_rows(
  est_all_icetotanner %>% mutate(Causal_Pathway = "Sea Ice → Tanner"),
  est_all_icetosnow %>% mutate(Causal_Pathway = "Sea Ice → Snow"),
  est_all_snowtotanner %>% mutate(Causal_Pathway = "Snow → Tanner"))

#Write as output so we don't have to rerun models to create lag figure
est_all %>%
  select(-1) %>%
  mutate(Causal_Pathway = Causal_Pathway %>%
           str_replace_all(" → ", " to ")) %>%
  write.csv("./output/tanner_lags.csv")

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

#---------------------------------------------------#
# Fit final Tanner crab model ----
#---------------------------------------------------#

#define SEM based on best lag structure from above
sem_final <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_tanner_abundance -> log_tanner_abundance, 1, ar_tanner
  
  # Causal pathways
  sea_ice -> log_tanner_abundance, 1, icetotanner
  sea_ice -> log_snow_abundance, 1, icetosnow
  log_snow_abundance -> log_tanner_abundance, 3, snowtotanner"

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

# fix process variance SD = 0.1
pars_final$lnsigma_j <- rep(log(0.1), ncol(data))

# prevent estimation of SD
map_final$lnsigma_j <- factor(rep(NA, ncol(data)))

#run final model fit with Delta0 and fixed SD
  #this final model estimates AR1 coefficients, lagged effects, 
  #and delta0 while keeping process SD fixed
  
fit_dsem <- dsem(sem=sem_final, tsdata=data, 
                   family=family,
                   estimate_delta0=TRUE,
                   control=dsem_control(parameters = pars_final,
                                        map = map_final,
                                        quiet = TRUE,
                                        getsd = TRUE))

summary(fit_dsem)

#save as output
write.csv(summary(fit_dsem), "./output/tanner_final_dsem_summary.csv", row.names = FALSE)

#----------------------------------------#
# Compute direct/indirect effects ----
#----------------------------------------#

#extract estimate and p-value for causal pathways 
paths_of_interest <- subset(summary(fit_dsem),
  (first == "sea_ice" & second == "log_tanner_abundance") |
    (first == "sea_ice" & second == "log_snow_abundance") |
    (first == "log_snow_abundance" & second == "log_tanner_abundance"))

paths_of_interest[, c("first", "second", "lag", "Estimate", "Std_Error", "p_value")]

#indirect effect of sea ice
icetosnow <- paths_of_interest %>%
  filter((first == "sea_ice" & second == "log_snow_abundance")) %>%
  pull(Estimate)

snowtotanner <- paths_of_interest %>%
  filter((first == "log_snow_abundance" & second == "log_tanner_abundance")) %>%
  pull(Estimate)

indirect = icetosnow * snowtotanner #-0.11

#Direct effect of sea ice
icetotanner <- paths_of_interest %>%
  filter((first == "sea_ice" & second == "log_tanner_abundance")) %>%
  pull(Estimate)
direct = icetotanner #0.04

#total effect of sea ice
total = indirect + icetotanner #-0.07

#relative importance of direct vrs indirect effect
prop_indirect = abs(indirect) / (abs(indirect) + abs(direct))
#overall, sea ice has a negative effect on tanner crab, and ~73% of the total 
  #effect is driven by the indirect pathway

#-------------------------------------------#
# Final Tanner crab model diagnostics ----
#-------------------------------------------#
#Convergence diagnostics:

#Hessian/SE - should be no hessian warnings or NA SE estimates
fit_dsem$sdrep

#check maximum final gradient- should be < 0.001
max(fit_dsem$sdrep$gradient.fixed)

#Identifiability/overfitting
fit_dsem$sdrep$pdHess # TRUE = good here
summary(fit_dsem$sdrep, "fixed")[, "Std. Error"] #shouldn't be any NA/NaN

#Residual diagnostics:

res <- residuals(fit_dsem)

#autocorrelation of residuals
acf(na.omit(res[, 1]))
acf(na.omit(res[, 2]))
acf(na.omit(res[, 3]))

#normality of residuals
par(mfrow = c(1,3))

for(i in 1:ncol(res)) {
  qqnorm(res[, i], main = colnames(res)[i])
  qqline(res[, i], col = "red")}
#approximately normal, but some deviations in tails 

#Latent State Diagnostics:

#extract all latent states
states_vec <- fit_dsem$sdrep$value[names(fit_dsem$sdrep$value) == "z_tj"]
  #each element = one state (time x variable)

#reshape into matrix 
n_time <- nrow(data)
n_vars <- 3
states_mat <- matrix(states_vec, nrow = n_time, ncol = n_vars, byrow = FALSE)
colnames(states_mat) <- c("sea_ice", "log_tanner_abundance", "log_snow_abundance")

#extract variables
obs_ice   <- data[, "sea_ice"]
obs_tanner <- data[, "log_tanner_abundance"]
obs_snow  <- data[, "log_snow_abundance"]

#latent vrs observed: tanner abundance
df_plot <- data.frame(
  time = 1:nrow(data),
  observed = obs_tanner,
  latent = states_mat[, "log_tanner_abundance"])

ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = observed), color = "black", linewidth = 0.8) +
  geom_line(aes(y = latent), color = "blue", linewidth = 0.8) +
  theme_bw() +
  labs(y = "Log Tanner abundance",
       title = "Observed vs Latent State")
#this is probably expected since we fixed SD?
  #ie model is estimating trend vrs latent process noise 

#residual vrs fitted: tanner abundance
res_tanner <- obs_tanner - states_mat[, "log_tanner_abundance"]

df_res <- data.frame(
  fitted = states_mat[, "log_tanner_abundance"],
  res = res_tanner)

ggplot(df_res, aes(x = fitted, y = res)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw()
#looks good 

#Would be worth following up on these, there are very few published diagnostics
  #in the dsem package currently 

#-------------------------------------------------#
#Plot final Tanner crab model fit ----
#-------------------------------------------------#

#build function to plot each fitted variable
  #function revised from C. Monnahan plot_fit function
plot_fit <- function(Y, fit, start_year = NULL){
  
  # Extract latent states (model estimates of underlying system state)
  ParHat <- fit$obj$env$parList()
  pred <- ParHat$x_tj
  
  # Extract standard errors around latent states
  sd_list <- as.list(fit$sdrep, what = "Std.")
  SD <- sd_list$x_tj
  
  #create time variable for missing data
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
      title = "Tanner crab model DSEM fit: observed vs estimated latent state") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"))
}

#plot
plot_fit(data, fit_dsem, start_year = start_year)

#and now plot fitted DAG
plot(as_fitted_DAG(fit_dsem, what = "path_coefficient"))
plot(as_fitted_DAG(fit_dsem_auto, what = c("Estimate", "Std_Error", "p_value")), type = "width")
#need to follow up on this- dsem output doesn't retain path names that 
  #as_fitted_DAG() function needs since we manually specified map and par

#-------------------------------------------#
# Compute cumulative effects ----
#-------------------------------------------#

#lag/cumulative effects plots:
  #lag: defines when the effect enters the system
  #variance partitioning over time shows how long a lag's influence persists 

# Calculate cumulative lagged effects
  #ie if variable X changes at time t, what is the cumulative effect of variable Y 
  #after lag k, accounting for direct + indirect pathways
effect = total_effect(fit_dsem, n_lags = 6) 

# Plot total effect
ggplot( effect) + 
  geom_bar(aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )
#snow->tanner strong negative delayed effect, and biological responses appear to be
  #delayed and cumulative 

#relative importance of variables as predictors of tanner abundance
  #i.e. at each time step, where proportion of tanner crab variability is coming from
partition_variance(fit_dsem,
                   which_response = "log_tanner_abundance",
                   n_times = 10 )
#tanner crab variance primarily driven by internal population dynamics at short lags, 
  #but indirect ice->snow->tanner pathway has a strong influence, and accumulates over time

#plot
var_df <- as.data.frame(partition_variance(fit_dsem,
                                           which_response = "log_tanner_abundance",
                                           n_times = 10)$proportion_variance_explained)

var_df$time <- 1:nrow(var_df)

# Pivot longer
var_long <- var_df %>%
  pivot_longer(cols = c(sea_ice, log_snow_abundance, log_tanner_abundance),
               names_to = "Component",
               values_to = "Proportion")

# Plot
ggplot(var_long, aes(x = time, y = Proportion, fill = Component)) +
  geom_area(alpha = 0.8, color = "black") +
  labs(x = "Time step",
       y = "Proportion of variance explained",
       fill = "Component",
       title = "Variance partitioning of Tanner crab dynamics") +
  theme_minimal()

#second view
ggplot(var_long, aes(x = time, y = Proportion, color = Component)) +
  geom_line(size = 1.2) +
  geom_point() +
  labs(x = "Time step",
       y = "Proportion of Tanner crab variance explained",
       title = "Tanner crab variance contributions over time") +
  theme_minimal()

#plot at final time step 10 
var_df %>%
  filter(time == 10) %>%
  pivot_longer(cols = c(sea_ice, log_snow_abundance, log_tanner_abundance),
               names_to = "Component",
               values_to = "Proportion") %>%
  ggplot(aes(x = Component, y = Proportion, fill = Component)) +
  geom_col() +
  ylim(0,1) +
  theme_minimal() +
  labs(title = "Variance partitioning at equilibrium",
       y = "Proportion of variance explained")

#extract at final time step
var_df %>%
  filter(time == max(time))
#at time step 10, AR explains 30%, snow abundance explains 54% and sea ice
  #explains 16% of tanner crab abundance

#note that for DSEM we don't have a single global Rsq like a regression model because
  #goodness of fit is multi-dimensional 

#-----------------------------------------------------------#
# Evaluate sensitivity to fixed observation error ----
#-----------------------------------------------------------#

#test plausible range of fixed values
obs_sd_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)

results_list <- list()

#loop through values and refit final model
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

#Both ice->snow and snow->tanner pathways are robust to changes in process variance

#-------------------------------------------------#
# Monte Carlo Simulation-based validation ----
#-------------------------------------------------#

#Here, we'll simulate data from final Tanner dsem model, refit the model to
  #each simulated dataset, and compare parameter estimates to true values to check
  #whether the fitted model can recover its own parameters when new data is simulated

simTestDSEM <- function(fitDSEM, sem, n_sim, start, family) {
    
    DSEMlist <- list()
    
    for (i in 1:n_sim) {
      
      # --- 1. Simulate data ---
      simDSEM <- simulate(fitDSEM, resimulate_gmrf = TRUE, fill_missing = TRUE)
      
      try({
        
        sim_data <- simDSEM[[1]]
        
        # --- 2. Build model (no run yet) ---
        fit_build <- dsem(
          sem = sem,
          tsdata = ts(sim_data, start = start),
          family = family,
          estimate_delta0 = TRUE,
          control = dsem_control(run_model = FALSE))
        
        pars <- fit_build$tmb_inputs$parameters
        map  <- fit_build$tmb_inputs$map
        
        # --- 3. Fix observation SD ---
        pars$lnsigma_j <- rep(log(0.1), ncol(sim_data))
        map$lnsigma_j  <- factor(rep(NA, ncol(sim_data)))
        
        # --- 4. Refit model ---
        tempFit <- dsem(
          sem = sem,
          tsdata = ts(sim_data, start = start),
          family = family,
          estimate_delta0 = TRUE,
          control = dsem_control(
            parameters = pars,
            map = map,
            quiet = TRUE,
            getsd = FALSE))
        
        # --- 5. Extract parameters, fitted values and observed values ---
        DSEMlist[[i]] <- list(
          beta = tempFit$opt$par[names(tempFit$opt$par) == "beta_z"],
          pred = predict(tempFit),
          obs  = sim_data)
      }, silent = TRUE)
    }
    
    DSEMlist <- DSEMlist[!sapply(DSEMlist, is.null)]
    
    return(DSEMlist)
  }
  
#now run 1000 simulations on our model 
selfSimExp <- simTestDSEM(
  fitDSEM = fit_dsem,
  sem = sem_final,
  n_sim = 1000,   
  start = 2000,
  family = family)

#Hessian warnings are expected from fit_build b/c parameters are at
  #default/initial value

#parameter recovery:
beta_list <- lapply(selfSimExp, function(x) x$beta) #extract beta parameter vectors
beta_mat <- do.call(cbind, beta_list) #rows = parameters, cols=simulations
df <- as.data.frame(beta_mat)

#add parameter labels from simulation
df <- df %>%
  mutate(Path  = fit_dsem$sem_full$path[fit_dsem$sem_full$parameter != 0],
    Param = fit_dsem$sem_full$name[fit_dsem$sem_full$parameter != 0]) %>%
  pivot_longer(cols = starts_with("V"), values_to = "Beta")

#true parameters from fitted model
df_true <- data.frame(
  Param = fit_dsem$sem_full$name[fit_dsem$sem_full$parameter != 0],
  Path  = fit_dsem$sem_full$path[fit_dsem$sem_full$parameter != 0],
  Beta_true = fit_dsem$opt$par[names(fit_dsem$opt$par) == "beta_z"])

#and join
param_df <- df %>%
  inner_join(df_true, by = c("Path", "Param"))
#the v[] params are our process variance for each variable

#plot parameter recovery results
ggplot(param_df %>% filter(!str_starts(Param, "V")),
       aes(x = Beta, y = Param)) +
  geom_violin(fill = "lightblue", scale = "width") +
  geom_point(aes(x = Beta_true), color = "red") +
  geom_vline(xintercept = 0) +
  theme_minimal() +
  ggtitle("Parameter Recovery")

#violin plots represent the distribution of estimated beta parameters across all 
  #simulations, red dot is true parameter value from original fitted model, 0 line
  #is "no effect"
#Looks good! Model is recovering parameters reliably

#compute bias and RMSE
param_df %>%
  filter(!str_starts(Param, "V")) %>%
  group_by(Param) %>%
  summarise(bias = mean(Beta - Beta_true), #average error
             rmse = sqrt(mean((Beta - Beta_true)^2)), #bias + variance
             sd = sd(Beta)) #precision- variability across simulations
#negligible bias, RMSE ~= SD, which tells us that error is random, not structural

#check extreme failures, ie instability/identifiability issues
param_df %>%
  filter(!str_starts(Param, "V")) %>%
  group_by(Param) %>%
  summarise(min = min(Beta),
            max = max(Beta))
#low to moderate uncertainty (expected in ecological timeseries)

#Fit + residual diagnostics:
  #i.e. can our simulations reproduce residual structure seen in our final 
  #fitted model 

#build residual dataset from simulations 
sim_res_df <- selfSimExp %>%
  imap(~ tibble(
    sim_id = as.character(.y),
    fitted = as.numeric(.x$pred),
    resid  = as.numeric(.x$obs - .x$pred))) %>%
  bind_rows()

#extract residuals from final Tanner model fit
fit_res_df <- tibble(
  sim_id = "original",
  fitted = as.numeric(predict(fit_dsem)),
  resid  = as.numeric(residuals(fit_dsem)))

#combine
all_res_df <- bind_rows(
  sim_res_df %>% mutate(type = "Simulated"),
  fit_res_df %>% mutate(type = "Original"))

#plot to compare residual distributions
all_res_df %>%
  ggplot(aes(x = resid, fill = type)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Residual Distribution: Simulated vs Original")

#plot autocorrelation of residuals
acf_df <- all_res_df %>%
  group_by(type, sim_id) %>%
  summarise(acf1 = acf(resid, plot = FALSE, na.action = na.pass)$acf[2])

orig_acf1 <- acf_df %>%
  filter(type == "Original") %>%
  pull(acf1)

acf_df %>%
  filter(type == "Simulated") %>%
  ggplot(aes(x = acf1)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = orig_acf1, color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Lag-1 residual autocorrelation")

#residuals vrs fit plot
ggplot() +
  geom_point(data = all_res_df %>% filter(type == "Simulated"),
    aes(fitted, resid),alpha = 0.05, color = "grey") +
  geom_point(data = all_res_df %>% filter(type == "Original"),
    aes(fitted, resid), color = "red", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Residual vs Fitted: Simulated (grey) vs Original (red)")


#calculate direct, indirect and total effect + 95% CI from simulation output
  #b/c we couldn't propagate uncertainty from multiplicative indirect effects 
  #when we did this earlier, so we can use our simulation to do so

#compute effects
effects <- param_df %>%
  select(Param, Beta) %>%
  group_by(Param) %>%
  mutate(sim = row_number()) %>%   # create simulation ID
  ungroup() %>%
  pivot_wider(names_from = Param, values_from = Beta) %>%
  mutate(direct   = icetotanner,
          indirect = icetosnow * snowtotanner,
          total    = direct + indirect)

#calculate uncertainty
effects %>%
  summarise(
    direct_mean = mean(direct),
    direct_2.5perc  = quantile(direct, 0.025),
    direct_97.5perc = quantile(direct, 0.975),
    
    indirect_mean = mean(indirect),
    indirect_2.5perc  = quantile(indirect, 0.025),
    indirect_97.5perc = quantile(indirect, 0.975),
    
    total_mean = mean(total),
    total_2.5perc  = quantile(total, 0.025),
    total_97.5perc = quantile(total, 0.975))
#very similar to estimates from final model output!

#-------------------------------------------#
# Hindcasting validation ----
#-------------------------------------------#

#If we remove the MHW event (i.e. use 1988-2017 only), are inferred relationships
  #still observed? i.e. are pathways persistent, or driven by MHW

#filter dataset
pre_mhw <- sea_ice %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
  filter(year >= start_year) %>% 
  rename(sea_ice = Mar_Apr_ice_EBS_NBS) %>%
  full_join(crab_abund %>%
              filter(category == "population_50mm_plus" & species != "hybrid") %>%
              select(year, abundance, species) %>%
              pivot_wider(names_from = "species", values_from = "abundance") %>%
              rename(snow_abundance=snow, tanner_abundance=tanner)) %>%
  filter(year <= 2017)

#Scale all variables for dsem
vars <- c("sea_ice",
          "snow_abundance", "log_snow_abundance",
          "tanner_abundance", "log_tanner_abundance")

pre_mhw_data <- pre_mhw %>%
  mutate(log_tanner_abundance = log(tanner_abundance),
         log_snow_abundance = log(snow_abundance)) %>%
  mutate(across(all_of(c(vars)), ~ as.numeric(scale(.)))) %>%
  select(-tanner_abundance, -snow_abundance)

#prep data for dsem model
pre_mhw_data <- pre_mhw_data %>%
  select(-year) %>%
  ts()

#run model using same specifications/lags from above 

#two-stage modeling fitting procedure: 
##initial first model run without delta0 (to improve starting values)
fit0_pre <- dsem(sem = sem_final, tsdata = pre_mhw_data,
             family = family, estimate_delta0 = FALSE,
             control = dsem_control(quiet = FALSE, getsd = FALSE))

# extract starting parameters
parameters_pre <- fit0_pre$obj$env$parList()

# add starting values for delta0
parameters_pre$delta0_j <- rep(0, ncol(pre_mhw_data))

#build model without running it
#(needed so we can modify TMB inputs)
fit_build_final_pre <- dsem(sem = sem_final, tsdata = pre_mhw_data,
                        family = family,
                        estimate_delta0 = TRUE,
                        control = dsem_control(
                          run_model = FALSE,
                          parameters = parameters_pre))

pars_final_pre <- fit_build_final_pre$tmb_inputs$parameters
map_final_pre  <- fit_build_final_pre$tmb_inputs$map

# fix observation SD = 0.1
pars_final_pre$lnsigma_j <- rep(log(0.1), ncol(pre_mhw_data))

# prevent estimation of observation SD
map_final_pre$lnsigma_j <- factor(rep(NA, ncol(pre_mhw_data)))

#run final model fit with Delta0 and fixed observation error
fit_dsem_pre_mhw <- dsem(sem=sem_final, tsdata=pre_mhw_data, 
                 family=family,
                 estimate_delta0=TRUE,
                 control=dsem_control(parameters = pars_final_pre,
                                      map = map_final_pre,
                                      quiet = TRUE,
                                      getsd = TRUE))

summary(fit_dsem_pre_mhw)
#Interesting! So our snow-tanner causal relationship holds even without the MHW-
  #Evidence for consistent, negative effect of snow crab on tanner crab with 3 yr lag 
  #Mortality mechanism for snow crab is dependent on the MHW (though positive effect remains)
























