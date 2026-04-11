#Goals ----
#Fit a dsem model to examine causal linkages between the
#MHW and hybrid abundance increase

#Specifically, we'll test the following causal pathways:
#sea ice -> snow crab population abundance 
  #Mechanism testing: sea ice loss/mass mortality events 
  #Hypothesized lags to test: 1-2 years
#sea ice -> hybrid population abundance 
  #Mechanism testing: sea ice effects on survival/productivity
  #Hypothesized lags to test: 1-2 years
#snow crab population abundance -> hybrid population abundance
  #Mechanism: snow crab collapse results in competitive release/niche & habitat expansion
  #Hypothesized lags to test: 1-3 years 
#sea ice -> snow/tanner spatial overlap
  #Mechanism: temperature effects on spatial distributions
  #Hypothesized lags to test: 0-1 years
#snow/tanner spatial overlap -> hybrid abundance
  #Mechanism: increased interspecies mating opportunities
  #Hypothesized lags to test: 5-6 years

#NOTES: Because hybridization results from bidirectional parental crosses, 
  #we'll test both female snow/male tanner and male tanner/female snow overlap 
  #metrics in separate final models to ensure that overlap results are robust to 
  #both, but will only test lags and report results for female snow/tanner overlap
  #for parsimony and simplicity

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

#-----------------------------------------#
# Define covariates and standardize ---- 
#-----------------------------------------#

#wrangle data
dat_hybrid <- sea_ice %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
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
  patchwork::wrap_plots(plots)}

plot_histo(dat_hybrid %>% select(-year))

#Scale all variables and transform abundance variables 
vars <- c("sea_ice", "bhatta_snowfem_tanmale", "bhatta_tannerfem_snowmale",
          "snow_abundance", "log_snow_abundance",
          "hybrid_abundance", "log_hybrid_abundance")

hybrid_data <- dat_hybrid %>%
  mutate(log_hybrid_abundance = log(hybrid_abundance),
         log_snow_abundance = log(snow_abundance)) %>%
  mutate(across(all_of(c(vars)), ~ as.numeric(scale(.)))) %>%
  select(-hybrid_abundance, -snow_abundance)

#check distributions now
plot_histo(hybrid_data %>% select(-year))

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

#-------------------------------#
#Prep data for dsem models ----
#-------------------------------#

data <- hybrid_data %>%
  select(-year, -bhatta_tannerfem_snowmale) %>%
  ts()

#And we'll use this dataset for testing our second overlap metric
data2 <- hybrid_data %>%
  select(-year, -bhatta_snowfem_tanmale) %>%
  ts()

family <- rep("normal", ncol(data))

#---------------------------------------#
# Test for DAG-data consistency ----
#---------------------------------------#

#download specified DAG from dagitty.net - conditional independencies are 
#identified based on structure of DAG drawn in dagitty
hybrid_dag <- dagitty('dag {
"Hybrid abundance" [outcome,pos="-0.958,0.004"]
"Sea Ice" [exposure,pos="-1.800,0.078"]
"Snow and Tanner spatial overlap" [pos="-1.342,0.380"]
"Snow crab abundance" [pos="-1.382,-0.416"]
"Sea Ice" -> "Hybrid abundance"
"Sea Ice" -> "Snow and Tanner spatial overlap"
"Sea Ice" -> "Snow crab abundance" [pos="-1.453,-0.356"]
"Snow and Tanner spatial overlap" -> "Hybrid abundance"
"Snow crab abundance" -> "Hybrid abundance"
}
')

#plot DAG
ggdag(hybrid_dag, layout = "nicely") +
  theme_dag()

plot(hybrid_dag) 

ggdag_status(hybrid_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(hybrid_dag)
#3 open causal pathways 

#and plot paths
ggdag_paths(hybrid_dag, text = FALSE, use_labels = "name") +
  theme_dag()

#find adjustment sets for response variable 
adjustmentSets(hybrid_dag, exposure="Sea Ice", outcome="Hybrid abundance")
#This tells us that no covariate adjustment is necessary to identify the causal effect
#ie there are no open backdoor paths 

#and visualize adjustment sets, if there are any
ggdag_adjustment_set(hybrid_dag, shadow = TRUE) +
  theme_dag()

#find conditional independencies- i.e. two variables that are implied to 
# be independent and not correlated shouldn't be connected by a node
impliedConditionalIndependencies(hybrid_dag)
#DAG assumes that spatial overlap is conditionally independent of snow crab 
  #abundance given sea ice (fork example), i.e. sea ice is the only meaningful 
  #driver linking overlap and snow crab abundance and we don't need a node b/w
  #snow crab abundance and overlap

#d-sep test: does this conditional independency hold in our data?

#Assess for snow female/tanner male overlap metric first:
  #pre-lag dataset using final lag structure from below
lagged_dat <- hybrid_data %>% 
  select(-bhatta_tannerfem_snowmale) %>%
  arrange(year) %>%
  mutate(hybrid_t   = log_hybrid_abundance,
    snow_t   = log_snow_abundance,
    SeaIce_t = sea_ice,
    SeaIce_t1  = lag(sea_ice, 1),
    snow_t3    = lag(log_snow_abundance, 3), 
    overlap_t = bhatta_snowfem_tanmale,
    overlap_t5 = lag(bhatta_snowfem_tanmale, 5)) %>%
  select(hybrid_t, snow_t, SeaIce_t1, SeaIce_t, snow_t3, overlap_t, overlap_t5) %>%
  na.omit()

#build a lagged dag- we'll need to add missing temporal links for snow/overlap
  #variables being lagged at different times
lagged_dag <- dagitty('dag {
SeaIce_t1 -> snow_t
SeaIce_t  -> overlap_t
SeaIce_t1 -> hybrid_t

snow_t -> snow_t3
overlap_t -> overlap_t5

snow_t3 -> hybrid_t
overlap_t5 -> hybrid_t
}')

localTests(lagged_dag, lagged_dat)
#All ov_||_ sn_ terms are > 0.5 so implied snow crab/overlap independency is supported 
  #and we've confirmed DAG-data consistency

#now do the same for tanner female/snow male overlap metric using the same DAG:
lagged_dat2 <- hybrid_data %>% 
  select(-bhatta_snowfem_tanmale) %>%
  arrange(year) %>%
  mutate(hybrid_t   = log_hybrid_abundance,
         snow_t   = log_snow_abundance,
         SeaIce_t = sea_ice,
         SeaIce_t1  = lag(sea_ice, 1),
         snow_t3    = lag(log_snow_abundance, 3), 
         overlap_t = bhatta_tannerfem_snowmale,
         overlap_t5 = lag(bhatta_tannerfem_snowmale, 5)) %>%
  select(hybrid_t, snow_t, SeaIce_t1, SeaIce_t, snow_t3, overlap_t, overlap_t5) %>%
  na.omit()

localTests(lagged_dag, lagged_dat2)
#Implied independancy using other overlap metric also supported by data

#----------------------------------------------------#
#Lag testing: sea ice -> hybrid causal pathway ----
#----------------------------------------------------#

#Note that because the hybrid models have more causal pathways to estimate compared
  #to the tanner models, we won't estimate delta0 and will fix observation error at 0.1

#Define SEMs for each sea ice -> hybrid lag
sem_lag1_hybrid <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

# Causal pathways
  #starting with minimum lags to be tested for all other causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

sem_lag2_hybrid <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

  # Causal pathways
  sea_ice -> log_hybrid_abundance, 2, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

# ----------------------------#
# Fit Lag 1 ----
#-----------------------------#
#build model without running it (needed so we can modify TMB inputs)
fit_build1_hybrid <- dsem(
  sem = sem_lag1_hybrid,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

#Extract parameters & map from full model
pars1_hybrid <- fit_build1_hybrid$tmb_inputs$parameters
map1_hybrid  <- fit_build1_hybrid$tmb_inputs$map

#Set reasonable starting values
n_vars <- ncol(data) #used to index parameters
pars1_hybrid$lnsigma_j <- rep(log(0.1), n_vars) #small observation error to improve stability
map1_hybrid$lnsigma_j <- factor(rep(NA, n_vars)) #tells TMB to fix observation error
n_beta1_hybrid <- length(pars1_hybrid$beta_z) #total # of coefficients in the SEM (beta_z)
pars1_hybrid$beta_z[(n_vars+1):n_beta1_hybrid] <- 0.05  # small lag starting values

fit_lag1_hybrid <- dsem(
  sem = sem_lag1_hybrid,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars1_hybrid,
    map = map1_hybrid,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 2 ----
#-----------------------------#
fit_build2_hybrid <- dsem(
  sem = sem_lag2_hybrid,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars2_hybrid <- fit_build2_hybrid$tmb_inputs$parameters
map2_hybrid  <- fit_build2_hybrid$tmb_inputs$map

#set starting parameters
pars2_hybrid$lnsigma_j <- rep(log(0.1), n_vars)
map2_hybrid$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2_hybrid <- length(pars2_hybrid$beta_z)
pars2_hybrid$beta_z[(n_vars+1):n_beta2_hybrid] <- 0.05

fit_lag2_hybrid <- dsem(
  sem = sem_lag2_hybrid,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars2_hybrid,
    map = map2_hybrid,
    quiet = TRUE,
    getsd = TRUE))

#---------------------------------#
# Compare ice/hybrid models ----
#---------------------------------#
aic_values_icetohybrid <- c(
  lag1 = AIC(fit_lag1_hybrid),
  lag2 = AIC(fit_lag2_hybrid))

loglik_values_icetohybrid <- c(
  lag1 = logLik(fit_lag1_hybrid),
  lag2 = logLik(fit_lag2_hybrid))

aic_values_icetohybrid
loglik_values_icetohybrid

summary(fit_lag1_hybrid)
summary(fit_lag2_hybrid)

#extract icetohybrid estimates from each model
sm_lag1_hybrid <- as.data.frame(summary(fit_lag1_hybrid))
sm_lag2_hybrid <- as.data.frame(summary(fit_lag2_hybrid))

est_lag1_hybrid <- sm_lag1_hybrid %>% filter(name == "icetohybrid") %>% mutate(lag = "Lag 1")
est_lag2_hybrid <- sm_lag2_hybrid %>% filter(name == "icetohybrid") %>% mutate(lag = "Lag 2")

# Combine into one data frame
est_all_icetohybrid <- bind_rows(est_lag1_hybrid, est_lag2_hybrid)

# Plot effect sizes with error bars
ggplot(est_all_icetohybrid, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 1" = "#1f78b4", 
                               "Lag 2" = "#33a02c")) +
  labs(title = "Sea Ice → hybrid Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#Very small insignificant effect size and AIC can't distinguish between lag structures- we'll stick
  #with lag 1 for parsimony/consistency with snow crab lag

#------------------------------------------------#
#Lag testing: ice -> snow causal pathway ----
#------------------------------------------------#

#Define SEMs for each sea ice -> snow crab lag
sem_lag1_snow <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

  # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

sem_lag2_snow <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

  # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 2, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

#-----------------------------#
# Fit Lag 1 ----
#-----------------------------#
fit_build1_snow <- dsem(
  sem = sem_lag1_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars1_snow <- fit_build1_snow$tmb_inputs$parameters
map1_snow  <- fit_build1_snow$tmb_inputs$map

# Safe starting values
n_vars <- ncol(data)
pars1_snow$lnsigma_j <- rep(log(0.1), n_vars) #observation error
map1_snow$lnsigma_j <- factor(rep(NA, n_vars))
n_beta1_snow <- length(pars1_snow$beta_z)
pars1_snow$beta_z[(n_vars+1):n_beta1_snow] <- 0.05  # small lag starting values

fit_lag1_snow <- dsem(
  sem = sem_lag1_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
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
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars2_snow <- fit_build2_snow$tmb_inputs$parameters
map2_snow  <- fit_build2_snow$tmb_inputs$map

pars2_snow$lnsigma_j <- rep(log(0.1), n_vars)
map2_snow$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2_snow <- length(pars2_snow$beta_z)
pars2_snow$beta_z[(n_vars+1):n_beta2_snow] <- 0.05

fit_lag2_snow <- dsem(
  sem = sem_lag2_snow,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
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

#AIC supports lag 1

#--------------------------------------------------#
#Lag testing: snow -> hybrid causal pathway ----
#--------------------------------------------------#

#Define SEMs for each snow to hybrid (sth) lag being tested
sem_lag1_sth <- "
# AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

 # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 1, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

sem_lag2_sth <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

   # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

sem_lag3_sth <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

   # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 3, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

# -----------------------------#
# Fit Lag 1 ----
#-----------------------------#
fit_build1_sth <- dsem(
  sem = sem_lag1_sth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars1_sth <- fit_build1_sth$tmb_inputs$parameters
map1_sth  <- fit_build1_sth$tmb_inputs$map

# Safe starting values
n_vars <- ncol(data)
pars1_sth$lnsigma_j <- rep(log(0.1), n_vars) 
map1_sth$lnsigma_j <- factor(rep(NA, n_vars))
n_beta1_sth <- length(pars1_sth$beta_z)
pars1_sth$beta_z[(n_vars+1):n_beta1_sth] <- 0.05  # small lag starting values

fit_lag1_sth <- dsem(
  sem = sem_lag1_sth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars1_sth,
    map = map1_sth,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 2 ----
#-----------------------------#
fit_build2_sth <- dsem(
  sem = sem_lag2_sth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars2_sth <- fit_build2_sth$tmb_inputs$parameters
map2_sth  <- fit_build2_sth$tmb_inputs$map

pars2_sth$lnsigma_j <- rep(log(0.1), n_vars)
map2_sth$lnsigma_j <- factor(rep(NA, n_vars))
n_beta2 <- length(pars2_sth$beta_z)
pars2_sth$beta_z[(n_vars+1):n_beta2] <- 0.05

fit_lag2_sth <- dsem(
  sem = sem_lag2_sth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars2_sth,
    map = map2_sth,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 3 ----
#-----------------------------#
fit_build3_sth <- dsem(
  sem = sem_lag3_sth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars3_sth <- fit_build3_sth$tmb_inputs$parameters
map3_sth  <- fit_build3_sth$tmb_inputs$map

pars3_sth$lnsigma_j <- rep(log(0.1), n_vars)
map3_sth$lnsigma_j <- factor(rep(NA, n_vars))
n_beta3 <- length(pars3_sth$beta_z)
pars3_sth$beta_z[(n_vars+1):n_beta3] <- 0.05

fit_lag3_sth <- dsem(
  sem = sem_lag3_sth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars3_sth,
    map = map3_sth,
    quiet = TRUE,
    getsd = TRUE))

#---------------------------------#
# Compare snow/hybrid models ----
#---------------------------------#
aic_values_snowtohybrid <- c(
  lag1 = AIC(fit_lag1_sth),
  lag2 = AIC(fit_lag2_sth),
  lag3 = AIC(fit_lag3_sth))

loglik_values_snowtohybrid <- c(
  lag1 = logLik(fit_lag1_sth),
  lag2 = logLik(fit_lag2_sth),
  lag3 = logLik(fit_lag3_sth))

aic_values_snowtohybrid
loglik_values_snowtohybrid

summary(fit_lag1_sth)
summary(fit_lag2_sth)
summary(fit_lag3_sth)

#extract snow to hybrid estimates from each model
sm_lag1_sth <- as.data.frame(summary(fit_lag1_sth))
sm_lag2_sth <- as.data.frame(summary(fit_lag2_sth))
sm_lag3_sth <- as.data.frame(summary(fit_lag3_sth))

est_lag1_sth <- sm_lag1_sth %>% filter(name == "snowtohybrid") %>% mutate(lag = "Lag 1")
est_lag2_sth <- sm_lag2_sth %>% filter(name == "snowtohybrid") %>% mutate(lag = "Lag 2")
est_lag3_sth <- sm_lag3_sth %>% filter(name == "snowtohybrid") %>% mutate(lag = "Lag 3")

# Combine into one data frame
est_all_snowtohybrid <- bind_rows(est_lag1_sth, est_lag2_sth, est_lag3_sth)

# Plot effect sizes with error bars
ggplot(est_all_snowtohybrid, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 1" = "#1f78b4", 
                               "Lag 2" = "#33a02c",
                               "Lag 3" = "#e31a1c")) +
  labs(title = "Snow Abundance → Hybrid Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#Best support for lag 2

#----------------------------------------------------------------#
#Lag testing: spatial overlap -> hybrid causal pathway ----
#-----------------------------------------------------------------#

#Define SEMs for each overlap to hybrid (oth) lag being tested
sem_lag5_oth <- "
# AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

 # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

sem_lag6_oth <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

   # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
 bhatta_snowfem_tanmale -> log_hybrid_abundance, 6, overlaptohybrid"

# -----------------------------#
# Fit Lag 5 ----
#-----------------------------#
fit_build5_oth <- dsem(
  sem = sem_lag5_oth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars5_oth <- fit_build5_oth$tmb_inputs$parameters
map5_oth  <- fit_build5_oth$tmb_inputs$map

# Safe starting values
n_vars <- ncol(data)
pars5_oth$lnsigma_j <- rep(log(0.1), n_vars) 
map5_oth$lnsigma_j <- factor(rep(NA, n_vars))
n_beta5_oth <- length(pars5_oth$beta_z)
pars5_oth$beta_z[(n_vars+1):n_beta5_oth] <- 0.05  # small lag starting values

fit_lag5_oth <- dsem(
  sem = sem_lag5_oth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars5_oth,
    map = map5_oth,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 6 ----
#-----------------------------#
fit_build6_oth <- dsem(
  sem = sem_lag6_oth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars6_oth <- fit_build6_oth$tmb_inputs$parameters
map6_oth  <- fit_build6_oth$tmb_inputs$map

pars6_oth$lnsigma_j <- rep(log(0.1), n_vars)
map6_oth$lnsigma_j <- factor(rep(NA, n_vars))
n_beta6 <- length(pars6_oth$beta_z)
pars6_oth$beta_z[(n_vars+1):n_beta6] <- 0.05

fit_lag6_oth <- dsem(
  sem = sem_lag6_oth,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars6_oth,
    map = map6_oth,
    quiet = TRUE,
    getsd = TRUE))

#---------------------------------#
# Compare overlap/hybrid models ----
#---------------------------------#
aic_values_overlaptohybrid <- c(
  lag5 = AIC(fit_lag5_oth),
  lag6 = AIC(fit_lag6_oth))

loglik_values_overlaptohybrid <- c(
  lag5 = logLik(fit_lag5_oth),
  lag6 = logLik(fit_lag6_oth))

aic_values_overlaptohybrid
loglik_values_overlaptohybrid

summary(fit_lag5_oth)
summary(fit_lag6_oth)

#extract overlap to hybrid estimates from each model
sm_lag5_oth <- as.data.frame(summary(fit_lag5_oth))
sm_lag6_oth <- as.data.frame(summary(fit_lag6_oth))

est_lag5_oth <- sm_lag5_oth %>% filter(name == "overlaptohybrid") %>% mutate(lag = "Lag 5")
est_lag6_oth <- sm_lag6_oth %>% filter(name == "overlaptohybrid") %>% mutate(lag = "Lag 6")

# Combine into one data frame
est_all_overlaptohybrid <- bind_rows(est_lag5_oth, est_lag6_oth)

# Plot effect sizes with error bars
ggplot(est_all_overlaptohybrid, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 5" = "#1f78b4", 
                               "Lag 6" = "#33a02c")) +
  labs(title = "Spatial Overlap → Hybrid Abundance Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#AIC scores are similar, effect sizes are small and insignificant. We'll 
  #go with lowest AIC/largest effect size and use lag 6 

#----------------------------------------------------#
#Lag testing: sea ice -> overlap causal pathway ----
#----------------------------------------------------#

#Define SEMs for each sea ice -> overlap (ito) lag
sem_lag0_ito <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

# Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

sem_lag1_ito <- "
  # AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

  # Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 1, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 5, overlaptohybrid"

# ----------------------------#
# Fit Lag 0 ----
#-----------------------------#
#build model without running it (needed so we can modify TMB inputs)
fit_build0_ito <- dsem(
  sem = sem_lag0_ito,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

#Extract parameters & map from full model
pars0_ito <- fit_build0_ito$tmb_inputs$parameters
map0_ito  <- fit_build0_ito$tmb_inputs$map

#Set reasonable starting values
n_vars <- ncol(data) #used to index parameters
pars0_ito$lnsigma_j <- rep(log(0.1), n_vars) #small observation error to improve stability
map0_ito$lnsigma_j <- factor(rep(NA, n_vars)) #tells TMB to fix observation error
n_beta0_ito <- length(pars0_ito$beta_z) #total # of coefficients in the SEM (beta_z)
pars0_ito$beta_z[(n_vars+1):n_beta0_ito] <- 0.05  # small lag starting values

fit_lag0_ito <- dsem(
  sem = sem_lag0_ito,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars0_ito,
    map = map0_ito,
    quiet = TRUE,
    getsd = TRUE))

#-----------------------------#
# Fit Lag 1 ----
#-----------------------------#
fit_build1_ito <- dsem(
  sem = sem_lag1_ito,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(run_model = FALSE))

pars1_ito <- fit_build1_ito$tmb_inputs$parameters
map1_ito  <- fit_build1_ito$tmb_inputs$map

#set starting parameters
pars1_ito$lnsigma_j <- rep(log(0.1), n_vars)
map1_ito$lnsigma_j <- factor(rep(NA, n_vars))
n_beta1_ito <- length(pars1_ito$beta_z)
pars1_ito$beta_z[(n_vars+1):n_beta1_ito] <- 0.05

fit_lag1_ito <- dsem(
  sem = sem_lag1_ito,
  tsdata = data,
  family = family,
  estimate_delta0 = FALSE,
  control = dsem_control(
    parameters = pars1_ito,
    map = map1_ito,
    quiet = TRUE,
    getsd = TRUE))

#---------------------------------#
# Compare ice/overlap models ----
#---------------------------------#
aic_values_icetooverlap <- c(
  lag0 = AIC(fit_lag0_ito),
  lag1 = AIC(fit_lag1_ito))

loglik_values_icetooverlap <- c(
  lag0 = logLik(fit_lag0_ito),
  lag1 = logLik(fit_lag1_ito))

aic_values_icetooverlap
loglik_values_icetooverlap

summary(fit_lag0_ito)
summary(fit_lag1_ito)

#extract icetooverlap estimates from each model
sm_lag0_ito <- as.data.frame(summary(fit_lag0_ito))
sm_lag1_ito <- as.data.frame(summary(fit_lag1_ito))

est_lag0_ito <- sm_lag0_ito %>% filter(name == "icetooverlap") %>% mutate(lag = "Lag 0")
est_lag1_ito <- sm_lag1_ito %>% filter(name == "icetooverlap") %>% mutate(lag = "Lag 1")

# Combine into one data frame
est_all_icetooverlap <- bind_rows(est_lag0_ito, est_lag1_ito)

# Plot effect sizes with error bars
ggplot(est_all_icetooverlap, aes(x = lag, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Lag 0" = "#1f78b4", 
                               "Lag 1" = "#33a02c")) +
  labs(title = "Sea Ice → Spatial Overlap Effect Sizes by Lag",
       x = "",
       y = "Effect Size (Estimate ± Std Error)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#Again, AIC can't distinguish between lag structures- we'll stick
#with lag 0 for parsimony

#---------------------------------------------------#
# Model comparison for causal pathways ----
#---------------------------------------------------#

#create dataframe with all AIC scores
aic_long <- bind_rows(
  tibble(Pathway = "Sea Ice → Hybrid", Lag = c(1,2), AIC = aic_values_icetohybrid),
  tibble(Pathway = "Sea Ice → Snow",   Lag = c(1,2),   AIC = aic_values_icetosnow),
  tibble(Pathway = "Snow → Hybrid",    Lag = c(1,2,3), AIC = aic_values_snowtohybrid),
  tibble(Pathway = "Overlap → Hybrid",    Lag = c(5,6), AIC = aic_values_overlaptohybrid),
  tibble(Pathway = "Ice → Overlap",    Lag = c(0,1), AIC = aic_values_icetooverlap))

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
  row_spec(4, extra_css = "border-bottom: 2px solid black;") %>%
  row_spec(6, extra_css = "border-bottom: 2px solid black;") %>%
  row_spec(8, extra_css = "border-bottom: 2px solid black;")

# Combine effect sizes
est_all <- bind_rows(
  est_all_icetohybrid %>% mutate(Causal_Pathway = "Sea Ice → Hybrid"),
  est_all_icetosnow %>% mutate(Causal_Pathway = "Sea Ice → Snow"),
  est_all_snowtohybrid %>% mutate(Causal_Pathway = "Snow → Hybrid"),
  est_all_overlaptohybrid %>% mutate(Causal_Pathway = "Overlap → Hybrid"),
  est_all_icetooverlap %>% mutate(Causal_Pathway = "Ice → Overlap"))

#Write as output so we don't have to rerun models to create lag figure
est_all %>%
  select(-1) %>%
  mutate(Causal_Pathway = Causal_Pathway %>%
           str_replace_all(" → ", " to ")) %>%
  write.csv("./output/hybrid_lags.csv")

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
# Fit final hybrid model ----
#---------------------------------------------------#

#define SEM based on best lag structure from above
sem_final <- "
  #AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
  bhatta_snowfem_tanmale -> bhatta_snowfem_tanmale, 1, ar_overlap

  #Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_snowfem_tanmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_snowfem_tanmale -> log_hybrid_abundance, 6, overlaptohybrid"

#build model without running it
#(needed so we can modify TMB inputs)
fit_build_final <- dsem(sem = sem_final, tsdata = data,
                        family = family,
                        estimate_delta0 = FALSE,
                        control = dsem_control(
                          run_model = FALSE))

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
                 estimate_delta0=FALSE,
                 control=dsem_control(parameters = pars_final,
                                      map = map_final,
                                      quiet = TRUE,
                                      getsd = TRUE))

summary(fit_dsem)

#save as output
write.csv(summary(fit_dsem), "./output/hybrid_final_dsem_summary.csv", row.names = FALSE)

#-------------------------------------------------------#
# Fit final hybrid model using second overlap index ----
#-------------------------------------------------------#

#define SEM based on best lag structure from above, but using female Tanner/male
  #overlap index instead to verify if overlap results are robust to either parental
  #cross measure of overlap
sem_final_2 <- "
  #AR terms
  sea_ice -> sea_ice, 1, ar_seaice
  log_snow_abundance -> log_snow_abundance, 1, ar_snow
  log_hybrid_abundance -> log_hybrid_abundance, 1, ar_hybrid
 bhatta_tannerfem_snowmale -> bhatta_tannerfem_snowmale, 1, ar_overlap

  #Causal pathways
  sea_ice -> log_hybrid_abundance, 1, icetohybrid
  sea_ice -> log_snow_abundance, 1, icetosnow
  sea_ice -> bhatta_tannerfem_snowmale, 0, icetooverlap
  log_snow_abundance -> log_hybrid_abundance, 2, snowtohybrid
  bhatta_tannerfem_snowmale -> log_hybrid_abundance, 6, overlaptohybrid"

#build model without running it
fit_build_final2 <- dsem(sem = sem_final_2, tsdata = data2,
                        family = family,
                        estimate_delta0 = FALSE,
                        control = dsem_control(
                          run_model = FALSE))

pars_final_b <- fit_build_final2$tmb_inputs$parameters
map_final_b  <- fit_build_final2$tmb_inputs$map

# fix process variance SD = 0.1
pars_final_b$lnsigma_j <- rep(log(0.1), ncol(data2))

# prevent estimation of SD
map_final_b$lnsigma_j <- factor(rep(NA, ncol(data2)))

#run final model fit with Delta0 and fixed SD
fit_dsem2 <- dsem(sem=sem_final_2, tsdata=data2, 
                 family=family,
                 estimate_delta0=FALSE,
                 control=dsem_control(parameters = pars_final_b,
                                      map = map_final_b,
                                      quiet = TRUE,
                                      getsd = TRUE,
                                      newton_loops = 3))

summary(fit_dsem2)

#confirming that our results are robust regardless which overlap metric we use-
  #sea ice->overlap and overlap->hybrid effects are the same direction and insignificant
  #both both parental cross metrics 

#----------------------------------------#
# Compute direct/indirect effects ----
#----------------------------------------#

#extract estimate and p-value for causal pathways 
paths_of_interest <- subset(summary(fit_dsem),
                            (first == "sea_ice" & second == "log_hybrid_abundance") |
                              (first == "sea_ice" & second == "log_snow_abundance") |
                              (first == "log_snow_abundance" & second == "log_hybrid_abundance") |
                              (first == "bhatta_snowfem_tanmale" & second == "log_hybrid_abundance") |
                              (first == "sea_ice" & second == "bhatta_snowfem_tanmale"))

paths_of_interest[, c("first", "second", "lag", "Estimate", "Std_Error", "p_value")]

#indirect effect of sea ice
icetosnow <- paths_of_interest %>%
  filter((first == "sea_ice" & second == "log_snow_abundance")) %>%
  pull(Estimate)

snowtohybrid <- paths_of_interest %>%
  filter((first == "log_snow_abundance" & second == "log_hybrid_abundance")) %>%
  pull(Estimate)

icetooverlap <- paths_of_interest %>%
  filter((first == "sea_ice" & second == "bhatta_snowfem_tanmale")) %>%
  pull(Estimate)

overlaptohybrid <- paths_of_interest %>%
  filter((first == "bhatta_snowfem_tanmale" & second == "log_hybrid_abundance")) %>%
  pull(Estimate)

indirect_snow = icetosnow * snowtohybrid
indirect_overlap = icetooverlap * overlaptohybrid

#indirect effect = sum of indirect pathways
indirect = indirect_snow + indirect_overlap # -0.09

#Direct effect of sea ice
icetohybrid <- paths_of_interest %>%
  filter((first == "sea_ice" & second == "log_hybrid_abundance")) %>%
  pull(Estimate)

direct = icetohybrid # -0.01

#total effect of sea ice
total = indirect + direct # -0.10

#relative importance of direct vrs indirect effect

#identify most important indirect effect
most_important = max(abs(indirect_snow), abs(indirect_overlap)) #indirect snow, 0.09 estimate

prop_indirect = most_important / 
  (abs(direct) + abs(indirect_snow) + abs(indirect_overlap))
#overall, sea ice has a negative effect on hybrids, and the dominant indirect effect, 
  #ice->snow crab-> hybrids, accounts for ~76% of the total 
  #effect of sea ice on hybrid abundance

#-------------------------------------------#
# Final hybrid model diagnostics ----
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

res_df <- as.data.frame(res) %>%
  mutate(time = 1:n()) %>%
  pivot_longer(-time, names_to = "variable", values_to = "residual")

#autocorrelation of residuals
acf_df <- res_df %>%
  group_by(variable) %>%
  summarise(acf = list(acf(na.omit(residual), plot = FALSE))) %>%
  mutate(lag = map(acf, ~ .x$lag),
         acf_val = map(acf, ~ .x$acf)) %>%
  unnest(c(lag, acf_val))

ggplot(acf_df, aes(x = lag, y = acf_val)) +
  geom_hline(yintercept = 0) +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "Lag",
       y = "ACF",
       title = "Autocorrelation of Residuals") +
  theme_minimal()

#normality of residuals
ggplot(res_df, aes(sample = residual)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "QQ Plots of Residuals",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()
#approximately normal, but some deviations in tails 

#Latent State Diagnostics:

#extract all latent states
states_vec <- fit_dsem$sdrep$value[names(fit_dsem$sdrep$value) == "z_tj"]
#each element = one state (time x variable)

#reshape into matrix 
n_time <- nrow(data)
n_vars <- 4
states_mat <- matrix(states_vec, nrow = n_time, ncol = n_vars, byrow = FALSE)
colnames(states_mat) <- c("sea_ice", "bhatta_snowfem_tanmale", 
                          "log_hybrid_abundance", "log_snow_abundance")

#extract variables
obs_ice   <- data[, "sea_ice"]
obs_overlap <- data[, "bhatta_snowfem_tanmale"]
obs_hybrid <- data[, "log_hybrid_abundance"]
obs_snow  <- data[, "log_snow_abundance"]

#latent vrs observed: hybrid abundance
df_plot <- data.frame(
  time = 1:nrow(data),
  observed = obs_hybrid,
  latent = states_mat[, "log_hybrid_abundance"])

ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = observed), color = "black", linewidth = 0.8) +
  geom_line(aes(y = latent), color = "blue", linewidth = 0.8) +
  theme_bw() +
  labs(y = "Log hybrid abundance",
       title = "Observed vs Latent State")
#this is probably expected since we fixed SD?
#ie model is estimating trend vrs latent process noise 

#residual vrs fitted: hybrid abundance
res_hybrid <- obs_hybrid - states_mat[, "log_hybrid_abundance"]

df_res <- data.frame(
  fitted = states_mat[, "log_hybrid_abundance"],
  res = res_hybrid)

ggplot(df_res, aes(x = fitted, y = res)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw()
#looks good 

#-------------------------------------------------#
#Plot final hybrid model fit ----
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
plot_fit(data, fit_dsem, start_year = 1980)

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
effect = total_effect(fit_dsem, n_lags = 7) 

# Plot total effect
ggplot( effect) + 
  geom_bar(aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )
#snow->hybrid strong negative delayed effect, and biological responses appear to be
#delayed and cumulative 

#relative importance of variables as predictors of hybrid abundance
#i.e. at each time step, where proportion of hybrid crab variability is coming from
partition_variance(fit_dsem,
                   which_response = "log_hybrid_abundance",
                   n_times = 10 )
#hybrid crab variance primarily driven by internal population dynamics at short lags, 
#but indirect ice->snow->hybrid pathway has a strong influence, and accumulates over time

#plot
var_df <- as.data.frame(partition_variance(fit_dsem,
                                           which_response = "log_hybrid_abundance",
                                           n_times = 10)$proportion_variance_explained)

var_df$time <- 1:nrow(var_df)

# Pivot longer
var_long <- var_df %>%
  pivot_longer(cols = c(sea_ice, log_snow_abundance, log_hybrid_abundance),
               names_to = "Component",
               values_to = "Proportion")

# Plot
ggplot(var_long, aes(x = time, y = Proportion, fill = Component)) +
  geom_area(alpha = 0.8, color = "black") +
  labs(x = "Time step",
       y = "Proportion of variance explained",
       fill = "Component",
       title = "Variance partitioning of hybrid crab dynamics") +
  theme_minimal()

#second view
ggplot(var_long, aes(x = time, y = Proportion, color = Component)) +
  geom_line(size = 1.2) +
  geom_point() +
  labs(x = "Time step",
       y = "Proportion of hybrid crab variance explained",
       title = "hybrid crab variance contributions over time") +
  theme_minimal()

#plot at final time step 10 
var_df %>%
  filter(time == 10) %>%
  pivot_longer(cols = c(sea_ice, log_snow_abundance, log_hybrid_abundance),
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
#at time step 10, AR explains 36%, snow abundance explains 49% and sea ice
#explains 13% of hybrid crab abundance

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
    estimate_delta0 = FALSE,
    control = dsem_control(
      parameters = pars_tmp,
      map = map_tmp,
      quiet = TRUE,
      getsd = TRUE))
  
  sm <- as.data.frame(summary(fit_tmp))
  
  # extract key pathways
  effects <- sm %>%
    filter(name %in% c("icetohybrid", "icetosnow", "snowtohybrid")) %>%
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

#Both ice->snow and snow->hybrid pathways are robust to changes in process variance.
  #BUT we're getting convergence warnings with higher fixed values, suggested that
  #model is more unstable when more noise is included 

#-------------------------------------------------#
# Monte Carlo Simulation-based validation ----
#-------------------------------------------------#

#Here, we'll simulate data from final hybrid dsem model, refit the model to
#each simulated dataset, and compare parameter estimates to true values to check
#whether the fitted model can recover its own parameters when new data is simulated

#function modified from J. Bigman ATF ESP 

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

#extract residuals from final hybrid model fit
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
#suggest that simulation residuals have a higher autocorrelation than the model 
  #predicts and there's still autocorrelation left in residuals - not as concerning
  #since we're not using simulations for inference, but suggests that model is 
  #not capturing full temporal dependence of the system (to be expected)

#residuals vrs fit plot
ggplot() +
  geom_point(data = all_res_df %>% filter(type == "Simulated"),
             aes(fitted, resid),alpha = 0.05, color = "grey") +
  geom_point(data = all_res_df %>% filter(type == "Original"),
             aes(fitted, resid), color = "red", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Residual vs Fitted: Simulated (grey) vs Original (red)")





















