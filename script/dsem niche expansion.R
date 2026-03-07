#Testing niche expansion as a causal mechanism to explain hybrid increase
  #Here we're characterizing niche expansion using realized thermal niche
  #and hybrid population spatial extent 

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

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
thermal <- read.csv("./output/thermal_niche.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")
area <- read.csv("./output/area_occupied_output.csv")

#########################################################
#Join covariates and response
sea_ice %>%
  rename_with(tolower) %>%
  rename(sea_ice = ice_avg) %>%
  filter(year >= 1988) %>%
  full_join(thermal %>%
              rename_with(tolower) %>%
              rename(thermal_niche = temp_occ)) %>%
  full_join(crab_abund %>%
              filter(category == "population_subset") %>%
              select(year, abundance) %>%
              rename(hybrid_abundance=abundance)) %>%
  full_join(area %>%
              rename_with(tolower) %>%
              select(year, hybrid_area)) %>%
  drop_na() -> dat

#plot 
dat %>%
  pivot_longer(2:5, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#Lets look at distributions of covariates
plot_histo <- function(data) {
  # Use imap to iterate over columns (values) and their names (keys)
  plots <- imap(data, ~ggplot(data, aes(x = .x)) +
                  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
                  labs(title = paste("Histogram of", .y), x = .y, y = "Frequency") +
                  theme_minimal())
  # Return the list of plots
  return(plots)
}

plot_histo(dat)

#Normalize and standardize covariates
dat %>%
  mutate(log_hybrid_abundance = log(hybrid_abundance)) %>%
  select(-hybrid_abundance) %>%
  mutate(across(c(2:5), ~ (.-mean(.,na.rm=T))/sd(.,na.rm=T), .names = "z_{.col}")) %>%
  select(-sea_ice, -thermal_niche, -hybrid_area, -log_hybrid_abundance) %>%
  rename(sea_ice = z_sea_ice, thermal_niche = z_thermal_niche, 
         hybrid_area = z_hybrid_area,
         hybrid_abundance = z_log_hybrid_abundance) -> hybrid_data


#Assess collinearity b/w variables
#if correlated, we should have a node in our DAG below
hybrid_data %>% 
  select(-year) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(method="number")

#Now look at distributions of normalized/standardized covariates
plot_histo(hybrid_data)
#looks much better

#plot z-scored variables
hybrid_data %>%
  pivot_longer(2:5, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#########################################################
#First let's use a DAG to test for conditional independencies and validate our
#model structure 

#download specified DAG from dagitty.net - conditional independencies are 
#identified based on structure of DAG drawn in dagitty
niche_dag <- dagitty('dag {
hybrid_abundance [outcome,pos="0.722,-0.101"]
hybrid_spatial_extent [pos="-0.346,-1.003"]
realized_thermal_niche [pos="-0.337,0.257"]
sea_ice [exposure,pos="-1.541,-0.241"]
hybrid_spatial_extent -> hybrid_abundance
realized_thermal_niche -> hybrid_abundance
sea_ice -> hybrid_spatial_extent
sea_ice -> realized_thermal_niche
}
')

#plot DAG
ggdag(niche_dag, layout = "nicely") +
  theme_dag()

plot(niche_dag) 

ggdag_status(niche_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(niche_dag)
#2 open causal pathways 

#and plot paths
ggdag_paths(niche_dag, text = FALSE, use_labels = "name") +
  theme_dag()

#find adjustment sets for response variable 
adjustmentSets(niche_dag, exposure="sea_ice", outcome="hybrid_abundance")
#This tells us that no covariate adjustment is necessary to identify the causal effect
#of marine heatwave on tanner abundance, ie there are no open backdoor paths 

#and visualize adjustment sets, if there are any
ggdag_adjustment_set(niche_dag, shadow = TRUE) +
  theme_dag()

#find conditional independencies- i.e. two variables that are implied to 
# be independent and not correlated shouldn't be connected by a node
impliedConditionalIndependencies(niche_dag)

############################################
#Now since our DAG looks good, let's specify a SEM

#remove year
hybrid_data <- hybrid_data %>% select(-year)
data = ts(hybrid_data)

#define formulation: link, lag, param_name, start_value (see ?make_dsem_ram())
sem ="
#AR1 for each variable (change 1 to 0 for iid)
 #Estimate in output shows correlation (if very small, might not need AR1)
sea_ice -> sea_ice, 1, ar1 
thermal_niche -> thermal_niche, 1, ar2
hybrid_area -> hybrid_area, 1, ar3
hybrid_abundance -> hybrid_abundance, 1, ar4 

#causal links with mechanistic lags
sea_ice -> hybrid_area, 0, icetoarea
sea_ice -> thermal_niche, 0, icetothermal
hybrid_area -> hybrid_abundance, 0, areatohybrid
thermal_niche -> hybrid_abundance, 0, thermaltohybrid
"
#define family
family <- rep('fixed', ncol(hybrid_data))
## family <- c('normal', 'normal', 'normal', 'fixed', 'fixed', 'fixed')

#control section
control <- dsem_control(use_REML=FALSE, run_model=TRUE,
                        getsd=TRUE, trace=100, newton_loops=0)


## initial first model run without delta0 (to improve starting values)
fit0 <- dsem(sem=sem, 
             tsdata=data, 
             family=family,
             estimate_delta0=FALSE, 
             control = dsem_control(
               quiet = FALSE,
               getsd = FALSE))

parameters = fit0$obj$env$parList()
parameters$delta0_j = rep( 0, ncol(data) )


## refit model with delta0
fit_dsem <- dsem(sem=sem, tsdata=data, 
                 family=family,
                 estimate_delta0=TRUE,
                 control=dsem_control(quiet=TRUE,
                                      parameters = parameters))

summary(fit_dsem)
fit_dsem$opt$par
ParHat = fit_dsem$obj$env$parList()
knitr::kable( summary(fit_dsem), digits=3 )

#plot dsem output
plot_fit(data, fit_dsem)
#Need to revise function: also see https://james-thorson-noaa.github.io/dsem/articles/features.html
#for timeseries plot and dag with effects 

# sample-based quantile residuals
samples = loo_residuals(fit_dsem, what="samples", track_progress=FALSE)
which_use = which(!is.na(hybrid_data))
fitResp = loo_residuals( fit_dsem, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]

# Build and display DHARMa object
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(hybrid_data)[which_use],
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

#D-separation test
test_dsep(fit_dsem) #null hypothesis- the data came from the model- 
#hmmm, need to investigate this, we don't have missing data  
