#GOAL: Identify causal links between the marine heatwave, collapse and tanner
  #crab abundance

#Author: EJF

#load 
library(dsem)
library(ggplot2)
library(dplyr)
library(dagitty)
library(ggdag)
library(knitr)
library(corrplot)

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")
area <- read.csv("./output/area_occupied_output.csv")

#########################################################
#Join covariates and response, Z score variables 
sea_ice %>%
  rename_with(tolower) %>%
  rename(sea_ice = ice_avg) %>%
  filter(year >= 1988) %>%
  full_join(area %>%
              rename_with(tolower) %>%
              rename(tanner_spatial_extent = tanner_area) %>%
              select(-snow_area)) %>%
  full_join(crab_abund %>%
              filter(category == "population" & species %in% c("snow", "tanner")) %>%
              select(year, abundance, species) %>%
              pivot_wider(names_from = "species", values_from = "abundance") %>%
              rename(snow_abundance=snow, tanner_abundance=tanner)) %>%
  drop_na() %>%
  mutate(across(c(2:5), ~ (.-mean(.,na.rm=T))/sd(.,na.rm=T), .names = "z_{.col}")) %>%
  select(-sea_ice, -tanner_spatial_extent, -snow_abundance, -tanner_abundance) %>%
  rename(sea_ice = z_sea_ice, tanner_spatial_extent = z_tanner_spatial_extent, 
         snow_abundance = z_snow_abundance,
         tanner_abundance = z_tanner_abundance) -> tanner_data

#Assess collinearity b/w variables
#if correlated, we should have a node in our DAG below
tanner_data %>% 
  select(-year) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(method="number")

#plot z-scored variables
tanner_data %>%
  pivot_longer(2:5, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#visually explore lags
ccf(tanner_data$snow_abundance, tanner_data$tanner_abundance)
ccf(tanner_data$snow_abundance, tanner_data$tanner_spatial_extent)
ccf(tanner_data$sea_ice, tanner_data$tanner_spatial_extent)
ccf(tanner_data$sea_ice, tanner_data$tanner_abundance)
ccf(tanner_data$tanner_spatial_extent, tanner_data$tanner_abundance)

##########################################
#First let's use a DAG to test for conditional independencies and validate our
#model structure 

#download specified DAG from dagitty.net - conditional independencies are 
#identified based on structure of DAG drawn in dagitty
tanner_dag <- dagitty('dag {
sea_ice [exposure,pos="-0.185,1.150"]
snow_abundance [pos="-0.008,1.032"]
tanner_abundance [outcome,pos="0.145,1.184"]
tanner_spatial_extent [pos="-0.012,1.338"]
sea_ice -> snow_abundance
sea_ice -> tanner_abundance
sea_ice -> tanner_spatial_extent
snow_abundance -> tanner_abundance
snow_abundance -> tanner_spatial_extent
tanner_spatial_extent -> tanner_abundance
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
#5 open causal pathways 

#and plot paths
ggdag_paths(tanner_dag, text = FALSE, use_labels = "name") +
  theme_dag()

#find adjustment sets
adjustmentSets(tanner_dag, exposure="sea_ice", outcome="tanner_abundance")
#This tells us that no covariate adjustment is necessary to identify the causal effect
#of marine heatwave on tanner abundance, ie there are no open backdoor paths 

#and visualize adjustment sets, if there are any
ggdag_adjustment_set(tanner_dag, shadow = TRUE) +
  theme_dag()

#find conditional independencies- i.e. two variables that are implied to 
# be independent and not correlated shouldn't be connected by a node
impliedConditionalIndependencies(tanner_dag)
#because this is empty (ie no independencies exist, localTests()
#below will produce no results)

# evaluate the d-separation implications of our DAG with our simulated dataset
#i.e. test for correlation between conditional independencies that are assumed
#based on DAG
test <- localTests(tanner_dag, tanner_data, type="cis")
#Note that this is probably not a useful tool for ecology since we have MUCH less
#data, and probably working with time series that have missing data

# perform Holm-Bonferrino correction to mitigate problems around multiple testing 
test$p.value <- p.adjust(test$p.value) 

test # should show all p values above 0.05, suggesting DAG-data consistency- aka you don't need 
#to reject DAG based on results of test
#If you fail data consistency test, you need to go back and rethink DAG - there is either a 
#common cause that you're not accounting for (or can't account for if latent), 
#or a spurious correlation with no mechanism- but we shouldn't be held to this, data shouldn't be informing DAG!

############################################
#Now since our DAG looks good, let's specify a SEM

#remove year
tanner_data <- tanner_data %>% select(-year)
data = ts(tanner_data)

#define formulation: link, lag, param_name, start_value (see ?make_dsem_ram())

sem ="
#AR1 for each variable
 #Estimate in output shows correlation (if very small, might not need AR1)
sea_ice -> sea_ice, 1, ar1 
tanner_abundance -> tanner_abundance, 1, ar2
snow_abundance -> snow_abundance, 1, ar3
tanner_spatial_extent -> tanner_spatial_extent, 1, ar4

#causal links with mechanistic lags
sea_ice -> snow_abundance, 1, icetosnow
sea_ice -> tanner_abundance, 1, icetotanner
sea_ice -> tanner_spatial_extent, 1, icetoarea
snow_abundance -> tanner_abundance, 3, snowtotanner
snow_abundance -> tanner_spatial_extent, 2, snowtoarea
tanner_spatial_extent -> tanner_abundance, 2, areatotanner
"
#define family
family <- rep('fixed', ncol(tanner_data))
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
fit_dsem_tanner <- dsem(sem=sem, tsdata=data, 
                     family=family,
                     estimate_delta0=TRUE,
                     control=dsem_control(quiet=TRUE,
                                          parameters = parameters))

summary(fit_dsem_tanner)
fit_dsem_tanner$opt$par
ParHat = fit_dsem_tanner$obj$env$parList()
knitr::kable( summary(fit_dsem_tanner), digits=3 )

#plot dsem output
plot_fit(data, fit_dsem_tanner)
#Need to revise function: also see https://james-thorson-noaa.github.io/dsem/articles/features.html
#for timeseries plot and dag with effects 

# sample-based quantile residuals
samples = loo_residuals(fit_dsem_tanner, what="samples", track_progress=FALSE)
which_use = which(!is.na(tanner_data))
fitResp = loo_residuals( fit_dsem_tanner, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]

# Build and display DHARMa object
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(tanner_data)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)
#looks good

# Calculate total effects
effect = total_effect(fit_dsem_tanner)

# Plot total effect
ggplot( effect) + 
  geom_bar(aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )

#relative importance of variables as predictors of female size structure
partition_variance(fit_dsem_tanner,
                   which_response = "female_size_structure",
                   n_times = 10 )
#hmm maybe this is in development version of dsem only? 

#D-separation test
test_dsep(fit_dsem_tanner) #null hypothesis- the data came from the model- 
#hmmm, need to investigate this, we don't have missing data  


