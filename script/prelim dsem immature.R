#GOAL: Test a preliminary DSEM model using snow/tanner crab recruitment 

#moving average of ice, ccf plot function 

#Author: EJF

#load 
library(dsem)
library(ggplot2)
library(dplyr)
library(dagitty)
library(ggdag)
library(knitr)
library(corrplot)
library(zoo)

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")

#########################################################
#Join covariates and response, Z score variables 
sea_ice %>%
  rename_with(tolower) %>%
  rename(sea_ice = ice_avg) %>%
  filter(year >= 1985) %>%
  #calculate a two year moving average
  mutate(sea_ice_2yr = rollmean(sea_ice, k = 2, fill = NA, align = "right")) %>%
  full_join(crab_abund %>%
              filter(category == "immature") %>%
              select(year, abundance, species) %>%
              pivot_wider(names_from = species, values_from = abundance) %>%
              rename(snow_immature_abundance=snow,
                     tanner_immature_abundance=tanner)) %>%
  drop_na() %>%
  mutate(across(c(2:5), ~ (.-mean(.,na.rm=T))/sd(.,na.rm=T), .names = "z_{.col}")) %>%
  select(-sea_ice, -sea_ice_2yr, -snow_immature_abundance,-tanner_immature_abundance) %>%
  rename(sea_ice = z_sea_ice, sea_ice_2yr = z_sea_ice_2yr, 
         snow_immature_abundance = z_snow_immature_abundance,
         tanner_immature_abundance = z_tanner_immature_abundance) -> imm_data

#Assess collinearity b/w variables
#if correlated, we should have a node in our DAG below
imm_data %>% 
  select(-year) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(method="number")

#plot z-scored variables
imm_data %>%
  pivot_longer(2:5, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#visually explore lags
ccf(imm_data$snow_immature_abundance, imm_data$sea_ice)
ccf(imm_data$tanner_immature_abundance, imm_data$sea_ice) #4-5yrs
ccf(imm_data$snow_immature_abundance, imm_data$tanner_immature_abundance)

##########################################
#First let's use a DAG to test for conditional independencies and validate our
#model structure 

#download specified DAG from dagitty.net - conditional independencies are 
#identified based on structure of DAG drawn in dagitty
imm_dag <- dagitty('dag {
sea_ice [exposure,pos="-2.144,0.239"]
snow_immature_abundance [pos="-0.565,-0.773"]
tanner_immature_abundance [outcome,pos="-0.669,0.930"]
sea_ice -> snow_immature_abundance
sea_ice -> tanner_immature_abundance
snow_immature_abundance -> tanner_immature_abundance
}
')

#plot DAG
ggdag(imm_dag, layout = "nicely") +
  theme_dag()
plot(imm_dag)

ggdag_status(imm_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(imm_dag)
#2 open causal pathways 

############################################
#Now since our DAG looks good, let's specify a SEM
#Need to use AIC to compare full, AR1 and iid models 

#remove year
imm_data <- imm_data %>% select(-year)
data = ts(imm_data)

#define formulation: link, lag, param_name, start_value (see ?make_dsem_ram())

sem ="
#AR1 for each variable
 #Estimate in output shows correlation (if very small, might not need AR1)
sea_ice -> sea_ice, 1, ar1 
snow_immature_abundance -> snow_immature_abundance, 1, ar2
tanner_immature_abundance -> tanner_immature_abundance, 1, ar3

#causal links with mechanistic lags
sea_ice -> snow_immature_abundance, 3, icesnow
sea_ice -> tanner_immature_abundance,  4, icetanner
snow_immature_abundance -> tanner_immature_abundance,  4, snowtanner
"
#define family
family <- rep('fixed', ncol(imm_data))
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
fit_dsem_imm <- dsem(sem=sem, tsdata=data, 
                     family=family,
                     estimate_delta0=TRUE,
                     control=dsem_control(quiet=TRUE,
                                          parameters = parameters))

summary(fit_dsem_imm)

#plot dsem output
plot_fit(data, fit_dsem_fem)
#Need to revise function: also see https://james-thorson-noaa.github.io/dsem/articles/features.html
#for timeseries plot and dag with effects 

# sample-based quantile residuals
samples = loo_residuals(fit_dsem_fem, what="samples", track_progress=FALSE)
which_use = which(!is.na(imm_data))
fitResp = loo_residuals( fit_dsem_fem, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]

# Build and display DHARMa object
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(imm_data)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)
#looks good

# Calculate total effects
effect = total_effect(fit_dsem_fem)

# Plot total effect
ggplot( effect) + 
  geom_bar(aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )

#relative importance of variables as predictors of female size structure
partition_variance(fit_dsem_fem,
                   which_response = "female_size_structure",
                   n_times = 10 )
#hmm maybe this is in development version of dsem only? 

#D-separation test
test_dsep(fit_dsem_fem) #null hypothesis- the data came from the model- 
#hmmm, need to investigate this, we don't have missing data  






