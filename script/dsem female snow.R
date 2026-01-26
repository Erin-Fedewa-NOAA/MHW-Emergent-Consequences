#GOAL: Identify causal links between the marine heatwave, collapse and female 
  #snow crab size structure

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
female_size <- read.csv("./output/female_size.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")

#########################################################
# quick fucntion to plot dsem output
  #Author: Cole M. 
plot_fit <- function(data, fit_dsem_fem){
  ParHat <- fit_dsem_fem$obj$env$parList()
  out <- lapply(1:ncol(data), function(i) {
    tmp <- data.frame(year=1988+1:nrow(data), variable=colnames(data)[i], 
                      obs=data[,i], pred=ParHat$x_tj[,i] )
    SD = as.list(fit_dsem_fem$sdrep,what="Std.")$x_tj[,i]
    cbind( tmp, "lower"=tmp$pred - ifelse(is.na(SD),0,SD),
           "upper"=tmp$pred + ifelse(is.na(SD),0,SD) )
  }) |> bind_rows() 
  out$variable <- factor(out$variable, levels=colnames(data))
  g <- ggplot(out, aes(x=year, y=pred, ymin=lower, ymax=upper)) + geom_line() + 
    facet_wrap('variable', scales='free', dir='v') + geom_ribbon(alpha=.2) 
  g + geom_point(data=na.omit(out), mapping=aes(y=obs), col='red')
}


######################################################

#Join covariates and response, Z score variables 
sea_ice %>%
  rename_with(tolower) %>%
  rename(sea_ice = ice_avg) %>%
  filter(year >= 1988) %>%
  full_join(female_size %>%
              rename_with(tolower) %>%
              rename(female_size_structure = proportion_large)) %>%
  full_join(crab_abund %>%
              filter(category == "large_male") %>%
              select(year, abundance) %>%
              rename(large_male_abundance = abundance)) %>%
  drop_na() %>%
  mutate(across(c(2:4), ~ (.-mean(.,na.rm=T))/sd(.,na.rm=T), .names = "z_{.col}")) %>%
  select(-sea_ice, -female_size_structure,-large_male_abundance) %>%
  rename(sea_ice = z_sea_ice, female_size_structure = z_female_size_structure, 
          large_male_abundance = z_large_male_abundance) -> fem_data

#Assess collinearity b/w variables
  #if correlated, we should have a node in our DAG below
fem_data %>% 
  select(-year) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(method="number")

#plot z-scored variables
fem_data %>%
  pivot_longer(2:4, names_to="variable", values_to="value") %>%
  ggplot(aes(year, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow=3) +
  theme_bw()

#visually explore lags
ccf(fem_data$female_size_structure, fem_data$sea_ice)
ccf(fem_data$female_size_structure, fem_data$large_male_abundance)
ccf(fem_data$sea_ice, fem_data$large_male_abundance)

##########################################
#First let's use a DAG to test for conditional independencies and validate our
  #model structure 

#download specified DAG from dagitty.net - conditional independencies are 
  #identified based on structure of DAG drawn in dagitty
fem_dag <- dagitty('dag {
female_size_structure [outcome,pos="0.335,1.306"]
large_male_abundance [pos="0.020,0.913"]
sea_ice [exposure,pos="-0.292,1.378"]
large_male_abundance -> female_size_structure
sea_ice -> female_size_structure
sea_ice -> large_male_abundance
}
')

#plot DAG
ggdag(fem_dag, layout = "nicely") +
  theme_dag()

ggdag_status(fem_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(fem_dag)
#2 open causal pathways 

#and plot paths
ggdag_paths(fem_dag, text = FALSE, use_labels = "name") +
  theme_dag()

#find adjustment sets
adjustmentSets(fem_dag, exposure="Marine heatwave", outcome="Female size structure")
#This tells us that no covariate adjustment is necessary to identify the causal effect
  #of marine heatwave on female size structure, ie there are no open backdoor paths 

#and visualize adjustment sets, if there are any
ggdag_adjustment_set(fem_dag, shadow = TRUE) +
  theme_dag()

#find conditional independencies- i.e. two variables that are implied to 
  # be independent and not correlated shouldn't be connected by a node
impliedConditionalIndependencies(fem_dag)
#because this is empty (ie no independencies exist, localTests()
  #below will produce no results)

# evaluate the d-separation implications of our DAG with our simulated dataset
#i.e. test for correlation between conditional independencies that are assumed
#based on DAG
test <- localTests(fem_dag, fem_data, type="cis")
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
#Need to use AIC to compare full, AR1 and iid models 

#remove year
fem_data <- fem_data %>% select(-year)
data = ts(fem_data)

#define formulation: link, lag, param_name, start_value (see ?make_dsem_ram())

sem ="
#AR1 for each variable
 #Estimate in output shows correlation (if very small, might not need AR1)
sea_ice -> sea_ice, 1, ar1 
female_size_structure -> female_size_structure, 1, ar2
large_male_abundance -> large_male_abundance, 1, ar3

#causal links with mechanistic lags
sea_ice -> female_size_structure, 2, icefem
sea_ice -> large_male_abundance,  1, icelgmale
large_male_abundance -> female_size_structure,  1, lgmalefem
"
#define family
family <- rep('fixed', ncol(fem_data))
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
fit_dsem_fem <- dsem(sem=sem, tsdata=data, 
                     family=family,
                     estimate_delta0=TRUE,
                     control=dsem_control(quiet=TRUE,
                                          parameters = parameters))

summary(fit_dsem_fem)
fit_dsem_fem$opt$par
ParHat = fit_dsem_fem$obj$env$parList()
knitr::kable( summary(fit_dsem_fem), digits=3 )

#plot dsem output
plot_fit(data, fit_dsem_fem)
#Need to revise function: also see https://james-thorson-noaa.github.io/dsem/articles/features.html
  #for timeseries plot and dag with effects 

# sample-based quantile residuals
samples = loo_residuals(fit_dsem_fem, what="samples", track_progress=FALSE)
which_use = which(!is.na(fem_data))
fitResp = loo_residuals( fit_dsem_fem, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]

# Build and display DHARMa object
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(fem_data)[which_use],
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






