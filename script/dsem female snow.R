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
  ParHat <- fit$obj$env$parList()
  out <- lapply(1:ncol(data), function(i) {
    tmp <- data.frame(year=1988+1:nrow(data), variable=colnames(data)[i], 
                      obs=data[,i], pred=ParHat$x_tj[,i] )
    SD = as.list(fit$sdrep,what="Std.")$x_tj[,i]
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
  facet_wrap(~variable, scales = "free_y") +
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
ggdag(fem_dag, layout = "nicely")

ggdag_status(fem_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(fem_dag)
#2 open causal pathways 

#and plot paths
ggdag_paths(fem_dag) +
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

#define formulation: link, lag, param_name, start_value

sem ="
sea_ice -> female_size_structure, 1, icefem
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
#quantile deviations detected, though qq plot looks good 

# Calculate total effects
effect = total_effect(fit_dsem_fem)

# Plot total effect
ggplot( effect) + 
  geom_bar(aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )











Z = ts( bering_sea[,c('log_CP', 'log_Cfall','log_Esummer'), drop=FALSE] )
# ar1 - should model each variable as an ar1 process in DSEM! 
sem <- '
 # AR1 for each variable
log_CP -> log_CP,1,ar1 #variable vrs using residuals to get at causal effect 
  #Estimate in output shows correlation (if very small, might not need AR1)
log_Cfall -> log_Cfall,1,ar2
log_Esummer -> log_Esummer,1,ar3

 # causal links - this is a fork example
log_CP -> log_Cfall, 0, CP_to_Cfall
log_CP -> log_Esummer, 0, CP_to_E
'
fit <- dsem(sem=sem, tsdata=Z) 
fit$opt$par #3 means and 8 betas (variance and correlation for 
# each variable + 2 causal effects- lines 4 and 5)
summary(fit)
plot_fit(Z, fit) #in this case, log_cfall is driven by strength of causal rxn to log_CP (even 
#though no data for the variable)

#d-seperation test 
test_dsep(fit) #null hypothesis- the data came from the model- here we pass our 
#data-consistency test







# quick fucntion to plot dsem output
plot_fit <- function(Y, fit){
  ParHat <- fit$obj$env$parList()
  out <- lapply(1:ncol(Y), function(i) {
    tmp <- data.frame(year=1959+1:nrow(Y), variable=colnames(Y)[i], 
                      obs=Y[,i], pred=ParHat$x_tj[,i] )
    SD = as.list(fit$sdrep,what="Std.")$x_tj[,i]
    cbind( tmp, "lower"=tmp$pred - ifelse(is.na(SD),0,SD),
           "upper"=tmp$pred + ifelse(is.na(SD),0,SD) )
  }) |> bind_rows() 
  out$variable <- factor(out$variable, levels=colnames(Y))
  g <- ggplot(out, aes(x=year, y=pred, ymin=lower, ymax=upper)) + geom_line() + 
    facet_wrap('variable', scales='free', dir='v') + geom_ribbon(alpha=.2) 
  g + geom_point(data=na.omit(out), mapping=aes(y=obs), col='red')
}











# start with just a single variable, log cold pool extent
data(bering_sea) #y is the data in the paper
Z <- ts(bering_sea[,'log_CP', drop=FALSE])

# iid
sem <- 'log_CP <-> log_CP,0,sd' #DAG is single variable here (arrow and lag notation)
fit <- dsem(sem=sem, tsdata=Z)
fit$opt$par
summary(fit)
# note that the means of the variables are hidden using the default REML approach
fit$obj$env$parList()$mu_j
plot_fit(Z, fit)

# ar1 (new dimension including time, cold pool causes change in cold pool the next yr)
sem <- 'log_CP -> log_CP,1,ar1'
fit <- dsem(sem=sem, tsdata=Z)
fit$opt$par
summary(fit)
plot_fit(Z, fit)
fit$obj$env$parList()$mu_j 

# random walk, force correlation of AR1 to be 1 by naming it "NA"
sem <- 'log_CP -> log_CP,1,NA,1'
fit <- dsem(sem=sem, tsdata=Z)
fit$opt$par
summary(fit)
plot_fit(Z, fit)

# ar2
sem <- '
log_CP -> log_CP,1,ar1
log_CP -> log_CP,2,ar2
'
fit <- dsem(sem=sem, tsdata=Z[,'log_CP', drop=FALSE])
fit$opt$par
summary(fit)
plot_fit(Z, fit)

# ar1 with observation error estimated - family = normal, family = fixed means no obsv. error
#ln_sigma_j is observation error in log space
sem <- 'log_CP -> log_CP,1,ar1'
fit <- dsem(sem=sem, tsdata=Z, family='normal')
fit$opt$par
summary(fit)
plot_fit(Z, fit)

# ar1 with fixed observation error
# first build model without estimating it
fit <- dsem(sem=sem, tsdata=Z, family='normal', control=dsem_control(run_model=FALSE))
# now modify TMB inputs to specify a fixed lnsigma
fit$tmb_inputs |> str()
pars <- fit$tmb_inputs$parameters
pars$lnsigma_j <- log(.1)
map <- fit$tmb_inputs$map
map$lnsigma_j <- factor(NA)
fit <- dsem(sem=sem, tsdata=Z, family='normal', 
            control=dsem_control(map=map, parameters = pars,
                                 getJointPrecision = TRUE))
fit$opt$par
summary(fit)
plot_fit(Z, fit) #process error is zero here if .1 is used


# now some simple causal models
Z = ts( bering_sea[,c('log_CP', 'log_Cfall','log_Esummer'), drop=FALSE] )
# ar1 - should model each variable as an ar1 process in DSEM! 
sem <- '
 # AR1 for each variable
log_CP -> log_CP,1,ar1 #variable vrs using residuals to get at causal effect 
  #Estimate in output shows correlation (if very small, might not need AR1)
log_Cfall -> log_Cfall,1,ar2
log_Esummer -> log_Esummer,1,ar3

 # causal links - this is a fork example
log_CP -> log_Cfall, 0, CP_to_Cfall
log_CP -> log_Esummer, 0, CP_to_E
'
fit <- dsem(sem=sem, tsdata=Z) 
fit$opt$par #3 means and 8 betas (variance and correlation for 
# each variable + 2 causal effects- lines 4 and 5)
summary(fit)
plot_fit(Z, fit) #in this case, log_cfall is driven by strength of causal rxn to log_CP (even 
#though no data for the variable)

# Specify model
sem = "
  # Link, lag, param_name
  log_seaice -> log_CP, 0, seaice_to_CP
  log_CP -> log_Cfall, 0, CP_to_Cfall
  log_CP -> log_Esummer, 0, CP_to_E
  log_PercentEuph -> log_RperS, 0, Seuph_to_RperS
  log_PercentCop -> log_RperS, 0, Scop_to_RperS
  log_Esummer -> log_PercentEuph, 0, Esummer_to_Suph
  log_Cfall -> log_PercentCop, 0, Cfall_to_Scop
  SSB -> log_RperS, 0, SSB_to_RperS

  log_seaice -> log_seaice, 1, AR1, 0.001
  log_CP -> log_CP, 1,  AR2, 0.001
  log_Cfall -> log_Cfall, 1, AR4, 0.001
  log_Esummer -> log_Esummer, 1, AR5, 0.001
  SSB -> SSB, 1, AR6, 0.001
  log_RperS ->  log_RperS, 1, AR7, 0.001
  log_PercentEuph -> log_PercentEuph, 1, AR8, 0.001
  log_PercentCop -> log_PercentCop, 1, AR9, 0.001
"
Z <- ts(bering_sea)
fit = dsem( sem = sem,
            tsdata = Z,
            control = dsem_control(use_REML=FALSE, quiet=TRUE) )
summary(fit)
plot_fit(Z, fit)

# d-separation tests, see: 
# https://dagitty.net/dags.html?id=SmvjV5p8
?test_dsep
example('test_dsep')
test_dsep(fit) #null hypothesis- the data came from the model- here we pass our 
#data-consistency test 


#Matrix::image(fit$obj$report()$Q)
