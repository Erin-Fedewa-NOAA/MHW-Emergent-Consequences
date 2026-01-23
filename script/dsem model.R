




#first off, combine all timeseries


library(dsem)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)

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

# explore simulation a bit
?simulate.dsem
Ysim <- simulate(fit, nsim=10) |> bind_cols()
matplot(Ysim)
Ysim <- simulate(fit, nsim=10, fill_missing = TRUE) |> bind_cols()
matplot(Ysim)
Ysim <- simulate(fit, nsim=10, variance = 'both') |> bind_cols() #random draw from NAs in 
#early timeseries 
matplot(Ysim)
Ysim <- simulate(fit, nsim=10, variance = 'both', fill_missing = TRUE) |> bind_cols()
matplot(Ysim)
#would use this simulation to generate time series data given a DAG to test 
Ysim <- simulate(fit, nsim=10, fill_missing = TRUE, resimulate_gmrf = TRUE) |> bind_cols()
matplot(Ysim)

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
