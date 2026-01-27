#Now that we've built complexity with female, tanner and hybrid DSEM models, 
  #let's fit a full model integrating all responses of interest 

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
overlap <- read.csv("./output/Emily output/bhatt_snowTanner.csv")

#########################################################

full_dag <- dagitty('dag {
female_size_structure [outcome,pos="-0.068,1.088"]
hybrid_abundance [outcome,pos="0.018,1.201"]
sea_ice [exposure,pos="-0.076,1.195"]
snow_abundance [pos="-0.026,1.106"]
snow_tanner_overlap [pos="-0.021,1.179"]
tanner_abundance [outcome,pos="-0.020,1.293"]
sea_ice -> female_size_structure
sea_ice -> hybrid_abundance [pos="0.004,1.254"]
sea_ice -> snow_abundance
sea_ice -> snow_tanner_overlap
sea_ice -> tanner_abundance
snow_abundance -> female_size_structure
snow_abundance -> hybrid_abundance [pos="0.037,1.112"]
snow_abundance -> snow_tanner_overlap
snow_abundance -> tanner_abundance [pos="-0.065,1.097"]
snow_tanner_overlap -> hybrid_abundance
tanner_abundance -> hybrid_abundance [pos="0.026,1.255"]
}
')

#plot DAG
ggdag(full_dag, layout = "nicely") +
  theme_dag()

plot(full_dag) 
