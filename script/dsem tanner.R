#GOAL: Identify causal links between the marine heatwave, collapse and tanner
  #crab abundance

#Author: EJF

#load 
library(dsem)
library(ggplot2)
library(dplyr)
library(dagitty)
library(ggdag)

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")

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
  select(-year) %>%
  mutate(across(c(1:3), ~ (.-mean(.,na.rm=T))/sd(.,na.rm=T), .names = "z_{.col}")) %>%
  select(-sea_ice, -female_size_structure,-large_male_abundance) %>%
  rename(sea_ice = z_sea_ice, female_size_structure = z_female_size_structure, 
         large_male_abundance = z_large_male_abundance) -> fem_data

##########################################
#First let's use a DAG to test for conditional independencies and validate our
#model structure 

#download specified DAG from dagitty.net - conditional independencies are 
#identified based on structure of DAG drawn in dagitty
tanner_dag <- dagitty('dag {
sea_ice [exposure,pos="-0.268,1.129"]
snow_abundance [pos="0.020,0.913"]
tanner_abundance [outcome,pos="0.178,1.194"]
tanner_spatial_extent [pos="-0.033,1.342"]
sea_ice -> snow_abundance
sea_ice -> tanner_abundance
sea_ice -> tanner_spatial_extent
snow_abundance -> tanner_abundance
tanner_spatial_extent -> tanner_abundance
}
')

#plot DAG
ggdag(tanner_dag, layout = "nicely")

ggdag_status(tanner_dag, text = FALSE, use_labels = "name") +
  #guides(color = "none") +  # Turn off legend
  theme_dag()

#identify paths
paths(tanner_dag)
#3 open causal pathways 

#and plot paths
ggdag_paths(tanner_dag) +
  theme_dag()

#find adjustment sets
adjustmentSets(tanner_dag)
#This tells us that no covariate adjustment is necessary to identify the causal effect
#of marine heatwave on tanner crab abundance, ie there are no open backdoor paths 

#and visualize adjustment sets, if there are any
ggdag_adjustment_set(tanner_dag, shadow = TRUE) +
  theme_dag()

#find conditional independencies- i.e. two variables that are implied to 
# be independent and not correlated shouldn't be connected by a node
impliedConditionalIndependencies(fem_dag)
#because this is empty (ie no independencies exist, localTests()
#below will produce no results)

# evaluate the d-separation implications of our DAG with our simulated dataset
#i.e. test for correlation between conditional independencies that are assumed
#based on DAG
test <- localTests(tanner_dag, fem_data, type="cis")
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