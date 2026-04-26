#Assess cross-correlations between causal relationships for initial  
  #evaluation of biologically-constrained lags

#Script Authors: E. Fedewa, M.Litzow, E. Ryznar

library(tidyverse)
library(zoo)
library(janitor)
library(gt)
library(forecast)

#read in data
sea_ice <- read.csv("./output/seaice_output.csv")
crab_abund <- read.csv("./output/crab_abundance.csv")
overlap <- read.csv("./output/overlap_output.csv")

#-------------------------------------
# Define covariates and scale
#-------------------------------------

dat <- sea_ice %>%
  select(year, Mar_Apr_ice_EBS_NBS) %>%
  filter(year >= 1988) %>%
  rename(sea_ice = Mar_Apr_ice_EBS_NBS) %>%
  full_join(overlap %>%
              select(-bhatta_global)) %>%
  full_join(crab_abund %>%
              filter(category == "population_50mm_plus") %>%
              select(year, abundance, species) %>%
              pivot_wider(names_from = "species", values_from = "abundance") %>%
              rename(snow_abundance=snow, tanner_abundance=tanner,
                     hybrid_abundance=hybrid))  

#Normalize and standardize covariates

#First define covariates used in the model
vars <- c("sea_ice",
          "bhatta_snowfem_tanmale",
          "snow_abundance",
          "tanner_abundance",
          "bhatta_tannerfem_snowmale")

#And center covariates
hybrid_data <- dat %>%
  mutate(log_hybrid_abundance = log(hybrid_abundance)) %>%
  mutate(across(all_of(c(vars)), ~ as.numeric(scale(.)))) %>%
  select(-year, -hybrid_abundance) 

#----------------------------------------------------------
# Define causal links from full DSEM model for assessing lags
#---------------------------------------------------------

links <- tibble::tribble(
  ~driver, ~response,
  "sea_ice","log_hybrid_abundance",
  "sea_ice","snow_abundance",
  "sea_ice","bhatta_snowfem_tanmale",
  "sea_ice","bhatta_tannerfem_snowmale",
  "sea_ice","tanner_abundance",
  "snow_abundance","log_hybrid_abundance",
  "snow_abundance","tanner_abundance",
  "bhatta_snowfem_tanmale","log_hybrid_abundance",
  "bhatta_tannerfem_snowmale","log_hybrid_abundance",
  "tanner_abundance","log_hybrid_abundance")

#----------------------------------------------------------
# Define biologically-constrained lags for each causal link
#---------------------------------------------------------

bio_lags <- tibble::tribble(
  ~driver, ~response, ~min_lag, ~max_lag,
  "sea_ice","log_hybrid_abundance", -5, -7,
  "sea_ice","snow_abundance", -1, -2,
  "sea_ice","bhatta_snowfem_tanmale",-5, -7, 
  "sea_ice","bhatta_tannerfem_snowmale",-5, -7,
  "sea_ice","tanner_abundance",-5, -7, 
  "snow_abundance","log_hybrid_abundance",-1, -3, 
  "snow_abundance","tanner_abundance",-1, -3, 
  "bhatta_snowfem_tanmale","log_hybrid_abundance",-5, -7, 
  "bhatta_tannerfem_snowmale","log_hybrid_abundance",-5, -7, 
  "tanner_abundance","log_hybrid_abundance", -1, -3) %>%
  mutate(link = paste(driver, response, sep = " → "))

#-----------------------------
# Functions: rolling means and cross-correlations
#-----------------------------
roll_mean <- function(x, k){
  zoo::rollmean(x, k = k, fill = NA, align = "right")
}

lag_cor <- function(x, y, lag){
  #negative lag shifts driver forward to time t + lag,
    #e.g. compare earlier values of driver with later values of response 
    if(lag < 0){
    cor(x[1:(length(x)+lag)],
        y[(1-lag):length(y)],
        use="complete.obs")
  #positive lag shifts driver backwards to time t - lag
  } else if(lag > 0){
    cor(x[(1+lag):length(x)],
        y[1:(length(y)-lag)],
        use="complete.obs")
    
  } else{
    cor(x,y,use="complete.obs")
  }
}

#-------------------------------------------------------------
# Compute cross correlations at 1, 2 and 3-year rolling means
#-------------------------------------------------------------
lags <- -7:0 #only negative lags are biologically meaningful here

results <- links %>%
  mutate(res = purrr::map2(driver,response,function(d,r){
    
    x <- hybrid_data[[d]]
    y <- hybrid_data[[r]]
    
    tibble(lag = lags,
      
      `One-year mean` =
        purrr::map_dbl(lags, ~lag_cor(x,y,.x)),
      
      `Two-year rolling mean` =
        purrr::map_dbl(lags, ~lag_cor(x, roll_mean(y,2), .x)),
      
      `Three-year rolling mean` =
        purrr::map_dbl(lags, ~lag_cor(x, roll_mean(y,3), .x))) %>%
      
      pivot_longer(-lag, names_to="mean_type", values_to="correlation")
  })) %>%
  tidyr::unnest(res) %>%
  mutate(link = paste(driver,"→",response))

#-----------------------------
# Plot
#-----------------------------
#order rolling means correctly
results$mean_type <- factor(results$mean_type,
                            levels = c("One-year mean",
                                       "Two-year rolling mean",
                                       "Three-year rolling mean"))
#and plot
results %>%
  ggplot(aes(x = lag, y = correlation, fill = mean_type)) +
  geom_col(width = .8, position = position_dodge(width = .8)) +
  #add lag constraints 
  geom_rect(data = bio_lags,
            aes(xmin = min_lag, xmax = max_lag, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.08,
            inherit.aes = FALSE) +
  facet_wrap(~link, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = -6:6) +
  scale_fill_manual(values = c(
    "One-year mean" = "#E69F00",
    "Two-year rolling mean" = "#56B4E9",
    "Three-year rolling mean" = "#009E73")) +
  labs(x = "Lag (years)",
       y = "Pearson correlation",
       fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        strip.background = element_rect(fill="grey90"),
        strip.text = element_text(face="bold"),
        panel.grid.minor = element_blank())

#-----------------------------------------------
# Select top negative lags for each causal effect
#-----------------------------------------------
#select top negative lag by one-year/rolling means 
  #(driver leads response by t minus lag years)
best_negative_lags <- results %>%
  filter(lag < 0) %>% # only negative lags
  group_by(link, mean_type) %>%
  slice_max(correlation, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(link, mean_type)

#table of results
best_negative_lags %>%
  select(link, mean_type, lag, correlation) %>%
  # Combine lag and correlation into a single string for neater table
  mutate(lag_corr = paste0(lag, " (", round(correlation, 2), ")")) %>%
  select(link, mean_type, lag_corr) %>%
  pivot_wider(names_from = mean_type, values_from = lag_corr) %>%
  gt() %>%
  tab_header(
    title = "Best Negative Lags by Rolling Mean") %>%
  cols_label(
    `One-year mean` = "1-Year Rolling Mean",
    `Two-year rolling mean` = "2-Year Rolling Mean",
    `Three-year rolling mean` = "3-Year Rolling Mean",
    link = "Link Type") %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    table.border.top.color = "black",
    table.border.bottom.color = "black")

#plot best 3 negative lags for one-year mean only (no rolling means)
best3_negative_lags <- results %>%
  filter(mean_type == "One-year mean",
         lag < 0) %>%
  group_by(link) %>%
  slice_max(correlation, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(link, mean_type) %>%
  group_by(link) %>%
  mutate(is_best = if_else(row_number(desc(correlation)) == 1, TRUE, FALSE)) %>%
  ungroup()

# Plot with different color for the best lag
best3_negative_lags %>%
  ggplot(aes(x = lag, y = correlation, fill = is_best)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#D55E00",   # best lag in dark orange/red
                               "FALSE" = "#E69F00")) +  # other bars in base orange
  facet_wrap(~link, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  labs(title = "Best 3 Negative Lags (One-Year Mean Only)",
    x = "Best Leading Lag (years)",
    y = "Pearson Correlation",
    fill = "Top Lag") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    legend.position = "top")

#-----------------------------------------------
# Pre-whitened cross-correlation function
#-----------------------------------------------
#Because autocorrelation in time series can create spurious lag peaks, we can
  #also evaluate lags after removing autocorrelation in the driver by fitting
  #an ARIMA model to the driver and computing ccfs on the residuals

prewhiten_ccf <- function(driver, response, data, lag_max = 7){
  
  x <- data[[driver]]
  y <- data[[response]]
  
  # Fit ARIMA to driver
  fit <- forecast::auto.arima(x)
  
  # Residuals of driver
  x_res <- residuals(fit)
  
  # Extract phi (AR coefficient) to remove AR structure from response variable
  phi <- fit$model$phi
  
  # Handle case where there are NO AR terms, i.e. there's no AR structure
  if (length(phi) == 0) {
    y_filt <- as.numeric(y)
  } else {
    y_filt <- as.numeric(stats::filter(y,
                                       filter = phi,
                                       method = "recursive"))
  }
  
  # Remove NA pairs
  keep <- complete.cases(x_res, y_filt)
  x_res <- x_res[keep]
  y_filt <- y_filt[keep]
  
  if (length(x_res) < 3) return(NULL)
  
  # --- manual cross-correlation ---
  lags <- -lag_max:0
  
  cor_vals <- sapply(lags, function(l){
    if (l < 0) {
      cor(x_res[1:(length(x_res)+l)],
          y_filt[(1-l):length(y_filt)])
    } else if (l > 0) {
      cor(x_res[(1+l):length(x_res)],
          y_filt[1:(length(y_filt)-l)])
    } else {
      cor(x_res, y_filt)
    }
  })
  
  tibble(
    lag = lags,
    correlation = cor_vals
  )
}
#-----------------------------------------------
# Run function across causal links from above
#-----------------------------------------------

prewhite_results <- links %>%
  mutate(res = map2(driver,response,
                    ~prewhiten_ccf(.x,.y,hybrid_data))) %>%
  tidyr::unnest(res) %>%
  mutate(link = paste(driver,"→",response))

#------------------------------------------------------
# plot to identify lag relationships that persist after
  #removing autocorrelation
#------------------------------------------------------
ggplot(prewhite_results,
       aes(lag, correlation)) +
  #biologically constrained lag window
  geom_rect(data = bio_lags,
            aes(xmin = min_lag, xmax = max_lag, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.08,
            inherit.aes = FALSE) +
  geom_col(fill="steelblue") +
  facet_wrap(~link, scales="free_y") +
  theme_bw() +
  labs(
    x="Lag (years)",
    y="Prewhitened correlation")

