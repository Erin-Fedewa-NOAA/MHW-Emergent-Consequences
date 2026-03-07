#Assess cross-correlations between causal relationships for testing putative lags

#data contains slamon, pink, sockeye, coho, driver

library(tidyverse)
library(zoo)

#-----------------------------
# Example structure of dataset
#-----------------------------
# df should contain:
# year, pink, sockeye, coho, driver

# df <- read_csv("your_data.csv")

#-----------------------------
# Helper function: lagged correlation
#-----------------------------
lag_cor <- function(x, y, lags = -4:4) {
  tibble(
    lag = lags,
    cor = map_dbl(lags, ~cor(dplyr::lag(x, .x), y, use = "complete.obs"))
  )
}

#-----------------------------
# Create rolling means
#-----------------------------
df_roll <- df %>%
  arrange(year) %>%
  mutate(
    driver_1yr = driver,
    driver_2yr = rollmean(driver, 2, fill = NA, align = "right"),
    driver_3yr = rollmean(driver, 3, fill = NA, align = "right")
  )

#-----------------------------
# Pivot salmon species
#-----------------------------
salmon_long <- df_roll %>%
  pivot_longer(
    cols = c(pink, sockeye, coho),
    names_to = "species",
    values_to = "abundance"
  )

#-----------------------------
# Define time periods
#-----------------------------
salmon_long <- salmon_long %>%
  mutate(
    period = case_when(
      year >= 1965 & year <= 1988 ~ "1965–1988",
      year >= 1989 & year <= 2013 ~ "1989–2013"
    )
  )

#-----------------------------
# Compute cross correlations
#-----------------------------
ccf_results <- salmon_long %>%
  pivot_longer(
    cols = starts_with("driver_"),
    names_to = "window",
    values_to = "driver_value"
  ) %>%
  group_by(species, period, window) %>%
  group_modify(~lag_cor(.x$driver_value, .x$abundance)) %>%
  ungroup()

# Clean window labels
ccf_results <- ccf_results %>%
  mutate(
    window = recode(
      window,
      driver_1yr = "One-year mean",
      driver_2yr = "Two-year rolling mean",
      driver_3yr = "Three-year rolling mean"
    )
  )

#-----------------------------
# Plot
#-----------------------------
ggplot(ccf_results,
       aes(x = lag, y = cor, fill = window)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(species ~ period) +
  scale_fill_manual(
    values = c(
      "One-year mean" = "#E69F00",
      "Two-year rolling mean" = "#56B4E9",
      "Three-year rolling mean" = "#009E73"
    )
  ) +
  labs(
    x = "Lag (years)",
    y = "Pearson correlation",
    fill = ""
  ) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey85"),
    strip.background = element_rect(fill = "grey85"),
    legend.position = "top"
  )