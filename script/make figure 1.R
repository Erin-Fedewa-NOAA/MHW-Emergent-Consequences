#Goals ----
#Create Figure 1 for paper 

#need to revise crab abundance script to loop thru cpue calc function for error bars
#





So the code you want will:
  
  Filter to males
Bin or subset by size (~4 inches ≈ ~100 mm carapace width)
Identify hybrids vs non-hybrids
Compute:
  
  
  df %>%
  filter(sex == "male") %>%
  mutate(size_class = cut(width_mm, breaks = ...)) %>%
  group_by(size_class) %>%
  summarise(
    prop_hybrid = sum(hybrid == 1) / n()
  )

#--------------------------------#
#Lag effect size figures ----
#--------------------------------#

#Data
tanner_lags <- read.csv("./output/tanner_lags.csv")
hybrid_lags <- read.csv("./output/hybrid_lags.csv")

# Color palette (colorblind-friendly)
lag_colors <- c(
  "Lag 1" = "#E69F00",
  "Lag 2" = "#56B4E9",
  "Lag 3" = "#009E73",
  "Lag 4" = "#F0E442",
  "Lag 5" = "#0072B2",
  "Lag 6" = "#D55E00")


# Tanner Plot
ggplot(tanner_lags, aes(x = Causal_Pathway, y = Estimate, fill = lag)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error),
                position = position_dodge(width = 0.7),
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = lag_colors) +
  labs(x = "",
       y = "Effect Size (± Std Error)",
       fill = "Lag") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10))













#option 2 faceted plot:
tanner_lags_plot <- tanner_lags %>%
  mutate(Causal_Pathway = factor(Causal_Pathway,
                                 levels = unique(Causal_Pathway)))

ggplot(tanner_lags_plot, aes(x = lag, y = Estimate, fill = lag)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Estimate - Std_Error,
                    ymax = Estimate + Std_Error),
                width = 0.15, size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_fill_manual(values = lag_colors, name = "Lag") +
  facet_wrap(~Causal_Pathway, nrow = 1) +  # single row = cleaner comparison
  labs(x = "", y = "Effect size (± SE)") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(face = "bold", size = 14, hjust = 0))