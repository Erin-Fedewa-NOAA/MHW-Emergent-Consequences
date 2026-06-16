#Goals ----
#Create Figure 2 for manuscript 

#Author: EJF

## Load packages
library(RColorBrewer)
library(tidyverse)
library(crabpack)

#Data
tanner_lags <- read.csv("./output/tanner_lags.csv")
hybrid_lags <- read.csv("./output/hybrid_lags.csv")

#--------------------------------#
#Lag effect figures ----
#--------------------------------#

# Tanner Plot
tanner_lag <- ggplot(tanner_lags, aes(x = Causal_Pathway, y = Estimate, fill = lag)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, color = "black",linewidth = 0.35) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error),
                position = position_dodge(width = 0.65),
                width = 0.15, linewidth = 0.45) +
  geom_hline(yintercept = 0, color = "black",linewidth = 0.3) +
  scale_fill_manual(values = c("#F7FBFF","#DEEBF7","#C6DBEF",
                               "#9ECAE1","#6BAED6","#3182BD","#08519C")) +
  labs(x = "", y = "Estimate (± SE)", fill = "Lag") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size=9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank())

#hybrid plot
hybrid_lag <- ggplot(hybrid_lags, aes(x = Causal_Pathway, y = Estimate, fill = lag)) +
  geom_col(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.65, color = "black",linewidth = 0.35) +
  geom_errorbar(aes(ymin = Estimate - Std_Error, ymax = Estimate + Std_Error),
                position = position_dodge(width = 0.65),
                width = 0.15, linewidth = 0.45) +
  geom_hline(yintercept = 0, color = "black",linewidth = 0.3) +
  scale_fill_manual(values = c("#F7FBFF","#DEEBF7","#C6DBEF",
                               "#9ECAE1","#6BAED6","#3182BD","#08519C")) +
  labs(x = "", y = "Estimate (± SE)", fill = "Lag") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size=9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank())

#--------------------------------#
#Combine and save figures ----
#--------------------------------#

blank_panel <- ggplot() + theme_void() +
  theme(panel.border = element_rect(
      color = NA, fill = NA))

fig2 <- ((blank_panel + tanner_lag) /
  (blank_panel + hybrid_lag)) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0, 1.04))

ggsave("./figures/fig2_tiff.tiff", fig2,
       width = 6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")

ggsave("./figures/fig2_pdf.pdf", fig2,
       width = 6.5, height = 4.5, units = "in", dpi = 600)
