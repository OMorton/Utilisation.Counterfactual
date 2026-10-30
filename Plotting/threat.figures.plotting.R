library(tidyverse)
library(gratia)
library(ggpubr)

source("functions.R")
# Load master dataset
load(file = "Data/data_models.RData")
#draws_att<-read_rds( "Outputs/Draws/counterfac.draws.inter.rds")

## Threat coverage summary -----------------------------------------------------
threat_vars <- c("threat_habitat_loss_change_f", "threat_pollution_f",
  "threat_disease_f",  "threat_climate_change_f",   "threat_invasive_spp_genes_f")
threat_cols <- c("used", "threat_habitat_loss_change", "threat_pollution", 
                 "threat_disease", "threat_climate_change", "threat_invasive_spp_genes")

mod_dat_full <- mod_dat_full %>%
  mutate(series = factor(series), fyear = factor(fyear), Site = factor(Site),
         across(all_of(threat_cols), ~ factor(.x), .names = "{.col}_f"),
         across(all_of(threat_cols), ~ ordered(.x, levels = c("0", "1")), .names = "{.col}_ord"),
         used_num = ifelse(used_f == "1", 1, 0),
         inter_luc_num = used_num * ifelse(threat_habitat_loss_change_f == "1", 1, 0),
         inter_pol_num = used_num * ifelse(threat_pollution_f == "1", 1, 0),
         inter_dis_num = used_num * ifelse(threat_disease_f == "1", 1, 0),
         inter_cli_num = used_num * ifelse(threat_climate_change_f == "1", 1, 0),
         inter_inv_num = used_num * ifelse(threat_invasive_spp_genes_f == "1", 1, 0)) %>%
  arrange(Class, series, year.z)

no_threats_count <- mod_dat_full %>%
  filter(used ==1,
    threat_habitat_loss_change_f == "0" & threat_pollution_f == "0" &
      threat_disease_f == "0" & threat_climate_change_f == "0" &
      threat_invasive_spp_genes_f == "0") %>%
  summarise(threat_clean = "No Threats", n_unique_series = n_distinct(series))

has_threats_counts <- mod_dat_full %>% filter(used ==1) %>%
  pivot_longer(cols = all_of(threat_vars),  names_to  = "threat",
    values_to = "present") %>%
  filter(present == "1") %>%
  mutate(threat_clean = case_when(
      threat == "threat_habitat_loss_change_f" ~ "Land-Use Change (LUC)",
      threat == "threat_pollution_f" ~ "Pollution (Pol)",
      threat == "threat_disease_f"  ~ "Disease (Dis)",
      threat == "threat_climate_change_f" ~ "Climate Change (CC)",
      threat == "threat_invasive_spp_genes_f"  ~ "Invasive Species (Inv)")) %>%
  group_by(threat_clean) %>%
  summarise(n_unique_series = n_distinct(series),.groups = "drop")

all_series_counts <- bind_rows(no_threats_count, has_threats_counts) %>%
  arrange(desc(n_unique_series))

## not used
no_threats_count_use0 <- mod_dat_full %>%
  filter(used ==0,
         threat_habitat_loss_change_f == "0" & threat_pollution_f == "0" &
           threat_disease_f == "0" & threat_climate_change_f == "0" &
           threat_invasive_spp_genes_f == "0") %>%
  summarise(threat_clean = "No Threats", n_unique_series = n_distinct(series))

has_threats_counts_use0 <- mod_dat_full %>% filter(used ==0) %>%
  pivot_longer(cols = all_of(threat_vars),  names_to  = "threat",
               values_to = "present") %>%
  filter(present == "1") %>%
  mutate(threat_clean = case_when(
    threat == "threat_habitat_loss_change_f" ~ "Land-Use Change (LUC)",
    threat == "threat_pollution_f" ~ "Pollution (Pol)",
    threat == "threat_disease_f"  ~ "Disease (Dis)",
    threat == "threat_climate_change_f" ~ "Climate Change (CC)",
    threat == "threat_invasive_spp_genes_f"  ~ "Invasive Species (Inv)")) %>%
  group_by(threat_clean) %>%
  summarise(n_unique_series = n_distinct(series),.groups = "drop")

all_series_counts_use0 <- bind_rows(no_threats_count_use0, has_threats_counts_use0) %>%
  arrange(desc(n_unique_series))

max((mod_dat_full %>% filter(threat_invasive_spp_genes == 1, used == 1))$Year) #2024
max((mod_dat_full %>% filter(threat_disease == 1, used == 1))$Year) # 2014
max((mod_dat_full %>% filter(threat_climate_change == 1, used == 1))$Year) # 2023
max((mod_dat_full %>% filter(threat_pollution == 1, used == 1))$Year) # 2024
max((mod_dat_full %>% filter(threat_habitat_loss_change == 1, used == 1))$Year) #2024

## Marginal effects of use accross co-occuring threats -------------------------
draws_inter <- read_rds("Outputs/Draws/ave.draws.inter.rds")

contrast_draws <- draws_inter %>%
  select(scenario, .draw, year, used_f, rate_rel, rate_rel_index) %>%
  pivot_wider(
    names_from = used_f,
    values_from = c(rate_rel, rate_rel_index),
    names_prefix = "used_") %>%
  mutate(
    # Log-Scale Relative Difference: (r_t,used - r_1970,used) - (r_t,notused - r_1970,notused)
    diff_rate_rel = rate_rel_used_1 - rate_rel_used_0,
    # Ratio of Relative Indices: Index_used_t / Index_notused_t = exp(diff_rate_rel)
    ratio_rel_index = rate_rel_index_used_1 / rate_rel_index_used_0,
    # Percentage Difference relative to Not Used performance
    pct_diff = (ratio_rel_index - 1) * 100) %>%
  ungroup() %>%
  filter(scenario %in% c("Pol", "LUC", "Inv", "No Threats")|(scenario == "CC" & year < 2024)|(scenario == "Dis" & year < 2015))

contrast_sum <- contrast_draws %>% group_by(scenario, year) %>%
  summarise(rel_ratio_mn = mean(ratio_rel_index),
            rel_ratio_lo  = quantile(ratio_rel_index, probs = 0.025, na.rm = TRUE),
            rel_ratio_hi  = quantile(ratio_rel_index, probs = 0.975, na.rm = TRUE))


# NO THREATS
nt_plt <- ggplot(filter(contrast_draws, scenario == "No Threats"), aes(x = year, y = ratio_rel_index)) +
  geom_line(aes(group = .draw), alpha = 0.08, color = "#4393c3") +
  geom_ribbon(data = filter(contrast_sum, scenario == "No Threats"), 
              aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = filter(contrast_sum, scenario == "No Threats"), aes(year, rel_ratio_mn), colour = "#4393c3", linewidth = .7) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
  coord_cartesian(xlim= c(1970, 2024), ylim = c(0, 3), expand = F) +
  scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  labs(x = "Year", y = "Index") +
  theme_bw() +
  theme(strip.background =  element_blank(), legend.title = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

# LUC
luc_plt <- ggplot(filter(contrast_draws, scenario == "LUC"), aes(x = year, y = ratio_rel_index)) +
  geom_line(aes(group = .draw), alpha = 0.08, color = "#8c510a") +
  geom_ribbon(data = filter(contrast_sum, scenario == "LUC"), 
              aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = filter(contrast_sum, scenario == "LUC"), aes(year, rel_ratio_mn), colour = "#8c510a", linewidth = .5) +
  geom_line(data = filter(contrast_sum, scenario == "No Threats"), aes(year, rel_ratio_mn), colour = "#4393c3", linewidth = .7) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
  coord_cartesian(xlim= c(1970, 2024), ylim = c(0, 3), expand = F) +
   scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  labs(x = "Year", y = "Index") +
  theme_bw() +
  theme(strip.background =  element_blank(), legend.title = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))
 
 # CC
cc_plt <- ggplot(filter(contrast_draws, scenario == "CC"), aes(x = year, y = ratio_rel_index)) +
   geom_line(aes(group = .draw), alpha = 0.08, color = "#8c510a") +
   geom_ribbon(data = filter(contrast_sum, scenario == "CC"), 
               aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
               fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "CC"), aes(year, rel_ratio_mn), colour = "#8c510a", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "No Threats"), aes(year, rel_ratio_mn), colour = "#4393c3", linewidth = .7) +
   geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
   coord_cartesian(xlim= c(1970, 2024), ylim = c(0, 3), expand = F) +
   scale_x_continuous(breaks = c(1980, 2000, 2020)) +
   labs(x = "Year", y = "Index") +
   theme_bw() +
   theme(strip.background =  element_blank(), legend.title = element_blank(),
         legend.position = "bottom", strip.text = element_text(face = "bold"))
 
 # Pol
pol_plt <- ggplot(filter(contrast_draws, scenario == "Pol"), aes(x = year, y = ratio_rel_index)) +
   geom_line(aes(group = .draw), alpha = 0.08, color = "#8c510a") +
   geom_ribbon(data = filter(contrast_sum, scenario == "Pol"), 
               aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
               fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "Pol"), aes(year, rel_ratio_mn), colour = "#8c510a", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "No Threats"), aes(year, rel_ratio_mn), colour = "#4393c3", linewidth = .7) +
   geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
   coord_cartesian(xlim= c(1970, 2024), ylim = c(0, 3), expand = F) +
   scale_x_continuous(breaks = c(1980, 2000, 2020)) +
   labs(x = "Year", y = "Index") +
   theme_bw() +
   theme(strip.background =  element_blank(), legend.title = element_blank(),
         legend.position = "bottom", strip.text = element_text(face = "bold"))
 
 ## Dis
dis_plt <- ggplot(filter(contrast_draws, scenario == "Dis"), aes(x = year, y = ratio_rel_index)) +
   geom_line(aes(group = .draw), alpha = 0.08, color = "#8c510a") +
   geom_ribbon(data = filter(contrast_sum, scenario == "Dis"), 
               aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
               fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "Dis"), aes(year, rel_ratio_mn), colour = "#8c510a", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "No Threats"), aes(year, rel_ratio_mn), colour = "#4393c3", linewidth = .7) +
   geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
   coord_cartesian(xlim= c(1970, 2024), ylim = c(0, 3), expand = F) +
   scale_x_continuous(breaks = c(1980, 2000, 2020)) +
   labs(x = "Year", y = "Index") +
   theme_bw() +
   theme(strip.background =  element_blank(), legend.title = element_blank(),
         legend.position = "bottom", strip.text = element_text(face = "bold"))
 
 ## Inv
inv_plt <- ggplot(filter(contrast_draws, scenario == "Inv"), aes(x = year, y = ratio_rel_index)) +
   geom_line(aes(group = .draw), alpha = 0.08, color = "#8c510a") +
   geom_ribbon(data = filter(contrast_sum, scenario == "Inv"), 
               aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
               fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "Inv"), aes(year, rel_ratio_mn), colour = "#8c510a", linewidth = .5) +
   geom_line(data = filter(contrast_sum, scenario == "No Threats"), aes(year, rel_ratio_mn), colour = "#4393c3", linewidth = .7) +
   geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
   coord_cartesian(xlim= c(1970, 2024), ylim = c(0, 45), expand = F) +
   scale_x_continuous(breaks = c(1980, 2000, 2020)) +
   labs(x = "Year", y = "Index") +
   theme_bw() +
   theme(strip.background =  element_blank(), legend.title = element_blank(),
         legend.position = "bottom", strip.text = element_text(face = "bold"))

thr.arr <- ggarrange(
  ggarrange(NULL, nt_plt, NULL, widths = c(1, 3, 1), nrow = 1, ncol = 3,
            labels = c("", "A.", "")),
  ggarrange(luc_plt, pol_plt,  dis_plt,cc_plt, inv_plt, nrow = 1, ncol = 5, 
            labels = c("B.", "C.", "D.", "E.", "F."), label.x = -0.02),
  nrow = 2, heights = c(1, 1))

library(cowplot)
thr.arr2 <- ggdraw(thr.arr) +
  draw_plot(thr.arr) +
  draw_image("Outputs/Misc/digger.png",
             x = 0.04, y = 0.40, width = 0.09, height = 0.09) + 
  draw_image("Outputs/Misc/pollution.png",
             x = 0.24, y = 0.41, width = 0.06, height = 0.06) + 
  draw_image("Outputs/Misc/disease.png",
             x = 0.44, y = 0.40, width = 0.08, height = 0.08) + 
  draw_image("Outputs/Misc/thermometer.png",
             x = 0.63, y = 0.40, width = 0.08, height = 0.08) + 
  draw_image("Outputs/Misc/ant.png",
             x = 0.86, y = 0.41, width = 0.05, height = 0.05) +
  draw_label("No threats", x = 0.36, y = 0.94, size = 14, fontface = "bold")

ggsave(path = "Outputs/Figures/Tidy", filename = "Threat.Est.v1.png",
       plot = thr.arr2, width = 8, height = 6, units = "in",bg = "white")


## v2
thr.arr2 <- ggdraw(thr.arr) +
  draw_plot(thr.arr) +
  draw_label("LUC",
             x = 0.05, y = 0.45, size = 9, fontface = "bold", hjust = 0, colour = "#8c510a") + 
  draw_label("Climate change",
             x = 0.25, y = 0.45, size = 9, fontface = "bold", hjust = 0, colour = "#8c510a") + 
  draw_label("Pollution",
             x = 0.45, y = 0.45, size = 9, fontface = "bold", hjust = 0, colour = "#8c510a") + 
  draw_label("Disease",
             x = 0.65, y = 0.45, size = 9, fontface = "bold", hjust = 0, colour = "#8c510a") + 
  draw_label("Invasion",
             x = 0.86, y = 0.45, size = 9, fontface = "bold", hjust = 0, colour = "#8c510a") +
  draw_label("No threats", x = 0.36, y = 0.94, size = 14, fontface = "bold", colour = "#4393c3")

ggsave(path = "Outputs/Figures/Tidy", filename = "Threat.Est.v2.png",
       plot = thr.arr2, width = 8, height = 6, units = "in",bg = "white")

## Alternate format
thr.arr <- ggarrange(ggarrange(nt_plt,
                                # right side 2 rows
                                ggarrange(luc_plt, pol_plt,
                                  ncol = 1, labels = c("B.", "C."), label.x = -0.02),
                                ncol = 2, widths = c(2, 1), labels = c("A."), label.x = 0),
                      # bottom row
                      ggarrange(dis_plt,cc_plt, inv_plt,
                        ncol = 3, nrow = 1, labels = c("D.", "E.", "F."), label.x = -0.02),
                      nrow = 2, heights = c(2, 1))


library(cowplot)
thr.arr2 <- ggdraw(thr.arr) +
  draw_plot(thr.arr) +
  draw_label("No threats", x = 0.15, y = 0.96, size = 14, fontface = "bold") + 
  draw_image("Outputs/Misc/digger.png",
             x = 0.72, y = 0.92, width = 0.07, height = 0.07) +
  draw_image("Outputs/Misc/pollution.png",
             x = 0.72, y = 0.59, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/disease.png",
             x = 0.05, y = 0.25, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/thermometer.png",
             x = 0.37, y = 0.25, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/ant.png",
             x = 0.73, y = 0.25, width = 0.05, height = 0.05) # tweak place and size due to orientation

ggsave(path = "Outputs/Figures/Tidy", filename = "Threat.Est.v3.png",
       plot = thr.arr2, width = 8, height = 6, units = "in",bg = "white")

# alternate log
## Alternate format
thr.arr.log <- ggarrange(ggarrange(nt_plt+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F),
                               # right side 2 rows
                               ggarrange(luc_plt+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F),
                                         pol_plt+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F),
                                         ncol = 1, labels = c("B.", "C."), label.x = -0.02),
                               ncol = 2, widths = c(2, 1), labels = c("A."), label.x = 0),
                     # bottom row
                     ggarrange(dis_plt+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F),
                               cc_plt+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F),
                               inv_plt+scale_y_log10() + coord_cartesian(ylim = c(0.01, 80),xlim = c(1970, 2024), expand = F),
                               ncol = 3, nrow = 1, labels = c("D.", "E.", "F."), label.x = -0.02),
                     nrow = 2, heights = c(2, 1))


library(cowplot)
thr.arr.log2 <- ggdraw(thr.arr.log) +
  draw_plot(thr.arr.log) +
  draw_label("No threats", x = 0.15, y = 0.96, size = 14, fontface = "bold") + 
  draw_image("Outputs/Misc/digger.png",
             x = 0.73, y = 0.93, width = 0.07, height = 0.07) +
  draw_image("Outputs/Misc/pollution.png",
             x = 0.73, y = 0.61, width = 0.04, height = 0.04) +
  draw_image("Outputs/Misc/disease.png",
             x = 0.06, y = 0.26, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/thermometer.png",
             x = 0.39, y = 0.27, width = 0.05, height = 0.05) +
  draw_image("Outputs/Misc/ant.png",
             x = 0.75, y = 0.27, width = 0.04, height = 0.04) # tweak place and size due to orientation

ggsave(path = "Outputs/Figures/Tidy", filename = "Threat.Est.v3.log.png",
       plot = thr.arr.log2, width = 8, height = 6, units = "in",bg = "white")

## Explicit contrast -----------------------------------------------------------
threat_contrast_draws <- contrast_draws %>%
  # Separate "No Threats" baseline from active threat scenarios
  filter(scenario != "No Threats") %>%
  left_join(contrast_draws %>% filter(scenario == "No Threats") %>% 
      select(year, .draw, ratio_no_threats = ratio_rel_index, rate_no_threats = diff_rate_rel),
    by = c("year", ".draw")) %>%
  mutate(
    # Multiplicative ratio contrast (1.0 = threat has no additional effect beyond utilisation alone)
    # < 1.0 means utilisation causes MORE severe decline when this threat is present
    threat_synergy_ratio = ratio_rel_index / ratio_no_threats,
    # Log-scale additive difference (0.0 = no synergistic difference)
    # < 0.0 means additive exacerbation of decline
    threat_synergy_diff = diff_rate_rel - rate_no_threats)

# Summarise posterior distribution of the contrast per scenario per year
threat_contrast_sum <- threat_contrast_draws %>%
  group_by(scenario, year) %>%
  summarise(
    # Ratio contrast summary
    synergy_ratio_mn = mean(threat_synergy_ratio, na.rm = TRUE),
    synergy_ratio_lo = quantile(threat_synergy_ratio, probs = 0.025, na.rm = TRUE),
    synergy_ratio_hi = quantile(threat_synergy_ratio, probs = 0.975, na.rm = TRUE),
    prob_worse = mean(threat_synergy_ratio < 1.0, na.rm = TRUE),
    .groups = "drop")

scenarios <- unique(threat_contrast_sum$scenario)
plot_list <- lapply(scenarios, function(scen) {
  
  sum_sub   <- threat_contrast_sum %>% filter(scenario == scen)
  draws_sub <- threat_contrast_draws %>% filter(scenario == scen)
  
  ggplot(sum_sub, aes(x = year, y = synergy_ratio_mn)) +
    geom_line(data = draws_sub, aes(x = year, y = threat_synergy_ratio, group = .draw),
              colour = "#fcbba1", alpha = 0.2) +
    geom_line(colour = "#fb6a4a", linewidth = 0.8) +
    geom_ribbon(aes(ymin = synergy_ratio_lo, ymax = synergy_ratio_hi), 
                fill = NA, colour = "black", linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "solid", colour = "black", alpha = 0.7) +
    labs(x = "Year",y = "Proportionate change in \nutilisation impact") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
})

names(plot_list) <- scenarios

thr.contrast.plt <- ggarrange(plot_list[["LUC"]] + coord_cartesian(ylim = c(0,4)),
          plot_list[["Pol"]] + coord_cartesian(ylim = c(0,4)),
          plot_list[["Dis"]] + coord_cartesian(ylim = c(0,4)),
          plot_list[["CC"]] + coord_cartesian(ylim = c(0,4)), 
          plot_list[["Inv"]],
          nrow = 3, ncol = 2, labels = c("A.", "B.", "C.", "D.", "E."))


ggsave(path = "Outputs/Figures/SM", filename = "Threat.contrast.png",
       plot = thr.contrast.plt, width = 6, height = 7, units = "in",bg = "white")

## Number for text -------------------------------------------------------------

contrast_sum %>% filter(scenario == "No Threats", year == max(year)) %>%
  mutate(rel_ratio_mn = 100- rel_ratio_mn*100,
         rel_ratio_lo = 100- rel_ratio_lo*100,
         rel_ratio_hi = 100- rel_ratio_hi*100)

contrast_sum %>% filter(scenario == "LUC", year == max(year)) %>%
  mutate(rel_ratio_mn = 100- rel_ratio_mn*100,
         rel_ratio_lo = 100- rel_ratio_lo*100,
         rel_ratio_hi = 100- rel_ratio_hi*100)

contrast_sum %>% filter(scenario == "Pol", year == max(year)) %>%
  mutate(rel_ratio_mn = 100- rel_ratio_mn*100,
         rel_ratio_lo = 100- rel_ratio_lo*100,
         rel_ratio_hi = 100- rel_ratio_hi*100)

contrast_sum %>% filter(scenario == "Dis", year == max(year)) %>%
  mutate(rel_ratio_mn = 100- rel_ratio_mn*100,
         rel_ratio_lo = 100- rel_ratio_lo*100,
         rel_ratio_hi = 100- rel_ratio_hi*100)
