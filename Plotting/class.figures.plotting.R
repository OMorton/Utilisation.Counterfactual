library(tidyverse)
library(gratia)
library(ggpubr)

source("functions.R")
# Load master dataset
load(file = "Data/data_models.RData")

## Data Prep -------------------------------------------------------------------
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



## Plotting --------------------------------------------------------------------

actin_att_chained <- readRDS(paste0("Outputs/Draws/Class/Actinopteri_draws.v2.csv"))
actin_chained_sum <- read.csv(paste0("Outputs/Tables/Class/Actinopteri_summary_contrast.v2.csv"))
aves_att_chained <- readRDS(paste0("Outputs/Draws/Class/Aves_draws.v2.csv"))
aves_chained_sum <- read.csv(paste0("Outputs/Tables/Class/Aves_summary_contrast.v2.csv"))
elas_att_chained <- readRDS(paste0("Outputs/Draws/Class/Elasmobranchii_draws.v2.csv"))
elas_chained_sum <- read.csv(paste0("Outputs/Tables/Class/Elasmobranchii_summary_contrast.v2.csv"))
mam_att_chained <- readRDS(paste0("Outputs/Draws/Class/Mammalia_draws.v2.csv"))
mam_chained_sum <- read.csv(paste0("Outputs/Tables/Class/Mammalia_summary_contrast.v2.csv"))
rept_att_chained <- readRDS(paste0("Outputs/Draws/Class/Reptilia_draws.v2.csv"))
rept_chained_sum <- read.csv(paste0("Outputs/Tables/Class/Reptilia_summary_contrast.v2.csv"))

actin_plt <- class_att_plt(draws = actin_att_chained, summary = actin_chained_sum, common.y = T, taxa.cols = T)
aves_plt <- class_att_plt(draws = aves_att_chained, summary = aves_chained_sum, common.y = T, taxa.cols = T)
elas_plt <- class_att_plt(draws = elas_att_chained, summary = elas_chained_sum, common.y = T, taxa.cols = T)
mam_plt <- class_att_plt(draws = mam_att_chained, summary = mam_chained_sum, common.y = T, taxa.cols = T)
rept_plt <- class_att_plt(draws = rept_att_chained, summary = rept_chained_sum, common.y = T, taxa.cols = T)

## overall
draws_att<-read_rds( "Outputs/Draws/counterfac.draws.inter.rds")

att_chained <- draws_att %>%
  mutate(year = as.numeric(as.character(fyear))) %>%
  select(series, year, .draw, arm, .fitted) %>%
  pivot_wider(names_from = arm, values_from = .fitted) %>%
  # Log treatment difference per series per year (REs cancel)
  mutate(log_effect = factual - counterfactual) %>%
  # Calculate year-on-year incremental change (delta) per series
  group_by(series, .draw) %>%
  arrange(year) %>%
  mutate(
    # Annual change in treatment effect: delta_t = effect_t - effect_{t-1}
    inc_log_effect = log_effect - lag(log_effect)) %>%
  ungroup() %>%
  # drop NA increments (first year of each series)
  filter(!is.na(inc_log_effect)) %>%
  # Average annual increments across ALL series active in year t
  group_by(year, .draw) %>%
  summarise(mean_inc_att = mean(inc_log_effect, na.rm = TRUE), .groups = "drop") %>%
  # Chain-link by cumulatively summing annual increments from 1970
  group_by(.draw) %>%
  arrange(year) %>%
  mutate(cum_log_att = cumsum(ifelse(year == 1970, 0, mean_inc_att)),
         # Exponentiate to get the global chained ATT ratio (1970 = 1.0)
         att_ratio = exp(cum_log_att)) %>%
  ungroup() %>%
  # Filter boundary years
  filter(year >= 1970)

att_chained_sum <- att_chained %>% group_by(year) %>%
  summarise(att_ratio_mn = mean(att_ratio),
            att_ratio_lo  = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
            att_ratio_hi  = quantile(att_ratio, probs = 0.975, na.rm = TRUE))

write.csv(att_chained_sum, "Outputs/Tables/overall.counterfac.sum.csv")

chained_plt<- ggplot(att_chained, aes(x = year, y = att_ratio)) +
  geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = .3) +
  geom_ribbon(data = att_chained_sum, 
              aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = att_chained_sum, aes(year, att_ratio_mn), colour = "black", linewidth = .7) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
  coord_cartesian(ylim = c(0, 1.3), expand = F) +
  labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  theme_bw()

taxa.arr <- ggarrange(ggarrange(chained_plt,
        # right side 2 cols
          ggarrange(
            ggarrange(actin_plt$bar, actin_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
            ggarrange(aves_plt$bar, aves_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
            ncol = 1, labels = c("B.", "C."), label.x = -0.02),
          ncol = 2, widths = c(2, 1), labels = c("A."), label.x = 0),
        # bottom row
        ggarrange(
          ggarrange(elas_plt$bar, elas_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
          ggarrange(mam_plt$bar, mam_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
          ggarrange(rept_plt$bar, rept_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
          ncol = 3, nrow = 1, labels = c("D.", "E.", "F."), label.x = -0.02),
      nrow = 2, heights = c(2, 1))


library(cowplot)
taxa.arr2 <- ggdraw(taxa.arr) +
  draw_plot(taxa.arr) +
  draw_image("Outputs/Misc/Gadus_morhua.png",
             x = 0.93, y = 0.88, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Phalacrocorax_capensis.png",
             x = 0.93, y = 0.54, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Prionace_glauca.png",
             x = 0.26, y = 0.21, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Canis_lupus.png",
             x = 0.59, y = 0.21, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Chelonia_mydas.png",
             x = 0.94, y = 0.21, width = 0.05, height = 0.05) # tweak place and size due to orientation

ggsave(path = "Outputs/Figures/Tidy", filename = "Class.Est.fyear.v2.png",
       plot = taxa.arr2, width = 8, height = 6, units = "in",bg = "white")

## Alternate log scale 
taxa.arr.log <- ggarrange(ggarrange(chained_plt +scale_y_log10() + coord_cartesian(ylim = c(0.03, 5), expand = F),
                                # right side 2 cols
                                ggarrange(
                                  ggarrange(actin_plt$bar, actin_plt$att +scale_y_log10() + coord_cartesian(ylim = c(0.03, 5), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
                                  ggarrange(aves_plt$bar, aves_plt$att +scale_y_log10() + coord_cartesian(ylim = c(0.03, 5), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
                                  ncol = 1, labels = c("B.", "C."), label.x = -0.02),
                                ncol = 2, widths = c(2, 1), labels = c("A."), label.x = 0),
                      # bottom row
                      ggarrange(
                        ggarrange(elas_plt$bar, elas_plt$att +scale_y_log10() + coord_cartesian(ylim = c(0.03, 5), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
                        ggarrange(mam_plt$bar, mam_plt$att +scale_y_log10() + coord_cartesian(ylim = c(0.03, 5), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
                        ggarrange(rept_plt$bar, rept_plt$att +scale_y_log10() + coord_cartesian(ylim = c(0.03, 5), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
                        ncol = 3, nrow = 1, labels = c("D.", "E.", "F."), label.x = -0.02),
                      nrow = 2, heights = c(2, 1))

taxa.arr.log2 <- ggdraw(taxa.arr.log) +
  draw_plot(taxa.arr.log) +
  draw_image("Outputs/Misc/Gadus_morhua.png",
             x = 0.93, y = 0.88, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Phalacrocorax_capensis.png",
             x = 0.93, y = 0.54, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Prionace_glauca.png",
             x = 0.26, y = 0.21, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Canis_lupus.png",
             x = 0.59, y = 0.21, width = 0.06, height = 0.06) +
  draw_image("Outputs/Misc/Chelonia_mydas.png",
             x = 0.94, y = 0.21, width = 0.05, height = 0.05) # tweak place and size due to orientation

ggsave(path = "Outputs/Figures/Tidy", filename = "Class.Est.fyear.log.v2.png",
       plot = taxa.arr.log2, width = 8, height = 6, units = "in",bg = "white")

## Numbers reported in the main text --------------------------------------------

att_chained_sum %>% filter(year == 2024) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

mod_dat_full %>% filter(used == 1) %>% group_by(Class) %>% summarise(n_distinct(series))
mod_dat_full %>% filter(used == 0) %>% group_by(Class) %>% summarise(n_distinct(series))
mod_dat_full %>% group_by(Class) %>% summarise(n_distinct(series))

actin_chained_sum %>% filter(att_ratio_hi < 1) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

aves_chained_sum %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

elas_chained_sum %>% filter(year == 2023) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

mam_chained_sum  %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)
