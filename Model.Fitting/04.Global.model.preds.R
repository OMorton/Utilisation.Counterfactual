
library(tidyverse)
library(mgcv)
library(gratia)
library(ggpubr)

load(file = "Data/data_models.RData")
glob.mod.gam.rhotuned <- readRDS("Outputs/Models/GLOB.use.threat.syn.full.v3.rds")

## Data Preparation -----------------------------------------------------------
threat_cols <- c("used", "threat_habitat_loss_change", "threat_pollution", 
                 "threat_disease", "threat_climate_change", "threat_invasive_spp_genes")

mod_dat_full <- mod_dat_full %>%
  mutate(series = factor(series), fyear = factor(fyear), Site = factor(Site),
         Class_f = as.factor(Class),
         # ordered threat factors
         across(all_of(threat_cols), ~ factor(.x), .names = "{.col}_f"),
         across(all_of(threat_cols), ~ ordered(.x, levels = c("0", "1")), .names = "{.col}_ord"),
         used_num = ifelse(used_f == "1", 1, 0),
         # Synergistic interaction multipliers (Used = 1 AND Threat = 1)
         inter_luc_num = used_num * ifelse(threat_habitat_loss_change_f == "1", 1, 0),
         inter_pol_num = used_num * ifelse(threat_pollution_f == "1", 1, 0),
         inter_dis_num = used_num * ifelse(threat_disease_f == "1", 1, 0),
         inter_cli_num = used_num * ifelse(threat_climate_change_f == "1", 1, 0),
         inter_inv_num = used_num * ifelse(threat_invasive_spp_genes_f == "1", 1, 0)) %>%
  arrange(series, year.z) %>%
  group_by(series) %>%
  mutate(AR_start = (row_number() == 1)) %>% # Marks beginning of each time-series
  ungroup()

## Chain linked counterfactual -------------------------------------------------

# Aims to combines series-level causal cancellation without data loss or survivor
# bias. First, the net treatment effect is calculated within each series per draw,
# cancelling out series-specific REs. Second, the annual step-by-step change 
# in treatment effect between consecutive years is computed for each series.
# These annual incremental changes are then averaged across all populations
# in a given year, and the annual averages are cumulatively summed from 1970
# to build a continuous log-scale index. Exponentiating this cumulative val gives
# the final geometric mean ratio index.

# Two alternate counterfactual approaches were considered and discarded in favour of the
# above approach. Methods considered are detailed below.

# (1) Annual avverage approach
# Calculates global mean trajectories by averaging model predictions 
# across all treated populations for each year, under both factual (harvested) and 
# counterfactual (unharvested) conditions. Each global trajectory is anchored to 
# its 1970 value to produce relative growth indices, and the treatment effect is 
# measured as the ratio between the two resulting curves.
# vulnerable to sample composition bias across years. 
# Because averaging occurs across all raw predictions prior to isolating the 
# treatment effect, temporal shifts in dataset composition (such as new series 
# entering in later years) can distort the global trajectory.

# (2) Alternate Series level approach
# This method computes the log-scale treatment effect for each population 
# series individually before aggregating. By subtracting the counterfactual prediction 
# from the factual prediction within each series, all shared random effects (such 
# as location, baseline density, and series-specific slopes) cancel out 
# Each series-level effect is anchored to its own 1970 baseline value, and the 1970-anchored 
# log effects are averaged across active populations per year to form a geometric 
# mean index ratio.
# Clean causal isolation at the series level. Requires every included series 
# to have an observation in 1970. Thus substantial data loss and  
# survivor bias toward older monitoring programs.

dat_used_obs <- mod_dat_full %>%
  filter(used_f == "1", as.numeric(as.character(fyear)) >= 1969)

# Make observed and counterfactual data (inc zeroing use interactions)
grid_factual <- dat_used_obs %>%
  mutate(arm = "factual")

grid_counterfactual <- dat_used_obs %>%
  mutate(arm = "counterfactual",
         used_f = factor("0", levels = c("0", "1")), used_num = 0,
         inter_luc_num = 0, inter_pol_num = 0, inter_dis_num = 0, 
         inter_cli_num = 0, inter_inv_num = 0)

# Combine
pred_grid_att <- bind_rows(grid_factual, grid_counterfactual) %>%
  mutate(.row = row_number())

draws_att <- fitted_samples(glob.mod.gam.rhotuned,  n = 100,
                            data = pred_grid_att, exclude = NULL) %>%
  left_join(pred_grid_att, by = ".row")

saveRDS(draws_att, "Outputs/Draws/counterfac.draws.inter.rds", compress = "xz")
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

chained_plt<- ggplot(att_chained, aes(x = year, y = att_ratio)) +
  geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = .3) +
  geom_ribbon(data = att_chained_sum, 
              aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = att_chained_sum, aes(year, att_ratio_mn), colour = "black", linewidth = .5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  coord_cartesian(ylim = c(0, 1.3), expand = F) +
  labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  theme_bw()







