
library(tidyverse)
library(mgcv)
library(gratia)
library(ggpubr)
detach("package:MASS", unload = TRUE)

load(file = "Data/data_models.RData")
glob.mod.gam.rhotuned <- readRDS("Outputs/Models/GLOB.use.threat.syn.full.v2.rds")

## Data Preparation -----------------------------------------------------------
threat_cols <- c("used", "threat_habitat_loss_change", "threat_pollution", 
                 "threat_disease", "threat_climate_change", "threat_invasive_spp_genes")

mod_dat_full <- mod_dat_full %>%
  mutate(series = factor(series), fyear = factor(fyear), Site = factor(Site),
         # ordered threat factors
         across(all_of(threat_cols), ~ factor(.x), .names = "{.col}_f"),
         across(all_of(threat_cols), ~ ordered(.x, levels = c("0", "1")), .names = "{.col}_ord"),
         used_num = ifelse(used_f == "1", 1, 0),
         # Synergistic interaction multipliers (Used = 1 AND Threat = 1)
         inter_luc_num = used_num * ifelse(threat_habitat_loss_change_f == "1", 1, 0),
         inter_pol_num = used_num * ifelse(threat_pollution_f == "1", 1, 0),
         inter_dis_num = used_num * ifelse(threat_disease_f == "1", 1, 0),
         inter_cli_num = used_num * ifelse(threat_climate_change_f == "1", 1, 0),
         inter_inv_num = used_num * ifelse(threat_invasive_spp_genes_f == "1", 1, 0),
  ) %>%
  arrange(series, year.z) %>%
  group_by(series) %>%
  mutate(AR_start = (row_number() == 1)) %>% # Marks beginning of each time-series
  ungroup()

## Average populations predictions ---------------------------------------------

rank <- mod_dat_full %>% group_by(threat_habitat_loss_change_f, threat_pollution_f, threat_disease_f, threat_climate_change_f, threat_invasive_spp_genes_f) %>% tally()

# map fyear to year.z
year_map <- mod_dat_full %>%
  distinct(fyear, year.z) %>%
  arrange(year.z)

# Scenario helper function incorporating exact fyear-year.z pairs
make_scenario_grid <- function(scenario_name, luc = 0, pol = 0, dis = 0, cli = 0, inv = 0) {
  expand.grid(fyear    = year_map$fyear, 
              used_f   = factor(c("0", "1"), levels = c("0", "1")),
              threat_habitat_loss_change_ord = ordered(luc, levels = c("0", "1")),
              threat_pollution_ord           = ordered(pol, levels = c("0", "1")),
              threat_disease_ord             = ordered(dis, levels = c("0", "1")),
              threat_climate_change_ord     = ordered(cli, levels = c("0", "1")),
              threat_invasive_spp_genes_ord  = ordered(inv, levels = c("0", "1")),
    # Baseline dummy series and median coordinates excluded from prediction
    series   = mod_dat_full$series[500], 
    time_dec = 0, time_scaled = 0, scaled_series_year = 0,
    Latitude = median(mod_dat_full$Latitude, na.rm = TRUE),
    Longitude = median(mod_dat_full$Longitude, na.rm = TRUE)) %>%
    # Join the exact numeric year.z corresponding to each fyear level
    left_join(year_map, by = "fyear") %>%
    mutate(scenario = scenario_name,
          used_num = ifelse(used_f == "1", 1, 0),
          inter_luc_num = used_num * luc,
          inter_pol_num = used_num * pol,
          inter_dis_num = used_num * dis,
          inter_cli_num = used_num * cli,
          inter_inv_num = used_num * inv)
}

# build multi-scenario grid
pred_grid_inter <- bind_rows(
  make_scenario_grid("No Threats"),
  make_scenario_grid("LUC", luc = 1),
  make_scenario_grid("CC", cli = 1),
  make_scenario_grid("Pol", pol = 1),
  make_scenario_grid("LUC + Pol", luc = 1, pol = 1),
  make_scenario_grid("LUC + CC", luc = 1, cli = 1)) %>% 
  mutate(.row = row_number())

all_smooth_labels <- sapply(glob.mod.gam.rhotuned$smooth, function(x) x$label)
rand_effects <- all_smooth_labels[grep("series|Latitude|Longitude", all_smooth_labels)]
rand_effects_to_exclude <- c(
  "s(time_dec,series)", "s(time_scaled,series)" ,"s(series)","s(scaled_series_year,series)",
  "s(Latitude,Longitude)", "ti(year.z,Latitude,Longitude)")

draws_inter <- fitted_samples(glob.mod.gam.rhotuned,  n = 100, 
                              data = pred_grid_inter, exclude = rand_effects_to_exclude) %>% 
  left_join(pred_grid_inter, by = ".row") %>%
  group_by(used_f, scenario, .draw) %>%
  arrange(year.z) %>%
  mutate(
    # Raw Model Prediction: Predicted log-scale state (popchange) in year t.
    rate = .fitted,
    # Baseline Log Difference: Cumulative log-scale change relative to 1970 (r_t - r_1970).
    # A value of 0 represents no change from the 1970 baseline value.
    rate_rel = .fitted - .fitted[fyear == "1970"][1],
    # Relative LPI style Index: Back-transformed multiplicative index anchored to 1.0 in 1970.
    # (e.g., 0.80 = 20% decline relative to 1970 level).
    rate_rel_index = exp(.fitted - .fitted[fyear == "1970"][1]),
    # Absolute Exponentiated Scale: Back-transformed multiplicative state exp(r_t).
    # Exponentiated raw prediction without reference to a baseline year.
    rate_index = exp(.fitted),
    year = as.numeric(as.character(fyear))) %>%
  ungroup() %>% filter(year >= 1970)

saveRDS(draws_inter, "Outputs/Draws/ave.draws.inter.rds", compress = "xz")

ave_pred_plt <- ggplot(draws_inter, aes(x = year, y = rate_rel_index, group = interaction(.draw, used_f), color = used_f)) +
  geom_line(alpha = 0.12, linewidth = .3) +
  stat_summary(aes(group = used_f), fun = mean, geom = "line", linewidth = .5) +
  geom_hline(yintercept =  1, linewidth = 0.5) +
  facet_wrap(~factor(scenario, levels=c("No Threats","CC","LUC","Pol", 
                                        "LUC + CC", "LUC + Pol"))) +
  coord_cartesian(xlim = c(1970, 2025)) +
  scale_color_manual(values = c("0" = "grey25", "1" = "#bf812d"), labels = c("Not Used", "Used")) +
  xlab("Year") + ylab("Relative population index") +
  theme_bw() +
  theme(strip.background =  element_blank(), legend.title = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(path = "Outputs/Figures", filename = "ave.preds.v1.png",
       plot = ave_pred_plt, width = 8, height = 6, units = "in",bg = "white")

## Average population contrast -------------------------------------------------

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
    pct_diff = (ratio_rel_index - 1) * 100)

contrast_sum <- contrast_draws %>% group_by(scenario, year) %>%
  summarise(rel_ratio_mn = mean(ratio_rel_index),
            rel_ratio_lo  = quantile(ratio_rel_index, probs = 0.025, na.rm = TRUE),
            rel_ratio_hi  = quantile(ratio_rel_index, probs = 0.975, na.rm = TRUE))

ave_contrast_plt <- ggplot(contrast_draws, aes(x = year, y = ratio_rel_index)) +
  geom_line(aes(group = .draw), alpha = 0.08, color = "darkred") +
  geom_ribbon(data = contrast_sum, 
              aes(y = rel_ratio_mn, ymin = rel_ratio_lo, ymax = rel_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = contrast_sum, aes(year, rel_ratio_mn), colour = "black", linewidth = .5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~factor(scenario, levels=c("No Threats","CC","LUC","Pol", 
                                        "LUC + CC", "LUC + Pol"))) +
  labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  theme_bw() +
  theme(strip.background =  element_blank(), legend.title = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(path = "Outputs/Figures", filename = "ave.contrast.v1.png",
       plot = ave_contrast_plt, width = 8, height = 6, units = "in",bg = "white")
## Annual average counterfactual -----------------------------------------------

# Calculates global mean trajectories by averaging model predictions 
# across all treated populations for each year, under both factual (harvested) and 
# counterfactual (unharvested) conditions. Each global trajectory is anchored to 
# its 1970 value to produce relative growth indices, and the treatment effect is 
# measured as the ratio between the two resulting curves.

# vulnerable to sample composition bias across years. 
# Because averaging occurs across all raw predictions prior to isolating the 
# treatment effect, temporal shifts in dataset composition (such as new series 
# entering in later years) can distort the global trajectory.

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


# Aggregate across treated series, anchor to 1970, and compute ATT
att_annual_average <- draws_att %>%
  mutate(year = as.numeric(as.character(fyear))) %>%
  group_by(arm, year, .draw) %>%
  summarise(mean_fitted = mean(.fitted, na.rm = TRUE), .groups = "drop") %>%
  group_by(arm, .draw) %>%
  arrange(year) %>%
  mutate(
    # Anchor each arm's trajectory to its 1970 baseline value
    rate_rel = mean_fitted - mean_fitted[year == 1970][1],
    rate_rel_index = exp(rate_rel)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(year, .draw),
    names_from = arm, values_from = c(mean_fitted, rate_rel, rate_rel_index)) %>%
  mutate(
    # Log-scale ATT: Net log divergence attributable strictly to harvesting
    att_log = rate_rel_factual - rate_rel_counterfactual,
    # Relative ATT Index Ratio: (1970 Relative Index Factual) / (1970 Relative Index Counterfactual)
    att_ratio = rate_rel_index_factual / rate_rel_index_counterfactual) %>% 
  filter(year >= 1970)

att_annual_ave_sum <- att_annual_average %>% group_by(year) %>%
  summarise(att_ratio_mn = mean(att_ratio),
            att_ratio_lo  = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
            att_ratio_hi  = quantile(att_ratio, probs = 0.975, na.rm = TRUE))

ann_ave_plt <- ggplot(att_annual_average, aes(x = year, y = att_ratio)) +
  geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
  geom_ribbon(data = att_annual_ave_sum, 
              aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = att_annual_ave_sum, aes(year, att_ratio_mn), colour = "black", linewidth = .5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  coord_cartesian(ylim = c(0, 1.1), expand = F) +
  labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  theme_bw()

## Series-level counterfactual -------------------------------------------------

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

att_series_level <- draws_att %>%
  mutate(year = as.numeric(as.character(fyear))) %>%
  select(series, year, .draw, arm, .fitted) %>%
  pivot_wider(names_from = arm, values_from = .fitted) %>%
  # Log difference per series (random effects cancel out mathematically)
  mutate(log_effect = factual - counterfactual) %>%
  group_by(series, .draw) %>%
  arrange(year) %>%
  mutate( log_att_series = log_effect - log_effect[year == 1970][1]) %>%
  ungroup() %>%
  # Filter out series that were not monitored in 1970 (missing 1970 anchor)
  filter(!is.na(log_att_series)) %>%
  # Average log effects across active treated series per year per draw
  group_by(year, .draw) %>%
  summarise(mean_log_att = mean(log_att_series, na.rm = TRUE), .groups = "drop") %>%
  # exponentiate averaged log difference (Geometric Mean ATT Ratio)
  mutate(att_ratio = exp(mean_log_att)) %>% filter(year >= 1970)

att_series_sum <- att_series_level %>% group_by(year) %>%
  summarise(att_ratio_mn = mean(att_ratio),
            att_ratio_lo  = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
            att_ratio_hi  = quantile(att_ratio, probs = 0.975, na.rm = TRUE))

series_plt <- ggplot(att_series_level, aes(x = year, y = att_ratio)) +
  geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
  geom_ribbon(data = att_series_sum, 
              aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
              fill = NA, colour = "black", linetype = "dashed", linewidth = .5) +
  geom_line(data = att_series_sum, aes(year, att_ratio_mn), colour = "black", linewidth = .5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  coord_cartesian(ylim = c(0, 1.1), expand = F) +
  labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  theme_bw()

## Chain linked counterfactual -------------------------------------------------

# Aims to combines series-level causal cancellation without data loss or survivor
# bias. First, the net treatment effect is calculated within each series per draw,
# cancelling out series-specific REs. Second, the annual step-by-step change 
# in treatment effect between consecutive years is computed for each series.
# These annual incremental changes are then averaged across all populations
# in a given year, and the annual averages are cumulatively summed from 1970
# to build a continuous log-scale index. Exponentiating this cumulative val gives
# the final geometric mean ratio index.

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
  coord_cartesian(ylim = c(0, 1.1), expand = F) +
  labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  theme_bw()


counter_arr <- ggarrange(ann_ave_plt, series_plt, chained_plt, 
                         font.label = list(size = 10),
                         nrow = 1, ncol = 3, labels = c("A.", "B.", "C."))

ggsave(path = "Outputs/Figures", filename = "counterfac.est.v1.png",
       plot = counter_arr, width = 8, height = 3, units = "in",bg = "white")

