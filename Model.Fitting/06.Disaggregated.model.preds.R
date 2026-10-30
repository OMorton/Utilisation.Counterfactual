library(tidyverse)
library(gratia)
library(ggpubr)

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


## Class ATT loop --------------------------------------------------------------

taxa_list <- c("Aves", "Actinopteri", "Elasmobranchii", "Mammalia", "Reptilia")

taxa_plt_ls <- list()
taxa <- "Aves"

for (taxa in taxa_list) {
  
  cat(" Processing Taxa Group:", taxa, "\n")

  # Dynamic file paths
  model_path <- paste0("Outputs/Models/Class/", taxa, ".use.threat.syn.full.v2.rdsxz")
  
  # Load model and filter taxa subset
  class.mod <- readRDS(model_path)
  dat_sub   <- mod_dat_full %>% filter(Class == taxa)
  
  # Filter for observed used series >= 1969
  dat_used_obs <- dat_sub %>%
    filter(used_f == "1", as.numeric(as.character(fyear)) >= 1969)
  
  if (nrow(dat_used_obs) == 0) {
    warning("No 'used' observations >= 1969 for ", taxa, ". Skipping...")
    next
  }
  
  # Build factual & counterfactual grids
  grid_factual <- dat_used_obs %>% 
    mutate(arm = "factual")
  
  grid_counterfactual <- dat_used_obs %>%
    mutate(arm = "counterfactual",
           used_f = factor("0", levels = c("0", "1")), used_num = 0,
           inter_luc_num = 0, inter_pol_num = 0, inter_dis_num = 0, 
           inter_cli_num = 0, inter_inv_num = 0)
  
  pred_grid_att <- bind_rows(grid_factual, grid_counterfactual) %>%
    mutate(.row = row_number())
  
  # Posterior draws
  draws_att <- fitted_samples(class.mod, n = 100, data = pred_grid_att, exclude = NULL) %>%
    left_join(pred_grid_att, by = ".row")

  ## A. ATT contrasts -----------------------------------------------
  att_chained <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    pivot_wider(names_from = arm, values_from = .fitted) %>%
    mutate(log_effect = factual - counterfactual) %>%
    group_by(series, .draw) %>%
    arrange(year) %>%
    mutate(inc_log_effect = log_effect - lag(log_effect)) %>%
    ungroup() %>%
    filter(!is.na(inc_log_effect)) %>%
    group_by(year, .draw) %>%
    summarise(mean_inc_att = mean(inc_log_effect, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw) %>%
    arrange(year) %>%
    mutate(cum_log_att = cumsum(ifelse(year == 1970, 0, mean_inc_att)),
           att_ratio = exp(cum_log_att)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  write_rds(att_chained, 
             paste0("Outputs/Draws/Class/", taxa, "_draws.v2.csv"), compress = "xz")

  
  att_chained_sum <- att_chained %>% 
    group_by(year) %>%
    summarise(att_ratio_mn = mean(att_ratio),
              att_ratio_lo = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
              att_ratio_hi = quantile(att_ratio, probs = 0.975, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(taxa = taxa) %>%
    select(taxa, year, everything())
  
  # Save Summary Contrast Table
  write_csv(att_chained_sum, paste0("Outputs/Tables/Class/", taxa, "_summary_contrast.v2.csv"))
  
  ## Reference summary figures
  # p_att <- ggplot(att_chained, aes(x = year, y = att_ratio)) +
  #   geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
  #   geom_ribbon(data = att_chained_sum, 
  #               aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
  #               fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
  #   geom_line(data = att_chained_sum, aes(year, att_ratio_mn), colour = "black", linewidth = 0.5) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   coord_cartesian(ylim = c(0, max(att_chained$att_ratio)), expand = FALSE) +
  #   labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  #   theme_bw()
  # 
  # dat_counts <- dat_sub %>%
  #   mutate(year = as.numeric(as.character(fyear))) %>%
  #   filter(year >= 1970) %>%
  #   group_by(year, used_f) %>%
  #   summarise(n_series = n_distinct(series), .groups = "drop") %>%
  #   mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
  #                              levels = c("Used", "Not used")))
  # 
  # max_series_count <- dat_counts %>%
  #   group_by(year) %>%
  #   summarise(total_n = sum(n_series), .groups = "drop") %>%
  #   pull(total_n) %>%
  #   max()
  # 
  # p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
  #   geom_col(width = 0.8) +
  #   scale_fill_manual(values = c("Used" = "darkred", "Not used" = "grey60")) +
  #   coord_cartesian(xlim = c(1970, max(att_chained$year)), expand = FALSE) +
  #   scale_y_continuous(breaks = c(0, max_series_count), limits = c(0, max_series_count)) +
  #   labs(y = "Populations") +
  #   theme_bw() +
  #   theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #     legend.position = "none",
  #     panel.grid.minor = element_blank())
  # 
  # p_arr <- ggarrange(p_bar, p_att, nrow = 2, ncol = 1,
  #           heights = c(0.25,1))
  # 
  # ggsave(paste0("Outputs/Figures/Class/", taxa, "_fig_contrast.v2.png"), 
  #        plot = p_arr, width = 7, height = 5, dpi = 300)
  
  ## B. counter preds  ------------------------------------------
  chained_predictions <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    group_by(series, .draw, arm) %>%
    arrange(year) %>%
    mutate(inc_log = .fitted - lag(.fitted)) %>%
    ungroup() %>%
    filter(!is.na(inc_log)) %>%
    group_by(year, .draw, arm) %>%
    summarise(mean_inc = mean(inc_log, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw, arm) %>%
    arrange(year) %>%
    mutate(cum_log  = cumsum(ifelse(year == 1970, 0, mean_inc)),
           pop_index = exp(cum_log)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  chained_pred_sum <- chained_predictions %>%
    group_by(year, arm) %>%
    summarise(
      pop_index_mn = mean(pop_index),
      pop_index_lo = quantile(pop_index, probs = 0.025, na.rm = TRUE),
      pop_index_hi = quantile(pop_index, probs = 0.975, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(taxa = taxa) %>%
    select(taxa, year, arm, everything())
  
  # Save Summary Predictions Table
  write_csv(chained_pred_sum, paste0("Outputs/Tables/", taxa, "_summary_predictions.v2.csv"))
  
  ## Reference predictions plots
  # p_pred <- ggplot(chained_predictions, aes(x = year, y = pop_index, color = arm, fill = arm)) +
  #   geom_line(aes(group = interaction(.draw, arm)), alpha = 0.08, linewidth = 0.3) +
  #   geom_line(data = chained_pred_sum, aes(y = pop_index_mn, group = arm), linewidth = 0.8) +
  #   geom_ribbon(data = chained_pred_sum,
  #               aes(y = pop_index_mn, ymin = pop_index_lo, ymax = pop_index_hi, group = arm),
  #               alpha = 0.2, color = NA) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   scale_color_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                      labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   scale_fill_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                     labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   coord_cartesian(xlim = c(1970, max(chained_predictions$year)), expand = FALSE) +
  #   labs(x = "Year", y = "Relative Population Index\n(1970 = 1.0)") +
  #   theme_bw() +
  #   theme(legend.title = element_blank(), legend.position = "bottom")
  # 
  # ggsave(paste0("Outputs/Figures/Class/", taxa, "_fig_predictions.v2.png"), 
  #        plot = p_pred, width = 7, height = 5, dpi = 300)
  
  taxa_plt_ls[[taxa]] <- list(contrast_arr = p_arr, contrast = p_att,
    predictions = p_pred)
}

## System ATT loop -------------------------------------------------------------
system_list <- c("Freshwater", "Marine", "Terrestrial")

system_plt_ls <- list()

system <- "Freshwater"
for (system in system_list) {
  
  cat(" Processing system Group:", system, "\n")
  
  # Dynamic file paths
  model_path <- paste0("Outputs/Models/", system, ".use.threat.syn.full.v2.rdsxz")
  
  # Load model and filter system subset
  class.mod <- readRDS(model_path)
  dat_sub   <- mod_dat_full %>% filter(System == system)
  
  # Filter for observed used series >= 1969
  dat_used_obs <- dat_sub %>%
    filter(used_f == "1", as.numeric(as.character(fyear)) >= 1969)
  
  if (nrow(dat_used_obs) == 0) {
    warning("No 'used' observations >= 1969 for ", system, ". Skipping...")
    next
  }
  
  # Build factual & counterfactual grids
  grid_factual <- dat_used_obs %>% 
    mutate(arm = "factual")
  
  grid_counterfactual <- dat_used_obs %>%
    mutate(arm = "counterfactual",
           used_f = factor("0", levels = c("0", "1")), used_num = 0,
           inter_luc_num = 0, inter_pol_num = 0, inter_dis_num = 0, 
           inter_cli_num = 0, inter_inv_num = 0)
  
  pred_grid_att <- bind_rows(grid_factual, grid_counterfactual) %>%
    mutate(.row = row_number())
  
  # Posterior draws
  draws_att <- fitted_samples(class.mod, n = 100, data = pred_grid_att, exclude = NULL) %>%
    left_join(pred_grid_att, by = ".row")
  
  
  
  ## A. ATT contrasts -----------------------------------------------
  att_chained <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    pivot_wider(names_from = arm, values_from = .fitted) %>%
    mutate(log_effect = factual - counterfactual) %>%
    group_by(series, .draw) %>%
    arrange(year) %>%
    mutate(inc_log_effect = log_effect - lag(log_effect)) %>%
    ungroup() %>%
    filter(!is.na(inc_log_effect)) %>%
    group_by(year, .draw) %>%
    summarise(mean_inc_att = mean(inc_log_effect, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw) %>%
    arrange(year) %>%
    mutate(cum_log_att = cumsum(ifelse(year == 1970, 0, mean_inc_att)),
           att_ratio = exp(cum_log_att)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  write_rds(att_chained, 
            paste0("Outputs/Draws/System/", system, "_draws.v2.csv"), compress = "xz")
  
  
  att_chained_sum <- att_chained %>% 
    group_by(year) %>%
    summarise(att_ratio_mn = mean(att_ratio),
              att_ratio_lo = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
              att_ratio_hi = quantile(att_ratio, probs = 0.975, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(system = system) %>%
    select(system, year, everything())
  
  # Save Summary Contrast Table
  write_csv(att_chained_sum, paste0("Outputs/Tables/System/", system, "_summary_contrast.v2.csv"))
  
  # p_att <- ggplot(att_chained, aes(x = year, y = att_ratio)) +
  #   geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
  #   geom_ribbon(data = att_chained_sum, 
  #               aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
  #               fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
  #   geom_line(data = att_chained_sum, aes(year, att_ratio_mn), colour = "black", linewidth = 0.5) +
  #   geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
  #   coord_cartesian(ylim = c(0, max(att_chained$att_ratio)), expand = FALSE) +
  #   labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  #   theme_bw(base_size = 9)
  #   
  # dat_counts <- dat_sub %>%
  #   mutate(year = as.numeric(as.character(fyear))) %>%
  #   filter(year >= 1970) %>%
  #   group_by(year, used_f) %>%
  #   summarise(n_series = n_distinct(series), .groups = "drop") %>%
  #   mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
  #                              levels = c("Used", "Not used")))
  # 
  # max_series_count <- dat_counts %>%
  #   group_by(year) %>%
  #   summarise(total_n = sum(n_series), .groups = "drop") %>%
  #   pull(total_n) %>%
  #   max()
  # 
  # p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
  #   geom_col(width = 0.8) +
  #   scale_fill_manual(values = c("Used" = "darkred", "Not used" = "grey60")) +
  #   coord_cartesian(xlim = c(1970, max(att_chained$year)), expand = FALSE) +
  #   scale_y_continuous(breaks = c(0, max_series_count), limits = c(0, max_series_count)) +
  #   #labs(y = "Populations") +
  #   labs(y = "") +
  #   theme_bw(base_size = 9) +
  #   theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         legend.position = "none",
  #         panel.grid.minor = element_blank())
  # 
  # p_arr <- ggarrange(p_bar, p_att, nrow = 2, ncol = 1,
  #                    heights = c(0.25,1), align = "v")
  # 
  # ggsave(paste0("Outputs/Figures/System/", system, "_fig_contrast.v2.png"), 
  #        plot = p_arr, width = 7, height = 5, dpi = 300)
  
  ## B. counter preds  ------------------------------------------
  chained_predictions <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    group_by(series, .draw, arm) %>%
    arrange(year) %>%
    mutate(inc_log = .fitted - lag(.fitted)) %>%
    ungroup() %>%
    filter(!is.na(inc_log)) %>%
    group_by(year, .draw, arm) %>%
    summarise(mean_inc = mean(inc_log, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw, arm) %>%
    arrange(year) %>%
    mutate(cum_log  = cumsum(ifelse(year == 1970, 0, mean_inc)),
           pop_index = exp(cum_log)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  chained_pred_sum <- chained_predictions %>%
    group_by(year, arm) %>%
    summarise(
      pop_index_mn = mean(pop_index),
      pop_index_lo = quantile(pop_index, probs = 0.025, na.rm = TRUE),
      pop_index_hi = quantile(pop_index, probs = 0.975, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(system = system) %>%
    select(system, year, arm, everything())
  
  # Save Summary Predictions Table
  write_csv(chained_pred_sum, paste0("Outputs/Tables/System/", system, "_summary_predictions.v2.csv"))
  
  # p_pred <- ggplot(chained_predictions, aes(x = year, y = pop_index, color = arm, fill = arm)) +
  #   geom_line(aes(group = interaction(.draw, arm)), alpha = 0.08, linewidth = 0.3) +
  #   geom_line(data = chained_pred_sum, aes(y = pop_index_mn, group = arm), linewidth = 0.8) +
  #   geom_ribbon(data = chained_pred_sum,
  #               aes(y = pop_index_mn, ymin = pop_index_lo, ymax = pop_index_hi, group = arm),
  #               alpha = 0.2, color = NA) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   scale_color_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                      labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   scale_fill_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                     labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   coord_cartesian(xlim = c(1970, max(chained_predictions$year)), expand = FALSE) +
  #   labs(x = "Year", y = "Relative Population Index\n(1970 = 1.0)") +
  #   theme_bw() +
  #   theme(legend.title = element_blank(), legend.position = "bottom")
  # 
  # ggsave(paste0("Outputs/Figures/System/", system, "_fig_predictions.v2.png"), 
  #        plot = p_pred, width = 7, height = 5, dpi = 300)
  
  system_plt_ls[[system]] <- list(contrast_arr = p_arr, contrast = p_att,
                              predictions = p_pred)
}

# sys.plt <- ggarrange(system_plt_ls$Freshwater$contrast_arr, 
#                      system_plt_ls$Marine$contrast_arr,
#           system_plt_ls$Terrestrial$contrast_arr, nrow = 3,
#           labels = c("A.", "B.", "C."))
# 
# ggsave(paste0("Outputs/Figures/SM/sysste.counterfac.v1.png"), 
#        plot = sys.plt, width = 6, height = 9, dpi = 300)

## Country ATT loop ------------------------------------------------------------

model_files <- list.files("Outputs/Models/Country", pattern = "\\.use\\.threat\\.syn\\.full\\.v1\\.rdsxz$", full.names = TRUE)

country_plt_ls <- list()

mod_dat_full <- mod_dat_full %>%
  mutate(Country = iconv(Country, to = "UTF-8", sub = ""),
    Country_clean = gsub("[^A-Za-z0-9_-]", "_", Country))

for (model_path in model_files) {
  table_dir  <- "Outputs/Tables/Country"
  figure_dir <- "Outputs/Figures/Country"
  
  # Extract country name directly from the file name
  file_name <- basename(model_path)
  country_clean <- sub("\\.use\\.threat\\.syn\\.full\\.v1\\.rdsxz$", "", file_name)
  
  cat("\n Processing Country Model:", country_clean, "\n")
  
  # Load country model
  country.mod <- readRDS(model_path)
  
  # Match country in dataset (handling sanitized underscores if necessary)
  dat_sub <- mod_dat_full %>% 
    filter(Country_clean == country_clean)
  
  # Fallback filter if original country string is present
  if (nrow(dat_sub) == 0) {
    dat_sub <- mod_dat_full %>% filter(Country == country_clean)
  }
  
  # Filter for observed used series >= 1969
  dat_used_obs <- dat_sub %>%
    filter(used_f == "1", as.numeric(as.character(fyear)) >= 1969)
  
  # Skip if no valid 'used' observations exist for this country
  if (nrow(dat_used_obs) == 0) {
    warning("No 'used' observations >= 1969 for ", country_clean, ". Skipping...")
    next
  }
  
  # Build factual & counterfactual grids
  grid_factual <- dat_used_obs %>% 
    mutate(arm = "factual")
  
  grid_counterfactual <- dat_used_obs %>%
    mutate(arm = "counterfactual",
           used_f = factor("0", levels = c("0", "1")), used_num = 0,
           inter_luc_num = 0, inter_pol_num = 0, inter_dis_num = 0, 
           inter_cli_num = 0, inter_inv_num = 0)
  
  pred_grid_att <- bind_rows(grid_factual, grid_counterfactual) %>%
    mutate(.row = row_number())
  
  # Posterior draws
  draws_att <- tryCatch({
    fitted_samples(country.mod, n = 100, data = pred_grid_att, exclude = NULL) %>%
      left_join(pred_grid_att, by = ".row")
  }, error = function(e) {
    warning("Failed to calculate fitted samples for ", country_clean, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(draws_att)) next
  
  ##  ATT contrasts
  att_chained <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    pivot_wider(names_from = arm, values_from = .fitted) %>%
    mutate(log_effect = factual - counterfactual) %>%
    group_by(series, .draw) %>%
    arrange(year) %>%
    mutate(inc_log_effect = log_effect - lag(log_effect)) %>%
    ungroup() %>%
    filter(!is.na(inc_log_effect)) %>%
    group_by(year, .draw) %>%
    summarise(mean_inc_att = mean(inc_log_effect, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw) %>%
    arrange(year) %>%
    mutate(cum_log_att = cumsum(ifelse(year == 1970, 0, mean_inc_att)),
           att_ratio = exp(cum_log_att)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  if (nrow(att_chained) == 0) {
    warning("Insufficient temporal overlap post-1970 for ", country_clean, ". Skipping...")
    next
  }
  write_rds(att_chained, 
            paste0("Outputs/Draws/Country/", country_clean, "_draws.csv"), compress = "xz")
  
  att_chained_sum <- att_chained %>% 
    group_by(year) %>%
    summarise(att_ratio_mn = mean(att_ratio),
              att_ratio_lo = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
              att_ratio_hi = quantile(att_ratio, probs = 0.975, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(country = country_clean) %>%
    select(country, year, everything())
  
  # Save Summary Contrast Table directly in Outputs/Tables/
  write_csv(att_chained_sum, file.path(table_dir, paste0(country_clean, "_summary_contrast.csv")))
  
  # p_att <- ggplot(att_chained, aes(x = year, y = att_ratio)) +
  #   geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
  #   geom_ribbon(data = att_chained_sum, 
  #               aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
  #               fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
  #   geom_line(data = att_chained_sum, aes(year, att_ratio_mn), colour = "black", linewidth = 0.5) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   coord_cartesian(ylim = c(0, max(att_chained$att_ratio, na.rm = TRUE)), expand = FALSE) +
  #   labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  #   theme_bw()
  # 
  # dat_counts <- dat_sub %>%
  #   mutate(year = as.numeric(as.character(fyear))) %>%
  #   filter(year >= 1970) %>%
  #   group_by(year, used_f) %>%
  #   summarise(n_series = n_distinct(series), .groups = "drop") %>%
  #   mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
  #                              levels = c("Used", "Not used")))
  # 
  # max_series_count <- dat_counts %>%
  #   group_by(year) %>%
  #   summarise(total_n = sum(n_series), .groups = "drop") %>%
  #   pull(total_n) %>%
  #   max(na.rm = TRUE)
  # 
  # p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
  #   geom_col(width = 0.8) +
  #   scale_fill_manual(values = c("Used" = "darkred", "Not used" = "grey60")) +
  #   coord_cartesian(xlim = c(1970, max(att_chained$year)), expand = FALSE) +
  #   scale_y_continuous(breaks = c(0, max_series_count), limits = c(0, max_series_count)) +
  #   labs(y = "Populations") +
  #   theme_bw() +
  #   theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         legend.position = "none",
  #         panel.grid.minor = element_blank())
  # 
  # p_arr <- ggarrange(p_bar, p_att, nrow = 2, ncol = 1, heights = c(0.25, 1), 
  #                    align = "v")
  # 
  # # Save Contrast Figure directly in Outputs/Figures/
  # ggsave(file.path(figure_dir, paste0(country_clean, "_fig_contrast.png")), 
  #        plot = p_arr, width = 7, height = 5, dpi = 300)
  
  ## Counterfactual predictions
  chained_predictions <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    group_by(series, .draw, arm) %>%
    arrange(year) %>%
    mutate(inc_log = .fitted - lag(.fitted)) %>%
    ungroup() %>%
    filter(!is.na(inc_log)) %>%
    group_by(year, .draw, arm) %>%
    summarise(mean_inc = mean(inc_log, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw, arm) %>%
    arrange(year) %>%
    mutate(cum_log  = cumsum(ifelse(year == 1970, 0, mean_inc)),
           pop_index = exp(cum_log)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  chained_pred_sum <- chained_predictions %>%
    group_by(year, arm) %>%
    summarise(
      pop_index_mn = mean(pop_index),
      pop_index_lo = quantile(pop_index, probs = 0.025, na.rm = TRUE),
      pop_index_hi = quantile(pop_index, probs = 0.975, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(country = country_clean) %>%
    select(country, year, arm, everything())
  
  # Save Summary Predictions Table directly in Outputs/Tables/
  write_csv(chained_pred_sum, file.path(table_dir, paste0(country_clean, "_summary_predictions.csv")))
  
  # p_pred <- ggplot(chained_predictions, aes(x = year, y = pop_index, color = arm, fill = arm)) +
  #   geom_line(aes(group = interaction(.draw, arm)), alpha = 0.08, linewidth = 0.3) +
  #   geom_line(data = chained_pred_sum, aes(y = pop_index_mn, group = arm), linewidth = 0.8) +
  #   geom_ribbon(data = chained_pred_sum,
  #               aes(y = pop_index_mn, ymin = pop_index_lo, ymax = pop_index_hi, group = arm),
  #               alpha = 0.2, color = NA) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   scale_color_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                      labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   scale_fill_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                     labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   coord_cartesian(xlim = c(1970, max(chained_predictions$year)), expand = FALSE) +
  #   labs(x = "Year", y = "Relative Population Index\n(1970 = 1.0)") +
  #   theme_bw() +
  #   theme(legend.title = element_blank(), legend.position = "bottom")
  # 
  # # Save Predictions Figure directly in Outputs/Figures/
  # ggsave(file.path(figure_dir, paste0(country_clean, "_fig_predictions.png")), 
  #        plot = p_pred, width = 7, height = 5, dpi = 300)
  
  # Store plot outputs in list
  country_plt_ls[[country_clean]] <- list(contrast_arr = p_arr, contrast = p_att,
    predictions  = p_pred)
}

## Region ATT loop -------------------------------------------------------------

model_files <- list.files("Outputs/Models/Region", pattern = "\\.use\\.threat\\.syn\\.full\\.v2\\.rdsxz$", full.names = TRUE)
model_path <- "Outputs/Models/Region/North_America.use.threat.syn.full.v2.rdsxz"
region_plt_ls <- list()

mod_dat_full <- mod_dat_full %>%
  mutate(Country = iconv(Country, to = "UTF-8", sub = ""),
         Country_clean = gsub("[^A-Za-z0-9_-]", "_", Country))

for (model_path in model_files) {
  table_dir  <- "Outputs/Tables/Region"
  figure_dir <- "Outputs/Figures/Region"
  
  # Extract country name directly from the file name
  file_name <- basename(model_path)
  Region_clean <- gsub("_", " ", sub("\\.use\\.threat\\.syn\\.full\\.v2\\.rdsxz$", "", file_name))
  
  cat("\n Processing Region Model:", Region_clean, "\n")
  
  # Load country model
  region.mod <- readRDS(model_path)
  
  # Match country in dataset (handling sanitized underscores if necessary)
  dat_sub <- mod_dat_full %>% 
    filter(Region == Region_clean)
  
  # Filter for observed used series >= 1969
  dat_used_obs <- dat_sub %>%
    filter(used_f == "1", as.numeric(as.character(fyear)) >= 1969)
  
  # Skip if no valid 'used' observations exist for this country
  if (nrow(dat_used_obs) == 0) {
    warning("No 'used' observations >= 1969 for ", Region_clean, ". Skipping...")
    next
  }
  
  # Build factual & counterfactual grids
  grid_factual <- dat_used_obs %>% 
    mutate(arm = "factual")
  
  grid_counterfactual <- dat_used_obs %>%
    mutate(arm = "counterfactual",
           used_f = factor("0", levels = c("0", "1")), used_num = 0,
           inter_luc_num = 0, inter_pol_num = 0, inter_dis_num = 0, 
           inter_cli_num = 0, inter_inv_num = 0)
  
  pred_grid_att <- bind_rows(grid_factual, grid_counterfactual) %>%
    mutate(.row = row_number())
  
  # Posterior draws
  draws_att <- tryCatch({
    fitted_samples(region.mod, n = 100, data = pred_grid_att, exclude = NULL) %>%
      left_join(pred_grid_att, by = ".row")
  }, error = function(e) {
    warning("Failed to calculate fitted samples for ", Region_clean, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(draws_att)) next
  
  ##  ATT contrasts
  att_chained <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    pivot_wider(names_from = arm, values_from = .fitted) %>%
    mutate(log_effect = factual - counterfactual) %>%
    group_by(series, .draw) %>%
    arrange(year) %>%
    mutate(inc_log_effect = log_effect - lag(log_effect)) %>%
    ungroup() %>%
    filter(!is.na(inc_log_effect)) %>%
    group_by(year, .draw) %>%
    summarise(mean_inc_att = mean(inc_log_effect, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw) %>%
    arrange(year) %>%
    mutate(cum_log_att = cumsum(ifelse(year == 1970, 0, mean_inc_att)),
           att_ratio = exp(cum_log_att)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  write_rds(att_chained, 
            paste0("Outputs/Draws/Region/", Region_clean, "_draws.csv"), compress = "xz")
  
  if (nrow(att_chained) == 0) {
    warning("Insufficient temporal overlap post-1970 for ", Region_clean, ". Skipping...")
    next
  }
  
  att_chained_sum <- att_chained %>% 
    group_by(year) %>%
    summarise(att_ratio_mn = mean(att_ratio),
              att_ratio_lo = quantile(att_ratio, probs = 0.025, na.rm = TRUE),
              att_ratio_hi = quantile(att_ratio, probs = 0.975, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(country = Region_clean) %>%
    select(country, year, everything())
  
  # Save Summary Contrast Table directly in Outputs/Tables/
  write_csv(att_chained_sum, file.path(table_dir, paste0(Region_clean, "_summary_contrast.csv")))
  
  # p_att <- ggplot(att_chained, aes(x = year, y = att_ratio)) +
  #   geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
  #   geom_ribbon(data = att_chained_sum, 
  #               aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
  #               fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
  #   geom_line(data = att_chained_sum, aes(year, att_ratio_mn), colour = "black", linewidth = 0.5) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   coord_cartesian(ylim = c(0, max(att_chained$att_ratio, na.rm = TRUE)), expand = FALSE) +
  #   labs(x = "Year", y = "Use impact index\n(1.0 = no net effect)") +
  #   theme_bw()
  # 
  # dat_counts <- dat_sub %>%
  #   mutate(year = as.numeric(as.character(fyear))) %>%
  #   filter(year >= 1970) %>%
  #   group_by(year, used_f) %>%
  #   summarise(n_series = n_distinct(series), .groups = "drop") %>%
  #   mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
  #                              levels = c("Used", "Not used")))
  # 
  # max_series_count <- dat_counts %>%
  #   group_by(year) %>%
  #   summarise(total_n = sum(n_series), .groups = "drop") %>%
  #   pull(total_n) %>%
  #   max(na.rm = TRUE)
  # 
  # p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
  #   geom_col(width = 0.8) +
  #   scale_fill_manual(values = c("Used" = "darkred", "Not used" = "grey60")) +
  #   coord_cartesian(xlim = c(1970, max(att_chained$year)), expand = FALSE) +
  #   scale_y_continuous(breaks = c(0, max_series_count), limits = c(0, max_series_count)) +
  #   labs(y = "Populations") +
  #   theme_bw() +
  #   theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         legend.position = "none",
  #         panel.grid.minor = element_blank())
  # 
  # p_arr <- ggarrange(p_bar, p_att, nrow = 2, ncol = 1, heights = c(0.25, 1), 
  #                    align = "v")
  # 
  # # Save Contrast Figure directly in Outputs/Figures/
  # ggsave(file.path(figure_dir, paste0(Region_clean, "_fig_contrast.png")), 
  #        plot = p_arr, width = 7, height = 5, dpi = 300)
  
  ## Counterfactual predictions
  chained_predictions <- draws_att %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    select(series, year, .draw, arm, .fitted) %>%
    group_by(series, .draw, arm) %>%
    arrange(year) %>%
    mutate(inc_log = .fitted - lag(.fitted)) %>%
    ungroup() %>%
    filter(!is.na(inc_log)) %>%
    group_by(year, .draw, arm) %>%
    summarise(mean_inc = mean(inc_log, na.rm = TRUE), .groups = "drop") %>%
    group_by(.draw, arm) %>%
    arrange(year) %>%
    mutate(cum_log  = cumsum(ifelse(year == 1970, 0, mean_inc)),
           pop_index = exp(cum_log)) %>%
    ungroup() %>%
    filter(year >= 1970)
  
  chained_pred_sum <- chained_predictions %>%
    group_by(year, arm) %>%
    summarise(
      pop_index_mn = mean(pop_index),
      pop_index_lo = quantile(pop_index, probs = 0.025, na.rm = TRUE),
      pop_index_hi = quantile(pop_index, probs = 0.975, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(Region = Region_clean) %>%
    select(Region, year, arm, everything())
  
  # Save Summary Predictions Table directly in Outputs/Tables/
  write_csv(chained_pred_sum, file.path(table_dir, paste0(Region_clean, "_summary_predictions.csv")))
  
  # p_pred <- ggplot(chained_predictions, aes(x = year, y = pop_index, color = arm, fill = arm)) +
  #   geom_line(aes(group = interaction(.draw, arm)), alpha = 0.08, linewidth = 0.3) +
  #   geom_line(data = chained_pred_sum, aes(y = pop_index_mn, group = arm), linewidth = 0.8) +
  #   geom_ribbon(data = chained_pred_sum,
  #               aes(y = pop_index_mn, ymin = pop_index_lo, ymax = pop_index_hi, group = arm),
  #               alpha = 0.2, color = NA) +
  #   geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #   scale_color_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                      labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   scale_fill_manual(values = c("factual" = "#bf812d", "counterfactual" = "#35978f"),
  #                     labels = c("Counterfactual (No Use)", "Factual (Observed Use)")) +
  #   coord_cartesian(xlim = c(1970, max(chained_predictions$year)), expand = FALSE) +
  #   labs(x = "Year", y = "Relative Population Index\n(1970 = 1.0)") +
  #   theme_bw() +
  #   theme(legend.title = element_blank(), legend.position = "bottom")
  # 
  # # Save Predictions Figure directly in Outputs/Figures/
  # ggsave(file.path(figure_dir, paste0(Region_clean, "_fig_predictions.png")), 
  #        plot = p_pred, width = 7, height = 5, dpi = 300)
  
  # Store plot outputs in list
  region_plt_ls[[Region_clean]] <- list(contrast_arr = p_arr, contrast = p_att,
                                          predictions  = p_pred)
}
