library(tidyverse)
library(mgcv)
library(MASS)
library(gratia)

load(file = "Data/data_models.RData")

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

## Make precision matrix for mgcv (REDUNDANT CURRENTLY)-------------------------
# distances
locs <- distinct(mod_dat_full[,c("Longitude","Latitude")]) %>%
  mutate(Site = paste(Latitude,Longitude,sep = "_"))
dist_mat <- geosphere::distm(locs[, c("Longitude", "Latitude")], 
                             fun = geosphere::distHaversine) / 1000

# converts  distance into correlation between sites using exponential decay
phi <- median(dist_mat[dist_mat > 0]) 
cov_mat <- exp(-dist_mat / phi)

# add a 2% Nugget (0.02) to cap maximum eigenvalues
nugget <- 0.02 
cov_mat_nugget <- cov_mat + diag(nugget, nrow(cov_mat))

colnames(cov_mat_nugget) <- locs$Site
rownames(cov_mat_nugget) <- locs$Site

# precision matrix
spa_precision <- solve(cov_mat_nugget)
rm(cov_mat_nugget)
rm(cov_mat)
rm(dist_mat)

## Fit AR Naive model ----------------------------------------------------------
glob.mod.gam.rhonaive <- bam(
  popchange ~ 
    # Focal Use baseline offset and main trajectory
    used_f + 
    s(year.z, by = used_f, k = 15, bs = "cr") + 
    # Threat main effects (parametric offsets + difference smooths)
    threat_habitat_loss_change_ord + threat_pollution_ord + 
    threat_disease_ord + threat_climate_change_ord + threat_invasive_spp_genes_ord +
    s(year.z, by = threat_habitat_loss_change_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_pollution_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_disease_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_climate_change_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_invasive_spp_genes_ord, k = 15, bs = "cr") +
    # Synergistic Use x Threat interaction smooths
    s(year.z, by = inter_luc_num, k = 10, bs = "cr") +
    s(year.z, by = inter_pol_num, k = 10, bs = "cr") +
    s(year.z, by = inter_dis_num, k = 10, bs = "cr") +
    s(year.z, by = inter_cli_num, k = 10, bs = "cr") +
    s(year.z, by = inter_inv_num, k = 10, bs = "cr") +
    # Spatial field & Spatiotemporal dynamics
    s(Latitude, Longitude, bs = "tp", k = 20) +
    ti(Latitude, Longitude, year.z, d = c(2, 1), bs = c("tp", "cr")) +
    # Random slopes, intercepts for series & annual macro shocks
    s(series, bs = "re") +
    s(series, time_scaled, bs = "re") + 
    s(fyear, used_f, bs = "re"),
  #rho = 0.2462075, AR.start = mod_dat_full$AR_start,
  data = mod_dat_full, discrete = TRUE, nthreads = 4)

rho_est <- acf(residuals(glob.mod.gam.rhonaive), plot = FALSE)$acf[2] # 0.2422293

## Fit AR Model ----------------------------------------------------------------
glob.mod.gam.rhotuned <- bam(
  popchange ~ 
    # Focal Use baseline offset and main trajectory
    used_f + 
    s(year.z, by = used_f, k = 15, bs = "cr") + 
    # Threat main effects (parametric offsets + difference smooths)
    threat_habitat_loss_change_ord + threat_pollution_ord + 
    threat_disease_ord + threat_climate_change_ord + threat_invasive_spp_genes_ord +
    s(year.z, by = threat_habitat_loss_change_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_pollution_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_disease_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_climate_change_ord, k = 15, bs = "cr") +
    s(year.z, by = threat_invasive_spp_genes_ord, k = 15, bs = "cr") +
    # Synergistic Use x Threat interaction smooths
    s(year.z, by = inter_luc_num, k = 10, bs = "cr") +
    s(year.z, by = inter_pol_num, k = 10, bs = "cr") +
    s(year.z, by = inter_dis_num, k = 10, bs = "cr") +
    s(year.z, by = inter_cli_num, k = 10, bs = "cr") +
    s(year.z, by = inter_inv_num, k = 10, bs = "cr") +
    # Spatial field & Spatiotemporal dynamics
    s(Latitude, Longitude, bs = "tp", k = 20) +
    ti(Latitude, Longitude, year.z, d = c(2, 1), bs = c("tp", "cr")) +
    # Random slopes, intercepts for series & annual macro shocks
    s(series, bs = "re") +
    s(series, time_scaled, bs = "re") + 
    s(fyear, used_f, bs = "re"),
  rho = 0.2422293, AR.start = mod_dat_full$AR_start,
  data = mod_dat_full, discrete = TRUE, nthreads = 4)

# Check residual autocorr (~0.00)
itsadug::acf_resid(glob.mod.gam.rhotuned, plot = FALSE)[2]

saveRDS(glob.mod.gam.rhotuned, "GLOB.use.threat.syn.full.v1.rds", compress = "xz")