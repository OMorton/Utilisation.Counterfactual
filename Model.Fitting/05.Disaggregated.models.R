library(tidyverse)
library(mgcv)
library(gratia)

load(file = "Data/data_models.RData")
source("functions.R")
min.series <- 10

## Data Summaries --------------------------------------------------------------
mod_dat_full %>% group_by(Class, used) %>% summarise(n = length(unique(series))) %>%
  group_by(Class) %>% filter(n() == 2, all(n >= min.series))
mod_dat_full %>% group_by(System, used) %>% summarise(length(unique(series)))
mod_dat_full %>% group_by(Region, used) %>% summarise(n = length(unique(series)))
mod_dat_full %>% group_by(Country, used) %>% summarise(n = length(unique(series))) %>%
  group_by(Country) %>% filter(n() == 2, all(n >= min.series))

## Data Preparation -----------------------------------------------------------
threat_cols <- c("used", "threat_habitat_loss_change", "threat_pollution", 
                 "threat_disease", "threat_climate_change", "threat_invasive_spp_genes")

mod_dat_full <- mod_dat_full %>%
  mutate(series = factor(series), fyear = factor(fyear),  Site = factor(Site),
    across(all_of(threat_cols), ~ factor(.x), .names = "{.col}_f"),
    across(all_of(threat_cols), ~ ordered(.x, levels = c("0", "1")), .names = "{.col}_ord"),
    used_num = ifelse(used_f == "1", 1, 0),
    inter_luc_num = used_num * ifelse(threat_habitat_loss_change_f == "1", 1, 0),
    inter_pol_num = used_num * ifelse(threat_pollution_f == "1", 1, 0),
    inter_dis_num = used_num * ifelse(threat_disease_f == "1", 1, 0),
    inter_cli_num = used_num * ifelse(threat_climate_change_f == "1", 1, 0),
    inter_inv_num = used_num * ifelse(threat_invasive_spp_genes_f == "1", 1, 0)) %>%
  arrange(Class, series, year.z)

# Threat mapping metadata for dynamic formula construction
threat_lookup <- list(
  luc = list(ord = "threat_habitat_loss_change_ord", inter = "inter_luc_num"),
  pol = list(ord = "threat_pollution_ord", inter = "inter_pol_num"),
  dis = list(ord = "threat_disease_ord",inter = "inter_dis_num"),
  cli = list(ord = "threat_climate_change_ord", inter = "inter_cli_num"),
  inv = list(ord = "threat_invasive_spp_genes_ord", inter = "inter_inv_num"))

## Modelling: Taxonomic Class --------------------------------------------------
class.n <- mod_dat_full %>% 
  group_by(Class, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= min.series))

class.ls <- unique(class.n$Class)
model_results <- list()

i <- "Reptilia"
for (i in class.ls) {
  message(paste("\n--- Processing Class:", i, "---"))
  
  safe_class_name <- gsub("[^A-Za-z0-9_]", "_", i)
  
  # Subsetting, AR_start assignment, fyear factor creation, and year z-scoring
  dat_sub <- mod_dat_full %>% 
    filter(Class == i) %>% 
    droplevels() %>% 
    mutate(AR_start = !duplicated(series), fyear = factor(Year)) %>% 
    ungroup() %>%
    mutate(year.z = (Year - mean(Year)) / sd(Year))
  
  # Assemble base candidate terms
  candidate_terms <- c("- 1", "used_f", "s(year.z, by = used_f, k = 15, bs = 'cr')",
    "s(Latitude, Longitude, bs = 'tp', k = 20)",
    "ti(Latitude, Longitude, year.z, d = c(2, 1), bs = c('tp', 'cr'))",
    "s(scaled_series_year, series, k = 3, bs = 'fs')",
    "s(fyear, used_f, bs = 're')")
  
  # Enforce Model Marginality during candidate term assembly
  for (t_name in names(threat_lookup)) {
    ord_var   <- threat_lookup[[t_name]]$ord
    inter_var <- threat_lookup[[t_name]]$inter
    
    has_ord <- ord_var %in% names(dat_sub) && length(unique(na.omit(dat_sub[[ord_var]]))) > 1
    
    if (has_ord) {
      # Add Main Parametric & Smooth Effects
      candidate_terms <- c(candidate_terms, ord_var, 
        paste0("s(year.z, by = ", ord_var, ", k = 15, bs = 'cr')"))
      
      # Add Interaction Smooth ONLY if Main Effect exists AND interaction has >= 10 non-zero values
      has_inter <- inter_var %in% names(dat_sub) && sum(dat_sub[[inter_var]] != 0, na.rm = TRUE) >= 10
      if (has_inter) {
        candidate_terms <- c(candidate_terms, 
          paste0("s(year.z, by = ", inter_var, ", k = 10, bs = 'cr')"))
      }
    }
  }
  
  # Dynamic Pruning with strict marginality enforcement via threat_lookup
  final_terms     <- get_clean_formula_terms(candidate_terms, dat_sub, threat_lookup = threat_lookup)
  dynamic_formula <- as.formula(paste("popchange ~", paste(final_terms, collapse = " + ")))
  
  # Fit Naive Model
  mod_naive <- tryCatch({
    bam(formula = dynamic_formula, data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Naive GAM]:", e$message))
    return(NULL)
  })
  
  if (is.null(mod_naive)) next
  
  # Extract Rho & fit tuned AR(1) model
  rho_est <- acf(residuals(mod_naive, type = "response"), plot = FALSE)$acf[2]
  
  mod_tuned <- tryCatch({
    bam(formula = dynamic_formula, rho = rho_est, AR.start = dat_sub$AR_start,
        data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Tuned GAM]:", e$message))
    return(NULL)
  })
  
  if (is.null(mod_tuned)) next
  
  # Save result summary & fitted model object
  model_results[[i]] <- list(formula = dynamic_formula, rho_est = rho_est)
  
  saveRDS(mod_tuned, paste0("Outputs/Models/", safe_class_name, ".use.threat.syn.full.v2.rds", 
                            compress = "xz"))
}

save(model_results, file = "Outputs/Models/taxonomic.models.v2.Rdata")

## Modelling: System  --------------------------------------------------
system.n <- mod_dat_full %>% 
  group_by(System, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= min.series))

system.ls <- unique(system.n$System)
model_results <- list()

i <- "Marine"
for (i in system.ls) {
  message(paste("\n--- Processing Class:", i, "---"))
  
  safe_system_name <- gsub("[^A-Za-z0-9_]", "_", i)
  
  # Subsetting, AR_start assignment, fyear factor creation, and year z-scoring
  dat_sub <- mod_dat_full %>% 
    filter(System == i) %>% 
    droplevels() %>% 
    mutate(AR_start = !duplicated(series), fyear = factor(Year)) %>% 
    ungroup() %>%
    mutate(year.z = (Year - mean(Year)) / sd(Year))
  
  # Assemble base candidate terms
  candidate_terms <- c("- 1", "used_f", "s(year.z, by = used_f, k = 15, bs = 'cr')",
                       "s(Latitude, Longitude, bs = 'tp', k = 20)",
                       "ti(Latitude, Longitude, year.z, d = c(2, 1), bs = c('tp', 'cr'))",
                       "s(scaled_series_year, series, k = 3, bs = 'fs')",
                       "s(fyear, used_f, bs = 're')")
  
  # Enforce Model Marginality during candidate term assembly
  for (t_name in names(threat_lookup)) {
    ord_var   <- threat_lookup[[t_name]]$ord
    inter_var <- threat_lookup[[t_name]]$inter
    
    has_ord <- ord_var %in% names(dat_sub) && length(unique(na.omit(dat_sub[[ord_var]]))) > 1
    
    if (has_ord) {
      # Add Main Parametric & Smooth Effects
      candidate_terms <- c(candidate_terms, ord_var, 
                           paste0("s(year.z, by = ", ord_var, ", k = 15, bs = 'cr')"))
      
      # Add Interaction Smooth ONLY if Main Effect exists AND interaction has >= 10 non-zero values
      has_inter <- inter_var %in% names(dat_sub) && sum(dat_sub[[inter_var]] != 0, na.rm = TRUE) >= 10
      if (has_inter) {
        candidate_terms <- c(candidate_terms, 
                             paste0("s(year.z, by = ", inter_var, ", k = 10, bs = 'cr')"))
      }
    }
  }
  
  # Dynamic Pruning with strict marginality enforcement via threat_lookup
  final_terms     <- get_clean_formula_terms(candidate_terms, dat_sub, threat_lookup = threat_lookup)
  dynamic_formula <- as.formula(paste("popchange ~", paste(final_terms, collapse = " + ")))
  
  # Fit Naive Model
  mod_naive <- tryCatch({
    bam(formula = dynamic_formula, data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Naive GAM]:", e$message))
    return(NULL)
  })
  
  if (is.null(mod_naive)) next
  
  # Extract Rho & fit tuned AR(1) model
  rho_est <- acf(residuals(mod_naive, type = "response"), plot = FALSE)$acf[2]
  
  mod_tuned <- tryCatch({
    bam(formula = dynamic_formula, rho = rho_est, AR.start = dat_sub$AR_start,
        data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Tuned GAM]:", e$message))
    return(NULL)
  })
  
  if (is.null(mod_tuned)) next
  
  # Save result summary & fitted model object
  model_results[[i]] <- list(formula = dynamic_formula, rho_est = rho_est)
  
  saveRDS(mod_tuned, paste0("Outputs/Models/", safe_system_name, ".use.threat.syn.full.v2.rds", 
                            compress = "xz"))
}

save(model_results, file = "Outputs/Models/system.models.v2.Rdata")


## Modelling: Region -----------------------------------------------------------

region.n <- mod_dat_full %>% 
  group_by(Region, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= min.series))

region.ls <- unique(region.n$Region) # 7 regions
region_model_results <- list()


for (i in region.ls) {
  message(paste("\n--- Processing Region:", i, "---"))
  
  safe_region_name <- gsub("[^A-Za-z0-9_]", "_", i)
  
  dat_sub <- mod_dat_full %>% 
    filter(Region == i) %>% 
    droplevels() %>% 
    mutate(
      AR_start = !duplicated(series),
      fyear = factor(Year)) %>% 
    ungroup() %>%
    mutate(year.z = (Year - mean(Year)) / sd(Year))
  
  # Base candidate terms
  candidate_terms <- c( "- 1",
    "used_f","s(year.z, by = used_f, k = 15, bs = 'cr')",
    "s(Latitude, Longitude, bs = 'tp', k = 20)",
    "ti(Latitude, Longitude, year.z, d = c(2, 1), bs = c('tp', 'cr'))",
    "s(scaled_series_year, series, k = 3, bs = 'fs')",
    "s(fyear, used_f, bs = 're')")
  
  # Enforce Model Marginality during term assembly
  for (t_name in names(threat_lookup)) {
    ord_var   <- threat_lookup[[t_name]]$ord
    inter_var <- threat_lookup[[t_name]]$inter
    
    has_ord <- ord_var %in% names(dat_sub) && length(unique(na.omit(dat_sub[[ord_var]]))) > 1
    
    if (has_ord) {
      # 1. Add Main Parametric & Smooth Effects
      candidate_terms <- c(
        candidate_terms, 
        ord_var, 
        paste0("s(year.z, by = ", ord_var, ", k = 15, bs = 'cr')"))
      
      # 2. Add Interaction Smooth ONLY if Main Effect exists AND interaction has >= 10 non-zero values
      has_inter <- inter_var %in% names(dat_sub) && sum(dat_sub[[inter_var]] != 0, na.rm = TRUE) >= 10
      if (has_inter) {
        candidate_terms <- c(
          candidate_terms, 
          paste0("s(year.z, by = ", inter_var, ", k = 10, bs = 'cr')"))
      }
    }
  }
  
  # Dynamic Pruning
  final_terms     <- get_clean_formula_terms(candidate_terms, dat_sub, threat_lookup = threat_lookup)
  dynamic_formula <- as.formula(paste("popchange ~", paste(final_terms, collapse = " + ")))
  
  # Naive model
  mod_naive <- tryCatch({
    bam(formula = dynamic_formula, data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Naive GAM]:", e$message))
    return(NULL)
  })
  
  if (is.null(mod_naive)) next
  
  # Extract Rho & fit tuned model
  rho_est <- acf(residuals(mod_naive, type = "response"), plot = FALSE)$acf[2]
  
  mod_tuned <- tryCatch({
    bam(formula = dynamic_formula, rho = rho_est, AR.start = dat_sub$AR_start,
        data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Tuned GAM]:", e$message))
    return(NULL)
  })
  
  if (is.null(mod_tuned)) next
  
  region_model_results[[i]] <- list(formula = dynamic_formula, rho_est = rho_est)
  
  saveRDS(mod_tuned, paste0("Outputs/Models/Region/", safe_region_name, ".use.threat.syn.full.v2.rdsxz"), 
          compress = "xz")
}

save(region_model_results, file = "Outputs/Models/region.models.v2.Rdata")
## Modelling: National Level ---------------------------------------------------

country.n <- mod_dat_full %>% 
  group_by(Country, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= min.series))

country.ls <- unique(country.n$Country) # 18 countries
country_model_results <- list()


for (i in country.ls) {
  message(paste("\n--- Processing Country:", i, "---"))
  
  safe_country_name <- gsub("[^A-Za-z0-9_]", "_", i)
  
  # subsetting, AR_start assignment, re-scale global year smooth
  dat_sub <- mod_dat_full %>% 
    filter(Country == i) %>% 
    droplevels() %>% 
    mutate(AR_start = !duplicated(series)) %>% 
    ungroup() %>%
    mutate(year.z = (Year - mean(Year)) / sd(Year))
  
  # assemble base candidate terms
  candidate_terms <- c("- 1",
    "used_f", "s(year.z, by = used_f, k = 15, bs = 'cr')",
    "s(Latitude, Longitude, bs = 'tp', k = 20)",
    "ti(Latitude, Longitude, year.z, d = c(2, 1), bs = c('tp', 'cr'))",
    "s(scaled_series_year, series, k = 3, bs = 'fs')") # due to small datasizes fyear excluded
  
  # compile threat and threat interactions where combinations exist
  for (t_name in names(threat_lookup)) {
    ord_var   <- threat_lookup[[t_name]]$ord
    inter_var <- threat_lookup[[t_name]]$inter
    
    if (ord_var %in% names(dat_sub) && length(unique(na.omit(dat_sub[[ord_var]]))) > 1) {
      candidate_terms <- c(
        candidate_terms, 
        ord_var, 
        paste0("s(year.z, by = ", ord_var, ", k = 15, bs = 'cr')"))
    }
    
    if (inter_var %in% names(dat_sub) && length(unique(na.omit(dat_sub[[inter_var]]))) > 1) {
      candidate_terms <- c(
        candidate_terms, 
        paste0("s(year.z, by = ", inter_var, ", k = 10, bs = 'cr')"))
    }
  }
  
  # Dynamic Pruning to ensure full-rank parametric matrix
  final_terms     <- get_clean_formula_terms(candidate_terms, dat_sub)
  dynamic_formula <- as.formula(paste("popchange ~", paste(final_terms, collapse = " + ")))
  
  # Naive model
  mod_naive <- tryCatch({
    bam(formula = dynamic_formula,
        data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Naive GAM]:", e$message))
    return(NULL)})
  
  if (is.null(mod_naive)) next
  
  # Extract Rho & fit tuned model
  rho_est <- acf(residuals(mod_naive, type = "response"), plot = FALSE)$acf[2]
  
  mod_tuned <- tryCatch({
    bam(formula = dynamic_formula,
        rho = rho_est, AR.start = dat_sub$AR_start,
        data = dat_sub, discrete = TRUE, nthreads = 4)
  }, error = function(e) {
    message(paste("  [Error in Tuned GAM]:", e$message))
    return(NULL)})
  
  if (is.null(mod_tuned)) next
  
  country_model_results[[i]] <- list(formula = dynamic_formula,
                                     rho_est = rho_est)
  
  # Save tuned country GAM
  saveRDS(mod_tuned, paste0("Outputs/Models/Country/", safe_country_name, ".use.threat.syn.full.v1.rdsxz"), 
          compress = "xz")
}

save(country_model_results, file = "Outputs/Models/country.models.Rdata")
