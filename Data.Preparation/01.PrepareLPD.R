
#############################
##    DATA PREPARATION     ##
#############################

# Script to process the LPD as provided by the ZSL LPD team (NOTE - this includes
# the additional threat data not available in the initial public download).
# Processing filters the data to only include populations where threat and utilisation
# state is knownn (present or absent). Subsequently filtering focuses on series
# of at least 7 years length, and 4 obs (main text) longer cut offs tested in the SM.

# Initial code adapted from https://github.com/duncanobrien/multiple-threats/tree/v1.0.2/Code


library(tidyverse)

# Load the living planet data
lpd_raw <- data.table::fread("Data/LPD - use and trade/LPD_20260424.csv")
length(unique(lpd_raw$ID))
length(unique(lpd_raw$Binomial))
min()
# Not all time series have been checked for threats - remove unknowns
lpd_raw %>% group_by(Threat_status) %>% tally()
lpd_raw <- lpd_raw %>% filter(!grepl("Unknown", Threat_status), Threat_status != "Not yet coded")
length(unique(lpd_raw$ID))
# Not all the time series have infomation on use
lpd_raw %>% group_by(Utilised) %>%tally()
lpd_raw <- lpd_raw %>% filter(!grepl("Unknown", Utilised), Utilised != "Not yet coded")
length(unique(lpd_raw$ID))
## This leaves 6307 populations series.

# Make long format (479,332 obs)
lpd_long <- lpd_raw %>% 
  pivot_longer(cols = `1950`:`2025`, names_to = "Year", values_to = "Count") %>%
  mutate(Count = na_if(Count, "NULL"), Count = as.numeric(Count),
         Year = as.integer(Year))

lpd_long <- lpd_long %>%
  group_by(ID) %>%  
  # Calculate population change
  # (max(Count, na.rm=T) /100 this get onehundreth of the max population size
  # and adds this as an offset to remove log(0) events.
  # then logs the value
  mutate(popchange=log(Count+(max(Count, na.rm=T) /100))) %>%
         #Protected_status=gsub(" .*", "", Protected_status)) 
  ungroup()

# Filter for 7 or more years of data with at least 50% of the years with records
##  Note: Sci Adv paper kept only 10 or more years we're trying with 7.
lpd_final <- lpd_long %>%
  group_by(ID) %>%
  drop_na(popchange) %>%
  mutate(n = n(),
         Duration=(max(Year)-min(Year))+1,
         ratio=n/Duration) %>% 
  ungroup() %>% 
  filter(Duration>6,
         ratio>=0.5)
length(unique(lpd_final$ID)) # 2792 populations

lpd_final %>% group_by(ID, Year) %>% filter(n()>1)
lpd_final %>% group_by(ID) %>% tally() %>% group_by(n) %>% tally() %>% arrange(n)
threat_dummies <- lpd_final %>%
  select(ID, Year, Primary_threat, Secondary_threat, Tertiary_threat) %>%
  pivot_longer(cols = ends_with("_threat"),values_to = "threat") %>%
  #filter(!is.na(threat), !threat %in% c("NULL", "None", "")) %>%
  distinct(ID, Year, threat) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = threat,values_from = presence,values_fill = 0,
    names_prefix = "threat_", names_repair = janitor::make_clean_names) %>%
  rename("Year" = "year", "ID" = "id") %>%
  mutate(threat_habitat_loss_change = ifelse(threat_habitat_loss ==1|threat_habitat_degradation_change == 1, 1, 0)) %>%
  select(!threat_null)

# Further wrangling in to long format and prepare explanatory variables
mod_dat_full <- lpd_final %>% ungroup() %>%
  left_join(threat_dummies) %>%
  mutate(Year = as.numeric(Year),
         series = as.character(ID)) %>%
  group_by(ID) %>%
  mutate(scaled_series_year = c(scale(Year,center = TRUE,scale = TRUE)), #z-score time to 0 for each timeseries
         time = seq_along(Year) - 1,
         time_dec = time/10,
         time_scaled = c(scale(time,center = TRUE,scale = TRUE)),
         popchange_centered = c(scale(popchange,center = TRUE,scale = FALSE))) %>% 
  group_by(ID) %>%
  arrange(Year) %>%
  mutate(popchange_t0 = popchange[1],
    popchange_0 = popchange - popchange_t0,
    first_year = min(Year)) %>%
  ungroup(ID) %>%
  mutate(Taxon = case_when(Class %in% c("Holocephali", "Elasmobranchii", "Myxini", 
                   "Cephalaspidomorphi", "Actinopteri", "Sarcopterygii", "Dipneusti") ~ "Fish",
                   Class == "Aves" ~ "Birds",
                   Class == "Mammalia" ~ "Mammals",
                   Class == "Amphibia" ~ "Amphibians",
                   Class == "Reptilia" ~ "Reptiles",
                   .default = NA_character_),
         # do not consider exploitation here and hab change and hab loss combined
         threats = threat_invasive_spp_genes + threat_climate_change + threat_disease +
           threat_pollution + threat_habitat_loss_change, 
         used = ifelse(Utilised == "Yes",1, 0),
         year.z = (Year - mean(Year))/sd(Year),
         fyear = as.character(Year),
         Site = paste(Latitude,Longitude,sep = "_"))

save(mod_dat_full, file = "Data/data_models.RData")
