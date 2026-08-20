library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)

load(file = "Data/data_models.RData")


mod_dat_full %>% group_by(Utilised) %>% summarise(length(unique(ID)))
mod_dat_full %>% group_by(Utilised) %>% summarise(length(unique(Binomial)))

mod_dat_full %>% group_by(threats) %>% summarise(length(unique(ID)))
mod_dat_full %>% group_by(threats) %>% summarise(length(unique(Binomial)))

d <- mod_dat_full %>% group_by(Country) %>% summarise(length(unique(ID)))
unique(mod_dat_full$Country)

unique(world$sovereignt)


mod_dat_full <- mod_dat_full %>%
  mutate(
    Country_clean = iconv(Country, from = "latin1", to = "UTF-8", sub = "?"),
    Country_clean = case_when(
      grepl("rkiye", Country_clean, ignore.case = TRUE) ~ "Turkey",
      grepl("land Islands", Country_clean, ignore.case = TRUE)~ "Aland Islands",
      grepl("Cura", Country_clean, ignore.case = TRUE)~ "Curacao",
      grepl("union", Country_clean, ignore.case = TRUE)~ "Reunion",
      grepl("Ivoire", Country_clean, ignore.case = TRUE) ~ "Cote d'Ivoire",
      TRUE ~ Country_clean))

mod_dat_full$iso3c <- countrycode(
  mod_dat_full$Country_clean, 
  origin = "country.name", 
  destination = "iso3c",
  custom_match = c(
    "International Waters" = NA_character_,
    "Reunion"              = "REU",
    "R?Runion"             = "REU",
    "Curacao"              = "CUW",
    "Aland Islands"        = "ALA"))

national.sum <- mod_dat_full %>% group_by(iso3c) %>%
  rename("iso_a3_eh" = "iso3c") %>%
  summarise(n.series = length(unique(ID)),
            n.used.series = n_distinct(ID[used == 1])) %>%
  mutate(series.cat = case_when(n.series < 6 ~ "up to 5",
                                n.series >5 & n.series <= 50 ~ "6-50",
                                n.series >50 & n.series <= 100 ~ "51-100",
                                n.series >100 ~ ">100"),
         perc.used = n.used.series/n.series *100,
         series.cat = factor(series.cat, levels = c("up to 5", "6-50", "51-100", ">100")))

national.sum <- world %>% left_join(national.sum)


ggplot(national.sum) +
  geom_sf(aes(fill = series.cat), color = "grey30", size = 0.1) +
  scale_fill_manual(values = c("#fff7bc", "#fec44f", "#ec7014", "#993404"),
                    na.value = NA, na.translate = F, "Time series")+
  #scale_fill_viridis_d(option = "inferno") +
  coord_sf(crs = "+proj=moll") +
  theme_minimal() +
  theme( legend.position = "bottom",panel.grid = element_blank(),
    axis.text = element_blank(),axis.ticks = element_blank()) 

ggplot(national.sum) +
  geom_sf(aes(fill = perc.used), color = "grey30", size = 0.1) +
  scale_fill_gradientn(
    colors = c("#fff7bc", "#993404"),
    na.value = "grey90",                        # Light grey for NA countries
    limits = c(0, 100),                         # Pin scale to 0–100%
    labels = function(x) paste0(x, "%"),        # Adds % sign to legend breaks
    name = "Used series (%)"
  ) +
  #scale_fill_viridis_d(option = "inferno") +
  coord_sf(crs = "+proj=moll") +
  theme_minimal() +
  theme( legend.position = "bottom",panel.grid = element_blank(),
         axis.text = element_blank(),axis.ticks = element_blank()) 


ggplot(mod_dat_full, aes(Year, popchange_centered, group = ID)) +
  geom_line(, alpha = .1)

ggplot(mod_dat_full, aes(scaled_year, popchange_centered, group = ID)) +
  geom_line(, alpha = .1)

ggplot(mod_dat_full, aes(scaled_time, popchange_centered, group = ID)) +
  geom_line(, alpha = .1)

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "gray92", color = "gray75", linewidth = 0.2) +
  geom_point(
    data = mod_dat_full,
    aes(x = Longitude, y = Latitude),
    color = "firebrick",
    alpha = 0.4,
    size = 1.2) +
  coord_sf(expand = FALSE) +
  theme_minimal() 
