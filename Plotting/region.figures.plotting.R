library(tidyverse)
library(gratia)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)
library(sf)
library(ggrepel)
library(cowplot)

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

## Map plotting ----------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent != "Seven seas (open ocean)")

country.n <- mod_dat_full %>% 
  group_by(Country, Region, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= 10))
country.ls <- unique(country.n$Country) # 18 countries

iso.ls <- countrycode(
  country.ls, 
  origin = "country.name", 
  destination = "iso3c") 


world.regions <- world %>% mutate(country.mod = ifelse(iso_a3_eh %in% iso.ls, "Yes", "No"))

high_contrast_continents <- c("#C0392B","#F39C12", "#00B4D8", "#388E3C",
                              "#8E44AD", "#D81B60",  "#4E6E58")

map.col.plt <- ggplot() +
  geom_sf(data = world.regions, aes(fill = continent, geometry = geometry, 
                                    alpha = country.mod, colour = country.mod),
          size = 0.1) +
  scale_color_manual(values = c(NA, "black")) +
  scale_alpha_manual(values = c(0.25, 1)) +
  scale_fill_manual(values = high_contrast_continents) +
  coord_sf(crs = "+proj=moll") +
  theme_void() +
  theme(legend.position = "none")


## Region plotting -------------------------------------------------------------
afr_att_chained <- readRDS(paste0("Outputs/Draws/Region/Africa_draws.csv"))
afr_chained_sum <- read.csv(paste0("Outputs/Tables/Region/Africa_summary_contrast.csv"))
ant_att_chained <- readRDS(paste0("Outputs/Draws/Region/Antarctic_draws.csv"))
ant_chained_sum <- read.csv(paste0("Outputs/Tables/Region/Antarctic_summary_contrast.csv"))
asia_att_chained <- readRDS(paste0("Outputs/Draws/Region/Asia_draws.csv"))
asia_chained_sum <- read.csv(paste0("Outputs/Tables/Region/Asia_summary_contrast.csv"))
eur_att_chained <- readRDS(paste0("Outputs/Draws/Region/Europe_draws.csv"))
eur_chained_sum <- read.csv(paste0("Outputs/Tables/Region/Europe_summary_contrast.csv"))
lam_att_chained <- readRDS(paste0("Outputs/Draws/Region/Latin America and Caribbean_draws.csv"))
lam_chained_sum <- read.csv(paste0("Outputs/Tables/Region/Latin America and Caribbean_summary_contrast.csv"))
nam_att_chained <- readRDS(paste0("Outputs/Draws/Region/North America_draws.csv"))
nam_chained_sum <- read.csv(paste0("Outputs/Tables/Region/North America_summary_contrast.csv"))
oce_att_chained <- readRDS(paste0("Outputs/Draws/Region/Oceania_draws.csv"))
oce_chained_sum <- read.csv(paste0("Outputs/Tables/Region/Oceania_summary_contrast.csv"))

afr_plt <- region_att_plt(draws = afr_att_chained, summary = afr_chained_sum, common.y = T, reg.cols = T)
ant_plt <- region_att_plt(draws = ant_att_chained, summary = ant_chained_sum, common.y = T, reg.cols = T)
asia_plt <- region_att_plt(draws = asia_att_chained, summary = asia_chained_sum, common.y = T, reg.cols = T)
eur_plt <- region_att_plt(draws = eur_att_chained, summary = eur_chained_sum, common.y = T, reg.cols = T)
lam_plt <- region_att_plt(draws = lam_att_chained, summary = lam_chained_sum, common.y = T, reg.cols = T)
nam_plt <- region_att_plt(draws = nam_att_chained, summary = nam_chained_sum, common.y = T, reg.cols = T)
oce_plt <- region_att_plt(draws = oce_att_chained, summary = oce_chained_sum, common.y = T, reg.cols = T)


reg.arr <- ggarrange(
  # top row
  ggarrange(NULL,
    ggarrange(nam_plt$bar, nam_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(eur_plt$bar, eur_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(asia_plt$bar, asia_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    NULL, widths = c(0.5, 1, 1, 1, 0.5), nrow = 1, labels = c("", "A.", "B.", "C.", ""),
    label.x = -0.02),
  # mid row
  map.col.plt,
  # bottom row
  ggarrange(
    ggarrange(lam_plt$bar, lam_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(afr_plt$bar, afr_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(ant_plt$bar, ant_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(oce_plt$bar, oce_plt$att, ncol = 1, heights = c(0.2, 1), align = "v"),
    nrow = 1, labels = c("D.", "E.", "F.", "G."), label.x = -0.02), 
  nrow =3, heights = c(1, 1.5, 1))


arrow_style <- arrow(length = unit(0.08, "inches"), type = "closed")

reg.arr2 <- ggdraw(reg.arr) +
  # North America
  draw_line(x = c(0.28, 0.32), y = c(0.71, 0.62), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Europe
  draw_line(x = c(0.53, 0.60), y = c(0.71, 0.66), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Asia
  draw_line(x = c(0.78, 0.64), y = c(0.71, 0.60), color = "grey30", size = 0.6, arrow = arrow_style) +
    # Latin America
  draw_line(x = c(0.15, 0.39), y = c(0.29, 0.47), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Africa 
  draw_line(x = c(0.40, 0.55), y = c(0.29, 0.51), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Antarctica 
  draw_line(x = c(0.65, 0.52), y = c(0.29, 0.32), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Oceania
  draw_line(x = c(0.90, 0.78), y = c(0.29, 0.43), color = "grey30", size = 0.6, arrow = arrow_style)

ggsave(path = "Outputs/Figures/Tidy", filename = "Region.Est.v2.png",
       plot = reg.arr2, width = 8, height = 8, units = "in",bg = "white")

## log scale
reg.arr.log <- ggarrange(
  # top row
  ggarrange(NULL,
            ggarrange(nam_plt$bar, nam_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
            ggarrange(eur_plt$bar, eur_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
            ggarrange(asia_plt$bar, asia_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
            NULL, widths = c(0.5, 1, 1, 1, 0.5), nrow = 1, labels = c("", "A.", "B.", "C.", ""),
            label.x = -0.02),
  # mid row
  map.col.plt,
  # bottom row
  ggarrange(
    ggarrange(lam_plt$bar, lam_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(afr_plt$bar, afr_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(ant_plt$bar, ant_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
    ggarrange(oce_plt$bar, oce_plt$att+scale_y_log10() + coord_cartesian(ylim = c(0.01, 3),xlim = c(1970, 2024), expand = F), ncol = 1, heights = c(0.2, 1), align = "v"),
    nrow = 1, labels = c("D.", "E.", "F.", "G."), label.x = -0.02), 
  nrow =3, heights = c(1, 1.5, 1))


arrow_style <- arrow(length = unit(0.08, "inches"), type = "closed")

reg.arr.log2 <- ggdraw(reg.arr.log) +
  # North America
  draw_line(x = c(0.28, 0.32), y = c(0.71, 0.62), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Europe
  draw_line(x = c(0.53, 0.60), y = c(0.71, 0.66), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Asia
  draw_line(x = c(0.78, 0.64), y = c(0.71, 0.60), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Latin America
  draw_line(x = c(0.15, 0.39), y = c(0.29, 0.47), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Africa 
  draw_line(x = c(0.40, 0.55), y = c(0.29, 0.51), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Antarctica 
  draw_line(x = c(0.65, 0.52), y = c(0.29, 0.32), color = "grey30", size = 0.6, arrow = arrow_style) +
  # Oceania
  draw_line(x = c(0.90, 0.78), y = c(0.29, 0.43), color = "grey30", size = 0.6, arrow = arrow_style)

ggsave(path = "Outputs/Figures/Tidy", filename = "Region.Est.Log.v2.png",
       plot = reg.arr.log2, width = 8, height = 8, units = "in",bg = "white")

## Country level plotting ------------------------------------------------------

country_region_map <- mod_dat_full %>% 
  group_by(Country, Region, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= 10)) %>% select(Country, Region) %>% distinct() %>%
  ungroup() %>%
  mutate(iso = iso.ls)

mod_dat_full %>% 
  group_by(Country, Region, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= 10)) %>% group_by(Region) %>% distinct(Country)


country_plots <- set_names(iso.ls, country.ls) %>%
  imap(function(iso, country_name) {
    file_name <- gsub(" ", "_", country_name)
    
    # Extract the pre-assigned region for this country
    current_region <- country_region_map$Region[country_region_map$Country == country_name]
    
    draws_data <- readRDS(paste0("Outputs/Draws/Country/", file_name, "_draws.csv"))
    
    sum_data   <- read.csv(paste0("Outputs/Tables/Country/", file_name, "_summary_contrast.csv")) %>%
      mutate(
        iso = iso, 
        country = country_name,
        region = current_region
      )
    
    country_att_plt(
      draws = draws_data, 
      summary = sum_data, 
      common.y = TRUE, 
      reg.cols = TRUE
    )
  })

ggarrange(
  ggarrange(country_plots$Canada$bar, country_plots$Canada$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$`United States`$bar, country_plots$`United States`$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$`United Kingdom`$bar, country_plots$`United Kingdom`$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Spain$bar, country_plots$Spain$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$France$bar, country_plots$France$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Germany$bar, country_plots$Germany$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Italy$bar, country_plots$Italy$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Poland$bar, country_plots$Poland$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Norway$bar, country_plots$Norway$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Finland$bar, country_plots$Finland$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Sweden$bar, country_plots$Sweden$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$`Russian Federation`$bar, country_plots$`Russian Federation`$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Japan$bar, country_plots$Japan$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Brazil$bar, country_plots$Brazil$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$`South Africa`$bar, country_plots$`South Africa`$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Zimbabwe$bar, country_plots$Zimbabwe$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$Australia$bar, country_plots$Australia$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  ggarrange(country_plots$`New Zealand`$bar, country_plots$`New Zealand`$att, ncol = 1, heights = c(0.2, 1), align = "v"),
  nrow = 5, ncol = 4
)


## national alt--------------

country_region_map <- c("Australia" = "Oceania", "Brazil" = "Latin America and Caribbean",
  "Canada" = "North America", "Finland" = "Europe", "France" = "Europe",
  "Germany" = "Europe",  "Italy" = "Europe",  "Japan" = "Asia",
  "New Zealand" = "Oceania",  "Norway" = "Europe",
  "Poland" = "Europe",  "Russian Federation" = "Europe",
  "South Africa" = "Africa",  "Spain" = "Europe",  "Sweden" = "Europe",
  "United Kingdom" = "Europe",  "United States" = "North America",
  "Zimbabwe" = "Africa")

all_summaries <- set_names(iso.ls, country.ls) %>%
  imap_dfr(function(iso, country_name) {
    file_name <- gsub(" ", "_", country_name)
    file_path <- paste0("Outputs/Tables/Country/", file_name, "_summary_contrast.csv")
    read.csv(file_path) %>%
      mutate(
        iso = iso,
        country = country_name,
        region = country_region_map[country_name])
  })

national.sum <- all_summaries %>% filter(year %in% c(2010, 2020)) %>%
  mutate(region_order = case_when(
    region == "North America" ~ 1,
    region == "Europe" ~ 2,
    region == "Asia" ~ 3,
    region == "Latin America and Caribbean" ~ 4,
    region == "Africa" ~ 5,
    region == "Oceania" ~ 6,
    TRUE ~ 7),
  region = factor(region, 
    levels = c("North America", "Europe", "Asia", "Latin America and Caribbean", "Africa", "Oceania")),
    sig = if_else(att_ratio_mn < 1 & att_ratio_hi < 1 | att_ratio_mn > 1 & att_ratio_lo > 1,
                  "Yes", "No"))


eur.plt <- ggplot(filter(national.sum, region == "Europe"), aes(
  x = country, y = att_ratio_mn, colour = region, group = factor(year), alpha = sig)) +
  geom_hline(yintercept = 1.0, linetype = "solid", linewidth = 0.5) +
  geom_pointrange(aes(ymin = att_ratio_lo, ymax = att_ratio_hi),
                      position = position_dodge(width = 0.6, preserve = "single"), size = 0.65, 
                  linewidth = 0.6, shape = 16, colour = "#388E3C") +
  coord_cartesian(ylim = c(0,3), expand = F) +
  scale_alpha_manual(values = c(0.4, 1)) +
  xlab("Country") +ylab("Index") +
  theme_bw() +
    theme(legend.position = "none",
          plot.margin = margin(l = 5, r = 25, t = 5, unit = "pt"))

nat.plt <- ggplot(filter(national.sum, region != "Europe"), aes(
  x = reorder(country, region_order), y = att_ratio_mn, colour = region, group = factor(year), alpha = sig)) +
  geom_hline(yintercept = 1.0, linetype = "solid", linewidth = 0.5) +
  geom_pointrange(aes(ymin = att_ratio_lo, ymax = att_ratio_hi),
                  position = position_dodge(width = 0.6, preserve = "single"), size = 0.65, 
                  linewidth = 0.6, shape = 16) +
  coord_cartesian(ylim = c(0,3), expand = F) +
  scale_alpha_manual(values = c(0.4, 1)) +
  scale_color_manual(values = c("#8E44AD", "#00B4D8","#4E6E58", "#C0392B", "#D81B60")) +
  xlab("Country") +ylab("Index") +
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = margin(l = 5, r = 25, t = 5, unit = "pt"))

nat.arr <- ggarrange(eur.plt, nat.plt, nrow = 2, labels = c("A.", "B."), label.x = -0.01)

ggsave(path = "Outputs/Figures/Tidy", filename = "National.Est.v1.png",
       plot = nat.arr, width = 8, height = 4, units = "in",bg = "white")


## national alt 2--------------

country_region_map <- c("Australia" = "Oceania", "Brazil" = "Latin America and Caribbean",
                        "Canada" = "North America", "Finland" = "Europe", "France" = "Europe",
                        "Germany" = "Europe",  "Italy" = "Europe",  "Japan" = "Asia",
                        "New Zealand" = "Oceania",  "Norway" = "Europe",
                        "Poland" = "Europe",  "Russian Federation" = "Europe",
                        "South Africa" = "Africa",  "Spain" = "Europe",  "Sweden" = "Europe",
                        "United Kingdom" = "Europe",  "United States" = "North America",
                        "Zimbabwe" = "Africa")

all_summaries <- set_names(iso.ls, country.ls) %>%
  imap_dfr(function(iso, country_name) {
    file_name <- gsub(" ", "_", country_name)
    file_path <- paste0("Outputs/Tables/Country/", file_name, "_summary_contrast.csv")
    read.csv(file_path) %>%
      mutate(
        iso = iso,
        country = country_name,
        region = country_region_map[country_name])
  })


national.sum <- all_summaries %>% mutate(region_order = case_when(
    region == "North America" ~ 1,
    region == "Europe" ~ 2,
    region == "Asia" ~ 3,
    region == "Latin America and Caribbean" ~ 4,
    region == "Africa" ~ 5,
    region == "Oceania" ~ 6,
    TRUE ~ 7),
    region = factor(region, 
                    levels = c("North America", "Europe", "Asia", "Latin America and Caribbean", "Africa", "Oceania")),
    sig = if_else(att_ratio_mn < 1 & att_ratio_hi < 1 | att_ratio_mn > 1 & att_ratio_lo > 1,
                  TRUE, FALSE))

national_segments <- national.sum %>%
  mutate(year = as.numeric(as.character(year))) %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(xend = lead(year), yend = lead(att_ratio_mn),
    sig_seg = factor(
      ifelse(lead(sig) == TRUE | lead(sig) == "sig", "Significant", "Non-significant"),
      levels = c("Significant", "Non-significant"))) %>%
  filter(!is.na(xend)) %>%
  ungroup()

end_labels <- national_segments %>%
  group_by(country) %>%
  filter(xend == max(xend)) %>%
  ungroup()

eur.plt <- ggplot() +
  geom_hline(yintercept = 1.0, linetype = "solid", linewidth = 0.5) +
  geom_segment(data = filter(national_segments, region == "Europe"),
    aes( x = year, y = att_ratio_mn,xend = xend, yend = yend, 
      linetype = sig_seg), colour = "#388E3C", linewidth = 0.5) +
  geom_text_repel(data = filter(end_labels, region == "Europe"),
    aes(x = xend, y = yend, label = iso), colour = "black",
    hjust = -0.2, direction = "y",
    nudge_x = 0.5, segment.color = NA, 
    size = 2.5,  fontface = "bold") +
  scale_linetype_manual(name = "Significance",
    values = c("Significant" = "solid", "Non-significant" = "dotted")) +
  scale_y_log10(breaks = c(0.1, 0.3, 0.5, 1, 3, 5, 10), limits = c(NA, 12)) +
  labs(x = "Year", y = "Index") +
  theme_bw() +
  theme(legend.position = "none")

nat.plt <- ggplot() +
  geom_hline(yintercept = 1.0, linetype = "solid", linewidth = 0.5) +
  geom_segment(data = filter(national_segments, region != "Europe"),
               aes( x = year, y = att_ratio_mn,xend = xend, yend = yend, 
                    linetype = sig_seg, colour = region), linewidth = 0.5) +
  geom_text_repel(data = filter(end_labels, region != "Europe"),
                  aes(x = xend, y = yend, label = iso), colour = "black",
                  hjust = -0.2, direction = "y",
                  nudge_x = 0.5, segment.color = NA, 
                  size = 2.5,  fontface = "bold") +
  scale_linetype_manual(name = "Significance",
                        values = c("Significant" = "solid", "Non-significant" = "dotted")) +
  scale_color_manual(values = c("#8E44AD", "#00B4D8","#4E6E58", "#C0392B", "#D81B60")) +
  scale_y_log10(breaks = c(0.1, 0.3, 0.5, 1, 3, 5, 12), limits = c(NA, 10)) +
  labs(x = "Year", y = "Index") +
  theme_bw() +
  theme(legend.position = "none")


nat.arr <- ggarrange(eur.plt, nat.plt, nrow = 1, labels = c("A.", "B."), label.x = -0.01)

ggsave(path = "Outputs/Figures/Tidy", filename = "National.Est.ALT.v1.png",
       plot = nat.arr, width = 8, height = 3, units = "in",bg = "white")

## Numbers for the text ------------------------------------------------------
mod_dat_full %>% group_by(Region) %>% summarise(n_distinct(series))

nam_chained_sum %>% filter(year== max(year)) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

eur_chained_sum %>% filter(year== max(year)) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

asia_chained_sum %>% filter(year== max(year)) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

ant_chained_sum %>% filter(year== max(year)) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100)

nat.yr.sum <- all_summaries %>% group_by(country) %>% filter(year== max(year)) %>%
  mutate(att_ratio_mn = 100- att_ratio_mn*100,
         att_ratio_lo = 100- att_ratio_lo*100,
         att_ratio_hi = 100- att_ratio_hi*100) %>%
  ungroup() %>% arrange(region)

nat.yr.sum %>% filter(region == "Europe")

##
nrow(all_summaries)
nrow(all_summaries %>% filter(att_ratio_lo > 1))
nrow(all_summaries %>% filter(att_ratio_hi < 1))
914-19-340


region.n <- mod_dat_full %>% 
  group_by(Region, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= 10)) %>% pivot_wider(values_from = "n", names_from = "used") %>%
  rename("Region/Country" = "Region")

country.n <- mod_dat_full %>% 
  group_by(Country, used) %>% 
  summarise(n = n_distinct(series), .groups = "drop_last") %>% 
  filter(n() == 2, all(n >= 10)) %>% pivot_wider(values_from = "n", names_from = "used") %>%
  rename("Region/Country" = "Country")

write.csv(rbind(region.n, country.n), "Outputs/Tables/SM/Country.Region.n.summary.csv")
