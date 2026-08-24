library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)

load(file = "Data/data_models.RData")

## Summarise and tidy naming ---------------------------------------------------
mod_dat_full %>% group_by(Utilised) %>% summarise(length(unique(ID)))
mod_dat_full %>% group_by(Utilised) %>% summarise(length(unique(Binomial)))

mod_dat_full %>% group_by(threats) %>% summarise(length(unique(ID)))
mod_dat_full %>% group_by(threats) %>% summarise(length(unique(Binomial)))

mod_dat_full %>% group_by(Utilised) %>% summarise(mean(threats))

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
  custom_match = c("International Waters" = NA_character_,
    "Reunion" = "REU", "R?Runion" = "REU", "Curacao" = "CUW",
    "Aland Islands" = "ALA"))

## Mapping series coverage and use ---------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

# Summarize and create 4x4 classification
national.sum <- mod_dat_full %>% 
  group_by(iso3c) %>%
  summarise(n.series = n_distinct(ID),
    n.used.series = n_distinct(ID[used == 1])) %>%
  mutate(perc.used = (n.used.series / n.series) * 100,
         series.cat = case_when(n.series < 6 ~ "up to 5",
                                n.series > 5 & n.series <= 50 ~ "6-50",
                                n.series > 50 & n.series <= 100 ~ "51-100",
                                n.series > 100 ~ ">100"),
         series.cat = factor(series.cat, levels = c("up to 5", "6-50", "51-100", ">100")),
         var1_num = as.numeric(series.cat),
         perc.cat = cut(perc.used,breaks = c(-Inf, 25, 50, 75, 100),
                        labels = c("0-25%", "26-50%", "51-75%", "76-100%"),include.lowest = TRUE),
         var2_num = as.numeric(perc.cat),
         bi_class = paste0(var1_num, "-", var2_num))

d <- national.sum %>% left_join(world, by = c("iso3c" = "iso_a3_eh"))

# Build 4x4 Color Matrix by 2D Interpolation
create_2d_palette <- function(c11, c41, c14, c44, n = 4) {
  left_col <- colorRampPalette(c(c11, c14))(n)
  right_col <- colorRampPalette(c(c41, c44))(n)
  
  mat <- matrix(NA_character_, nrow = n, ncol = n)
  for (i in 1:n) {
    mat[i, ] <- colorRampPalette(c(left_col[i], right_col[i]))(n)
  }
  return(mat)
}

# Define 4 corner colors: [Low/Low, High Var1/Low Var2, Low Var1/High Var2, High/High]
pal_mat <- create_2d_palette(c11 = "#f0e4d7", # (Low Series, Low % - distinct from NA)
  c41 = "#d95f02", # (High Series, Low %)
  c14 = "#1b9e77", # (Low Series, High %)
  c44 = "#3b0040", # (High Series, High %)
  n = 4)

# Convert 4x4 matrix to lookup dataframe
pal_df <- expand.grid(var1_num = 1:4, var2_num = 1:4) %>%
  mutate(
    bi_class = paste0(var1_num, "-", var2_num),
    fill_color = apply(., 1, function(row) {
      x <- as.numeric(row["var1_num"])
      y <- as.numeric(row["var2_num"])
      pal_mat[y, x]
    })
  )

bi_pal <- setNames(pal_df$fill_color, pal_df$bi_class)

# Join with spatial world object
world_map <- world %>%
  left_join(national.sum, by = c("iso_a3_eh" = "iso3c"))

# Create Main Map
map_p <- ggplot(world_map) +
  geom_sf(aes(fill = bi_class), color = "grey30", size = 0.1) +
  scale_fill_manual(values = bi_pal, na.value = "grey90", na.translate = FALSE) +
  coord_sf(crs = "+proj=moll") +
  theme_minimal() +
  theme(legend.position = "none",  panel.grid = element_blank(),
    axis.text = element_blank(),  axis.ticks = element_blank())

# Create 4x4 Legend Plot
legend_p <- ggplot(pal_df, aes(x = factor(var1_num), y = factor(var2_num), fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  scale_x_discrete(labels = c("up to 5", "6-50", "51-100", ">100")) +
  scale_y_discrete(labels = c("0-25%", "26-50%", "51-75%", "76-100%")) +
  labs(x = "Total Series", y = "% Used") +
  theme_minimal() +
  theme(axis.title = element_text(size = 10, face = "plain"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(), aspect.ratio = 1)

# Combine map + legend
final_map <- map_p + inset_element(legend_p, left = 0.02, bottom = 0.05, right = 0.30, top = 0.40)
final_map

## Taxonomic diversification ---------------------------------------------------
mod_dat_full %>% 
  group_by(Taxon, used) %>%
  summarise(n.series = length(unique(ID)), .groups = "drop")

taxa.summary <- mod_dat_full %>% 
  group_by(Taxon, Class, used) %>%
  summarise(n.series = length(unique(ID)), .groups = "drop") %>%
  bind_rows(data.frame(Class = c("Dipneusti", "Holocephali", "Myxini"), 
                       Taxon = "Fish",used = c(0, 1, 0))) %>%
  mutate(n.series = ifelse(is.na(n.series), 0, n.series),
         used.cat = factor(ifelse(used == 1, "Used", "Not used"), 
                           levels = c("Used", "Not used"))) %>% 
  arrange(Taxon, Class) %>%
  mutate(Class = factor(Class, levels = unique(Class)))

muted_9 <- c( "#D95F02", "#D4A373", "#2B5C70", "#3B6982", "#4A8B99", 
              "#1B9E77", "#6A9360", "#7570B3", "#A08B76" )


taxa.plt <-  ggplot(taxa.summary, aes(x = Class, y = n.series, fill = Class, color = Class)) +
  geom_col(aes(alpha = used.cat),
           position = position_dodge(width = 0.8, preserve = "single"), width = 0.7) +
  geom_text(aes(label = n.series, group = used.cat),
            position = position_dodge(width = 0.8, preserve = "single"), # Matches bar dodge position
            vjust = -0.5,size = 3,color = "black",show.legend = FALSE) +
  scale_fill_manual(values = muted_9, guide = "none") +   
  scale_color_manual(values = muted_9, guide = "none") +  
  scale_alpha_manual( values = c("Used" = 1, "Not used" = 0), 
                      name = "Status",guide = guide_legend(override.aes = list(fill = "grey40", color = "black"))) +
  scale_y_continuous(breaks = c(0, 250, 500, 750)) +
  annotate( "segment",  x = 2.5, xend = 7.5, y = 900, yend = 900, 
            color = "black",  linewidth = 0.6) +
  annotate( "segment",  x = 2.5, xend = 2.5, y = 900, yend = 850, 
            color = "black",  linewidth = 0.6) +
  annotate( "segment",  x = 7.5, xend = 7.5, y = 900, yend = 850, 
            color = "black",  linewidth = 0.6) +
  annotate( "text", label = "Fish",  x = 5, y = 975, 
            color = "black", fontface = "bold", size = 3.5) +
  ylab("Total series") +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black",
                                   face = c("bold", "bold", "plain", "plain",
                                            "plain", "plain", "plain", "bold", "bold")),
        legend.position = "none")

library(ggpubr)

descr.arr <- ggarrange(final_map, taxa.plt, nrow = 2, ncol = 1,
                       heights = c(1, .7), labels = c("A.", "B."))

ggsave(path = "Outputs/Figures", filename = "rawdata.summary.v1.png",
       plot = descr.arr, width = 8, height = 8, units = "in",bg = "white")


