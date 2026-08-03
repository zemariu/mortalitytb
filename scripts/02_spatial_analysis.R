# Article: Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024:
#          a Bayesian spatiotemporal modelling study
# Repository: mortalitytb
# Script: Calculates municipal mortality rates, evaluates spatial autocorrelation, 
#         and produces descriptive mortality maps.
# Author: José Mário Nunes da Silva, PhD
# Last updated: 2026-08-01

# Packages ----
library(tidyverse)
library(sf)
library(spdep)
library(geobr)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(patchwork)
library(flextable)
library(officer)

# Output directories ----
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/spatial/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/spatial/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("objects", recursive = TRUE, showWarnings = FALSE)


# Map of Brazil----

# South American countries
south_america <- ne_countries(
  continent = "south america",
  returnclass = "sf"
)

# Brazilian states
states <- read_state(year = 2020) 

states <- states |>
  dplyr::mutate(
    region = factor(
      case_when(
        name_region == "Norte" ~ "North",
        name_region == "Nordeste" ~ "Northeast",
        name_region == "Sudeste" ~ "Southeast",
        name_region == "Sul" ~ "South",
        name_region == "Centro Oeste" ~ "Central-West"
      ),
      levels = c("North", "Northeast", "Southeast", 
                 "South", "Central-West")
    )
  )

# Brazilian municipalities
municipalities <- read_municipality(
  code_muni = "all",
  year = 2022
)

# Create state centroids 
state_centroids <- states |>
  st_centroid() |>
  st_coordinates() |>
  as.data.frame()

# Map of Brazil
gg_brazil <- ggplot() +
  
  # South America background
  geom_sf(
    data = south_america,
    aes(shape = "South America"),
    fill = "gray80",
    color = "gray10",
    linewidth = 0.2
  ) +
  
  # States (colored by region)
  geom_sf(
    data = states,
    aes(fill = region, shape = "Federative Units"),
    linewidth = 0.4
  ) +
  
  # Municipalities
  geom_sf(
    data = municipalities,
    aes(shape = "Municipalities"),
    color = "grey20",
    linewidth = 0.15
  ) +
  
  # State labels
  geom_sf_label(
    data = states,
    aes(label = name_state),
    size = 3,
    fill = "white",
    label.size = 0.1,
    fun.geometry = sf::st_point_on_surface
  ) +
  
  # Fill scale (regions)
  scale_fill_brewer(
    name = "Regions of Brazil:",
    palette = "Set2"
  ) +
  
  # Shape scale (map elements)
  scale_shape_manual(
    name = "",
    values = c(1, 2, 3),
    guide = guide_legend(
      override.aes = list(
        color = c("black", "grey20", "black"),
        fill = c("white", "white", "grey80")
      )
    )
  ) +
  
  # Ocean labels
  annotate("text", x = -35, y = -31, label = "Atlantic Ocean", size = 4) +
  annotate("text", x = -78, y = -31, label = "Pacific Ocean", size = 4) +
  
  # Map extent
  coord_sf(
    xlim = c(-80, -30),
    ylim = c(-35, 5)
  ) +
  
  ggspatial::annotation_scale(
    location = "br",
    text_cex = 0.8,
    height = unit(0.25, "cm"),
    pad_y = unit(0.6, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    style = ggspatial::north_arrow_fancy_orienteering(),
    location = "br",
    width = unit(1.2, "cm"),
    height = unit(1.5, "cm"),
    pad_x = unit(1.2, "cm"),
    pad_y = unit(1.2, "cm")
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "cm")
  )

# Export figure
ggsave(
  "results/spatial/figures/Brazil.tiff",
  plot = gg_brazil,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)

# Figure caption
figure_caption <- paste(
  "Fig S1. Geographic location of Brazil and its administrative divisions."
)

writeLines(
  figure_caption,
  "results/spatial/figures/Fig_S1_caption.txt"
)

# Input data ----
df_tb <- read_csv(
  "data/processed/tb_mortality.csv",
  show_col_types = FALSE
)

df_tb <- df_tb |>
  mutate(
    code_muni = str_sub(as.character(code_muni), 1, 6),
    age_group = as.character(age_group)
  )

# WHO standard population in three age groups ----
who_standard_3_age_groups <- tibble(
  age_group_3 = c("0–19", "20–59", "60+"),
  who_weight = c(0.3460789, 0.5344130, 0.1195082)
)

# Municipal age-standardised mortality rates by period ----
asmr_mun <- df_tb |>
  select(-any_of("pop_who")) |>
  mutate(
    period = case_when(
      between(year, 2010, 2014) ~ "2010-2014",
      between(year, 2015, 2019) ~ "2015-2019",
      between(year, 2020, 2024) ~ "2020-2024",
      TRUE ~ NA_character_
    ),
    age_group_3 = case_when(
      age_group %in% c("0-4", "5-9", "10-14", "15-19") ~ "0–19",
      age_group %in% c(
        "20-24", "25-29", "30-34", "35-39",
        "40-44", "45-49", "50-54", "55-59"
      ) ~ "20–59",
      age_group %in% c("60-64", "65-69", "70-74", "75-79", "80+") ~ "60+",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(period), !is.na(age_group_3)) |>
  group_by(period, code_muni, age_group_3) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    who_standard_3_age_groups,
    by = "age_group_3"
  ) |>
  mutate(
    weighted_rate = if_else(
      population > 0,
      (deaths / population) * 1e5 * who_weight,
      NA_real_
    )
  ) |>
  group_by(period, code_muni) |>
  summarise(
    rate = round(
      sum(weighted_rate, na.rm = TRUE) /
        sum(who_weight[!is.na(weighted_rate)]),
      1
    ),
    .groups = "drop"
  )

write_csv(
  asmr_mun,
  "data/processed/municipal_asmr_by_period.csv"
)

# Municipal and state geometries ----
brazil <- read_municipality(
  code_muni = "all",
  year = 2022,
  showProgress = FALSE
) |>
  dplyr::mutate(
    code_muni = str_sub(as.character(code_muni), 1, 6)
  )

states <- read_state(
  code_state = "all",
  year = 2020,
  showProgress = FALSE
)

period_levels <- c("2010-2014", "2015-2019", "2020-2024")

municipality_period_grid <- expand_grid(
  code_muni = brazil$code_muni,
  period = period_levels
)

# Spatial dataset for mapping ----
df_map <- brazil |>
  left_join(
    municipality_period_grid,
    by = "code_muni"
  ) |>
  left_join(
    asmr_mun,
    by = c("code_muni", "period")
  )

# Mortality-rate categories ----
rate_levels <- c(
  "0",
  "< 1.2",
  "1.2–1.9",
  "1.9–2.7",
  "2.7–4.2",
  "> 4.2"
)

rate_colours <- c(
  "0" = "#FFFFFF",
  "< 1.2" = "#FEF0D9",
  "1.2–1.9" = "#FDCC8A",
  "1.9–2.7" = "#FC8D59",
  "2.7–4.2" = "#E34A33",
  "> 4.2" = "#B30000"
)

df_map <- df_map |>
  mutate(
    rate_category = case_when(
      rate == 0 ~ "0",
      rate < 1.2 ~ "< 1.2",
      rate < 1.9 ~ "1.2–1.9",
      rate < 2.7 ~ "1.9–2.7",
      rate <= 4.2 ~ "2.7–4.2",
      rate > 4.2 ~ "> 4.2",
      TRUE ~ NA_character_
    ),
    rate_category = factor(
      rate_category,
      levels = rate_levels
    )
  )

rate_legend_counts <- df_map |>
  st_drop_geometry() |>
  dplyr::count(period, rate_category, .drop = FALSE) |>
  mutate(
    legend_label = paste0(
      rate_category,
      " [",
      scales::comma(n, big.mark = ",", decimal.mark = "."),
      "]"
    )
  )

legend_2010_2014 <- rate_legend_counts |>
  filter(period == "2010-2014") |>
  select(rate_category, legend_label) |>
  deframe()

legend_2015_2019 <- rate_legend_counts |>
  filter(period == "2015-2019") |>
  select(rate_category, legend_label) |>
  deframe()

legend_2020_2024 <- rate_legend_counts |>
  filter(period == "2020-2024") |>
  select(rate_category, legend_label) |>
  deframe()

# Common map theme ----
map_theme <- theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(
      face = "bold",
      size = 18,
      hjust = 0
    ),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(0.7, "cm")
  )

# Age-standardised mortality-rate map: 2010–2014 ----
gg_10_14 <- ggplot() +
  geom_sf(
    data = df_map |>
      filter(period == "2010-2014"),
    aes(fill = rate_category),
    colour = NA
  ) +
  geom_sf(
    data = states,
    aes(colour = "Federative units"),
    fill = NA,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    name = paste0(
      "Tuberculosis deaths per\n",
      "100,000 inhabitants\n",
      "[5,570 municipalities]"
    ),
    values = rate_colours,
    breaks = rate_levels,
    labels = legend_2010_2014[rate_levels],
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.5)
    )
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Federative units" = "black")
  ) +
  labs(
    title = "A) Age-standardized mortality rates",
    subtitle = "2010–2014"
  ) +
  map_theme

# Age-standardised mortality-rate map: 2015–2019 ----
gg_15_19 <- ggplot() +
  geom_sf(
    data = df_map |>
      filter(period == "2015-2019"),
    aes(fill = rate_category),
    colour = NA
  ) +
  geom_sf(
    data = states,
    aes(colour = "Federative units"),
    fill = NA,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    name = paste0(
      "Tuberculosis deaths per\n",
      "100,000 inhabitants\n",
      "[5,570 municipalities]"
    ),
    values = rate_colours,
    breaks = rate_levels,
    labels = legend_2015_2019[rate_levels],
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.5)
    )
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Federative units" = "black")
  ) +
  labs(subtitle = "2015–2019") +
  map_theme

# Age-standardised mortality-rate map: 2020–2024 ----
gg_20_24 <- ggplot() +
  geom_sf(
    data = df_map |>
      filter(period == "2020-2024"),
    aes(fill = rate_category),
    colour = NA
  ) +
  geom_sf(
    data = states,
    aes(colour = "Federative units"),
    fill = NA,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    name = paste0(
      "Tuberculosis deaths per\n",
      "100,000 inhabitants\n",
      "[5,570 municipalities]"
    ),
    values = rate_colours,
    breaks = rate_levels,
    labels = legend_2020_2024[rate_levels],
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(colour = "black", linewidth = 0.5)
    )
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Federative units" = "black")
  ) +
  labs(subtitle = "2020–2024") +
  map_theme

# Combined mortality-rate maps ----
map_rates_A <- gg_10_14 +
  gg_15_19 +
  gg_20_24 +
  plot_layout(ncol = 3)

map_rates_A

save(
  map_rates_A,
  file = "objects/spatial/map_rates_A.RData"
)

# Spatial-neighbourhood structure ----
brazil_projected <- st_transform(brazil, 5880)
municipality_points <- st_point_on_surface(brazil_projected)
municipality_coordinates <- st_coordinates(municipality_points)

neighbours_knn6 <- knearneigh(
  municipality_coordinates,
  k = 6
) |>
  knn2nb()

spatial_weights <- nb2listw(
  neighbours_knn6,
  style = "W",
  zero.policy = TRUE
)

# Overall municipal age-standardised mortality rates, 2010–2024 ----
asmr_mun_overall <- df_tb |>
  select(-any_of("pop_who")) |>
  mutate(
    age_group_3 = case_when(
      age_group %in% c("0-4", "5-9", "10-14", "15-19") ~ "0–19",
      age_group %in% c(
        "20-24", "25-29", "30-34", "35-39",
        "40-44", "45-49", "50-54", "55-59"
      ) ~ "20–59",
      age_group %in% c("60-64", "65-69", "70-74", "75-79", "80+") ~ "60+",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(age_group_3)) |>
  group_by(code_muni, age_group_3) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    who_standard_3_age_groups,
    by = "age_group_3"
  ) |>
  mutate(
    weighted_rate = if_else(
      population > 0,
      (deaths / population) * 1e5 * who_weight,
      NA_real_
    )
  ) |>
  group_by(code_muni) |>
  summarise(
    rate = round(
      sum(weighted_rate, na.rm = TRUE) /
        sum(who_weight[!is.na(weighted_rate)]),
      1
    ),
    .groups = "drop"
  )

# Overall spatial dataset ----
df_map_overall <- brazil |>
  left_join(
    asmr_mun_overall,
    by = "code_muni"
  )

# Global Moran's I helper ----
run_global_moran <- function(rate, period_label) {
  set.seed(2026)

  moran_result <- moran.mc(
    rate,
    spatial_weights,
    nsim = 999,
    zero.policy = TRUE
  )

  moran_i <- as.numeric(moran_result$statistic)
  z_score <- (moran_i - mean(moran_result$res)) / sd(moran_result$res)
  p_value <- moran_result$p.value

  distribution_pattern <- case_when(
    p_value >= 0.05 ~ "Random",
    moran_i > 0 ~ "Clustered",
    moran_i < 0 ~ "Dispersed",
    TRUE ~ "Random"
  )

  tibble(
    Period = period_label,
    `Moran's I` = round(moran_i, 3),
    `Z-score` = round(z_score, 2),
    `P-value` = p_value,
    `Distribution pattern` = distribution_pattern
  )
}

# Global spatial autocorrelation ----
moran_2010_2014 <- df_map |>
  filter(period == "2010–2014") |>
  pull(rate) |>
  run_global_moran("2010–2014")

moran_2015_2019 <- df_map |>
  filter(period == "2015–2019") |>
  pull(rate) |>
  run_global_moran("2015–2019")

moran_2020_2024 <- df_map |>
  filter(period == "2020–2024") |>
  pull(rate) |>
  run_global_moran("2020–2024")

moran_2010_2024 <- df_map_overall |>
  pull(rate) |>
  run_global_moran("2010–2024")

moran_results <- bind_rows(
  moran_2010_2014,
  moran_2015_2019,
  moran_2020_2024,
  moran_2010_2024
)

write_csv(
  moran_results,
  "results/spatial/tables/Table_S11_global_moran.csv"
)

# Global Moran's I table ----
table_s11 <- flextable(moran_results) |>
  set_caption(
      caption = as_paragraph(
        as_chunk("Table S11. ", props = fp_text(bold = TRUE,  font.size = 12)),
        as_chunk("Global spatial autocorrelation analysis of tuberculosis mortality in Brazil, 2010–2024.")
      )
  ) |>
  colformat_num(
    j = "Moran's I",
    digits = 3
  ) |>
  colformat_num(
    j = "Z-score",
    digits = 2
  ) |>
  colformat_num(
    j = "P-value",
    digits = 3
  ) |>
  bold(part = "header") |>
  theme_booktabs() |>
  bold(part = "header") |> 
  autofit()

table_s11

save_as_docx(
  `Table S11` = table_s11,
  path = "results/spatial/tables/Table_s11.docx"
)
