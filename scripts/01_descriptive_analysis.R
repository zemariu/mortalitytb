# Article: Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024:
# a Bayesian spatiotemporal modelling study
# Repository: mortalitytb
# Script: Prepares the mortality and population data and produces the descriptive 
#         epidemiological indicators used in subsequent analyses.
# Author: José Mário Nunes da Silva, PhD
# Last updated: 2026-08-01

# Packages ----
library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)
library(patchwork)

# Output directories ----
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/descriptive/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/descriptive/figures", recursive = TRUE, showWarnings = FALSE)

# Import mortality data ----
tb_deaths <- read_delim(
  "data/SIM_TB.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
)

# Recode mortality variables ----
age_group_levels <- c(
  "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
  "60-64", "65-69", "70-74", "75-79", "80+"
)

age_group_10_levels <- c(
  "0-9", "10-19", "20-29", "30-39", "40-49",
  "50-59", "60-69", "70-79", "80+"
)

region_levels <- c(
  "North", "Northeast", "Southeast", "South", "Central-West"
)

tb_deaths <- tb_deaths |>
  select(
    CODMUNRES, ANO, REGIAO, IDADE, FAIXA_ETARIA, SEXO,
    RACACOR, ESTCIV, ESC, LOCOCOR, CAUSABAS
  ) |>
  mutate(
    muni_code_original = sprintf("%06.0f", CODMUNRES),
    uf = substr(muni_code_original, 1, 2),
    code_muni = if_else(
      str_ends(muni_code_original, "0000"),
      NA_character_,
      muni_code_original
    ),
    year = as.integer(ANO),
    period = factor(
      case_when(
        year <= 2014 ~ "2010-2014",
        year <= 2019 ~ "2015-2019",
        year <= 2024 ~ "2020-2024",
        TRUE ~ NA_character_
      ),
      levels = c("2010-2014", "2015-2019", "2020-2024")
    ),
    region = factor(
      REGIAO,
      levels = c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste"),
      labels = region_levels
    ),
    age = suppressWarnings(as.numeric(IDADE)),
    age_group = case_when(
      FAIXA_ETARIA == "0 a 4 anos" ~ "0-4",
      FAIXA_ETARIA == "5 a 9 anos" ~ "5-9",
      FAIXA_ETARIA == "10 a 14 anos" ~ "10-14",
      FAIXA_ETARIA == "15 a 19 anos" ~ "15-19",
      FAIXA_ETARIA == "20 a 24 anos" ~ "20-24",
      FAIXA_ETARIA == "25 a 29 anos" ~ "25-29",
      FAIXA_ETARIA == "30 a 34 anos" ~ "30-34",
      FAIXA_ETARIA == "35 a 39 anos" ~ "35-39",
      FAIXA_ETARIA == "40 a 44 anos" ~ "40-44",
      FAIXA_ETARIA == "45 a 49 anos" ~ "45-49",
      FAIXA_ETARIA == "50 a 54 anos" ~ "50-54",
      FAIXA_ETARIA == "55 a 59 anos" ~ "55-59",
      FAIXA_ETARIA == "60 a 64 anos" ~ "60-64",
      FAIXA_ETARIA == "65 a 69 anos" ~ "65-69",
      FAIXA_ETARIA == "70 a 74 anos" ~ "70-74",
      FAIXA_ETARIA == "75 a 79 anos" ~ "75-79",
      FAIXA_ETARIA == "80 anos e mais" ~ "80+",
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = age_group_levels),
    age_group_10 = case_when(
      age_group %in% c("0-4", "5-9") ~ "0-9",
      age_group %in% c("10-14", "15-19") ~ "10-19",
      age_group %in% c("20-24", "25-29") ~ "20-29",
      age_group %in% c("30-34", "35-39") ~ "30-39",
      age_group %in% c("40-44", "45-49") ~ "40-49",
      age_group %in% c("50-54", "55-59") ~ "50-59",
      age_group %in% c("60-64", "65-69") ~ "60-69",
      age_group %in% c("70-74", "75-79") ~ "70-79",
      age_group == "80+" ~ "80+",
      TRUE ~ NA_character_
    ),
    age_group_10 = factor(age_group_10, levels = age_group_10_levels),
    sex = factor(
      case_when(
        SEXO == "Masculino" ~ "Male",
        SEXO == "Feminino" ~ "Female",
        TRUE ~ NA_character_
      ),
      levels = c("Male", "Female")
    ),
    race = factor(
      case_when(
        RACACOR == "Branca" ~ "White",
        RACACOR == "Preta" ~ "Black",
        RACACOR == "Parda" ~ "Brown",
        RACACOR == "Amarela" ~ "Yellow",
        RACACOR == "Indígena" ~ "Indigenous",
        TRUE ~ NA_character_
      ),
      levels = c("White", "Black", "Brown", "Yellow", "Indigenous")
    ),
    marital_status = factor(
      case_when(
        ESTCIV %in% c("Casado", "União consensual") ~ "Married",
        ESTCIV %in% c(
          "Outro", "Separado judicialmente", "Solteiro", "Viúvo"
        ) ~ "Unmarried",
        TRUE ~ NA_character_
      ),
      levels = c("Unmarried", "Married")
    ),
    education = factor(
      case_when(
        ESC == "Nenhuma" ~ "Illiterate",
        ESC %in% c("1 a 3 anos", "4 a 7 anos", "1o grau") ~ "0–7 years",
        ESC %in% c(
          "8 a 11 anos", "9 a 11 anos", "12 anos ou mais",
          "2o grau", "Superior"
        ) ~ "8+ years",
        TRUE ~ NA_character_
      ),
      levels = c("Illiterate", "0–7 years", "8+ years")
    ),
    place_of_death = factor(
      case_when(
        LOCOCOR == "Domicílio" ~ "Home",
        LOCOCOR %in% c(
          "Hospital", "Outro estabelecimento de saúde"
        ) ~ "Healthcare facility",
        LOCOCOR %in% c("Outros", "Via pública") ~ "Other",
        TRUE ~ NA_character_
      ),
      levels = c("Home", "Healthcare facility", "Other")
    ),
    cause_code = str_to_upper(
      str_remove_all(as.character(CAUSABAS), "[^A-Za-z0-9]")
    ),
    anatomical_classification = factor(
      case_when(
        str_detect(cause_code, "^(010|011|012|A15|A16)") ~ "Pulmonary",
        str_detect(cause_code, "^(013|A17)") ~ "Central nervous system",
        str_detect(cause_code, "^(018|A19)") ~ "Disseminated (miliary)",
        str_detect(cause_code, "^(014|015|016|017|A18)") ~ "Other organs",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Pulmonary",
        "Central nervous system",
        "Disseminated (miliary)",
        "Other organs"
      )
    )
  )

# Descriptive table ----
descriptive_data <- tb_deaths |>
  select(
    period,
    region,
    age,
    age_group_10,
    sex,
    race,
    marital_status,
    education,
    place_of_death,
    anatomical_classification
  )

descriptive_labels <- list(
  region ~ "Region of residence",
  age ~ "Age, years",
  age_group_10 ~ "Age group, years",
  sex ~ "Sex",
  race ~ "Race/ethnicity",
  marital_status ~ "Marital status",
  education ~ "Education, years",
  place_of_death ~ "Place of death",
  anatomical_classification ~ "Anatomical classification"
)

table_s9 <- descriptive_data |>
  tbl_summary(
    by = period,
    label = descriptive_labels,
    type = list(age ~ "continuous2"),
    statistic = list(
      all_categorical() ~ "{n} ({p})",
      age ~ c("{mean} ({sd})", "{median} ({p25}–{p75})")
    ),
    digits = list(
      all_categorical() ~ c(0, 1),
      age ~ c(1, 1, 1)
    ),
    missing = "ifany",
    missing_stat = "{N_miss} ({p_miss}%)",
    missing_text = "Missing"
  ) |>
  bold_labels() |>
  add_overall() |>
  add_stat_label(label = age ~ c("Mean (SD)", "Median (IQR)")) |>
  modify_header(label ~ "**Variables**") |>
  modify_spanning_header(all_stat_cols() ~ "**Period of death**") |>
  modify_caption(
    "**Table S9.** Epidemiological characteristics of tuberculosis-related deaths in Brazil, 2010–2024."
  )

table_s9

table_s9 |>
  as_flex_table() |>
  save_as_docx(path = "results/descriptive/tables/table_s9_descriptive_deaths.docx")

# Missing-data diagnostics ----
missing_before_imputation <- tibble(
  variable = c("Municipality code", "Age group", "Sex", "Race/ethnicity"),
  missing = c(
    sum(is.na(tb_deaths$code_muni)),
    sum(is.na(tb_deaths$age_group)),
    sum(is.na(tb_deaths$sex)),
    sum(is.na(tb_deaths$race))
  ),
  percentage = round(100 * missing / nrow(tb_deaths), 2)
)

print(missing_before_imputation)

# Hierarchical proportional imputation ----
impute_hierarchically <- function(data, variable, hierarchy) {
  original_levels <- if (is.factor(data[[variable]])) {
    levels(data[[variable]])
  } else {
    NULL
  }

  for (group_variables in hierarchy) {
    data <- data |>
      group_by(across(all_of(group_variables))) |>
      mutate(
        "{variable}" := {
          values <- .data[[variable]]
          missing_index <- which(is.na(values))
          observed_values <- values[!is.na(values)]

          if (length(missing_index) > 0 && length(observed_values) > 0) {
            values[missing_index] <- sample(
              observed_values,
              size = length(missing_index),
              replace = TRUE
            )
          }

          values
        }
      ) |>
      ungroup()
  }

  values <- data[[variable]]
  missing_index <- which(is.na(values))
  observed_values <- values[!is.na(values)]

  if (length(missing_index) > 0 && length(observed_values) > 0) {
    values[missing_index] <- sample(
      observed_values,
      size = length(missing_index),
      replace = TRUE
    )
  }

  if (!is.null(original_levels)) {
    values <- factor(values, levels = original_levels)
  }

  data[[variable]] <- values
  data
}

set.seed(123)

tb_deaths_imputed <- tb_deaths |>
  impute_hierarchically(
    variable = "code_muni",
    hierarchy = list(c("uf", "year"), "uf")
  ) |>
  impute_hierarchically(
    variable = "age_group",
    hierarchy = list(
      c("code_muni", "year"),
      "code_muni",
      c("uf", "year"),
      "uf"
    )
  ) |>
  impute_hierarchically(
    variable = "sex",
    hierarchy = list(
      c("code_muni", "year"),
      "code_muni",
      c("uf", "year"),
      "uf"
    )
  )

missing_after_imputation <- tibble(
  variable = c("Municipality code", "Age group", "Sex"),
  missing = c(
    sum(is.na(tb_deaths_imputed$code_muni)),
    sum(is.na(tb_deaths_imputed$age_group)),
    sum(is.na(tb_deaths_imputed$sex))
  )
)

print(missing_after_imputation)

# WHO standard population ----
who_standard <- tibble(
  age_group = age_group_levels,
  who_weight = c(
    0.0886, 0.0869, 0.0860, 0.0847, 0.0822, 0.0793,
    0.0761, 0.0715, 0.0659, 0.0604, 0.0537, 0.0455,
    0.0372, 0.0296, 0.0221, 0.0152, 0.01545
  )
)

# Municipality-level mortality dataset ----
population_data <- read_csv(
  "data/pop_br.csv",
  show_col_types = FALSE
) |>
  mutate(
    year = as.integer(year),
    region = as.character(region),
    uf = str_pad(as.character(uf), width = 2, pad = "0"),
    code_muni = str_pad(
      as.character(code_muni),
      width = 6,
      pad = "0"
    ),
    sex = as.character(sex),
    age_group = as.character(age_group)
  )

mortality_counts <- tb_deaths_imputed |>
  transmute(
    year,
    region = as.character(region),
    uf,
    code_muni,
    sex = as.character(sex),
    age_group = as.character(age_group)
  ) |>
  count(
    year,
    region,
    uf,
    code_muni,
    sex,
    age_group,
    name = "deaths"
  )

tb_mortality <- population_data |>
  left_join(
    mortality_counts,
    by = c(
      "year", "region", "uf", "code_muni", "sex", "age_group"
    )
  ) |>
  mutate(
    deaths = replace_na(deaths, 0L),
    age_group_10 = case_when(
      age_group %in% c("0-4", "5-9") ~ "0-9",
      age_group %in% c("10-14", "15-19") ~ "10-19",
      age_group %in% c("20-24", "25-29") ~ "20-29",
      age_group %in% c("30-34", "35-39") ~ "30-39",
      age_group %in% c("40-44", "45-49") ~ "40-49",
      age_group %in% c("50-54", "55-59") ~ "50-59",
      age_group %in% c("60-64", "65-69") ~ "60-69",
      age_group %in% c("70-74", "75-79") ~ "70-79",
      age_group == "80+" ~ "80+",
      TRUE ~ NA_character_
    ),
    age_group_10 = factor(age_group_10, levels = age_group_10_levels)
  ) |>
  left_join(who_standard, by = "age_group")

stopifnot(sum(tb_mortality$deaths) == nrow(tb_deaths_imputed))

write_csv(tb_mortality, "data/processed/tb_mortality.csv")

# Race-specific mortality dataset ----
race_population <- read_csv(
  "data/pop_race.csv",
  show_col_types = FALSE
  ) |>
  mutate(
    year = as.integer(year),
    sex = as.character(sex),
    race = as.character(race),
    age_group = as.character(age_group),
    population = as.numeric(population)
  ) |> 
  mutate(race = factor(race,
                       levels = c("White", "Black", "Brown",
                                  "Yellow", "Indigenous")
    )
 )

## Figure S2: race-specific population estimates

age_group_colours <- c(
  "#00468B", "#ED0000", "#42B540", "#FFC20A",
  "#0099B4", "#925E9F", "#FDAF91", "#AD002A", 
  "#ADB6B6", "#1B1919", "#597DBF", "#7F5112", 
  "#E76BF3", "#00A087", "#C1E168", "#DC0000",
  "#999933", "#00C2F2", "#B24745", "#676767"
)

figure_s2 <- ggplot(
  race_population,
  aes(
    x = year,
    y = population,
    colour = age_group
  )
) +
  geom_line(
    aes(
      linetype = "Estimated"),
    alpha = 0.4,
    linewidth = 0.7
  ) +
  geom_smooth(
    aes(
      linetype = "Smoothed (GAM)"),
    method = "gam",
    formula = y ~ s(x, k = 8),
    se = FALSE,
    linewidth = 1.1
  ) +
  geom_point(
    data = race_population |>
      filter(
        year %in% c(2010L, 2022L)
      ),
    aes(
      shape = "Population census"),
    colour = "black",
    fill = "grey",
    size = 1.5
  ) +
  facet_grid(
    rows = vars(race),
    cols = vars(sex),
    scales = "free_y",
    switch = "y"
  ) +
  scale_colour_manual(
    values = age_group_colours
  ) +
  scale_linetype_manual(
    values = c(
      "Estimated" = "dashed",
      "Smoothed (GAM)" = "solid")
  ) +
  scale_shape_manual(
    values = c(
      "Population census" = 21)
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ""),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  scale_x_continuous(
    limits = c(2010, 2024),
    breaks = seq(2010, 2024, by = 2)
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(
        colour = "black",
        linewidth = 1,
        alpha = 1)
    )
  ) +
  labs(
    x = "Year",
    y = "Race or ethnicity",
    colour = "Age group",
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(colour = "black"),
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.text.y.left = element_text(size = 14),
    strip.placement = "outside",
    strip.background.y = element_blank()
  )

ggsave(
  filename = "results/descriptive/figures/Fig_S2.tiff",
  plot = figure_s2,
  device = "tiff",
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  figure_caption,
  "results/descriptive/figures/Fig_S2_caption.txt"
)

## Population by race and age group
race_population_by_age <- race_population |>
  group_by(
    year,
    race,
    age_group
  ) |>
  summarise(
    population = sum(
      population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(
    year,
    race,
    age_group
  )

readr::write_csv(
  race_population_by_age,
  "data/processed/population_by_race.csv"
)

## Tuberculosis deaths by race and age group
race_mortality_counts <- tb_deaths_imputed |>
  filter(!is.na(race)) |>
  transmute(
    year,
    race = as.character(race),
    age_group = as.character(age_group)
  ) |>
  count(year, race, age_group, name = "deaths")

tb_mortality_race <- race_population_by_age |>
  left_join(
    race_mortality_counts,
    by = c("year", "race", "age_group")
  ) |>
  mutate(deaths = replace_na(deaths, 0L)) |>
  left_join(who_standard, by = "age_group")

stopifnot(
  sum(tb_mortality_race$deaths) == sum(!is.na(tb_deaths_imputed$race))
)

write_csv(tb_mortality_race, "data/processed/tb_mortality_race.csv")

# Overall age-standardised mortality rates ----
asmr_overall <- tb_mortality |>
  group_by(year, age_group) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    who_weight = first(who_weight),
    .groups = "drop"
  ) |>
  mutate(
    weighted_rate = (deaths / population) * 100000 * who_weight
  ) |>
  group_by(year) |>
  summarise(
    rate = round(sum(weighted_rate) / sum(who_weight), 1),
    .groups = "drop"
  )

# Sex-specific age-standardised mortality rates ----
asmr_sex <- tb_mortality |>
  group_by(year, sex, age_group) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    who_weight = first(who_weight),
    .groups = "drop"
  ) |>
  mutate(
    weighted_rate = (deaths / population) * 100000 * who_weight
  ) |>
  group_by(year, sex) |>
  summarise(
    rate = round(sum(weighted_rate) / sum(who_weight), 1),
    .groups = "drop"
  )

asmr_sex_wide <- asmr_sex |>
  mutate(series = paste0("sex_", str_to_lower(sex))) |>
  select(year, series, rate) |>
  pivot_wider(names_from = series, values_from = rate)

# Region-specific age-standardised mortality rates ----
asmr_region <- tb_mortality |>
  group_by(year, region, age_group) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    who_weight = first(who_weight),
    .groups = "drop"
  ) |>
  mutate(
    weighted_rate = (deaths / population) * 100000 * who_weight
  ) |>
  group_by(year, region) |>
  summarise(
    rate = round(sum(weighted_rate) / sum(who_weight), 1),
    .groups = "drop"
  )

asmr_region_wide <- asmr_region |>
  mutate(
    series = paste0(
      "region_",
      region |>
        str_to_lower() |>
        str_replace_all("-", "_")
    )
  ) |>
  select(year, series, rate) |>
  pivot_wider(names_from = series, values_from = rate)

# Federative-unit-specific age-standardised mortality rates ----

asmr_uf <- tb_mortality |>
  
  # Sort states by region
  mutate(
    region = factor(region, levels = region_levels)) |> 
  arrange(region, uf_name) |> 
  mutate(
    uf_name = factor(uf_name, levels = unique(uf_name))) |> 
  
  group_by(year, uf_name, age_group) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    who_weight = first(who_weight),
    .groups = "drop"
  ) |>
  mutate(
    weighted_rate = (deaths / population) * 100000 * who_weight
  ) |>
  group_by(year, uf_name) |>
  summarise(
    rate = round(sum(weighted_rate) / sum(who_weight), 1),
    .groups = "drop"
  )

asmr_uf_wide <- asmr_uf |>
  pivot_wider(names_from = uf_name, values_from = rate)

# Race-specific age-standardised mortality rates ----
asmr_race <- tb_mortality_race |>
  group_by(year, race, age_group) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    who_weight = first(who_weight),
    .groups = "drop"
  ) |>
  mutate(
    weighted_rate = (deaths / population) * 100000 * who_weight
  ) |>
  group_by(year, race) |>
  summarise(
    rate = round(sum(weighted_rate) / sum(who_weight), 1),
    .groups = "drop"
  )

asmr_race_wide <- asmr_race |>
  mutate(series = paste0("race_", str_to_lower(race))) |>
  select(year, series, rate) |>
  pivot_wider(names_from = series, values_from = rate)

# Age-specific mortality rates ----
age_specific_rates <- tb_mortality |>
  group_by(year, age_group_10) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    rate = round((deaths / population) * 100000, 1),
    rate = if_else(!is.na(rate) & rate == 0, 0.1, rate)
  )

age_specific_rates_wide <- age_specific_rates |>
  mutate(
    series = paste0(
      "age_",
      age_group_10 |>
        as.character() |>
        str_replace_all("-", "_") |>
        str_replace("\\+", "_plus")
    )
  ) |>
  select(year, series, rate) |>
  pivot_wider(names_from = series, values_from = rate)

# Joinpoint input dataset ----
joinpoint_rates <- asmr_overall |>
  rename(overall_asmr = rate) |>
  left_join(asmr_sex_wide, by = "year") |>
  left_join(age_specific_rates_wide, by = "year") |>
  left_join(asmr_race_wide, by = "year") |>
  left_join(asmr_region_wide, by = "year") |>
  left_join(asmr_uf_wide, by = "year") |>
  arrange(year)

write_csv(joinpoint_rates, "data/processed/joinpoint_mortality_rates.csv")

# Heatmap data ----
state_order <- tb_mortality |>
  distinct(region, uf_name) |>
  mutate(region = factor(region, levels = region_levels)) |>
  arrange(region, uf_name) |>
  pull(uf_name)

heatmap_top <- tb_mortality |>
  group_by(year) |>
  summarise(deaths = sum(deaths), .groups = "drop") |>
  left_join(asmr_overall, by = "year")

heatmap_central <- asmr_uf |>
  mutate(
    uf_name = factor(uf_name, levels = state_order)
  ) |>
  group_by(uf_name) |>
  mutate(rate_scaled = as.numeric(scale(rate))) |>
  ungroup()

heatmap_side <- tb_mortality |>
  mutate(
    region = factor(region, levels = region_levels),
    uf_name = factor(uf_name, levels = state_order)
  ) |>
  group_by(year, uf_name) |>
  summarise(deaths = sum(deaths), .groups = "drop")

# Heatmap figure ----
rate_multiplier <- max(heatmap_top$deaths) / max(heatmap_top$rate)

top_panel <- ggplot(heatmap_top, aes(x = year)) +
  geom_col(aes(y = deaths), fill = "#DE2D26", alpha = 0.85) +
  geom_line(aes(y = rate * rate_multiplier), linewidth = 0.9) +
  geom_point(aes(y = rate * rate_multiplier), size = 1.8) +
  geom_text(
    aes(y = rate * rate_multiplier, label = sprintf("%.1f", rate)),
    vjust = -0.9,
    size = 3,
    fontface = "bold"
  ) +
  geom_text(
    aes(y = deaths / 2, label = scales::label_comma()(deaths)),
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(
    name = "Deaths",
    sec.axis = sec_axis(
      ~ . / rate_multiplier,
      name = "ASMR per 100,000"
    ),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_x_continuous(
    breaks = 2010:2024,
    limits = c(2009.5, 2024.5),
    expand = c(0, 0)
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  )

heatmap_panel <- ggplot(
  heatmap_central,
  aes(x = year, y = fct_rev(uf_name), fill = rate_scaled)
) +
  geom_tile() +
  geom_text(
    aes(label = scales::label_number(accuracy = 0.1)(rate)),
    colour = "black",
    size = 3
  ) +
  scale_fill_gradientn(colours = rev(hcl.colors(10, "RdYlGn"))) +
  scale_x_continuous(
    breaks = 2010:2024,
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Year", y = "Federative units") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(colour = "black", size = 9),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
  )

side_panel <- ggplot(
  heatmap_side,
  aes(x = deaths, y = fct_rev(as.factor(year)))
) +
  geom_col(fill = "#DE2D26") +
  facet_grid(
    rows = vars(uf_name),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.22))
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Deaths", y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 0)
  )

figure_s6 <- (
  top_panel + plot_spacer()
) / (
  heatmap_panel + side_panel
) +
  plot_layout(
    widths = c(5, 1.4),
    heights = c(1.3, 6)
  )

figure_s6

ggsave(
  filename = "results/descriptive/figures/Fig_S6.tiff",
  plot = figure_s6,
  device = "tiff",
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

figure_caption <- paste(
  "Fig S6. Temporal distribution of the age-standardised mortality rate",
  "(ASMR) for tuberculosis across federative units in Brazil, 2010–2024.",
  "The top panel presents the annual ASMR and number of deaths. The main",
  "panel shows the variation in ASMR across federative units over the study",
  "period, while the right-hand panel displays the total number of deaths",
  "in each federative unit between 2010 and 2024."
)

writeLines(
  figure_caption,
  "results/descriptive/figures/Fig_S6_caption.txt"
)

# Reproducibility information ----
writeLines(
  capture.output(sessionInfo()),
  "results/descriptive/session_info.txt"
)
