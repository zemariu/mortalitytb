# Article: Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024:
#          a Bayesian spatiotemporal modelling study
# Repository: mortalitytb
# Script: Fits and compares Bayesian spatiotemporal models and estimates adjusted 
#         spatial risks and municipal temporal trends.
# Author: José Mário Nunes da Silva, PhD
# Last updated: 2026-08-01

# Packages ----
library(tidyverse)
library(sf)
library(spdep)
library(geobr)
library(INLA)
library(gt)
library(ggspatial)
library(patchwork)

# Output directories ----

# These directories are used by several analytical stages and are therefore
# defined once near the beginning of the script.
RESULTS_DATA_DIR <- "results/spatiotemporal/outputs"
RESULTS_TABLE_DIR <- "results/spatiotemporal/tables"
RESULTS_FIGURE_DIR <- "results/spatiotemporal/figures"

dir.create(RESULTS_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Input data ----
df_tb <- read_csv(
  "data/processed/tb_mortality.csv",
  show_col_types = FALSE
)

df_tb <- df_tb |>
  mutate(
    code_muni = stringr::str_pad(
      stringr::str_sub(as.character(code_muni), 1, 6),
      width = 6,
      side = "left",
      pad = "0"
    ),
    year = as.integer(year),
    deaths = as.numeric(deaths),
    population = as.numeric(population)
  ) |>
  filter(year >= 2010, year <= 2024)

# Distribution of municipal deaths ----

municipal_deaths <- df_tb |>
  group_by(code_muni) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(deaths)) |>
  mutate(rank = row_number())

gg_deaths <- ggplot(
  municipal_deaths,
  aes(x = rank, y = deaths)
) +
  geom_col(fill = "#b2182b", width = 1) +
  scale_x_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Municipalities ranked by number of deaths",
    y = "Number of deaths"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path(RESULTS_FIGURE_DIR, "Fig_S4.tiff"),
  plot = gg_deaths,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig S4. Distribution of tuberculosis deaths across municipalities",
    "in Brazil, 2010–2024."
  ),
  file.path(RESULTS_FIGURE_DIR, "Fig_S4_caption.txt")
)

# Expected deaths and standardised mortality ratios ----

national_rates <- df_tb |>
  group_by(year, sex, age_group) |>
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    national_rate = if_else(
      population > 0,
      deaths / population,
      NA_real_
    )
  )

df_tb_expected <- df_tb |>
  left_join(
    national_rates |>
      select(year, sex, age_group, national_rate),
    by = c("year", "sex", "age_group")
  ) |>
  mutate(
    expected_stratum = population * national_rate
  ) |>
  group_by(code_muni, year) |>
  summarise(
    O = sum(deaths, na.rm = TRUE),
    E = sum(expected_stratum, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    E = if_else(E > 0, E, NA_real_),
    log_expected = log(E),
    SMR = O / E
  )

readr::write_csv(
  df_tb_expected,
  file.path(RESULTS_DATA_DIR, "municipality_year_observed_expected.csv")
)

# Annual zero-death distribution ----

zero_summary <- df_tb_expected |>
  group_by(year) |>
  summarise(
    municipalities = n(),
    proportion_zero = mean(O == 0, na.rm = TRUE),
    mean_expected = mean(E, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  zero_summary,
  file.path(RESULTS_DATA_DIR, "annual_zero_death_summary.csv")
)

zero_plot_data <- zero_summary |>
  transmute(
    year,
    `Zero deaths` = proportion_zero,
    `Non-zero deaths` = 1 - proportion_zero
  ) |>
  pivot_longer(
    cols = c(`Zero deaths`, `Non-zero deaths`),
    names_to = "category",
    values_to = "proportion"
  ) |>
  mutate(
    category = factor(
      category,
      levels = c("Non-zero deaths", "Zero deaths")
    )
  )

gg_zero <- ggplot(
  zero_plot_data,
  aes(x = factor(year), y = proportion, fill = category)
) +
  geom_col(
    colour = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = scales::percent(proportion, accuracy = 0.1)),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold",
    size = 3.4
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c(
      "Zero deaths" = "#d73027",
      "Non-zero deaths" = "#1a9850"
    )
  ) +
  labs(
    x = "Year",
    y = "Proportion of municipalities",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  filename = file.path(RESULTS_FIGURE_DIR, "Fig_S5.tiff"),
  plot = gg_zero,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig S5. Annual proportion of municipalities with zero and non-zero",
    "tuberculosis deaths in Brazil, 2010–2024."
  ),
  file.path(RESULTS_FIGURE_DIR, "Fig_S5_caption.txt")
)

# Municipal spatial structure ----

# Parameters used only to construct the municipal geometries and adjacency.
BOUNDARY_YEAR <- 2022L
STATE_BOUNDARY_YEAR <- 2020L
PROJECTED_CRS <- 5880
K_NEIGHBOURS <- 6L

SPATIAL_OBJECT_DIR <- "objects/spatial"
dir.create(SPATIAL_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

shp_muni <- geobr::read_municipality(
  code_muni = "all",
  year = BOUNDARY_YEAR,
  showProgress = FALSE
) |>
  mutate(
    code_muni = stringr::str_sub(as.character(code_muni), 1, 6)
  ) |>
  arrange(code_muni)

data_codes <- sort(unique(df_tb_expected$code_muni))
geometry_codes <- sort(unique(shp_muni$code_muni))

if (!identical(data_codes, geometry_codes)) {
  missing_geometry <- setdiff(data_codes, geometry_codes)
  missing_data <- setdiff(geometry_codes, data_codes)

  stop(
    "Municipality codes differ between the mortality data and geobr geometry. ",
    "Missing geometry: ", length(missing_geometry),
    "; missing data: ", length(missing_data), "."
  )
}

shp_muni <- shp_muni |>
  filter(code_muni %in% data_codes) |>
  arrange(match(code_muni, data_codes)) |>
  mutate(id_muni = row_number())

shp_muni_projected <- sf::st_transform(shp_muni, PROJECTED_CRS)

nb_queen <- spdep::poly2nb(
  shp_muni_projected,
  queen = TRUE,
  snap = 1,
  row.names = shp_muni_projected$code_muni
)

isolated_municipalities <- which(spdep::card(nb_queen) == 0L)

if (length(isolated_municipalities) > 0L) {
  municipality_points <- sf::st_point_on_surface(shp_muni_projected)
  coordinates <- sf::st_coordinates(municipality_points)

  nb_knn6 <- spdep::knn2nb(
    spdep::knearneigh(coordinates, k = K_NEIGHBOURS),
    row.names = shp_muni_projected$code_muni
  )

  nb_mixed <- nb_queen

  for (i in isolated_municipalities) {
    neighbours_i <- nb_knn6[[i]]
    nb_mixed[[i]] <- neighbours_i

    for (j in neighbours_i) {
      nb_mixed[[j]] <- sort(unique(c(nb_mixed[[j]], i)))
    }
  }
} else {
  nb_mixed <- nb_queen
}

attr(nb_mixed, "region.id") <- shp_muni_projected$code_muni

if (any(spdep::card(nb_mixed) == 0L)) {
  stop("At least one municipality remained without neighbours.")
}

if (spdep::n.comp.nb(nb_mixed)$nc != 1L) {
  stop("The final municipal adjacency graph contains more than one connected component.")
}

adjacency_summary <- tibble(
  municipalities = length(nb_mixed),
  isolated_before_knn = length(isolated_municipalities),
  minimum_neighbours = min(spdep::card(nb_mixed)),
  mean_neighbours = mean(spdep::card(nb_mixed)),
  maximum_neighbours = max(spdep::card(nb_mixed)),
  connected_components = spdep::n.comp.nb(nb_mixed)$nc
)

readr::write_csv(
  adjacency_summary,
  file.path(RESULTS_DATA_DIR, "municipal_adjacency_summary.csv")
)

saveRDS(
  nb_mixed,
  file.path(SPATIAL_OBJECT_DIR, "municipal_neighbours.rds")
)

graph_file <- file.path(SPATIAL_OBJECT_DIR, "municipality_graph.adj")

spdep::nb2INLA(
  file = graph_file,
  nb = nb_mixed
)

graph_muni <- INLA::inla.read.graph(graph_file)

listw_muni <- spdep::nb2listw(
  nb_mixed,
  style = "W",
  zero.policy = TRUE
)

shp_states <- geobr::read_state(
  code_state = "all",
  year = STATE_BOUNDARY_YEAR,
  showProgress = FALSE
) |>
  sf::st_transform(PROJECTED_CRS)

# INLA analysis data ----

municipality_lookup <- shp_muni_projected |>
  st_drop_geometry() |>
  select(
    code_muni,
    id_muni,
    name_muni,
    name_state,
    name_region
  )

dados_inla <- df_tb_expected |>
  left_join(
    municipality_lookup,
    by = "code_muni"
  ) |>
  mutate(
    O = as.integer(O),
    id_muni = as.integer(id_muni),
    id_year = year - min(year) + 1L,
    time_centered = year - mean(year),
    id_muni_slope = id_muni,
    id_st = as.integer(interaction(id_muni, id_year, drop = TRUE))
  ) |>
  arrange(id_muni, id_year) |>
  as.data.frame()

if (anyNA(dados_inla$id_muni)) {
  stop("Municipal spatial identifiers could not be assigned to all observations.")
}

readr::write_csv(
  dados_inla,
  file.path(RESULTS_DATA_DIR, "bayesian_model_input.csv")
)

# Candidate spatiotemporal models ----

# Settings used only in the comparison of distributional assumptions and
# spatiotemporal specifications.
TEMPORAL_MODEL <- "rw1"
REFIT_DISTRIBUTION_MODELS <- FALSE # TRUE - Force model re-fitting (ignore existing .rds files)

MODEL_DIR <- "objects/bayesian_models"
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

formula_intercept <- O ~ 1

formula_spatial <- O ~ 1 +
  f(
    id_muni,
    model = "bym2",
    graph = graph_muni,
    scale.model = TRUE,
    constr = TRUE
  )

formula_spatiotemporal <- O ~ 1 +
  f(
    id_muni,
    model = "bym2",
    graph = graph_muni,
    scale.model = TRUE,
    constr = TRUE
  ) +
  f(
    id_year,
    model = TEMPORAL_MODEL,
    constr = TRUE
  )

formula_spatiotemporal_slope <- O ~ 1 +
  time_centered +
  f(
    id_muni,
    model = "bym2",
    graph = graph_muni,
    scale.model = TRUE,
    constr = TRUE
  ) +
  f(
    id_year,
    model = TEMPORAL_MODEL,
    constr = TRUE
  ) +
  f(
    id_muni_slope,
    time_centered,
    model = "besag",
    graph = graph_muni,
    scale.model = TRUE,
    constr = TRUE
  )

formula_spatiotemporal_slope_interaction <- update(
  formula_spatiotemporal_slope,
  . ~ . +
    f(
      id_st,
      model = "iid"
    )
)

candidate_formulas <- list(
  `0_intercept` = formula_intercept,
  `1_spatial` = formula_spatial,
  `2_spatiotemporal` = formula_spatiotemporal,
  `3_spatiotemporal_slope` = formula_spatiotemporal_slope,
  `4_spatiotemporal_slope_interaction` =
    formula_spatiotemporal_slope_interaction
)

model_specification_labels <- c(
  `0_intercept` = "Intercept only",
  `1_spatial` = "Spatial BYM2",
  `2_spatiotemporal` = paste0("Spatial BYM2 + temporal ", toupper(TEMPORAL_MODEL)),
  `3_spatiotemporal_slope` = paste0(
    "Spatial BYM2 + temporal ",
    toupper(TEMPORAL_MODEL),
    " + global and municipality-specific linear trends"
  ),
  `4_spatiotemporal_slope_interaction` = paste0(
    "Spatial BYM2 + temporal ",
    toupper(TEMPORAL_MODEL),
    " + global and municipality-specific linear trends + iid space–time interaction"
  )
)

family_labels <- c(
  poisson = "Poisson",
  nbinomial = "Negative binomial",
  zeroinflatedpoisson1 = "Zero-inflated Poisson",
  zeroinflatednbinomial1 = "Zero-inflated negative binomial"
)

fit_model_set <- function(family_code) {
  purrr::imap(
    candidate_formulas,
    function(formula_i, specification_i) {
      message("Fitting ", family_labels[[family_code]], ": ", specification_i)

      tryCatch(
        INLA::inla(
          formula = formula_i,
          family = family_code,
          data = dados_inla,
          offset = log_expected,
          control.predictor = list(
            compute = TRUE,
            link = 1
          ),
          control.compute = list(
            dic = TRUE,
            waic = TRUE,
            cpo = TRUE,
            config = TRUE,
            return.marginals.predictor = TRUE
          ),
          control.inla = list(
            strategy = "simplified.laplace"
          )
        ),
        error = function(e) {
          structure(
            list(message = conditionMessage(e)),
            class = "inla_fit_error"
          )
        }
      )
    }
  )
}

distribution_models_file <- file.path(
  MODEL_DIR,
  "distribution_and_specification_models.rds"
)

if (
  REFIT_DISTRIBUTION_MODELS ||
  !file.exists(distribution_models_file)
) {
  model_sets <- purrr::map(
    names(family_labels),
    fit_model_set
  )

  names(model_sets) <- names(family_labels)

  saveRDS(
    model_sets,
    distribution_models_file
  )
} else {
  model_sets <- readRDS(distribution_models_file)
}

# Table S5: model comparison ----

# First, extract the fit criteria for every successfully fitted candidate model.
# No family or model specification is selected before this comparison.
model_comparison <- purrr::imap_dfr(
  model_sets,
  function(models_i, family_i) {
    purrr::imap_dfr(
      models_i,
      function(model_i, specification_i) {
        failed <- inherits(model_i, "inla_fit_error")

        tibble(
          family_code = family_i,
          Distribution = family_labels[[family_i]],
          specification_code = specification_i,
          `Model specification` =
            model_specification_labels[[specification_i]],
          DIC = if (failed || is.null(model_i$dic$dic)) {
            NA_real_
          } else {
            model_i$dic$dic
          },
          WAIC = if (failed || is.null(model_i$waic$waic)) {
            NA_real_
          } else {
            model_i$waic$waic
          },
          Status = if (failed) {
            paste0("Failed: ", model_i$message)
          } else {
            "Completed"
          }
        )
      }
    )
  }
) |>
  mutate(
    `Delta DIC` = DIC - min(DIC, na.rm = TRUE),
    `Delta WAIC` = WAIC - min(WAIC, na.rm = TRUE)
  )

# Retain only models with valid DIC and WAIC values.
eligible_models <- model_comparison |>
  filter(
    Status == "Completed",
    is.finite(DIC),
    is.finite(WAIC)
  )


# Select the model only after all candidates have been compared.
# WAIC and DIC is used as a tie-breaker.
selected_model <- eligible_models |>
  arrange(WAIC, DIC) |>
  slice(1)

SELECTED_FAMILY <- selected_model$family_code[[1L]]
SELECTED_SPECIFICATION <- selected_model$specification_code[[1L]]

lowest_dic_model <- eligible_models |>
  arrange(DIC, WAIC) |>
  slice(1)

selection_criteria_agree <-
  SELECTED_FAMILY == lowest_dic_model$family_code[[1L]] &&
  SELECTED_SPECIFICATION == lowest_dic_model$specification_code[[1L]]

model_comparison <- model_comparison |>
  mutate(
    Selected = if_else(
      family_code == SELECTED_FAMILY &
        specification_code == SELECTED_SPECIFICATION,
      "Yes",
      "No"
    )
  ) |>
  arrange(WAIC, DIC)

m_best <- model_sets[[SELECTED_FAMILY]][[SELECTED_SPECIFICATION]]

selection_note <- if (selection_criteria_agree) {
  paste0(
    "The model was selected after fitting all candidates. The lowest WAIC and ",
    "lowest DIC identified the same model: **",
    family_labels[[SELECTED_FAMILY]],
    "** with the **",
    model_specification_labels[[SELECTED_SPECIFICATION]],
    "** specification."
  )
} else {
  paste0(
    "The model was selected after fitting all candidates using the lowest WAIC ",
    "as the primary criterion and DIC as the tie-breaker. The selected model was **",
    family_labels[[SELECTED_FAMILY]],
    "** with the **",
    model_specification_labels[[SELECTED_SPECIFICATION]],
    "** specification. The model with the lowest DIC differed from the model ",
    "with the lowest WAIC."
  )
}

message(
  "Selected unadjusted model: ",
  family_labels[[SELECTED_FAMILY]],
  " | ",
  model_specification_labels[[SELECTED_SPECIFICATION]]
)

writexl::write_xlsx(
  model_comparison |>
    select(
      Distribution,
      `Model specification`,
      DIC,
      `Delta DIC`,
      WAIC,
      `Delta WAIC`,
      Selected,
      Status
    ),
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S5_model_comparison.xlsx"
  )
)

table_s5 <- model_comparison |>
  select(
    Distribution,
    `Model specification`,
    DIC,
    `Delta DIC`,
    WAIC,
    `Delta WAIC`,
    Selected
  ) |>
  gt(groupname_col = "Distribution") |>
  fmt_number(
    columns = c(DIC, `Delta DIC`, WAIC, `Delta WAIC`),
    decimals = 2
  ) |>
  cols_label(
    Selected = "Selected model"
  ) |>
  tab_header(
    title = md(
      "**Table S5. Comparison of models for tuberculosis mortality in Brazil under different distributional assumptions.**"
    )
  ) |>
  tab_source_note(
    md(
      paste0(
        selection_note,
        " Lower DIC and WAIC values indicate better relative fit."
      )
    )
  )

gt::gtsave(
  table_s5,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S5_model_comparison.docx"
  )
)

# Unadjusted spatial effects and differential trends ----

# Posterior-probability thresholds used to classify spatial effects and
# municipality-specific temporal trends.
HOTSPOT_THRESHOLD <- 0.95
COLDSPOT_THRESHOLD <- 0.05
INCREASING_THRESHOLD <- 0.80
DECREASING_THRESHOLD <- 0.20

extract_spatiotemporal_classification <- function(
    model,
    model_data,
    adjustment_label
) {
  municipal_lookup_i <- model_data |>
    distinct(
      code_muni,
      id_muni,
      name_muni,
      name_state,
      name_region
    ) |>
    arrange(id_muni)

  n_municipalities <- nrow(municipal_lookup_i)

  spatial_marginals <- model$marginals.random$id_muni[
    seq_len(n_municipalities)
  ]

  if (length(spatial_marginals) != n_municipalities) {
    stop("The number of spatial marginals differs from the number of municipalities.")
  }

  spatial_results <- municipal_lookup_i |>
    mutate(
      RR_mean = vapply(
        spatial_marginals,
        function(marginal_i) {
          INLA::inla.emarginal(
            function(x) exp(x),
            marginal = marginal_i
          )
        },
        numeric(1)
      ),
      RR_lower_95 = vapply(
        spatial_marginals,
        function(marginal_i) {
          exp(INLA::inla.qmarginal(0.025, marginal = marginal_i))
        },
        numeric(1)
      ),
      RR_upper_95 = vapply(
        spatial_marginals,
        function(marginal_i) {
          exp(INLA::inla.qmarginal(0.975, marginal = marginal_i))
        },
        numeric(1)
      ),
      Pr_RR_gt1 = vapply(
        spatial_marginals,
        function(marginal_i) {
          1 - INLA::inla.pmarginal(0, marginal = marginal_i)
        },
        numeric(1)
      ),
      cluster = case_when(
        Pr_RR_gt1 >= HOTSPOT_THRESHOLD ~ "Hotspot",
        Pr_RR_gt1 <= COLDSPOT_THRESHOLD ~ "Coldspot",
        TRUE ~ "Neutralspot"
      )
    )

  slope_summary <- model$summary.random$id_muni_slope |>
    as_tibble() |>
    transmute(
      id_muni = as.integer(ID),
      slope_mean = mean,
      slope_sd = sd,
      slope_lower_95 = `0.025quant`,
      slope_upper_95 = `0.975quant`
    )

  slope_probability <- vapply(
    model$marginals.random$id_muni_slope,
    function(marginal_i) {
      1 - INLA::inla.pmarginal(0, marginal = marginal_i)
    },
    numeric(1)
  )

  slope_summary <- slope_summary |>
    mutate(
      Pr_slope_gt0 = slope_probability,
      trend = case_when(
        Pr_slope_gt0 >= INCREASING_THRESHOLD ~ "Increasing",
        Pr_slope_gt0 <= DECREASING_THRESHOLD ~ "Decreasing",
        TRUE ~ "Stable"
      )
    )

  spatial_results |>
    left_join(
      slope_summary,
      by = "id_muni"
    ) |>
    mutate(
      name_region = recode(
        name_region,
        "Norte" = "North",
        "Nordeste" = "Northeast",
        "Sudeste" = "Southeast",
        "Sul" = "South",
        "Centro-Oeste" = "Central-West",
        "Centro Oeste" = "Central-West",
        .default = name_region
      ),
      Adjustment = adjustment_label,
      profile = paste(cluster, trend, sep = " / ")
    )
}

classification_before <- extract_spatiotemporal_classification(
  model = m_best,
  model_data = dados_inla,
  adjustment_label = "Before"
)

readr::write_csv(
  classification_before,
  file.path(
    RESULTS_DATA_DIR,
    "spatial_temporal_classification_before_adjustment.csv"
  )
)

# Municipal covariates ----
indicators <- readr::read_delim("data/municipal_indicators.csv", 
                         delim = ";", 
                         escape_double = FALSE, 
                         trim_ws = TRUE)

indicators <- indicators |>
  mutate(
    code_muni = stringr::str_pad(
      stringr::str_sub(as.character(code_muni), 1, 6),
      width = 6,
      side = "left",
      pad = "0"
    ),
    year = as.integer(year)
  ) |>
  filter(year >= 2010, year <= 2024)

candidate_covariates <- c(
  "razao_MF",
  "prop_65mais",
  "prop_preta_parda",
  "densidade_d",
  "urbanizacao",
  "taxa_densidade_2_mais",
  "mun_prisao",
  "indice_gini",
  "ivs",
  "idhm",
  "renda_media",
  "prop_vulner_pobreza",
  "taxa_desocupacao",
  "taxa_analfabetismo_18_mais",
  "pib_cte_pc",
  "num_familias_bf",
  "tx_med",
  "tx_leito_sus",
  "desp_tot_saude_pc_mun",
  "cob_ab",
  "cob_vac_bcg",
  "Tx_HIV"
)

missing_covariates <- setdiff(candidate_covariates, names(indicators))

if (length(missing_covariates) > 0L) {
  stop(
    "Missing covariates in the indicator data: ",
    paste(missing_covariates, collapse = ", ")
  )
}

covariate_labels <- c(
  tx_leito_sus = "SUS hospital beds (per 100,000 inhabitants)",
  tx_med = "Physicians (per 1,000 inhabitants)",
  cob_vac_bcg = "BCG vaccination coverage (%)",
  cob_ab = "Primary healthcare coverage (%)",
  pib_cte_pc = "Gross domestic product per capita",
  num_familias_bf = "Families receiving Bolsa Família",
  desp_tot_saude_pc_mun = "Municipal health expenditure per capita",
  razao_MF = "Sex ratio (M/F)",
  prop_preta_parda = "Black and Brown population (%)",
  prop_65mais = "Population aged ≥65 years (%)",
  Tx_HIV = "AIDS incidence rate (per 100,000 inhabitants)",
  taxa_analfabetismo_18_mais = "Illiteracy rate among people aged ≥18 years (%)",
  indice_gini = "Gini index",
  prop_vulner_pobreza = "Population vulnerable to poverty (%)",
  taxa_desocupacao = "Unemployment rate (%)",
  taxa_densidade_2_mais = "Household crowding (%)",
  urbanizacao = "Urbanisation level (%)",
  densidade_d = "Population density (people/km²)",
  idhm = "Municipal Human Development Index",
  ivs = "Social Vulnerability Index",
  renda_media = "Average household income per capita",
  mun_prisao = "Prison in the municipality"
)


# Table S6: covariate correlations ----

correlation_data <- indicators |>
  select(all_of(candidate_covariates))

correlation_matrix <- cor(
  correlation_data,
  use = "pairwise.complete.obs",
  method = "spearman"
)

correlation_p_matrix <- matrix(
  NA_real_,
  nrow = length(candidate_covariates),
  ncol = length(candidate_covariates),
  dimnames = list(candidate_covariates, candidate_covariates)
)

format_p <- function(p) {
  ifelse(is.na(p),
    NA,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)
    )
  )
}

diag(correlation_p_matrix) <- 0

for (i in seq_len(length(candidate_covariates) - 1L)) {
  for (j in seq.int(i + 1L, length(candidate_covariates))) {
    complete_ij <- complete.cases(
      correlation_data[[i]],
      correlation_data[[j]]
    )

    correlation_test <- suppressWarnings(
      cor.test(
        correlation_data[[i]][complete_ij],
        correlation_data[[j]][complete_ij],
        method = "spearman",
        exact = FALSE
      )
    )

    correlation_p_matrix[i, j] <- correlation_test$p.value
    correlation_p_matrix[j, i] <- correlation_test$p.value
  }
}

correlation_matrix_fmt <- apply(
  correlation_matrix,
  c(1, 2),
  function(x) sprintf("%.2f", x)
)

correlation_p_matrix_fmt <- apply(
  correlation_p_matrix,
  c(1, 2),
  format_p
)

correlation_matrix_output <- as.data.frame(correlation_matrix_fmt) |>
  rownames_to_column("Variable") |>
  mutate(
    Variable = unname(covariate_labels[Variable])
  )

names(correlation_matrix_output)[-1] <- unname(
  covariate_labels[names(correlation_matrix_output)[-1]]
)

correlation_p_output <- as.data.frame(correlation_p_matrix_fmt) |>
  rownames_to_column("Variable") |>
  mutate(
    Variable = unname(covariate_labels[Variable])
  )

names(correlation_p_output)[-1] <- unname(
  covariate_labels[names(correlation_p_output)[-1]]
)

writexl::write_xlsx(
  correlation_matrix_output,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S6_spearman_correlations.xlsx"
  )
)

writexl::write_xlsx(
  correlation_p_output,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S6_spearman_p_values.xlsx"
  )
)

correlation_combined <- matrix(
  paste0(
    correlation_matrix_fmt,
    " (",
    correlation_p_matrix_fmt,
    ")"
  ),
  nrow = nrow(correlation_matrix_fmt),
  dimnames = dimnames(correlation_matrix_fmt)
)

correlation_output <- as.data.frame(correlation_combined) |>
  rownames_to_column("Variable") |>
  mutate(
    Variable = unname(covariate_labels[Variable])
  )

names(correlation_output)[-1] <- unname(
  covariate_labels[names(correlation_output)[-1]]
)

# Combine correlation coefficients (ρ) and p-values
writexl::write_xlsx(
  correlation_output,
  path = file.path(
    RESULTS_TABLE_DIR,
    "Table_S6_spearman_combined.xlsx"
  )
)

# Table S7: variance inflation factors and tolerance ----

VIF_THRESHOLD <- 5
VIF_TOLERANCE_THRESHOLD <- 0.20

calculate_vif <- function(data, variables, stage) {
  analysis_data <- data |>
    select(all_of(variables)) |>
    mutate(across(everything(), as.numeric)) |>
    filter(
      if_all(
        everything(),
        ~ !is.na(.x) & is.finite(.x)
      )
    )

  purrr::map_dfr(
    variables,
    function(variable_i) {
      other_variables <- setdiff(variables, variable_i)

      auxiliary_formula <- reformulate(
        other_variables,
        response = variable_i
      )

      auxiliary_model <- lm(
        auxiliary_formula,
        data = analysis_data
      )

      r_squared <- summary(auxiliary_model)$r.squared
      vif_value <- if (is.na(r_squared) || r_squared >= 1) {
        Inf
      } else {
        1 / (1 - r_squared)
      }

      tibble(
        Stage = stage,
        variable_code = variable_i,
        Variable = covariate_labels[[variable_i]],
        VIF = round(vif_value, 3),
        Tolerance = round(1 / vif_value, 3)
      )
    }
  )
}

vif_before <- calculate_vif(
  data = indicators,
  variables = candidate_covariates,
  stage = "Before exclusion"
) |>
  mutate(
    Decision = if_else(
      variable_code %in% collinear_exclusions,
      "Excluded",
      "Retained"
    )
  )

collinear_exclusions <- c(
  "taxa_analfabetismo_18_mais",
  "prop_vulner_pobreza",
  "renda_media",
  "idhm"
)

vif_after <- calculate_vif(
  data = indicators,
  variables = full_covariates,
  stage = "After exclusion"
) |>
  mutate(
    Decision = "Retained"
  )

post_collinearity_covariates <- setdiff(
  candidate_covariates,
  collinear_exclusions
)

table_s7_data <- bind_rows(
  vif_before,
  vif_after
) |>
  mutate(
    Stage = factor(
      Stage,
      levels = c("Before exclusion", "After exclusion")
    ),
    Flag = case_when(
      VIF >= VIF_THRESHOLD |
        Tolerance <= VIF_TOLERANCE_THRESHOLD ~ "Collinearity concern",
      TRUE ~ "Acceptable"
    )
  ) |>
  arrange(Stage, desc(VIF))

writexl::write_xlsx(
  table_s7_data |>
    select(
      Stage,
      Variable,
      VIF,
      Tolerance,
      Decision,
      Flag
    ),
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S7_vif_and_tolerance.xlsx"
  )
)

table_s7 <- table_s7_data |>
  select(
    Stage,
    Variable,
    VIF,
    Tolerance,
    Decision
  ) |>
  gt(groupname_col = "Stage") |>
  fmt_number(
    columns = c(VIF, Tolerance),
    decimals = 3
  ) |>
  tab_header(
    title = md(
      "**Table S7. Variance inflation factors and tolerance before and after exclusion of collinear variables.**"
    )
  ) |>
  tab_source_note(
    md(
      paste0(
        "Covariates with VIF ≥ ",
        VIF_THRESHOLD,
        " or tolerance ≤ ",
        format(VIF_TOLERANCE_THRESHOLD, nsmall = 2),
        " were considered to show potentially important collinearity."
      )
    )
  )

gt::gtsave(
  table_s7,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S7_vif_and_tolerance.docx"
  )
)

# Descriptive summary of municipal covariates ----

continuous_covariates <- setdiff(
  candidate_covariates,
  "mun_prisao"
)

continuous_summary <- purrr::map_dfr(
  continuous_covariates,
  function(variable_i) {
    values_i <- as.numeric(indicators[[variable_i]])

    tibble(
      Variable = covariate_labels[[variable_i]],
      Category = "",
      Summary = sprintf(
        "%.2f (%.2f)",
        mean(values_i, na.rm = TRUE),
        sd(values_i, na.rm = TRUE)
      )
    )
  }
)

total_mun <- indicators |>
  distinct(code_muni) |>
  nrow()

prison_summary <- indicators |>
  mutate(mun_prisao = as.integer(mun_prisao)) |>
  distinct(code_muni, mun_prisao) |>
  filter(!is.na(mun_prisao), mun_prisao == 1) |>
  count(mun_prisao, name = "n") |>
  mutate(
    Variable = covariate_labels[["mun_prisao"]],
    Category = "Yes",
    Summary = sprintf(
      "%s (%.1f%%)",
      scales::comma(n),
      100 * n / total_mun
    )
  ) |>
  select(Variable, Category, Summary)

descriptive_covariate_table <- bind_rows(
  continuous_summary,
  prison_summary
)

writexl::write_xlsx(
  descriptive_covariate_table,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_descriptive_covariates.xlsx"
  )
)

table_2 <- descriptive_covariate_table |>
  gt() |>
  cols_label(
    Variable = "Municipal characteristic",
    Category = "Category",
    Summary = "Mean (SD) or n (%)"
  ) |>
  tab_header(
    title = md(
      "**Table 2. Descriptive characteristics of municipal covariates included in the Bayesian spatiotemporal analysis, Brazil, 2010–2024.**"
    )
  )

gt::gtsave(
  table_2,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_descriptive_covariates.docx"
  )
)

# Standardised covariates and analytical panel ----

continuous_full_covariates <- setdiff(
  full_covariates,
  "mun_prisao"
)

indicators_standardised <- indicators |>
  mutate(
    across(
      all_of(continuous_full_covariates),
      ~ as.numeric(scale(as.numeric(.x)))
    ),
    mun_prisao = as.numeric(mun_prisao)
  ) |>
  select(
    code_muni,
    year,
    all_of(full_covariates)
  )

dados_inla_cov <- dados_inla |>
  left_join(
    indicators_standardised,
    by = c("code_muni", "year")
  )

covariate_missingness <- dados_inla_cov |>
  summarise(
    across(
      all_of(full_covariates),
      ~ sum(is.na(.x))
    )
  ) |>
  pivot_longer(
    everything(),
    names_to = "variable_code",
    values_to = "missing"
  ) |>
  mutate(
    Variable = unname(covariate_labels[variable_code])
  ) |>
  select(
    variable_code,
    Variable,
    missing
  )

readr::write_csv(
  covariate_missingness,
  file.path(
    RESULTS_DATA_DIR,
    "covariate_missingness.csv"
  )
)

dados_inla_cov <- dados_inla_cov |>
  filter(
    if_all(
      all_of(full_covariates),
      ~ !is.na(.x) & is.finite(.x)
    )
  ) |>
  arrange(id_muni, id_year) |>
  as.data.frame()

readr::write_csv(
  dados_inla_cov,
  file.path(
    RESULTS_DATA_DIR,
    "bayesian_model_input_with_covariates.csv"
  )
)

# Univariable covariate models ----

# Keep FALSE to reuse saved models. The models are fitted automatically when
# the saved file does not yet exist.
REFIT_UNIVARIABLE_MODELS <- FALSE

random_effect_terms <- paste0(
  "+ f(id_muni, model = 'bym2', graph = graph_muni, ",
  "scale.model = TRUE, constr = TRUE) ",
  "+ f(id_year, model = '", TEMPORAL_MODEL, "', constr = TRUE) ",
  "+ f(id_muni_slope, time_centered, model = 'besag', ",
  "graph = graph_muni, scale.model = TRUE, constr = TRUE) ",
  "+ f(id_st, model = 'iid')"
)

build_covariate_formula <- function(covariates) {
  fixed_terms <- c(
    "time_centered",
    covariates
  )
  
  as.formula(
    paste(
      "O ~ 1 +",
      paste(fixed_terms, collapse = " + "),
      random_effect_terms
    )
  )
}

fit_inla_model <- function(formula_i) {
  INLA::inla(
    formula = formula_i,
    family = SELECTED_FAMILY,
    data = dados_inla_cov,
    offset = log_expected,
    control.predictor = list(
      compute = TRUE,
      link = 1
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      config = TRUE,
      return.marginals.predictor = TRUE
    ),
    control.inla = list(
      strategy = "simplified.laplace"
    )
  )
}

null_formula <- build_covariate_formula(
  covariates = character(0)
)

univariable_model_file <- file.path(
  MODEL_DIR,
  "univariable_covariate_models.rds"
)

if (
  REFIT_UNIVARIABLE_MODELS ||
  !file.exists(univariable_model_file)
) {
  message("Fitting the covariate-free reference model...")
  
  null_model <- fit_inla_model(
    null_formula
  )
  
  univariable_models <- list()
  
  univariable_results <- purrr::map_dfr(
    post_collinearity_covariates,
    function(variable_i) {
      formula_i <- build_covariate_formula(
        covariates = variable_i
      )
      
      message("Fitting univariable model: ", variable_i)
      
      model_i <- tryCatch(
        fit_inla_model(formula_i),
        error = function(e) {
          structure(
            list(message = conditionMessage(e)),
            class = "inla_fit_error"
          )
        }
      )
      
      univariable_models[[variable_i]] <<- model_i
      
      if (inherits(model_i, "inla_fit_error")) {
        return(
          tibble(
            variable_code = variable_i,
            Variable = covariate_labels[[variable_i]],
            RR = NA_real_,
            CI_lower_95 = NA_real_,
            CI_upper_95 = NA_real_,
            DIC = NA_real_,
            Delta_DIC = NA_real_,
            WAIC = NA_real_,
            Delta_WAIC = NA_real_,
            Posterior_association = NA,
            Improves_DIC = NA,
            Improves_WAIC = NA,
            Status = paste0("Failed: ", model_i$message)
          )
        )
      }
      
      fixed_row <- model_i$summary.fixed[
        variable_i,
        ,
        drop = FALSE
      ]
      
      rr_i <- exp(fixed_row$mean)
      lower_i <- exp(fixed_row$`0.025quant`)
      upper_i <- exp(fixed_row$`0.975quant`)
      dic_i <- model_i$dic$dic
      waic_i <- model_i$waic$waic
      
      tibble(
        variable_code = variable_i,
        Variable = covariate_labels[[variable_i]],
        RR = rr_i,
        CI_lower_95 = lower_i,
        CI_upper_95 = upper_i,
        DIC = dic_i,
        Delta_DIC = dic_i - null_model$dic$dic,
        WAIC = waic_i,
        Delta_WAIC = waic_i - null_model$waic$waic,
        Posterior_association =
          lower_i > 1 | upper_i < 1,
        Improves_DIC = dic_i < null_model$dic$dic,
        Improves_WAIC = waic_i < null_model$waic$waic,
        Status = "Completed"
      )
    }
  ) |>
    arrange(WAIC, DIC)
  
  saveRDS(
    list(
      null_model = null_model,
      null_formula = null_formula,
      univariable_models = univariable_models,
      univariable_results = univariable_results
    ),
    univariable_model_file
  )
} else {
  univariable_objects <- readRDS(
    univariable_model_file
  )
  
  null_model <- univariable_objects$null_model
  null_formula <- univariable_objects$null_formula
  univariable_models <- univariable_objects$univariable_models
  univariable_results <- univariable_objects$univariable_results
}

writexl::write_xlsx(
  univariable_results,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_univariable_covariate_models.xlsx"
  )
)

univariable_table <- univariable_results |>
  select(
    Variable,
    RR,
    CI_lower_95,
    CI_upper_95,
    DIC,
    Delta_DIC,
    WAIC,
    Delta_WAIC,
    Posterior_association,
    Improves_DIC,
    Improves_WAIC,
    Status
  ) |>
  gt() |>
  fmt_number(
    columns = c(
      RR,
      CI_lower_95,
      CI_upper_95,
      DIC,
      Delta_DIC,
      WAIC,
      Delta_WAIC
    ),
    decimals = 2
  ) |>
  cols_label(
    CI_lower_95 = "Lower 95% CrI",
    CI_upper_95 = "Upper 95% CrI",
    Delta_DIC = "ΔDIC",
    Delta_WAIC = "ΔWAIC",
    Posterior_association = "95% CrI excludes 1",
    Improves_DIC = "DIC lower than null",
    Improves_WAIC = "WAIC lower than null"
  ) |>
  tab_header(
    title = md(
      "**Univariable Bayesian spatiotemporal models for municipal tuberculosis mortality, Brazil, 2010–2024.**"
    )
  ) |>
  tab_source_note(
    md(
      "ΔDIC and ΔWAIC are calculated relative to the covariate-free spatiotemporal model. Negative values indicate improvement."
    )
  )

gt::gtsave(
  univariable_table,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_univariable_covariate_models.docx"
  )
)

# Full multivariable model ----

# The full model is defined only after the univariable results have been
# calculated and reviewed. All non-collinear covariates with successfully
# fitted univariable models are included at this stage. 

full_covariates <- univariable_results |>
  filter(Status == "Completed") |>
  pull(variable_code)

if (length(full_covariates) == 0L) {
  stop("No covariates were available for the full multivariable model.")
}

full_formula <- build_covariate_formula(
  covariates = full_covariates
)

# Keep FALSE to reuse the saved full model. The model is fitted automatically
# when the saved file does not yet exist.
REFIT_FULL_MODEL <- FALSE

full_model_file <- file.path(
  MODEL_DIR,
  "full_covariate_model.rds"
)

if (
  REFIT_FULL_MODEL ||
  !file.exists(full_model_file)
) {
  message("Fitting the full multivariable model...")
  
  full_model <- fit_inla_model(
    full_formula
  )
  
  saveRDS(
    list(
      full_model = full_model,
      full_formula = full_formula,
      full_covariates = full_covariates
    ),
    full_model_file
  )
} else {
  full_model_objects <- readRDS(
    full_model_file
  )
  
  full_model <- full_model_objects$full_model
  full_formula <- full_model_objects$full_formula
  full_covariates <- full_model_objects$full_covariates
}

# Full-model results for manual review ----

extract_fixed_effects <- function(model, model_label) {
  model$summary.fixed |>
    as.data.frame() |>
    rownames_to_column("variable_code") |>
    as_tibble() |>
    transmute(
      Model = model_label,
      variable_code,
      Variable = case_when(
        variable_code == "(Intercept)" ~ "Intercept",
        variable_code == "time_centered" ~
          "Global linear temporal trend",
        variable_code %in% names(covariate_labels) ~
          unname(covariate_labels[variable_code]),
        TRUE ~ variable_code
      ),
      coefficient_mean = mean,
      coefficient_lower_95 = `0.025quant`,
      coefficient_upper_95 = `0.975quant`,
      RR = exp(mean),
      RR_lower_95 = exp(`0.025quant`),
      RR_upper_95 = exp(`0.975quant`)
    )
}

fixed_effects_full <- extract_fixed_effects(
  full_model,
  "Full model"
)

writexl::write_xlsx(
  fixed_effects_full,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_full_model_fixed_effects.xlsx"
  )
)

full_model_hyperparameters <- full_model$summary.hyperpar |>
  as.data.frame() |>
  rownames_to_column("Parameter") |>
  as_tibble() |>
  mutate(Model = "Full model", .before = 1)

writexl::write_xlsx(
  full_model_hyperparameters,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_full_model_hyperparameters.xlsx"
  )
)

full_model_fit <- tibble(
  Model = "Full model",
  Number_of_covariates = length(full_covariates),
  DIC = full_model$dic$dic,
  WAIC = full_model$waic$waic,
  effective_parameters = full_model$dic$p.eff
)

writexl::write_xlsx(
  full_model_fit,
  file.path(
    RESULTS_DATA_DIR,
    "full_model_fit_criteria.xlsx"
  )
)

readr::write_lines(
  c(
    "Covariates included in the full model:",
    paste(full_covariates, collapse = ", ")
  ),
  file.path(
    RESULTS_TABLE_DIR,
    "full_model_covariates_for_manual_review.txt"
  )
)

# Manual definition of the final covariate set ----

# IMPORTANT MANUAL STEP
# ---------------------
# Review the following outputs before completing `final_covariates`:
#
#   1. univariable_covariate_models.csv;
#   2. full_model_fixed_effects.csv;
#   3. full_model_fit_criteria.csv; and
#   4. the conceptual relevance of each covariate.
#
# Add below only the variable codes that should be tested in the candidate
# final model. The list must be a subset of `full_covariates` and must use the
# exact variable names from the analytical dataset.
#
# This is an intentionally manual and iterative step. To evaluate another
# candidate final model:
#
#   1. edit only `final_covariates`;
#   2. keep `REFIT_UNIVARIABLE_MODELS <- FALSE`;
#   3. keep `REFIT_FULL_MODEL <- FALSE`;
#   4. set `REFIT_FINAL_MODEL <- TRUE`; and
#   5. rerun this section and the subsequent sections.
#
# Compare the candidate model with the full model using posterior estimates,
# DIC, WAIC, effective parameters, residual diagnostics and conceptual
# plausibility. Repeat this manual process until an adequate final model is
# identified. 

final_covariates <- c(
  "urbanizacao",
  "taxa_densidade_2_mais",
  "taxa_desocupacao",
  "indice_gini",
  "desp_tot_saude_pc_mun",
  "ivs",
  "Tx_HIV",
  "mun_prisao"
)

final_formula <- build_covariate_formula(
  covariates = final_covariates
)

# Set TRUE whenever `final_covariates` is edited. After the selected candidate
# has been fitted and saved successfully, set FALSE to reuse that model.
REFIT_FINAL_MODEL <- TRUE

final_model_file <- file.path(
  MODEL_DIR,
  "final_covariate_model.rds"
)

if (
  REFIT_FINAL_MODEL ||
  !file.exists(final_model_file)
) {
  message(
    "Fitting the manually specified candidate final model with: ",
    paste(final_covariates, collapse = ", ")
  )
  
  final_model <- fit_inla_model(
    final_formula
  )
  
  saveRDS(
    list(
      final_model = final_model,
      final_formula = final_formula,
      final_covariates = final_covariates
    ),
    final_model_file
  )
} else {
  final_model_objects <- readRDS(
    final_model_file
  )
  
  saved_final_covariates <- final_model_objects$final_covariates
  
  if (!identical(final_covariates, saved_final_covariates)) {
    stop(
      paste0(
        "`final_covariates` differs from the covariate set stored with the saved ",
        "final model. Set `REFIT_FINAL_MODEL <- TRUE` to fit the new candidate."
      ),
      call. = FALSE
    )
  }
  
  final_model <- final_model_objects$final_model
  final_formula <- final_model_objects$final_formula
}

# Full and candidate final model results ----

fixed_effects_final <- extract_fixed_effects(
  final_model,
  "Candidate final model"
)

fixed_effects_combined <- bind_rows(
  fixed_effects_full,
  fixed_effects_final
)

writexl::write_xlsx(
  fixed_effects_combined,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_full_and_candidate_final_fixed_effects.xlsx"
  )
)

hyperparameters_combined <- bind_rows(
  full_model_hyperparameters,
  final_model$summary.hyperpar |>
    as.data.frame() |>
    rownames_to_column("Parameter") |>
    as_tibble() |>
    mutate(Model = "Candidate final model", .before = 1)
)

writexl::write_xlsx(
  hyperparameters_combined,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_full_and_candidate_final_hyperparameters.xlsx"
  )
)

adjusted_model_comparison <- tibble(
  Model = c("Full model", "Candidate final model"),
  Number_of_covariates = c(
    length(full_covariates),
    length(final_covariates)
  ),
  DIC = c(
    full_model$dic$dic,
    final_model$dic$dic
  ),
  WAIC = c(
    full_model$waic$waic,
    final_model$waic$waic
  ),
  effective_parameters = c(
    full_model$dic$p.eff,
    final_model$dic$p.eff
  )
) |>
  mutate(
    Delta_DIC_from_full = DIC - DIC[Model == "Full model"],
    Delta_WAIC_from_full = WAIC - WAIC[Model == "Full model"]
  )

writexl::write_xlsx(
  adjusted_model_comparison,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_full_and_candidate_final_model_comparison.xlsx"
  )
)

model_comparison_table <- adjusted_model_comparison |>
  gt() |>
  fmt_number(
    columns = c(
      DIC,
      WAIC,
      effective_parameters,
      Delta_DIC_from_full,
      Delta_WAIC_from_full
    ),
    decimals = 2
  ) |>
  cols_label(
    Number_of_covariates = "Covariates",
    effective_parameters = "Effective parameters",
    Delta_DIC_from_full = "ΔDIC from full",
    Delta_WAIC_from_full = "ΔWAIC from full"
  ) |>
  tab_header(
    title = md(
      "**Comparison of the full and manually specified candidate final Bayesian spatiotemporal models.**"
    )
  ) |>
  tab_source_note(
    md(
      "The candidate final model is manually specified in `final_covariates`. It should be regarded as final only after posterior estimates, information criteria and diagnostic results have been reviewed."
    )
  )

gt::gtsave(
  model_comparison_table,
  file.path(
    RESULTS_TABLE_DIR,
    "_Table_2_full_and_candidate_final_model_comparison.docx"
  )
)

readr::write_lines(
  c(
    "Full-model covariates:",
    paste(full_covariates, collapse = ", "),
    "",
    "Manually specified candidate final-model covariates:",
    paste(final_covariates, collapse = ", ")
  ),
  file.path(
    RESULTS_DATA_DIR,
    "full_and_candidate_final_covariate_sets.txt"
  )
)

# BYM2 precision components ----

summarise_bym2_precision <- function(model, model_label) {
  set.seed(123)

  hyperparameter_samples <- INLA::inla.hyperpar.sample(
    n = 10000,
    result = model
  )

  precision_name <- grep(
    "^Precision for id_muni$",
    colnames(hyperparameter_samples),
    value = TRUE
  )

  phi_name <- grep(
    "^Phi for id_muni$",
    colnames(hyperparameter_samples),
    value = TRUE
  )

  if (length(precision_name) != 1L || length(phi_name) != 1L) {
    stop(
      "The BYM2 precision or mixing hyperparameter could not be identified in ",
      model_label,
      "."
    )
  }

  tau_b <- hyperparameter_samples[, precision_name]
  phi <- hyperparameter_samples[, phi_name]

  precision_structured <- tau_b / phi
  precision_unstructured <- tau_b / (1 - phi)

  bind_rows(
    tibble(
      Model = model_label,
      Component = "Structured spatial component",
      Mean = mean(precision_structured, na.rm = TRUE),
      Median = median(precision_structured, na.rm = TRUE),
      Lower_95 = quantile(
        precision_structured,
        0.025,
        na.rm = TRUE
      ),
      Upper_95 = quantile(
        precision_structured,
        0.975,
        na.rm = TRUE
      )
    ),
    tibble(
      Model = model_label,
      Component = "Unstructured spatial component",
      Mean = mean(precision_unstructured, na.rm = TRUE),
      Median = median(precision_unstructured, na.rm = TRUE),
      Lower_95 = quantile(
        precision_unstructured,
        0.025,
        na.rm = TRUE
      ),
      Upper_95 = quantile(
        precision_unstructured,
        0.975,
        na.rm = TRUE
      )
    )
  )
}

bym2_precision_summary <- bind_rows(
  summarise_bym2_precision(
    full_model,
    "Full model"
  ),
  summarise_bym2_precision(
    final_model,
    "Final model"
  )
)

writexl::write_xlsx(
  bym2_precision_summary,
  file.path(
    RESULTS_DATA_DIR,
    "bym2_precision_components.xlsx"
  )
)

# Pearson-type dispersion ----

dispersion_summary <- purrr::imap_dfr(
  list(
    `Full model` = full_model,
    `Final model` = final_model
  ),
  function(model_i, model_label) {
    observed_i <- dados_inla_cov$O
    fitted_i <- model_i$summary.fitted.values$mean
    degrees_freedom <- nrow(dados_inla_cov) - model_i$dic$p.eff

    pearson_chisquare <- sum(
      (observed_i - fitted_i)^2 /
        pmax(fitted_i, 1e-12),
      na.rm = TRUE
    )

    dispersion <- pearson_chisquare / degrees_freedom
    alpha <- 0.05

    tibble(
      Model = model_label,
      Pearson_dispersion = dispersion,
      Lower_95 =
        degrees_freedom * dispersion /
        qchisq(1 - alpha / 2, degrees_freedom),
      Upper_95 =
        degrees_freedom * dispersion /
        qchisq(alpha / 2, degrees_freedom),
      degrees_freedom = degrees_freedom
    )
  }
)

writexl::write_xlsx(
  dispersion_summary,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_pearson_dispersion_summary.xlsx"
  )
)

# Goodness-of-fit diagnostics ----

observed <- dados_inla_cov$O
fitted_final <- final_model$summary.fitted.values$mean
fitted_full <- full_model$summary.fitted.values$mean

pearson_residual_final <- (
  observed - fitted_final
) / sqrt(pmax(fitted_final, 1e-12))

pearson_residual_full <- (
  observed - fitted_full
) / sqrt(pmax(fitted_full, 1e-12))

diagnostic_data <- tibble(
  observed,
  fitted_final,
  pearson_residual_final
)

diagnostic_plot_a <- ggplot(
  diagnostic_data,
  aes(x = fitted_final, y = observed)
) +
  geom_point(alpha = 0.30, size = 0.8) +
  geom_abline(
    intercept = 0,
    slope = 1,
    colour = "#b2182b",
    linewidth = 0.9
  ) +
  labs(
    title = "A)",
    x = "Fitted values",
    y = "Observed values"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold")
  )

diagnostic_plot_b <- ggplot(
  diagnostic_data,
  aes(x = fitted_final, y = pearson_residual_final)
) +
  geom_point(alpha = 0.30, size = 0.8) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    colour = "#b2182b",
    linetype = "dotted"
  ) +
  labs(
    title = "B)",
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold")
  )

diagnostic_figure <- diagnostic_plot_a + diagnostic_plot_b

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S7.tiff"
  ),
  plot = diagnostic_figure,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig S7. Goodness-of-fit diagnostics for the final Bayesian",
    "spatiotemporal model: (A) observed versus fitted values;",
    "(B) Pearson residuals versus fitted values."
  ),
  file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S7_caption.txt"
  )
)

# Spatial distribution of residuals ----

residual_data <- dados_inla_cov |>
  mutate(
    fitted_final = fitted_final,
    fitted_full = fitted_full,
    pearson_residual_final = pearson_residual_final,
    pearson_residual_full = pearson_residual_full
  )

municipal_residuals <- residual_data |>
  group_by(code_muni) |>
  summarise(
    residual_mean_final = mean(
      pearson_residual_final,
      na.rm = TRUE
    ),
    residual_mean_full = mean(
      pearson_residual_full,
      na.rm = TRUE
    ),
    residual_median_final = median(
      pearson_residual_final,
      na.rm = TRUE
    ),
    residual_sd_final = sd(
      pearson_residual_final,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

residual_map_data <- shp_muni_projected |>
  left_join(
    municipal_residuals,
    by = "code_muni"
  ) |>
  arrange(id_muni)

gg_residual <- ggplot() +
  geom_sf(
    data = residual_map_data,
    aes(fill = residual_mean_final),
    colour = NA
  ) +
  scale_fill_gradient2(
    low = "#2b83ba",
    mid = "#ffffbf",
    high = "#d7191c",
    midpoint = 0,
    name = "Mean Pearson\nresidual"
  ) +
  geom_sf(
    data = shp_states,
    fill = NA,
    colour = "black",
    linewidth = 0.5
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
  coord_sf(datum = NA) +
  theme_void(base_size = 12)

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S8.tiff"
  ),
  plot = gg_residual,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig S8. Spatial distribution of mean Pearson residuals from the final",
    "Bayesian spatiotemporal model, Brazil, 2010–2024."
  ),
  file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S8_caption.txt"
  )
)

# Residual spatial autocorrelation ----

MORAN_SIMULATIONS <- 999L

residual_moran_summary <- purrr::imap_dfr(
  list(
    `Full model` = residual_map_data$residual_mean_full,
    `Final model` = residual_map_data$residual_mean_final
  ),
  function(residuals_i, model_label) {
    set.seed(2026)

    moran_i <- spdep::moran.mc(
      residuals_i,
      listw = listw_muni,
      nsim = MORAN_SIMULATIONS,
      zero.policy = TRUE
    )

    tibble(
      Model = model_label,
      Moran_I = round(as.numeric(moran_i$statistic), 2),
      Z_score = (
        as.numeric(moran_i$statistic) -
          mean(moran_i$res)
      ) / sd(moran_i$res),
      P_value = moran_i$p.value,
      Simulations = MORAN_SIMULATIONS
    )
  }
)

writexl::write_xlsx(
  residual_moran_summary,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_2_residual_moran_diagnostics.xlsx"
  )
)

# Adjusted spatial effects and differential trends ----

classification_after <- extract_spatiotemporal_classification(
  model = final_model,
  model_data = dados_inla_cov,
  adjustment_label = "After"
)

writexl::write_xlsx(
  classification_after,
  file.path(
    RESULTS_DATA_DIR,
    "spatial_temporal_classification_after_adjustment.xlsx"
  )
)

classification_combined <- bind_rows(
  classification_before,
  classification_after
) |>
  mutate(
    Adjustment = factor(
      Adjustment,
      levels = c("Before", "After")
    ),
    cluster = factor(
      cluster,
      levels = c("Hotspot", "Coldspot", "Neutralspot")
    ),
    trend = factor(
      trend,
      levels = c("Increasing", "Decreasing", "Stable")
    ),
    name_region = factor(
      name_region,
      levels = c(
        "North",
        "Northeast",
        "Southeast",
        "South",
        "Central-West"
      )
    )
  )

# Table S12: national cross-classification before and after adjustment ----

table_s12_base <- classification_combined |>
  count(
    Adjustment,
    cluster,
    trend,
    name = "n"
  ) |>
  tidyr::complete(
    Adjustment,
    cluster,
    trend,
    fill = list(n = 0)
  ) |>
  pivot_wider(
    names_from = trend,
    values_from = n
  ) |>
  mutate(
    Total = Increasing + Decreasing + Stable
  ) |>
  arrange(
    Adjustment,
    cluster
  )

table_s12_totals <- table_s12_base |>
  group_by(Adjustment) |>
  summarise(
    cluster = "Total",
    Increasing = sum(Increasing),
    Decreasing = sum(Decreasing),
    Stable = sum(Stable),
    Total = sum(Total),
    .groups = "drop"
  )

table_s12_data <- bind_rows(
  table_s12_base |>
    mutate(cluster = as.character(cluster)),
  table_s12_totals
) |>
  mutate(
    Adjustment = factor(
      Adjustment,
      levels = c("Before", "After")
    ),
    cluster = factor(
      cluster,
      levels = c(
        "Hotspot",
        "Coldspot",
        "Neutralspot",
        "Total"
      )
    )
  ) |>
  arrange(
    Adjustment,
    cluster
  )

writexl::write_xlsx(
  table_s12_data |>
    rename(`Cluster type` = cluster),
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S12_national_cross_classification.xlsx"
  )
)

table_s12 <- table_s12_data |>
  rename(`Cluster type` = cluster) |>
  gt(groupname_col = "Adjustment") |>
  cols_label(
    Increasing = "Increasing",
    Decreasing = "Decreasing",
    Stable = "Stable",
    Total = "Total"
  ) |>
  tab_spanner(
    label = "Temporal trend",
    columns = c(
      Increasing,
      Decreasing,
      Stable
    )
  ) |>
  fmt_integer(
    columns = c(
      Increasing,
      Decreasing,
      Stable,
      Total
    ),
    use_seps = TRUE
  ) |>
  tab_header(
    title = md(
      "**Table S12. Number of municipalities according to the cross-classification of spatial clustering and temporal trends in tuberculosis mortality, before and after covariate adjustment, Brazil, 2010–2024.**"
    )
  )

gt::gtsave(
  table_s12,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S12_national_cross_classification.docx"
  )
)

# Table S13: regional and state cross-classification ----

geography_lookup <- classification_combined |>
  distinct(
    name_region,
    name_state
  ) |>
  arrange(
    name_region,
    name_state
  )

state_counts_observed <- classification_combined |>
  count(
    Adjustment,
    name_region,
    name_state,
    cluster,
    trend,
    name = "n"
  )

state_counts_complete <- tidyr::crossing(
  geography_lookup,
  Adjustment = factor(
    c("Before", "After"),
    levels = c("Before", "After")
  ),
  cluster = factor(
    c("Hotspot", "Coldspot", "Neutralspot"),
    levels = c("Hotspot", "Coldspot", "Neutralspot")
  ),
  trend = factor(
    c("Increasing", "Decreasing", "Stable"),
    levels = c("Increasing", "Decreasing", "Stable")
  )
) |>
  left_join(
    state_counts_observed,
    by = c(
      "Adjustment",
      "name_region",
      "name_state",
      "cluster",
      "trend"
    )
  ) |>
  mutate(
    n = replace_na(n, 0L)
  )

state_rows <- state_counts_complete |>
  pivot_wider(
    names_from = c(trend, Adjustment),
    values_from = n,
    names_glue = "{trend}_{Adjustment}"
  ) |>
  mutate(
    Total_Before =
      Increasing_Before +
      Decreasing_Before +
      Stable_Before,
    Total_After =
      Increasing_After +
      Decreasing_After +
      Stable_After,
    row_order = 1L
  ) |>
  rename(
    Region = name_region,
    State = name_state,
    `Cluster type` = cluster
  )

region_rows <- classification_combined |>
  count(
    Adjustment,
    name_region,
    trend,
    name = "n"
  ) |>
  tidyr::complete(
    Adjustment,
    name_region,
    trend,
    fill = list(n = 0)
  ) |>
  pivot_wider(
    names_from = c(trend, Adjustment),
    values_from = n,
    names_glue = "{trend}_{Adjustment}"
  ) |>
  mutate(
    Total_Before =
      Increasing_Before +
      Decreasing_Before +
      Stable_Before,
    Total_After =
      Increasing_After +
      Decreasing_After +
      Stable_After,
    State = "All federative units",
    `Cluster type` = "–",
    row_order = 0L
  ) |>
  rename(
    Region = name_region
  )

table_s13_data <- bind_rows(
  region_rows,
  state_rows
) |>
  mutate(
    Region = factor(
      Region,
      levels = c(
        "North",
        "Northeast",
        "Southeast",
        "South",
        "Central-West"
      )
    ),
    `Cluster type` = factor(
      as.character(`Cluster type`),
      levels = c(
        "–",
        "Hotspot",
        "Coldspot",
        "Neutralspot"
      )
    )
  ) |>
  arrange(
    Region,
    row_order,
    State,
    `Cluster type`
  ) |>
  select(
    -row_order
  )

writexl::write_xlsx(
  table_s13_data,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S13_region_state_cross_classification.xlsx"
  )
)

table_s13 <- table_s13_data |>
  gt(groupname_col = "Region") |>
  tab_spanner(
    label = "Increasing",
    columns = c(
      Increasing_Before,
      Increasing_After
    )
  ) |>
  tab_spanner(
    label = "Decreasing",
    columns = c(
      Decreasing_Before,
      Decreasing_After
    )
  ) |>
  tab_spanner(
    label = "Stable",
    columns = c(
      Stable_Before,
      Stable_After
    )
  ) |>
  tab_spanner(
    label = "Total",
    columns = c(
      Total_Before,
      Total_After
    )
  ) |>
  cols_label(
    State = "Region/state",
    `Cluster type` = "Cluster type",
    Increasing_Before = "Before",
    Increasing_After = "After",
    Decreasing_Before = "Before",
    Decreasing_After = "After",
    Stable_Before = "Before",
    Stable_After = "After",
    Total_Before = "Before",
    Total_After = "After"
  ) |>
  fmt_integer(
    columns = where(is.numeric),
    use_seps = TRUE
  ) |>
  tab_header(
    title = md(
      "**Table S13. Number of municipalities according to the cross-classification of spatial clusters and temporal trends in tuberculosis mortality by region and federative units, before and after covariate adjustment, Brazil, 2010–2024.**"
    )
  )

gt::gtsave(
  table_s13,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S13_region_state_cross_classification.docx"
  )
)

# Figure 3: spatial clusters and differential temporal trends ----

profile_levels <- c(
  "Hotspot / Decreasing",
  "Hotspot / Stable",
  "Hotspot / Increasing",
  "Coldspot / Decreasing",
  "Coldspot / Stable",
  "Coldspot / Increasing",
  "Neutralspot / Decreasing",
  "Neutralspot / Stable",
  "Neutralspot / Increasing"
)

profile_colours <- c(
  "#fcae91",
  "#fb6a4a",
  "#a50f15",
  "#bdd7e7",
  "#6baed6",
  "#08519c",
  "#c7e9c0",
  "#74c476",
  "#00441b"
)

before_map_classification <- classification_before |>
  mutate(
    profile = factor(
      paste(cluster, trend, sep = " / "),
      levels = profile_levels
    )
  ) |>
  group_by(profile) |>
  mutate(
    profile_count = n(),
    profile_label = paste0(
      profile,
      " [",
      scales::comma(profile_count),
      "]"
    )
  ) |>
  ungroup() |>
  mutate(
    profile_label = factor(
      profile_label,
      levels = unique(
        profile_label[order(profile)]
      )
    )
  )

after_map_classification <- classification_after |>
  mutate(
    profile = factor(
      paste(cluster, trend, sep = " / "),
      levels = profile_levels
    )
  ) |>
  group_by(profile) |>
  mutate(
    profile_count = n(),
    profile_label = paste0(
      profile,
      " [",
      scales::comma(profile_count),
      "]"
    )
  ) |>
  ungroup() |>
  mutate(
    profile_label = factor(
      profile_label,
      levels = unique(
        profile_label[order(profile)]
      )
    )
  )

map_before_data <- shp_muni_projected |>
  left_join(
    before_map_classification |>
      select(
        code_muni,
        profile,
        profile_label
      ),
    by = "code_muni"
  )

map_after_data <- shp_muni_projected |>
  left_join(
    after_map_classification |>
      select(
        code_muni,
        profile,
        profile_label
      ),
    by = "code_muni"
  )

before_colours <- profile_colours[
  match(
    levels(map_before_data$profile_label),
    before_map_classification |>
      distinct(profile, profile_label) |>
      arrange(profile) |>
      pull(profile_label)
  )
]

after_colours <- profile_colours[
  match(
    levels(map_after_data$profile_label),
    after_map_classification |>
      distinct(profile, profile_label) |>
      arrange(profile) |>
      pull(profile_label)
  )
]

names(before_colours) <- levels(map_before_data$profile_label)
names(after_colours) <- levels(map_after_data$profile_label)

map_before <- ggplot() +
  geom_sf(
    data = map_before_data,
    aes(fill = profile_label),
    colour = NA
  ) +
  geom_sf(
    data = shp_states,
    fill = NA,
    colour = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = before_colours,
    drop = FALSE,
    name = paste(
      "Spatial cluster / differential",
      "temporal effect"
    )
  ) +
  labs(
    title = "A) Before covariate adjustment"
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

map_after <- ggplot() +
  geom_sf(
    data = map_after_data,
    aes(fill = profile_label),
    colour = NA
  ) +
  geom_sf(
    data = shp_states,
    fill = NA,
    colour = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = after_colours,
    drop = FALSE,
    name = paste(
      "Spatial cluster / differential",
      "temporal effect"
    )
  ) +
  labs(
    title = "B) After covariate adjustment"
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )  +
  ggspatial::annotation_scale(
    location = "br",
    text_cex = 0.8,
    height = unit(0.25, "cm"),
    pad_y = unit(0.6, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    style = ggspatial::north_arrow_fancy_orienteering(),
    location = "br",
    width = unit(1.5, "cm"),
    height = unit(1.5, "cm"),
    pad_x = unit(1.2, "cm"),
    pad_y = unit(1.2, "cm")
  )

figure_3 <- map_before + map_after

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_3.tiff"
  ),
  plot = figure_3,
  width = 22,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig 3. Distribution of tuberculosis hotspots, coldspots and neutral areas",
    "according to differential temporal effects before (A) and after (B)",
    "covariate adjustment in Brazil, 2010–2024."
  ),
  file.path(
    RESULTS_FIGURE_DIR,
    "Fig_3_caption.txt"
  )
)

# Save workspace
save.image("results/data/bayesian_analysis.RData")

# Reproducibility information ----

writeLines(
  capture.output(sessionInfo()),
  file.path(
    RESULTS_DATA_DIR,
    "session_info_bayesian_analysis.txt"
  )
)
