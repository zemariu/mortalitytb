# Article: Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024:
#          a Bayesian spatiotemporal modelling study
# Repository: mortalitytb
# Script: Constructs the municipal priority index, incorporates uncertainty and 
#         robustness, and generates final priority classifications.
# Author: José Mário Nunes da Silva, PhD
# Last updated: 2026-08-02

# Packages and output directories ----

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(writexl)
library(gt)
library(ggplot2)
library(stringr)
library(scales)
library(geobr)
library(sf)
library(ggspatial)

RESULTS_DATA_DIR <- "results/priority_index/outputs"
RESULTS_TABLE_DIR <- "results/priority_index/tables"
RESULTS_FIGURE_DIR <- "results/priority_index/figures"

dir.create(
  RESULTS_DATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  RESULTS_TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  RESULTS_FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Import and validate classification-stability results ----

STABILITY_FILE <- file.path(
  "results/spatiotemporal/outputs",
  "municipal_classification_stability.rds"
)

classification_stability <- readRDS(
  STABILITY_FILE
)

priority_data <- classification_stability |>
  mutate(
    code_muni = sprintf(
      "%06d",
      as.integer(code_muni)
    )
  ) |>
  arrange(
    code_muni
  )


# Central priority score ----

rr_quantiles <- quantile(
  priority_data$RR_after,
  probs = c(
    0.50,
    0.75,
    0.90
  ),
  na.rm = TRUE,
  names = FALSE
)

names(rr_quantiles) <- c(
  "q50",
  "q75",
  "q90"
)

slope_quantiles <- quantile(
  priority_data$slope_after,
  probs = c(
    0.75,
    0.90
  ),
  na.rm = TRUE,
  names = FALSE
)

names(slope_quantiles) <- c(
  "q75",
  "q90"
)

priority_index <- priority_data |>
  mutate(
    spatial_score = case_when(
      cluster_after == "Hotspot" ~ 25,
      cluster_after == "Neutralspot" ~ 10,
      cluster_after == "Coldspot" ~ 0,
      TRUE ~ 0
    ),

    persistence_score = case_when(
      hotspot_status == "Persistent hotspot" ~ 20,
      hotspot_status == "Newly identified hotspot" ~ 15,
      TRUE ~ 0
    ),

    rr_score = case_when(
      is.na(RR_after) ~ 0,
      RR_after >= rr_quantiles[["q90"]] ~ 20,
      RR_after >= rr_quantiles[["q75"]] ~ 15,
      RR_after >= rr_quantiles[["q50"]] ~ 10,
      TRUE ~ 5
    ),

    rr_probability_score = case_when(
      is.na(Pr_RR_after) ~ 0,
      Pr_RR_after >= 0.999 ~ 20,
      Pr_RR_after >= 0.990 ~ 15,
      Pr_RR_after >= 0.975 ~ 10,
      Pr_RR_after >= 0.950 ~ 5,
      TRUE ~ 0
    ),

    trend_score = case_when(
      is.na(slope_after) ~ 0,
      slope_after >= slope_quantiles[["q90"]] ~ 15,
      slope_after >= slope_quantiles[["q75"]] ~ 10,
      slope_after > 0 ~ 5,
      TRUE ~ 0
    ),

    trend_probability_score = case_when(
      is.na(Pr_slope_after) ~ 0,
      Pr_slope_after >= 0.99 ~ 15,
      Pr_slope_after >= 0.95 ~ 12,
      Pr_slope_after >= 0.90 ~ 9,
      Pr_slope_after >= 0.80 ~ 5,
      TRUE ~ 0
    ),

    regime_score = case_when(
      is.na(regime_entropy_normalized) ~ 0,
      regime_entropy_normalized <= 0.15 ~ 10,
      regime_entropy_normalized <= 0.25 ~ 8,
      regime_entropy_normalized <= 0.35 ~ 6,
      TRUE ~ 4
    ),

    state_context_score = case_when(
      is.na(state_strict_persistence) ~ 0,
      state_strict_persistence >= 0.90 ~ 10,
      state_strict_persistence >= 0.80 ~ 8,
      state_strict_persistence >= 0.70 ~ 6,
      TRUE ~ 4
    ),

    central_raw_priority_score =
      spatial_score +
      persistence_score +
      rr_score +
      rr_probability_score +
      trend_score +
      trend_probability_score +
      regime_score +
      state_context_score
  )

central_score_max <- max(
  priority_index$central_raw_priority_score,
  na.rm = TRUE
)

priority_index <- priority_index |>
  mutate(
    central_priority_score =
      100 *
      central_raw_priority_score /
      central_score_max,

    priority_profile = case_when(
      cluster_after == "Hotspot" &
        cluster_persistent &
        Pr_RR_after >= 0.95 ~
        "Persistent critical hotspot",

      hotspot_status == "Newly identified hotspot" &
        Pr_RR_after >= 0.95 &
        Pr_slope_after >= 0.80 ~
        "Emerging hotspot",

      cluster_after != "Hotspot" &
        RR_after >= rr_quantiles[["q90"]] &
        Pr_RR_after >= 0.90 &
        Pr_slope_after >= 0.80 ~
        "Silent high-risk area",

      TRUE ~
        "Other areas"
    ),

    central_priority_class = case_when(
      central_priority_score >= 80 ~ "Very high",
      central_priority_score >= 60 ~ "High",
      central_priority_score >= 40 ~ "Moderate",
      TRUE ~ "Low"
    ),

    central_priority_class = factor(
      central_priority_class,
      levels = c(
        "Low",
        "Moderate",
        "High",
        "Very high"
      )
    ),

    spatial_evidence = case_when(
      cluster_after == "Hotspot" ~ "Current hotspot",
      cluster_after == "Neutralspot" ~ "Current neutral area",
      cluster_after == "Coldspot" ~ "Current coldspot",
      TRUE ~ NA_character_
    ),

    risk_evidence = case_when(
      RR_after > 1 ~ "Elevated risk",
      RR_after <= 1 ~ "No excess risk",
      TRUE ~ NA_character_
    ),

    trend_evidence = case_when(
      slope_after > 0 ~ "Increasing trend",
      slope_after <= 0 ~ "Non-increasing trend",
      TRUE ~ NA_character_
    ),

    regime_stability = if_else(
      strict_persistent,
      "Persistent regime",
      "Non-persistent regime",
      missing = NA_character_
    ),

    spatial_pattern_stability = if_else(
      cluster_persistent,
      "Persistent spatial pattern",
      "Non-persistent spatial pattern",
      missing = NA_character_
    ),

    trend_stability = if_else(
      trend_persistent,
      "Persistent trend",
      "Non-persistent trend",
      missing = NA_character_
    )
  )

priority_score_parameters <- tibble(
  Parameter = c(
    "RR 50th percentile",
    "RR 75th percentile",
    "RR 90th percentile",
    "Temporal coefficient 75th percentile",
    "Temporal coefficient 90th percentile",
    "Observed maximum raw priority score",
    "Theoretical maximum raw priority score"
  ),
  Value = c(
    rr_quantiles[["q50"]],
    rr_quantiles[["q75"]],
    rr_quantiles[["q90"]],
    slope_quantiles[["q75"]],
    slope_quantiles[["q90"]],
    central_score_max,
    135
  )
)

writexl::write_xlsx(
  priority_score_parameters,
  file.path(
    RESULTS_DATA_DIR,
    "priority_score_parameters.xlsx"
  )
)

write_csv(
  priority_index,
  file.path(
    RESULTS_DATA_DIR,
    "central_priority_index.csv"
  )
)

saveRDS(
  priority_index,
  file.path(
    RESULTS_DATA_DIR,
    "central_priority_index.rds"
  )
)

# Explicit uncertainty score ----

safe_zscore <- function(x) {
  mean_x <- mean(
    x,
    na.rm = TRUE
  )

  sd_x <- sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(sd_x) ||
      sd_x == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (
    x -
      mean_x
  ) /
    sd_x
}

safe_rescale_0_100 <- function(x) {
  minimum_x <- min(
    x,
    na.rm = TRUE
  )

  maximum_x <- max(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(minimum_x) ||
      !is.finite(maximum_x) ||
      minimum_x == maximum_x
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  100 *
    (
      x -
        minimum_x
    ) /
    (
      maximum_x -
        minimum_x
    )
}

priority_index <- priority_index |>
  mutate(
    rr_interval_width =
      RR_upper_95 -
      RR_lower_95,

    rr_relative_interval_width = case_when(
      !is.na(RR_after) &
        RR_after != 0 ~
        rr_interval_width /
        abs(RR_after),

      TRUE ~
        NA_real_
    ),

    slope_interval_width =
      slope_upper_95 -
      slope_lower_95,

    slope_relative_interval_width = case_when(
      !is.na(slope_after) &
        abs(slope_after) > 1e-8 ~
        slope_interval_width /
        abs(slope_after),

      TRUE ~
        NA_real_
    ),

    rr_probability_uncertainty = case_when(
      !is.na(Pr_RR_after) ~
        pmin(
          pmax(
            1 -
              2 *
              abs(
                Pr_RR_after -
                  0.5
              ),
            0
          ),
          1
        ),

      TRUE ~
        NA_real_
    ),

    slope_probability_uncertainty = case_when(
      !is.na(Pr_slope_after) ~
        pmin(
          pmax(
            1 -
              2 *
              abs(
                Pr_slope_after -
                  0.5
              ),
            0
          ),
          1
        ),

      TRUE ~
        NA_real_
    ),

    rr_uncertainty_component = rowMeans(
      cbind(
        safe_zscore(
          rr_interval_width
        ),
        safe_zscore(
          rr_relative_interval_width
        ),
        safe_zscore(
          rr_probability_uncertainty
        )
      ),
      na.rm = TRUE
    ),

    slope_uncertainty_component = rowMeans(
      cbind(
        safe_zscore(
          slope_interval_width
        ),
        safe_zscore(
          slope_relative_interval_width
        ),
        safe_zscore(
          slope_probability_uncertainty
        )
      ),
      na.rm = TRUE
    ),

    structural_uncertainty_component = rowMeans(
      cbind(
        safe_zscore(
          regime_entropy_normalized
        ),
        safe_zscore(
          state_cluster_entropy_normalized
        )
      ),
      na.rm = TRUE
    )
  ) |>
  mutate(
    rr_uncertainty_component = if_else(
      is.nan(
        rr_uncertainty_component
      ),
      NA_real_,
      rr_uncertainty_component
    ),

    slope_uncertainty_component = if_else(
      is.nan(
        slope_uncertainty_component
      ),
      NA_real_,
      slope_uncertainty_component
    ),

    structural_uncertainty_component = if_else(
      is.nan(
        structural_uncertainty_component
      ),
      NA_real_,
      structural_uncertainty_component
    ),

    raw_uncertainty_index = rowMeans(
      cbind(
        rr_uncertainty_component,
        slope_uncertainty_component,
        structural_uncertainty_component
      ),
      na.rm = TRUE
    ),

    raw_uncertainty_index = if_else(
      is.nan(
        raw_uncertainty_index
      ),
      NA_real_,
      raw_uncertainty_index
    ),

    uncertainty_score = safe_rescale_0_100(
      raw_uncertainty_index
    )
  )

uncertainty_quantiles <- quantile(
  priority_index$uncertainty_score,
  probs = c(
    0.25,
    0.75
  ),
  na.rm = TRUE,
  names = FALSE
)

priority_index <- priority_index |>
  mutate(
    confidence_class = case_when(
      uncertainty_score <= uncertainty_quantiles[1L] ~
        "High confidence",

      uncertainty_score <= uncertainty_quantiles[2L] ~
        "Moderate confidence",

      TRUE ~
        "Low confidence"
    ),

    confidence_class = factor(
      confidence_class,
      levels = c(
        "Low confidence",
        "Moderate confidence",
        "High confidence"
      )
    )
  )

# Uncertainty-adjusted priority score ----

UNCERTAINTY_PENALTY_WEIGHT <- 0.20

priority_index <- priority_index |>
  mutate(
    uncertainty_penalty =
      UNCERTAINTY_PENALTY_WEIGHT *
      uncertainty_score,

    adjusted_priority_score = pmax(
      central_priority_score -
        uncertainty_penalty,
      0
    ),

    absolute_score_loss =
      central_priority_score -
      adjusted_priority_score,

    relative_score_loss_pct = case_when(
      central_priority_score > 0 ~
        100 *
        absolute_score_loss /
        central_priority_score,

      TRUE ~
        NA_real_
    ),

    adjusted_priority_class = case_when(
      adjusted_priority_score >= 80 ~ "Very high",
      adjusted_priority_score >= 60 ~ "High",
      adjusted_priority_score >= 40 ~ "Moderate",
      TRUE ~ "Low"
    ),

    adjusted_priority_class = factor(
      adjusted_priority_class,
      levels = c(
        "Low",
        "Moderate",
        "High",
        "Very high"
      )
    ),

    class_changed_after_uncertainty =
      adjusted_priority_class !=
      central_priority_class,

    class_change_direction = case_when(
      as.integer(
        adjusted_priority_class
      ) <
        as.integer(
          central_priority_class
        ) ~
        "Downgraded",

      as.integer(
        adjusted_priority_class
      ) >
        as.integer(
          central_priority_class
        ) ~
        "Upgraded",

      TRUE ~
        "Unchanged"
    )
  )

write_csv(
  priority_index,
  file.path(
    RESULTS_DATA_DIR,
    "uncertainty_adjusted_priority_index.csv"
  )
)

# Posterior score simulation ----

SIMULATION_SEED <- 2026L
N_SIMULATIONS <- 1000L

set.seed(
  SIMULATION_SEED
)

simulation_input <- priority_index |>
  mutate(
    rr_sd_approx =
      (
        RR_upper_95 -
          RR_lower_95
      ) /
      3.92,

    slope_sd_approx =
      (
        slope_upper_95 -
          slope_lower_95
      ) /
      3.92,

    fixed_score_component =
      spatial_score +
      persistence_score +
      rr_probability_score +
      trend_probability_score +
      regime_score +
      state_context_score
  )

simulate_municipality_score <- function(
    rr_mean,
    rr_sd,
    slope_mean,
    slope_sd,
    fixed_score,
    n_simulations,
    rr_q50,
    rr_q75,
    rr_q90,
    slope_q75,
    slope_q90
) {
  rr_sd_used <- ifelse(
    is.na(rr_sd) |
      rr_sd <= 0,
    1e-6,
    rr_sd
  )

  slope_sd_used <- ifelse(
    is.na(slope_sd) |
      slope_sd <= 0,
    1e-6,
    slope_sd
  )

  rr_draw <- pmax(
    rnorm(
      n_simulations,
      mean = rr_mean,
      sd = rr_sd_used
    ),
    0
  )

  slope_draw <- rnorm(
    n_simulations,
    mean = slope_mean,
    sd = slope_sd_used
  )

  rr_score_draw <- case_when(
    rr_draw >= rr_q90 ~ 20,
    rr_draw >= rr_q75 ~ 15,
    rr_draw >= rr_q50 ~ 10,
    TRUE ~ 5
  )

  trend_score_draw <- case_when(
    slope_draw >= slope_q90 ~ 15,
    slope_draw >= slope_q75 ~ 10,
    slope_draw > 0 ~ 5,
    TRUE ~ 0
  )

  tibble(
    simulated_raw_priority_score =
      fixed_score +
      rr_score_draw +
      trend_score_draw
  )
}

simulated_scores_long <- simulation_input |>
  select(
    code_muni,
    municipality_name,
    state_name,
    RR_after,
    rr_sd_approx,
    slope_after,
    slope_sd_approx,
    fixed_score_component
  ) |>
  mutate(
    simulated_data = pmap(
      list(
        RR_after,
        rr_sd_approx,
        slope_after,
        slope_sd_approx,
        fixed_score_component
      ),
      ~ simulate_municipality_score(
        rr_mean = ..1,
        rr_sd = ..2,
        slope_mean = ..3,
        slope_sd = ..4,
        fixed_score = ..5,
        n_simulations = N_SIMULATIONS,
        rr_q50 = rr_quantiles[["q50"]],
        rr_q75 = rr_quantiles[["q75"]],
        rr_q90 = rr_quantiles[["q90"]],
        slope_q75 = slope_quantiles[["q75"]],
        slope_q90 = slope_quantiles[["q90"]]
      )
    )
  ) |>
  select(
    code_muni,
    municipality_name,
    state_name,
    simulated_data
  ) |>
  unnest(
    cols = simulated_data
  ) |>
  group_by(
    code_muni
  ) |>
  mutate(
    simulation_id = row_number()
  ) |>
  ungroup() |>
  mutate(
    simulated_priority_score = pmin(
      100,
      100 *
        simulated_raw_priority_score /
        central_score_max
    )
  )

simulation_summary <- simulated_scores_long |>
  group_by(
    code_muni,
    municipality_name,
    state_name
  ) |>
  summarise(
    simulated_score_mean = mean(
      simulated_priority_score,
      na.rm = TRUE
    ),

    simulated_score_median = median(
      simulated_priority_score,
      na.rm = TRUE
    ),

    simulated_score_sd = sd(
      simulated_priority_score,
      na.rm = TRUE
    ),

    simulated_score_lower_95 = quantile(
      simulated_priority_score,
      0.025,
      na.rm = TRUE
    ),

    simulated_score_q25 = quantile(
      simulated_priority_score,
      0.25,
      na.rm = TRUE
    ),

    simulated_score_q75 = quantile(
      simulated_priority_score,
      0.75,
      na.rm = TRUE
    ),

    simulated_score_upper_95 = quantile(
      simulated_priority_score,
      0.975,
      na.rm = TRUE
    ),

    Pr_score_very_high = mean(
      simulated_priority_score >= 80,
      na.rm = TRUE
    ),

    Pr_score_high_or_higher = mean(
      simulated_priority_score >= 60,
      na.rm = TRUE
    ),

    Pr_score_moderate_or_higher = mean(
      simulated_priority_score >= 40,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

write_csv(
  simulation_summary,
  file.path(
    RESULTS_DATA_DIR,
    "priority_score_simulation_summary.csv"
  )
)

priority_index <- priority_index |>
  left_join(
    simulation_summary,
    by = c(
      "code_muni",
      "municipality_name",
      "state_name"
    )
  ) |>
  mutate(
    simulated_score_interval_width =
      simulated_score_upper_95 -
      simulated_score_lower_95
  )

# Probabilistic ranking stability ----

ranking_long <- simulated_scores_long |>
  group_by(
    simulation_id
  ) |>
  mutate(
    simulated_rank = min_rank(
      desc(
        simulated_priority_score
      )
    )
  ) |>
  ungroup()

ranking_summary <- ranking_long |>
  group_by(
    code_muni,
    municipality_name,
    state_name
  ) |>
  summarise(
    mean_rank = mean(
      simulated_rank,
      na.rm = TRUE
    ),

    median_rank = median(
      simulated_rank,
      na.rm = TRUE
    ),

    sd_rank = sd(
      simulated_rank,
      na.rm = TRUE
    ),

    rank_lower_95 = quantile(
      simulated_rank,
      0.025,
      na.rm = TRUE
    ),

    rank_q25 = quantile(
      simulated_rank,
      0.25,
      na.rm = TRUE
    ),

    rank_q75 = quantile(
      simulated_rank,
      0.75,
      na.rm = TRUE
    ),

    rank_upper_95 = quantile(
      simulated_rank,
      0.975,
      na.rm = TRUE
    ),

    Pr_top10 = mean(
      simulated_rank <= 10,
      na.rm = TRUE
    ),

    Pr_top50 = mean(
      simulated_rank <= 50,
      na.rm = TRUE
    ),

    Pr_top100 = mean(
      simulated_rank <= 100,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>
  mutate(
    rank_interval_width =
      rank_upper_95 -
      rank_lower_95,

    ranking_robustness = case_when(
      Pr_top50 >= 0.90 ~ "Robust Top 50",
      Pr_top50 >= 0.50 ~ "Borderline Top 50",
      Pr_top50 >= 0.10 ~ "Occasional Top 50",
      TRUE ~ "Rarely Top 50"
    ),

    ranking_robustness = factor(
      ranking_robustness,
      levels = c(
        "Rarely Top 50",
        "Occasional Top 50",
        "Borderline Top 50",
        "Robust Top 50"
      )
    )
  )

write_csv(
  ranking_summary,
  file.path(
    RESULTS_DATA_DIR,
    "priority_ranking_simulation_summary.csv"
  )
)

priority_index <- priority_index |>
  left_join(
    ranking_summary,
    by = c(
      "code_muni",
      "municipality_name",
      "state_name"
    )
  )

rm(
  ranking_long,
  simulated_scores_long
)

gc(
  verbose = FALSE
)

# Sensitivity analysis ----

sensitivity_scenarios <- tribble(
  ~scenario, ~uncertainty_weight, ~w_spatial, ~w_persistence, ~w_rr, ~w_pr_rr, ~w_trend, ~w_pr_trend, ~w_regime, ~w_state, ~cut_very_high, ~cut_high, ~cut_moderate,
  "Base", 0.20, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 80, 60, 40,
  "Lower uncertainty penalty", 0.10, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 80, 60, 40,
  "Higher uncertainty penalty", 0.30, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 80, 60, 40,
  "Higher risk weight", 0.20, 1.00, 1.00, 1.20, 1.20, 1.00, 1.00, 1.00, 1.00, 80, 60, 40,
  "Higher persistence weight", 0.20, 1.00, 1.25, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 80, 60, 40,
  "Higher trend weight", 0.20, 1.00, 1.00, 1.00, 1.00, 1.20, 1.20, 1.00, 1.00, 80, 60, 40,
  "Lower structural-context weight", 0.20, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.75, 0.75, 80, 60, 40,
  "Stricter priority thresholds", 0.20, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 85, 65, 45,
  "Looser priority thresholds", 0.20, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 75, 55, 35
)

writexl::write_xlsx(
  sensitivity_scenarios,
  file.path(
    RESULTS_DATA_DIR,
    "priority_index_sensitivity_scenarios.xlsx"
  )
)

run_sensitivity_scenario <- function(
    data,
    scenario,
    uncertainty_weight,
    w_spatial,
    w_persistence,
    w_rr,
    w_pr_rr,
    w_trend,
    w_pr_trend,
    w_regime,
    w_state,
    cut_very_high,
    cut_high,
    cut_moderate
) {
  scenario_data <- data |>
    mutate(
      scenario_raw_score =
        w_spatial *
        spatial_score +
        w_persistence *
        persistence_score +
        w_rr *
        rr_score +
        w_pr_rr *
        rr_probability_score +
        w_trend *
        trend_score +
        w_pr_trend *
        trend_probability_score +
        w_regime *
        regime_score +
        w_state *
        state_context_score
    )

  scenario_maximum <- max(
    scenario_data$scenario_raw_score,
    na.rm = TRUE
  )

  scenario_data |>
    mutate(
      scenario_central_score =
        100 *
        scenario_raw_score /
        scenario_maximum,

      scenario_uncertainty_penalty =
        uncertainty_weight *
        uncertainty_score,

      scenario_adjusted_score = pmax(
        scenario_central_score -
          scenario_uncertainty_penalty,
        0
      ),

      scenario_priority_class = case_when(
        scenario_adjusted_score >= cut_very_high ~ "Very high",
        scenario_adjusted_score >= cut_high ~ "High",
        scenario_adjusted_score >= cut_moderate ~ "Moderate",
        TRUE ~ "Low"
      ),

      scenario_priority_class = factor(
        scenario_priority_class,
        levels = c(
          "Low",
          "Moderate",
          "High",
          "Very high"
        )
      ),

      scenario_rank = min_rank(
        desc(
          scenario_adjusted_score
        )
      ),

      scenario = scenario
    ) |>
    select(
      code_muni,
      municipality_name,
      state_name,
      scenario,
      scenario_central_score,
      scenario_uncertainty_penalty,
      scenario_adjusted_score,
      scenario_priority_class,
      scenario_rank
    )
}

sensitivity_results <- pmap_dfr(
  sensitivity_scenarios,
  function(
      scenario,
      uncertainty_weight,
      w_spatial,
      w_persistence,
      w_rr,
      w_pr_rr,
      w_trend,
      w_pr_trend,
      w_regime,
      w_state,
      cut_very_high,
      cut_high,
      cut_moderate
  ) {
    run_sensitivity_scenario(
      data = priority_index,
      scenario = scenario,
      uncertainty_weight = uncertainty_weight,
      w_spatial = w_spatial,
      w_persistence = w_persistence,
      w_rr = w_rr,
      w_pr_rr = w_pr_rr,
      w_trend = w_trend,
      w_pr_trend = w_pr_trend,
      w_regime = w_regime,
      w_state = w_state,
      cut_very_high = cut_very_high,
      cut_high = cut_high,
      cut_moderate = cut_moderate
    )
  }
)

base_scenario <- sensitivity_results |>
  filter(
    scenario == "Base"
  ) |>
  select(
    code_muni,
    base_rank = scenario_rank,
    base_class = scenario_priority_class,
    base_adjusted_score =
      scenario_adjusted_score
  )

base_score_difference <- priority_index |>
  select(
    code_muni,
    adjusted_priority_score
  ) |>
  left_join(
    base_scenario,
    by = "code_muni"
  ) |>
  summarise(
    maximum_absolute_difference = max(
      abs(
        adjusted_priority_score -
          base_adjusted_score
      ),
      na.rm = TRUE
    )
  ) |>
  pull(
    maximum_absolute_difference
  )

if (
  is.finite(base_score_difference) &&
    base_score_difference > 1e-8
) {
  warning(
    "The base sensitivity scenario does not exactly reproduce the uncertainty-adjusted score."
  )
}

base_top50 <- base_scenario |>
  arrange(
    base_rank
  ) |>
  slice_head(
    n = 50
  ) |>
  pull(
    code_muni
  )

base_very_high <- sensitivity_results |>
  filter(
    scenario == "Base",
    scenario_priority_class == "Very high"
  ) |>
  pull(
    code_muni
  )

sensitivity_comparison <- sensitivity_results |>
  left_join(
    base_scenario,
    by = "code_muni"
  ) |>
  group_by(
    scenario
  ) |>
  summarise(
    Number_of_municipalities = n(),

    Spearman_rank_correlation = suppressWarnings(
      cor(
        scenario_rank,
        base_rank,
        method = "spearman",
        use = "complete.obs"
      )
    ),

    Mean_absolute_rank_shift = mean(
      abs(
        scenario_rank -
          base_rank
      ),
      na.rm = TRUE
    ),

    Median_absolute_rank_shift = median(
      abs(
        scenario_rank -
          base_rank
      ),
      na.rm = TRUE
    ),

    Top_50_overlap_n = sum(
      code_muni %in%
        base_top50 &
        scenario_rank <= 50,
      na.rm = TRUE
    ),

    Very_high_overlap_n = sum(
      code_muni %in%
        base_very_high &
        scenario_priority_class == "Very high",
      na.rm = TRUE
    ),

    Current_very_high_n = sum(
      scenario_priority_class == "Very high",
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>
  mutate(
    Top_50_overlap_pct =
      100 *
      Top_50_overlap_n /
      50,

    Very_high_overlap_pct = case_when(
      length(
        base_very_high
      ) > 0 ~
        100 *
        Very_high_overlap_n /
        length(
          base_very_high
        ),

      TRUE ~
        NA_real_
    )
  )

municipality_sensitivity <- sensitivity_results |>
  group_by(
    code_muni,
    municipality_name,
    state_name
  ) |>
  summarise(
    Number_of_scenarios = n(),

    Mean_rank_across_scenarios = mean(
      scenario_rank,
      na.rm = TRUE
    ),

    SD_rank_across_scenarios = sd(
      scenario_rank,
      na.rm = TRUE
    ),

    Minimum_rank_across_scenarios = min(
      scenario_rank,
      na.rm = TRUE
    ),

    Maximum_rank_across_scenarios = max(
      scenario_rank,
      na.rm = TRUE
    ),

    Pr_top50_across_scenarios = mean(
      scenario_rank <= 50,
      na.rm = TRUE
    ),

    Pr_very_high_across_scenarios = mean(
      scenario_priority_class == "Very high",
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>
  mutate(
    scenario_robustness_class = case_when(
      Pr_top50_across_scenarios >= 0.90 ~
        "Robust across scenarios",

      Pr_top50_across_scenarios >= 0.50 ~
        "Moderately robust across scenarios",

      Pr_top50_across_scenarios > 0 ~
        "Sensitive to scenario assumptions",

      TRUE ~
        "Never in Top 50"
    ),

    scenario_robustness_class = factor(
      scenario_robustness_class,
      levels = c(
        "Never in Top 50",
        "Sensitive to scenario assumptions",
        "Moderately robust across scenarios",
        "Robust across scenarios"
      )
    )
  )

write_csv(
  sensitivity_results,
  file.path(
    RESULTS_DATA_DIR,
    "priority_index_sensitivity_results.csv"
  )
)

writexl::write_xlsx(
  sensitivity_comparison,
  file.path(
    RESULTS_DATA_DIR,
    "priority_index_sensitivity_summary.xlsx"
  )
)

write_csv(
  municipality_sensitivity,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_priority_sensitivity.csv"
  )
)

priority_index <- priority_index |>
  left_join(
    municipality_sensitivity,
    by = c(
      "code_muni",
      "municipality_name",
      "state_name"
    )
  )

# Final priority decision and operational classes ----

priority_index <- priority_index |>
  mutate(
    final_priority_class = case_when(
      adjusted_priority_score >= 80 &
        confidence_class == "High confidence" &
        Pr_top50 >= 0.80 &
        Pr_top50_across_scenarios >= 0.80 ~
        "Very high priority (robust)",

      adjusted_priority_score >= 80 ~
        "Very high priority (uncertain)",

      adjusted_priority_score >= 60 &
        (
          Pr_top50 >= 0.15 |
            Pr_top50_across_scenarios >= 0.15 |
            confidence_class %in%
            c(
              "High confidence",
              "Moderate confidence"
            )
        ) ~
        "High priority",

      adjusted_priority_score >= 40 ~
        "Moderate priority",

      adjusted_priority_score < 40 &
        (
          Pr_score_high_or_higher >= 0.05 |
            Pr_top100 >= 0.05 |
            Pr_slope_after >= 0.80 |
            hotspot_status ==
            "Newly identified hotspot"
        ) ~
        "Low priority (potential risk)",

      TRUE ~
        "Low priority"
    ),

    final_priority_class = factor(
      final_priority_class,
      levels = c(
        "Low priority",
        "Low priority (potential risk)",
        "Moderate priority",
        "High priority",
        "Very high priority (uncertain)",
        "Very high priority (robust)"
      )
    ),

    operational_priority = case_when(
      final_priority_class ==
        "Very high priority (robust)" ~
        "Immediate action required",

      final_priority_class ==
        "Very high priority (uncertain)" ~
        "Priority action with technical review",

      final_priority_class ==
        "High priority" ~
        "Priority planning",

      final_priority_class ==
        "Moderate priority" ~
        "Enhanced monitoring and review",

      final_priority_class ==
        "Low priority (potential risk)" ~
        "Monitor for early signals",

      TRUE ~
        "Routine monitoring"
    ),

    operational_priority = factor(
      operational_priority,
      levels = c(
        "Routine monitoring",
        "Monitor for early signals",
        "Enhanced monitoring and review",
        "Priority planning",
        "Priority action with technical review",
        "Immediate action required"
      )
    )
  ) |>
  arrange(
    desc(
      as.integer(
        operational_priority
      )
    ),
    desc(
      adjusted_priority_score
    ),
    desc(
      Pr_top50
    ),
    mean_rank
  )

final_decision_table <- priority_index |>
  select(
    code_muni,
    municipality_name,
    state_name,
    region_name,
    regime_after,
    cluster_after,
    trend_after,
    hotspot_status,
    priority_profile,
    central_raw_priority_score,
    central_priority_score,
    adjusted_priority_score,
    absolute_score_loss,
    uncertainty_score,
    confidence_class,
    simulated_score_median,
    simulated_score_lower_95,
    simulated_score_upper_95,
    mean_rank,
    median_rank,
    rank_lower_95,
    rank_upper_95,
    Pr_score_very_high,
    Pr_score_high_or_higher,
    Pr_top10,
    Pr_top50,
    Pr_top100,
    Pr_top50_across_scenarios,
    Pr_very_high_across_scenarios,
    scenario_robustness_class,
    ranking_robustness,
    final_priority_class,
    operational_priority
  )

priority_class_distribution <- priority_index |>
  dplyr::count(
    final_priority_class,
    name = "Number_of_municipalities",
    .drop = FALSE
  ) |>
  mutate(
    Percentage =
      100 *
      Number_of_municipalities /
      sum(
        Number_of_municipalities
      )
  )

operational_priority_distribution <- priority_index |>
  dplyr::count(
    operational_priority,
    name = "Number_of_municipalities",
    .drop = FALSE
  ) |>
  mutate(
    Percentage =
      100 *
      Number_of_municipalities /
      sum(
        Number_of_municipalities
      )
  )

report_summary <- final_decision_table |>
  transmute(
    code_muni,
    municipality_name,
    state_name,
    region_name,
    central_priority_score = round(
      central_priority_score,
      1
    ),
    adjusted_priority_score = round(
      adjusted_priority_score,
      1
    ),
    uncertainty_score = round(
      uncertainty_score,
      1
    ),
    simulated_score_95_interval = paste0(
      round(
        simulated_score_lower_95,
        1
      ),
      "–",
      round(
        simulated_score_upper_95,
        1
      )
    ),
    Pr_top50_pct = round(
      100 *
      Pr_top50,
      1
    ),
    Pr_top50_across_scenarios_pct = round(
      100 *
      Pr_top50_across_scenarios,
      1
    ),
    final_priority_class,
    operational_priority
  )

immediate_action_list <- final_decision_table |>
  filter(
    operational_priority ==
      "Immediate action required"
  )

technical_review_list <- final_decision_table |>
  filter(
    operational_priority ==
      "Priority action with technical review"
  )

priority_planning_list <- final_decision_table |>
  filter(
    operational_priority ==
      "Priority planning"
  )

write_csv(
  priority_index,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_priority_index.csv"
  )
)

saveRDS(
  priority_index,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_priority_index.rds"
  )
)

writexl::write_xlsx(
  final_decision_table,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_priority_decision_table.xlsx"
  )
)

write_csv(
  report_summary,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_priority_report_summary.csv"
  )
)

writexl::write_xlsx(
  priority_class_distribution,
  file.path(
    RESULTS_DATA_DIR,
    "priority_class_distribution.xlsx"
  )
)

writexl::write_xlsx(
  operational_priority_distribution,
  file.path(
    RESULTS_DATA_DIR,
    "operational_priority_distribution.xlsx"
  )
)

writexl::write_xlsx(
  immediate_action_list,
  file.path(
    RESULTS_DATA_DIR,
    "immediate_action_municipalities.xlsx"
  )
)

writexl::write_xlsx(
  technical_review_list,
  file.path(
    RESULTS_DATA_DIR,
    "technical_review_municipalities.xlsx"
  )
)

writexl::write_xlsx(
  priority_planning_list,
  file.path(
    RESULTS_DATA_DIR,
    "priority_planning_municipalities.xlsx"
  )
)

writexl::write_xlsx(
  list(
    Priority_decisions =
      final_decision_table,

    Report_summary =
      report_summary,

    Class_distribution =
      priority_class_distribution,

    Operational_distribution =
      operational_priority_distribution,

    Sensitivity_summary =
      sensitivity_comparison,

    Immediate_action =
      immediate_action_list,

    Technical_review =
      technical_review_list,

    Priority_planning =
      priority_planning_list
  ),
  file.path(
    RESULTS_DATA_DIR,
    "priority_index_results.xlsx"
  )
)

# Priority index summary national ----

br_summary_advanced <- final_decision_table |>
  summarise(
    region_name = "Brazil",
    
    state_name = "Brazil",
    
    total_municipalities = n(),
    
    n_very_high_robust = sum(
      final_priority_class == "Very high priority (robust)",
      na.rm = TRUE
    ),
    n_very_high_uncertain = sum(
      final_priority_class == "Very high priority (uncertain)",
      na.rm = TRUE
    ),
    n_high = sum(
      final_priority_class == "High priority",
      na.rm = TRUE
    ),
    n_moderate = sum(
      final_priority_class == "Moderate priority",
      na.rm = TRUE
    ),
    n_potential_risk = sum(
      final_priority_class == "Low priority (potential risk)",
      na.rm = TRUE
    ),
    n_low = sum(
      final_priority_class == "Low priority",
      na.rm = TRUE
    ),
    
    pct_very_high_robust = 100 * n_very_high_robust / total_municipalities,
    pct_very_high_uncertain = 100 * n_very_high_uncertain / total_municipalities,
    pct_high = 100 * n_high / total_municipalities,
    pct_moderate = 100 * n_moderate / total_municipalities,
    pct_potential_risk = 100 * n_potential_risk / total_municipalities,
    pct_low = 100 * n_low / total_municipalities,
    
    mean_score = mean(
      adjusted_priority_score,
      na.rm = TRUE
    ),
    sd_score = sd(
      adjusted_priority_score,
      na.rm = TRUE
    ),
    mean_sd_score = paste0(
      round(mean_score, 1),
      " (",
      round(sd_score, 1),
      ")"
    ),
    
    median_score = median(
      adjusted_priority_score,
      na.rm = TRUE
    ),
    q25_score = quantile(
      adjusted_priority_score,
      0.25,
      na.rm = TRUE
    ),
    q75_score = quantile(
      adjusted_priority_score,
      0.75,
      na.rm = TRUE
    ),
    median_iqr_score = paste0(
      round(median_score, 1),
      " (",
      round(q25_score, 1),
      "–",
      round(q75_score, 1),
      ")"
    ),
    
    median_uncertainty = median(
      uncertainty_score,
      na.rm = TRUE
    ),
    q25_uncertainty = quantile(
      uncertainty_score,
      0.25,
      na.rm = TRUE
    ),
    q75_uncertainty = quantile(
      uncertainty_score,
      0.75,
      na.rm = TRUE
    ),
    median_iqr_uncertainty = paste0(
      round(median_uncertainty, 1),
      " (",
      round(q25_uncertainty, 1),
      "–",
      round(q75_uncertainty, 1),
      ")"
    ),
    
    .groups = "drop"
  )


# Priority index summary by federative unit ----

uf_summary_advanced <- final_decision_table |>
  group_by(
    region_name,
    state_name
  ) |>
  summarise(
    total_municipalities = n(),

    n_very_high_robust = sum(
      final_priority_class == "Very high priority (robust)",
      na.rm = TRUE
    ),
    n_very_high_uncertain = sum(
      final_priority_class == "Very high priority (uncertain)",
      na.rm = TRUE
    ),
    n_high = sum(
      final_priority_class == "High priority",
      na.rm = TRUE
    ),
    n_moderate = sum(
      final_priority_class == "Moderate priority",
      na.rm = TRUE
    ),
    n_potential_risk = sum(
      final_priority_class == "Low priority (potential risk)",
      na.rm = TRUE
    ),
    n_low = sum(
      final_priority_class == "Low priority",
      na.rm = TRUE
    ),

    pct_very_high_robust = 100 * n_very_high_robust / total_municipalities,
    pct_very_high_uncertain = 100 * n_very_high_uncertain / total_municipalities,
    pct_high = 100 * n_high / total_municipalities,
    pct_moderate = 100 * n_moderate / total_municipalities,
    pct_potential_risk = 100 * n_potential_risk / total_municipalities,
    pct_low = 100 * n_low / total_municipalities,

    mean_score = mean(
      adjusted_priority_score,
      na.rm = TRUE
    ),
    sd_score = sd(
      adjusted_priority_score,
      na.rm = TRUE
    ),
    mean_sd_score = paste0(
      round(mean_score, 1),
      " (",
      round(sd_score, 1),
      ")"
    ),

    median_score = median(
      adjusted_priority_score,
      na.rm = TRUE
    ),
    q25_score = quantile(
      adjusted_priority_score,
      0.25,
      na.rm = TRUE
    ),
    q75_score = quantile(
      adjusted_priority_score,
      0.75,
      na.rm = TRUE
    ),
    median_iqr_score = paste0(
      round(median_score, 1),
      " (",
      round(q25_score, 1),
      "–",
      round(q75_score, 1),
      ")"
    ),

    median_uncertainty = median(
      uncertainty_score,
      na.rm = TRUE
    ),
    q25_uncertainty = quantile(
      uncertainty_score,
      0.25,
      na.rm = TRUE
    ),
    q75_uncertainty = quantile(
      uncertainty_score,
      0.75,
      na.rm = TRUE
    ),
    median_iqr_uncertainty = paste0(
      round(median_uncertainty, 1),
      " (",
      round(q25_uncertainty, 1),
      "–",
      round(q75_uncertainty, 1),
      ")"
    ),

    .groups = "drop"
  )

summary_advanced <- br_summary_advanced |> 
  bind_rows(uf_summary_advanced)

writexl::write_xlsx(
  summary_advanced,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S17_priority_index_by_state.xlsx"
  )
)

priority_index_by_state_table <- summary_advanced |>
  mutate(
    region_name = case_when(
      region_name %in% c("Brazil", "Brasil") ~ "Brazil",
      region_name %in% c("Norte", "North") ~ "North",
      region_name %in% c("Nordeste", "Northeast") ~ "Northeast",
      region_name %in% c("Sudeste", "Southeast") ~ "Southeast",
      region_name %in% c("Sul", "South") ~ "South",
      region_name %in% c("Centro Oeste", "Central-West") ~ "Central-West",
      TRUE ~ as.character(region_name)
    ),
    region_name = factor(
      region_name,
      levels = c(
        "Brazil",
        "North",
        "Northeast",
        "Southeast",
        "South",
        "Central-West"
      )
    ),
    very_high_robust = sprintf(
      "%d (%.1f%%)",
      n_very_high_robust,
      pct_very_high_robust
    ),
    very_high_uncertain = sprintf(
      "%d (%.1f%%)",
      n_very_high_uncertain,
      pct_very_high_uncertain
    ),
    high = sprintf(
      "%d (%.1f%%)",
      n_high,
      pct_high
    ),
    moderate = sprintf(
      "%d (%.1f%%)",
      n_moderate,
      pct_moderate
    ),
    potential_risk = sprintf(
      "%d (%.1f%%)",
      n_potential_risk,
      pct_potential_risk
    ),
    low = sprintf(
      "%d (%.1f%%)",
      n_low,
      pct_low
    )
  ) |>
  select(
    region_name,
    state_name,
    total_municipalities,
    mean_sd_score,
    median_iqr_score,
    median_iqr_uncertainty,
    very_high_robust,
    very_high_uncertain,
    high,
    moderate,
    potential_risk,
    low
  ) |>
  arrange(
    region_name,
    state_name
  ) |>
  gt(
    groupname_col = "region_name"
  ) |>
  cols_label(
    state_name = "Federative unit",
    total_municipalities = "Municipalities (n)",
    mean_sd_score = "Mean (SD)",
    median_iqr_score = "Median (IQR)",
    median_iqr_uncertainty = "Median (IQR)",
    very_high_robust = "Very high (robust)",
    very_high_uncertain = "Very high (uncertain)",
    high = "High",
    moderate = "Moderate",
    potential_risk = "Low (potential risk)",
    low = "Low"
  ) |>
  tab_spanner(
    label = "Priority index",
    columns = c(
      mean_sd_score,
      median_iqr_score
    )
  ) |>
  tab_spanner(
    label = "Uncertainty",
    columns = median_iqr_uncertainty
  ) |>
  tab_spanner(
    label = "Priority classification (n, %)",
    columns = c(
      very_high_robust,
      very_high_uncertain,
      high,
      moderate,
      potential_risk,
      low
    )
  ) |>
  cols_align(
    align = "center",
    columns = -state_name
  ) |>
  tab_header(
    title = "Table S17. Summary of the municipal tuberculosis mortality priority index by region and federative unit, Brazil, 2010–2024."
  ) |>
  tab_options(
    table.font.size = gt::px(12),
    data_row.padding = gt::px(3)
  ) |>
  tab_source_note(
    source_note = gt::md(
      "*Note:* The adjusted composite index ranges from 0 to 100 and combines spatial, temporal, and risk-related indicators, with a penalty applied for uncertainty. Priority categories are defined using the adjusted score together with probabilistic criteria and robustness across simulation scenarios.  \
*Abbreviations:* IQR, interquartile range (25th–75th percentiles); SD, standard deviation."
    )
  )

gtsave(
  priority_index_by_state_table,
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S17_priority_index_by_state.docx"
  )
)

# Priority-index map ----

MUNICIPAL_BOUNDARY_YEAR <- 2022L
STATE_BOUNDARY_YEAR <- 2020L

priority_levels <- c(
  "Low priority",
  "Low priority (potential risk)",
  "Moderate priority",
  "High priority",
  "Very high priority (uncertain)",
  "Very high priority (robust)"
)

priority_palette <- c(
  "Low priority" = "#f0f0f0",
  "Low priority (potential risk)" = "#d9f0ea",
  "Moderate priority" = "#80cdc1",
  "High priority" = "#fdae6b",
  "Very high priority (uncertain)" = "#e6550d",
  "Very high priority (robust)" = "#a50f15"
)

priority_map_counts <- final_decision_table |>
  mutate(
    final_priority_class = factor(
      final_priority_class,
      levels = priority_levels
    )
  ) |>
  dplyr::count(
    final_priority_class,
    .drop = FALSE,
    name = "n"
  ) |>
  mutate(
    priority_profile_cat = paste0(
      as.character(final_priority_class),
      " [",
      format(
        n,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      "]"
    )
  )

priority_map_data <- final_decision_table |>
  mutate(
    final_priority_class = factor(
      final_priority_class,
      levels = priority_levels
    )
  ) |>
  left_join(
    priority_map_counts,
    by = "final_priority_class"
  ) |>
  mutate(
    priority_profile_cat = factor(
      priority_profile_cat,
      levels = priority_map_counts$priority_profile_cat
    )
  ) |>
  select(
    code_muni,
    municipality_name,
    state_name,
    region_name,
    final_priority_class,
    priority_profile_cat,
    adjusted_priority_score,
    operational_priority
  )

write_csv(
  priority_map_data,
  file.path(
    RESULTS_DATA_DIR,
    "priority_index_map_data.csv"
  )
)

municipal_boundaries <- geobr::read_municipality(
  code_muni = "all",
  year = MUNICIPAL_BOUNDARY_YEAR,
  simplified = TRUE
) |>
  mutate(
    code_muni = substr(
      as.character(code_muni),
      1,
      6
    )
  )

state_boundaries <- geobr::read_state(
  code_state = "all",
  year = STATE_BOUNDARY_YEAR,
  simplified = TRUE
)

mapa_df <- municipal_boundaries |>
  left_join(
    priority_map_data,
    by = "code_muni"
  )

map_high <- mapa_df |>
  filter(
    final_priority_class == "High priority"
  )

map_very_high <- mapa_df |>
  filter(
    final_priority_class %in% c(
      "Very high priority (uncertain)",
      "Very high priority (robust)"
    )
  )

priority_profile_palette <- setNames(
  priority_palette[as.character(priority_map_counts$final_priority_class)],
  priority_map_counts$priority_profile_cat
)

gg_priority_map <- ggplot() +
  geom_sf(
    data = mapa_df,
    aes(fill = priority_profile_cat),
    color = NA
  ) +
  geom_sf(
    data = map_high,
    fill = NA,
    color = "#333333",
    linewidth = 0.15
  ) +
  geom_sf(
    data = map_very_high,
    fill = NA,
    color = "#000000",
    linewidth = 0.35
  ) +
  geom_sf(
    data = state_boundaries,
    fill = NA,
    color = "black",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    name = paste0(
      "Priority classification \n[",
      format(
        nrow(final_decision_table),
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      " municipalities]"
    ),
    values = priority_profile_palette,
    guide = guide_legend(
      override.aes = list(
        color = "black",
        linewidth = 0.3
      )
    )
  ) +
  coord_sf(datum = NA) +
  theme_void() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.key.size = grid::unit(0.8, "cm")
  ) +
  ggspatial::annotation_scale(
    location = "br",
    text_cex = 1,
    height = grid::unit(0.3, "cm"),
    pad_y = grid::unit(0.9, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    style = ggspatial::north_arrow_fancy_orienteering(),
    location = "br",
    width = grid::unit(1.5, "cm"),
    height = grid::unit(2, "cm"),
    pad_x = grid::unit(1.5, "cm"),
    pad_y = grid::unit(1.5, "cm")
  )

gg_priority_map

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_4_priority_index.tiff"
  ),
  plot = gg_priority_map,
  device = "tiff",
  width = 12,
  height = 8,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_4_priority_index.pdf"
  ),
  plot = gg_priority_map,
  device = "tiff",
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

writeLines(
  paste(
    "Fig 4. Classification of the municipal tuberculosis mortality priority index in Brazil, 2010–2024."
  ),
  file.path(
    RESULTS_FIGURE_DIR,
    "Fig_4_caption.txt"
  )
)

# Supplementary municipal table ----

MSTCAR_OVERALL_FILE <- file.path(
  "results",
  "mstcar",
  "combined_age_standardised_rates_2010_2024.csv"
)

mstcar_overall <- readr::read_csv(
  MSTCAR_OVERALL_FILE
)

required_mstcar_columns <- c(
    "code_muni",
    "posterior_median",
    "ci_lower_95",
    "ci_upper_95",
    "events",
    "average_annual_population"
  )

missing_mstcar_columns <- setdiff(
    required_mstcar_columns,
    names(mstcar_overall)
  )

if (length(missing_mstcar_columns) > 0L) {
    stop(
      "Columns missing from the MSTCAR overall-rate file: ",
      paste(
        missing_mstcar_columns,
        collapse = ", "
      )
    )
  }

mstcar_overall <- mstcar_overall |>
    mutate(
      code_muni = sprintf(
        "%06d",
        as.integer(code_muni)
      ),
      cmr =
        1e5 *
        events /
        population_person_years
    ) |>
    select(
      code_muni,
      observed_deaths = events,
      population = average_annual_population,
      cmr,
      asmr_estimated = posterior_median,
      asmr_p025 = ci_lower_95,
      asmr_p975 = ci_upper_95
    )

supplementary_table_s1 <- priority_index |>
    left_join(
      mstcar_overall,
      by = "code_muni"
    ) |>
    transmute(
      municipality_code=code_muni,
      municipality_name,
      state_name,
      region_name,

      observed_deaths,
      population,
      cmr,
      asmr_estimated,
      asmr_p025,
      asmr_p975,

      relative_risk = RR_after,
      relative_risk_p025 = RR_lower_95,
      relative_risk_p975 = RR_upper_95,
      Pr_RR_gt1 = Pr_RR_after,

      temporal_trend = slope_after,
      temporal_trend_p025 =
        slope_lower_95,
      temporal_trend_p975 =
        slope_upper_95,
      Pr_positive_trend = Pr_slope_after,

      regime_before,
      regime_after,
      cluster_before,
      cluster_after,
      trend_before,
      trend_after,
      hotspot_status,
      coldspot_status,
      cluster_transition,

      strict_regime_persistence =
        strict_persistent,
      spatial_classification_persistence =
        cluster_persistent,
      temporal_trend_persistence =
        trend_persistent,

      regime_entropy =
        regime_entropy_normalized,
      state_persistence =
        state_strict_persistence,
      state_cluster_entropy =
        state_cluster_entropy_normalized,

      spatial_score,
      persistence_score,
      rr_score,
      rr_probability_score,
      trend_score,
      trend_probability_score,
      regime_score,
      state_context_score,

      central_priority_score,
      adjusted_priority_score,
      uncertainty_score,
      confidence_class,

      simulated_score_median,
      simulated_score_p025=simulated_score_lower_95,
      simulated_score_p975=simulated_score_upper_95,

      median_rank,
      rank_p025=rank_lower_95,
      rank_p975=rank_upper_95,
      Pr_top10,
      Pr_top50,
      Pr_top100,

      Pr_score_very_high,
      Pr_score_high_or_higher,
      Pr_top50_across_scenarios,
      Pr_very_high_across_scenarios,
      ranking_robustness,
      scenario_robustness_class,

      priority_profile,
      final_priority_class,
      operational_priority
    ) |>
    mutate(
      population = round(population, 0),
      across(
        c(
          strict_regime_persistence,
          spatial_classification_persistence,
          temporal_trend_persistence
        ),
        ~ case_when(
          is.na(.x) ~ NA_character_,
          .x ~ "Yes",
          TRUE ~ "No"
        )
      ),

      across(
        c(
          cmr,
          asmr_estimated,
          asmr_p025,
          asmr_p975,
          relative_risk,
          relative_risk_p025,
          relative_risk_p975
        ),
        ~ round(
          .x,
          2
        )
      ),

      across(
        c(
          temporal_trend,
          temporal_trend_p025,
          temporal_trend_p975
        ),
        ~ round(
          .x,
          3
        )
      ),

      across(
        c(
          Pr_RR_gt1,
          Pr_positive_trend,
          Pr_top10,
          Pr_top50,
          Pr_top100,
          Pr_score_very_high,
          Pr_score_high_or_higher,
          Pr_top50_across_scenarios,
          Pr_very_high_across_scenarios
        ),
        ~ round(
          .x,
          3
        )
      ),

      across(
        c(
          regime_entropy,
          state_persistence,
          state_cluster_entropy,
          central_priority_score,
          adjusted_priority_score,
          uncertainty_score,
          simulated_score_median,
          simulated_score_p025,
          simulated_score_p975,
          median_rank,
          rank_p025,
          rank_p975
        ),
        ~ round(
          .x,
          2
        )
      )
    ) |>
    arrange(
      desc(
        adjusted_priority_score
      ),
      desc(
        Pr_top50
      )
    )


codebook <- tibble(
  variable = c(
    "municipality_code","municipality_name","state_name","region_name",
    "observed_deaths","population","cmr",
    "asmr_estimated","asmr_p025","asmr_p975",
    "relative_risk","relative_risk_p025","relative_risk_p975","Pr_RR_gt1",
    "temporal_trend","temporal_trend_p025","temporal_trend_p975","Pr_positive_trend",
    "regime_before","regime_after",
    "cluster_before","cluster_after",
    "trend_before","trend_after",
    "hotspot_status","coldspot_status","cluster_transition",
    "strict_regime_persistence",
    "spatial_classification_persistence",
    "temporal_trend_persistence",
    "regime_entropy","state_persistence","state_cluster_entropy",
    "spatial_score","persistence_score",
    "rr_score","rr_probability_score",
    "trend_score","trend_probability_score",
    "regime_score","state_context_score",
    "central_priority_score","adjusted_priority_score","uncertainty_score",
    "confidence_class",
    "simulated_score_median","simulated_score_p025","simulated_score_p975",
    "median_rank","rank_p025","rank_p975",
    "Pr_top10","Pr_top50","Pr_top100",
    "Pr_score_very_high","Pr_score_high_or_higher",
    "Pr_top50_across_scenarios","Pr_very_high_across_scenarios",
    "ranking_robustness","scenario_robustness_class",
    "priority_profile","final_priority_class","operational_priority"
  ),
  
  description = c(
    "Municipality code (6 digits)",
    "Municipality name",
    "State name",
    "Region of Brazil",
    
    "Observed number of deaths (2010–2024)",
    "Average annual population",
    "Crude mortality rate per 100,000 population",
    
    "Estimated age-standardised mortality rate (posterior median)",
    "Lower bound of 95% credible interval (ASMR)",
    "Upper bound of 95% credible interval (ASMR)",
    
    "Posterior median relative risk",
    "Lower bound of 95% credible interval (RR)",
    "Upper bound of 95% credible interval (RR)",
    "Posterior probability that RR > 1",
    
    "Estimated temporal trend (slope)",
    "Lower bound of 95% credible interval (trend)",
    "Upper bound of 95% credible interval (trend)",
    "Posterior probability of positive trend",
    
    "Spatiotemporal regime before adjustment",
    "Spatiotemporal regime after adjustment",
    
    "Spatial cluster classification before adjustment",
    "Spatial cluster classification after adjustment",
    
    "Temporal trend classification before adjustment",
    "Temporal trend classification after adjustment",
    
    "Hotspot transition classification",
    "Coldspot transition classification",
    "Overall cluster transition category",
    
    "Persistence of regime classification (Yes/No)",
    "Persistence of spatial classification (Yes/No)",
    "Persistence of temporal trend classification (Yes/No)",
    
    "Normalized entropy of regime classification",
    "State-level persistence measure",
    "State-level entropy of spatial clustering",
    
    "Spatial component score",
    "Persistence component score",
    
    "Relative risk score",
    "Probability-based RR score",
    
    "Temporal trend score",
    "Probability-based trend score",
    
    "Regime classification score",
    "State context score",
    
    "Combined central priority score",
    "Adjusted priority score (final ranking metric)",
    "Uncertainty score",
    
    "Confidence classification",
    
    "Median simulated priority score",
    "Lower bound of simulated score (95% CI)",
    "Upper bound of simulated score (95% CI)",
    
    "Median rank position",
    "Lower bound of rank (95% CI)",
    "Upper bound of rank (95% CI)",
    
    "Probability of being in top 10",
    "Probability of being in top 50",
    "Probability of being in top 100",
    
    "Probability of very high priority score",
    "Probability of high or higher priority score",
    
    "Probability of top 50 across scenarios",
    "Probability of very high priority across scenarios",
    
    "Ranking robustness metric",
    "Robustness classification across scenarios",
    
    "Priority profile classification",
    "Final priority class",
    "Operational priority classification"
  )
)


readr::write_csv(
  list(
    Data = supplementary_table_s1,
    Codebook = codebook
  ),
    file.path(
      RESULTS_TABLE_DIR,
      "Supplementary_Table_S1.csv"
    )
  )


writexl::write_xlsx(
  list(
    Data = supplementary_table_s1,
    Codebook = codebook
  ),
  file.path(
    RESULTS_TABLE_DIR,
    "Supplementary_Table_S1.xlsx"
  )
)

# Reproducibility information ----

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    RESULTS_DATA_DIR,
    "session_info_priority_index.txt"
  )
)
