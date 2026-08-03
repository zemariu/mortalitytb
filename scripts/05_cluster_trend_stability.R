# Article: Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024:
#          a Bayesian spatiotemporal modelling study
# Repository: mortalitytb
# Script: Compares classifications before and after covariate adjustment and evaluates 
#         cluster and trend stability.
# Author: José Mário Nunes da Silva, PhD
# Last updated: 2026-08-02

# Packages ----
library(tidyverse)
library(gt)
library(patchwork)

# Output directories ----

RESULTS_DATA_DIR <- "results/spatiotemporal/outputs"
RESULTS_TABLE_DIR <- "results/spatiotemporal/tables"
RESULTS_FIGURE_DIR <- "results/spatiotemporal/figures"

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

# Import municipal classifications ----

BEFORE_CLASSIFICATION_FILE <- file.path(
  RESULTS_DATA_DIR,
  "spatial_temporal_classification_before_adjustment.csv"
)

AFTER_CLASSIFICATION_FILE <- file.path(
  RESULTS_DATA_DIR,
  "spatial_temporal_classification_after_adjustment.csv"
)

if (!file.exists(BEFORE_CLASSIFICATION_FILE)) {
  stop(
    "The unadjusted municipal classification file was not found: ",
    BEFORE_CLASSIFICATION_FILE
  )
}

if (!file.exists(AFTER_CLASSIFICATION_FILE)) {
  stop(
    "The covariate-adjusted municipal classification file was not found: ",
    AFTER_CLASSIFICATION_FILE
  )
}

classification_before <- readr::read_csv(
  BEFORE_CLASSIFICATION_FILE,
  show_col_types = FALSE
)

classification_after <- readr::read_csv(
  AFTER_CLASSIFICATION_FILE,
  show_col_types = FALSE
)

required_classification_columns <- c(
  "code_muni",
  "name_muni",
  "name_state",
  "name_region",
  "profile",
  "cluster",
  "trend",
  "RR_mean",
  "RR_lower_95",
  "RR_upper_95",
  "Pr_RR_gt1",
  "slope_mean",
  "slope_lower_95",
  "slope_upper_95",
  "Pr_slope_gt0"
)

missing_before <- setdiff(
  required_classification_columns,
  names(classification_before)
)

missing_after <- setdiff(
  required_classification_columns,
  names(classification_after)
)

if (length(missing_before) > 0L) {
  stop(
    "Columns missing from the unadjusted classification file: ",
    paste(missing_before, collapse = ", ")
  )
}

if (length(missing_after) > 0L) {
  stop(
    "Columns missing from the adjusted classification file: ",
    paste(missing_after, collapse = ", ")
  )
}

classification_before <- classification_before |>
  mutate(
    code_muni = sprintf(
      "%06d",
      as.integer(code_muni)
    )
  )

classification_after <- classification_after |>
  mutate(
    code_muni = sprintf(
      "%06d",
      as.integer(code_muni)
    )
  )

if (anyDuplicated(classification_before$code_muni)) {
  stop(
    "The unadjusted classification contains duplicated municipality codes."
  )
}

if (anyDuplicated(classification_after$code_muni)) {
  stop(
    "The adjusted classification contains duplicated municipality codes."
  )
}

municipalities_missing_after <- setdiff(
  classification_before$code_muni,
  classification_after$code_muni
)

municipalities_missing_before <- setdiff(
  classification_after$code_muni,
  classification_before$code_muni
)

if (
  length(municipalities_missing_after) > 0L ||
  length(municipalities_missing_before) > 0L
) {
  stop(
    "The unadjusted and adjusted classifications do not contain the same municipalities."
  )
}

# Comparative municipal database ----

before_data <- classification_before |>
  transmute(
    code_muni,
    municipality_name = as.character(name_muni),
    state_name = as.character(name_state),
    region_name = recode(
      as.character(name_region),
      "Norte" = "North",
      "Nordeste" = "Northeast",
      "Sudeste" = "Southeast",
      "Sul" = "South",
      "Centro-Oeste" = "Central-West",
      "Centro Oeste" = "Central-West",
      .default = as.character(name_region)
    ),
    regime_before = as.character(profile),
    cluster_before = as.character(cluster),
    trend_before = as.character(trend),
    RR_before = as.numeric(RR_mean),
    Pr_RR_before = as.numeric(Pr_RR_gt1),
    slope_before = as.numeric(slope_mean),
    Pr_slope_before = as.numeric(Pr_slope_gt0)
  )

after_data <- classification_after |>
  transmute(
    code_muni,
    municipality_name_after = as.character(name_muni),
    state_name_after = as.character(name_state),
    region_name_after = recode(
      as.character(name_region),
      "Norte" = "North",
      "Nordeste" = "Northeast",
      "Sudeste" = "Southeast",
      "Sul" = "South",
      "Centro-Oeste" = "Central-West",
      "Centro Oeste" = "Central-West",
      .default = as.character(name_region)
    ),
    regime_after = as.character(profile),
    cluster_after = as.character(cluster),
    trend_after = as.character(trend),
    RR_after = as.numeric(RR_mean),
    RR_lower_95 = as.numeric(RR_lower_95),
    RR_upper_95 = as.numeric(RR_upper_95),
    Pr_RR_after = as.numeric(Pr_RR_gt1),
    slope_after = as.numeric(slope_mean),
    slope_lower_95 = as.numeric(slope_lower_95),
    slope_upper_95 = as.numeric(slope_upper_95),
    Pr_slope_after = as.numeric(Pr_slope_gt0)
  )

comp <- before_data |>
  inner_join(
    after_data,
    by = "code_muni"
  )

if (
  any(comp$municipality_name != comp$municipality_name_after) ||
  any(comp$state_name != comp$state_name_after) ||
  any(comp$region_name != comp$region_name_after)
) {
  stop(
    "Municipality, state, or region names differ between the two classification files."
  )
}

comp <- comp |>
  select(
    -municipality_name_after,
    -state_name_after,
    -region_name_after
  ) |>
  arrange(code_muni)

if (nrow(comp) != nrow(classification_before)) {
  stop(
    "The comparative database does not contain all municipalities."
  )
}

readr::write_csv(
  comp,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_classification_comparison.csv"
  )
)

saveRDS(
  comp,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_classification_comparison.rds"
  )
)

message(
  "Comparative database created with ",
  format(nrow(comp), big.mark = ","),
  " municipalities."
)

# Classification orders ----

REGION_ORDER <- c(
  "North",
  "Northeast",
  "Southeast",
  "South",
  "Central-West"
)

CLUSTER_ORDER <- c(
  "Hotspot",
  "Coldspot",
  "Neutralspot"
)

TREND_ORDER <- c(
  "Increasing",
  "Decreasing",
  "Stable"
)

REGIME_ORDER <- c(
  "Coldspot / Decreasing",
  "Coldspot / Stable",
  "Coldspot / Increasing",
  "Neutralspot / Decreasing",
  "Neutralspot / Stable",
  "Neutralspot / Increasing",
  "Hotspot / Decreasing",
  "Hotspot / Stable",
  "Hotspot / Increasing"
)

# Complete-regime transition matrix ----

regime_transition <- table(
  factor(
    comp$regime_before,
    levels = REGIME_ORDER
  ),
  factor(
    comp$regime_after,
    levels = REGIME_ORDER
  ),
  dnn = c(
    "Before adjustment",
    "After adjustment"
  )
)

regime_transition_long <- as.data.frame(
  regime_transition,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  rename(
    From = `Before.adjustment`,
    To = `After.adjustment`,
    Count = Freq
  ) |>
  mutate(
    From = factor(
      From,
      levels = REGIME_ORDER
    ),
    To = factor(
      To,
      levels = REGIME_ORDER
    ),
    diagonal = From == To
  )

readr::write_csv(
  regime_transition_long,
  file.path(
    RESULTS_DATA_DIR,
    "complete_regime_transition_matrix.csv"
  )
)

# Figure S9: complete-regime transition matrix ----

figure_s9 <- ggplot(
  regime_transition_long,
  aes(
    x = To,
    y = From
  )
) +
  geom_tile(
    aes(fill = Count),
    colour = "grey90",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  geom_tile(
    data = regime_transition_long |>
      filter(diagonal),
    fill = "#8DD3C7",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = Count),
    size = 3.6
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#D9D9D9"
  ) +
  scale_x_discrete(
    labels = \(x) stringr::str_wrap(x, width = 16)
  ) +
  scale_y_discrete(
    labels = \(x) stringr::str_wrap(x, width = 16)
  ) +
  labs(
    x = "After adjustment",
    y = "Before adjustment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title = element_text(
      colour = "black"
    )
  )

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S9.tiff"
  ),
  plot = figure_s9,
  width = 12,
  height = 9,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig S9. Cross-model classification matrix comparing municipal spatiotemporal",
    "regimes before and after covariate adjustment, Brazil, 2010–2024.",
    "Diagonal cells indicate classifications retained after adjustment, whereas",
    "off-diagonal cells indicate redistribution across model specifications."
  ),
  file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S9_caption.txt"
  )
)

# Spatial-cluster and temporal-trend transition matrices ----

cluster_transition <- table(
  factor(
    comp$cluster_before,
    levels = CLUSTER_ORDER
  ),
  factor(
    comp$cluster_after,
    levels = CLUSTER_ORDER
  ),
  dnn = c(
    "Before adjustment",
    "After adjustment"
  )
)

cluster_transition_long <- as.data.frame(
  cluster_transition,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  rename(
    From = `Before.adjustment`,
    To = `After.adjustment`,
    Count = Freq
  ) |>
  mutate(
    From = factor(
      From,
      levels = CLUSTER_ORDER
    ),
    To = factor(
      To,
      levels = CLUSTER_ORDER
    ),
    diagonal = From == To
  )

trend_transition <- table(
  factor(
    comp$trend_before,
    levels = TREND_ORDER
  ),
  factor(
    comp$trend_after,
    levels = TREND_ORDER
  ),
  dnn = c(
    "Before adjustment",
    "After adjustment"
  )
)

trend_transition_long <- as.data.frame(
  trend_transition,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  rename(
    From = `Before.adjustment`,
    To = `After.adjustment`,
    Count = Freq
  ) |>
  mutate(
    From = factor(
      From,
      levels = TREND_ORDER
    ),
    To = factor(
      To,
      levels = TREND_ORDER
    ),
    diagonal = From == To
  )

readr::write_csv(
  cluster_transition_long,
  file.path(
    RESULTS_DATA_DIR,
    "spatial_cluster_transition_matrix.csv"
  )
)

readr::write_csv(
  trend_transition_long,
  file.path(
    RESULTS_DATA_DIR,
    "differential_trend_transition_matrix.csv"
  )
)

# Figure S10: cluster and trend transition matrices ----

cluster_panel <- ggplot(
  cluster_transition_long,
  aes(
    x = To,
    y = From
  )
) +
  geom_tile(
    aes(fill = Count),
    colour = "grey90",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  geom_tile(
    data = cluster_transition_long |>
      filter(diagonal),
    fill = "#FB8072",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = Count),
    size = 4
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#D9D9D9"
  ) +
  labs(
    title = "A) Spatial cluster",
    x = NULL,
    y = "Before adjustment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      colour = "black"
    )
  )

trend_panel <- ggplot(
  trend_transition_long,
  aes(
    x = To,
    y = From
  )
) +
  geom_tile(
    aes(fill = Count),
    colour = "grey90",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  geom_tile(
    data = trend_transition_long |>
      filter(diagonal),
    fill = "#BEBADA",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = Count),
    size = 4
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#D9D9D9"
  ) +
  labs(
    title = "B) Differential trend",
    x = "After adjustment",
    y = "Before adjustment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      colour = "black"
    )
  )

figure_s10 <- cluster_panel / trend_panel &
  theme(
    plot.title = element_text(hjust = -0.1)
  )

ggsave(
  filename = file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S10.tiff"
  ),
  plot = figure_s10,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

writeLines(
  paste(
    "Fig S10. Cross-model classification matrices comparing municipal spatial-cluster",
    "classifications (A) and differential temporal-trend classifications (B) before and",
    "after covariate adjustment, Brazil, 2010–2024."
  ),
  file.path(
    RESULTS_FIGURE_DIR,
    "Fig_S10_caption.txt"
  )
)

# Global stability and entropy metrics ----

shannon_entropy <- function(tab) {
  probabilities <- as.numeric(
    prop.table(tab)
  )

  probabilities <- probabilities[
    is.finite(probabilities) &
      probabilities > 0
  ]

  -sum(
    probabilities * log(probabilities)
  )
}

row_entropy <- function(tab) {
  row_probabilities <- prop.table(
    tab,
    margin = 1
  )

  entropy_i <- apply(
    row_probabilities,
    1,
    function(probabilities) {
      probabilities <- probabilities[
        is.finite(probabilities) &
          probabilities > 0
      ]

      if (length(probabilities) == 0L) {
        return(NA_real_)
      }

      -sum(
        probabilities * log(probabilities)
      )
    }
  )

  as.numeric(entropy_i)
}

strict_persistence_global <-
  sum(diag(regime_transition)) /
  sum(regime_transition)

global_entropy <- shannon_entropy(
  regime_transition
)

global_entropy_normalized <-
  global_entropy /
  log(length(regime_transition))

regime_entropy_values <- row_entropy(
  regime_transition
)

regime_frequencies <- rowSums(
  regime_transition
)

regime_weights <-
  regime_frequencies /
  sum(regime_frequencies)

srce_global <- sum(
  regime_weights * regime_entropy_values,
  na.rm = TRUE
)

srce_global_normalized <-
  srce_global /
  log(ncol(regime_transition))

global_metrics <- tibble(
  Number_of_municipalities = nrow(comp),
  Classification_persistence_index = strict_persistence_global,
  Global_entropy = global_entropy,
  Normalized_global_entropy = global_entropy_normalized,
  Spatial_regime_classification_entropy = srce_global,
  Normalized_spatial_regime_classification_entropy =
    srce_global_normalized
)

writexl::write_xlsx(
  global_metrics,
  file.path(
    RESULTS_DATA_DIR,
    "global_classification_stability_metrics.xlsx"
  )
)

# Table S14: regime-specific entropy ----

table_s14_data <- tibble(
  regime = rownames(regime_transition),
  Number_of_municipalities = as.integer(regime_frequencies),
  Entropy = regime_entropy_values,
  Normalized_entropy =
    regime_entropy_values /
    log(ncol(regime_transition))
) |>
  tidyr::separate_wider_delim(
    regime,
    delim = " / ",
    names = c(
      "Spatial_cluster",
      "Temporal_trend"
    )
  ) |>
  mutate(
    Spatial_cluster = factor(
      Spatial_cluster,
      levels = CLUSTER_ORDER
    ),
    Temporal_trend = factor(
      Temporal_trend,
      levels = c(
        "Decreasing",
        "Increasing",
        "Stable"
      )
    )
  ) |>
  arrange(
    Spatial_cluster,
    Temporal_trend
  )

writexl::write_xlsx(
  table_s14_data |>
    mutate(
      Spatial_cluster = as.character(Spatial_cluster),
      Temporal_trend = as.character(Temporal_trend)
    ),
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S14_regime_specific_entropy.xlsx"
  )
)

table_s14 <- table_s14_data |>
  gt() |>
  cols_label(
    Spatial_cluster = "Spatial cluster",
    Temporal_trend = "Temporal trend",
    Number_of_municipalities = "Number of municipalities",
    Entropy = "Entropy",
    Normalized_entropy = "Normalized entropy"
  ) |>
  fmt_integer(
    columns = Number_of_municipalities,
    use_seps = TRUE
  ) |>
  fmt_number(
    columns = c(
      Entropy,
      Normalized_entropy
    ),
    decimals = 3
  ) |>
  cols_align(
    align = "center",
    columns = -c(
      Spatial_cluster,
      Temporal_trend
    )
  ) |>
  tab_header(
    title = md(
      "**Table S14.** Regime-specific entropy of spatiotemporal classifications for tuberculosis mortality, Brazil, 2010–2024."
    )
  ) |>
  tab_source_note(
    source_note = md(
      paste0(
        "*Note:* Entropy measures the heterogeneity of classification changes between ",
        "the unadjusted and covariate-adjusted models within each spatiotemporal regime. ",
        "Higher values indicate greater sensitivity to model specification, whereas lower ",
        "values indicate greater cross-model classification stability. Normalized entropy ",
        "ranges from 0 to 1, allowing comparison across regimes. Normalized global entropy = ",
        sprintf("%.3f", global_entropy_normalized),
        ". The normalized global Spatial Regime Classification Entropy (SRCE) was ",
        sprintf("%.3f", srce_global_normalized),
        "."
      )
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(3)
  )

gt::gtsave(
  table_s14,
  filename = file.path(
    RESULTS_TABLE_DIR,
    "Table_S14_regime_specific_entropy.docx"
  )
)

# Transition and persistence indicators ----

comp <- comp |>
  mutate(
    strict_persistent =
      regime_before == regime_after,
    cluster_persistent =
      cluster_before == cluster_after,
    trend_persistent =
      trend_before == trend_after,

    hotspot_status = case_when(
      cluster_before == "Hotspot" &
        cluster_after == "Hotspot" ~
        "Persistent hotspot",

      cluster_before == "Hotspot" &
        cluster_after != "Hotspot" ~
        "Attenuated hotspot",

      cluster_before != "Hotspot" &
        cluster_after == "Hotspot" ~
        "Newly identified hotspot",

      TRUE ~
        "No hotspot event"
    ),

    coldspot_status = case_when(
      cluster_before == "Coldspot" &
        cluster_after == "Coldspot" ~
        "Persistent coldspot",

      cluster_before == "Coldspot" &
        cluster_after != "Coldspot" ~
        "Attenuated coldspot",

      cluster_before != "Coldspot" &
        cluster_after == "Coldspot" ~
        "Newly identified coldspot",

      TRUE ~
        "No coldspot event"
    ),

    cluster_transition = paste(
      cluster_before,
      cluster_after,
      sep = " -> "
    )
  )

persistence_summary <- comp |>
  summarise(
    Number_of_municipalities = n(),
    Strict_regime_persistence_n = sum(
      strict_persistent,
      na.rm = TRUE
    ),
    Strict_regime_persistence = mean(
      strict_persistent,
      na.rm = TRUE
    ),
    Spatial_classification_persistence_n = sum(
      cluster_persistent,
      na.rm = TRUE
    ),
    Spatial_classification_persistence = mean(
      cluster_persistent,
      na.rm = TRUE
    ),
    Temporal_trend_persistence_n = sum(
      trend_persistent,
      na.rm = TRUE
    ),
    Temporal_trend_persistence = mean(
      trend_persistent,
      na.rm = TRUE
    )
  )

writexl::write_xlsx(
  persistence_summary,
  file.path(
    RESULTS_DATA_DIR,
    "national_classification_persistence_summary.xlsx"
  )
)

# Table S15: persistence by federative unit ----

persistence_by_state <- comp |>
  group_by(
    region_name,
    state_name
  ) |>
  summarise(
    Number_of_municipalities = n(),
    Strict_regime_persistence_n = sum(
      strict_persistent,
      na.rm = TRUE
    ),
    Strict_regime_persistence = mean(
      strict_persistent,
      na.rm = TRUE
    ),
    Spatial_classification_persistence_n = sum(
      cluster_persistent,
      na.rm = TRUE
    ),
    Spatial_classification_persistence = mean(
      cluster_persistent,
      na.rm = TRUE
    ),
    Temporal_trend_persistence_n = sum(
      trend_persistent,
      na.rm = TRUE
    ),
    Temporal_trend_persistence = mean(
      trend_persistent,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    region_name = factor(
      region_name,
      levels = REGION_ORDER
    )
  ) |>
  arrange(
    region_name,
    state_name
  )

writexl::write_xlsx(
  persistence_by_state |>
    mutate(
      region_name = as.character(region_name)
    ),
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S15_classification_stability_by_state.xlsx"
  )
)

table_s15 <- persistence_by_state |>
  gt(
    groupname_col = "region_name"
  ) |>
  cols_label(
    state_name = "Region/State",
    Number_of_municipalities = "Total",
    Strict_regime_persistence_n = "n",
    Strict_regime_persistence = "%",
    Spatial_classification_persistence_n = "n",
    Spatial_classification_persistence = "%",
    Temporal_trend_persistence_n = "n",
    Temporal_trend_persistence = "%"
  ) |>
  tab_spanner(
    label = "Strict regime persistence",
    columns = c(
      Strict_regime_persistence_n,
      Strict_regime_persistence
    )
  ) |>
  tab_spanner(
    label = "Spatial classification persistence",
    columns = c(
      Spatial_classification_persistence_n,
      Spatial_classification_persistence
    )
  ) |>
  tab_spanner(
    label = "Temporal trend persistence",
    columns = c(
      Temporal_trend_persistence_n,
      Temporal_trend_persistence
    )
  ) |>
  fmt_integer(
    columns = c(
      Number_of_municipalities,
      Strict_regime_persistence_n,
      Spatial_classification_persistence_n,
      Temporal_trend_persistence_n
    ),
    use_seps = TRUE
  ) |>
  fmt_percent(
    columns = c(
      Strict_regime_persistence,
      Spatial_classification_persistence,
      Temporal_trend_persistence
    ),
    decimals = 1
  ) |>
  cols_align(
    align = "center",
    columns = -state_name
  ) |>
  tab_header(
    title = md(
      "**Table S15.** Stability of spatiotemporal classifications across unadjusted and covariate-adjusted models by region and federative units in Brazil, 2010–2024."
    )
  ) |>
  tab_source_note(
    source_note = md(
      paste(
        "*Note:* Strict regime persistence refers to municipalities assigned to the same",
        "combination of spatial and temporal categories in both the unadjusted and",
        "covariate-adjusted models. Spatial classification persistence refers to retention",
        "of hotspot, coldspot, or neutralspot status. Temporal trend persistence refers to",
        "retention of increasing, decreasing, or stable trend classification."
      )
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) |>
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(3)
  )

gt::gtsave(
  table_s15,
  filename = file.path(
    RESULTS_TABLE_DIR,
    "Table_S15_classification_stability_by_state.docx"
  )
)

# State-level cluster-transition entropy ----

entropy_by_state <- comp |>
  group_by(
    region_name,
    state_name
  ) |>
  summarise(
    Number_of_municipalities = n(),
    Cluster_transition_entropy = {
      transition_i <- table(
        factor(
          cluster_before,
          levels = CLUSTER_ORDER
        ),
        factor(
          cluster_after,
          levels = CLUSTER_ORDER
        )
      )

      shannon_entropy(
        transition_i
      )
    },
    Normalized_cluster_transition_entropy = {
      transition_i <- table(
        factor(
          cluster_before,
          levels = CLUSTER_ORDER
        ),
        factor(
          cluster_after,
          levels = CLUSTER_ORDER
        )
      )

      shannon_entropy(
        transition_i
      ) / log(length(transition_i))
    },
    .groups = "drop"
  )

writexl::write_xlsx(
  entropy_by_state,
  file.path(
    RESULTS_DATA_DIR,
    "cluster_transition_entropy_by_state.xlsx"
  )
)

# Add stability context to the comparative database ----

regime_entropy_lookup <- table_s14_data |>
  mutate(
    regime_before = paste(
      Spatial_cluster,
      Temporal_trend,
      sep = " / "
    )
  ) |>
  transmute(
    regime_before,
    regime_n = Number_of_municipalities,
    regime_entropy = Entropy,
    regime_entropy_normalized = Normalized_entropy
  )

state_context_lookup <- persistence_by_state |>
  select(
    region_name,
    state_name,
    state_n = Number_of_municipalities,
    state_strict_persistence =
      Strict_regime_persistence,
    state_spatial_persistence =
      Spatial_classification_persistence,
    state_trend_persistence =
      Temporal_trend_persistence
  ) |>
  left_join(
    entropy_by_state |>
      select(
        region_name,
        state_name,
        state_cluster_entropy =
          Cluster_transition_entropy,
        state_cluster_entropy_normalized =
          Normalized_cluster_transition_entropy
      ),
    by = c(
      "region_name",
      "state_name"
    )
  )

comp <- comp |>
  left_join(
    regime_entropy_lookup,
    by = "regime_before"
  ) |>
  left_join(
    state_context_lookup,
    by = c(
      "region_name",
      "state_name"
    )
  ) |>
  mutate(
    global_strict_persistence =
      strict_persistence_global,
    normalized_global_entropy =
      global_entropy_normalized,
    normalized_global_srce =
      srce_global_normalized
  )

# Table S16: changes in hotspot and coldspot classification ----

format_n_over_N <- function(numerator, denominator) {
  ifelse(
    denominator > 0,
    paste0(
      numerator,
      "/",
      denominator,
      " (",
      sprintf(
        "%.1f",
        100 * numerator / denominator
      ),
      "%)"
    ),
    paste0(
      numerator,
      "/",
      denominator,
      " (–)"
    )
  )
}

hotspot_brazil <- comp |>
  summarise(
    region_name = "Brazil",
    state_name = "Brazil",
    resolved = sum(
      hotspot_status == "Attenuated hotspot",
      na.rm = TRUE
    ),
    persistent = sum(
      hotspot_status == "Persistent hotspot",
      na.rm = TRUE
    ),
    gained = sum(
      hotspot_status == "Newly identified hotspot",
      na.rm = TRUE
    ),
    baseline_denominator = sum(
      cluster_before == "Hotspot",
      na.rm = TRUE
    ),
    gained_denominator = sum(
      cluster_before != "Hotspot",
      na.rm = TRUE
    )
  )

hotspot_states <- comp |>
  group_by(
    region_name,
    state_name
  ) |>
  summarise(
    resolved = sum(
      hotspot_status == "Attenuated hotspot",
      na.rm = TRUE
    ),
    persistent = sum(
      hotspot_status == "Persistent hotspot",
      na.rm = TRUE
    ),
    gained = sum(
      hotspot_status == "Newly identified hotspot",
      na.rm = TRUE
    ),
    baseline_denominator = sum(
      cluster_before == "Hotspot",
      na.rm = TRUE
    ),
    gained_denominator = sum(
      cluster_before != "Hotspot",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

hotspot_summary <- bind_rows(
  hotspot_brazil,
  hotspot_states
) |>
  transmute(
    region_name,
    state_name,
    Hotspot_resolved = format_n_over_N(
      resolved,
      baseline_denominator
    ),
    Hotspot_persistent = format_n_over_N(
      persistent,
      baseline_denominator
    ),
    Hotspot_gained = format_n_over_N(
      gained,
      gained_denominator
    )
  )

coldspot_brazil <- comp |>
  summarise(
    region_name = "Brazil",
    state_name = "Brazil",
    resolved = sum(
      coldspot_status == "Attenuated coldspot",
      na.rm = TRUE
    ),
    persistent = sum(
      coldspot_status == "Persistent coldspot",
      na.rm = TRUE
    ),
    gained = sum(
      coldspot_status == "Newly identified coldspot",
      na.rm = TRUE
    ),
    baseline_denominator = sum(
      cluster_before == "Coldspot",
      na.rm = TRUE
    ),
    gained_denominator = sum(
      cluster_before != "Coldspot",
      na.rm = TRUE
    )
  )

coldspot_states <- comp |>
  group_by(
    region_name,
    state_name
  ) |>
  summarise(
    resolved = sum(
      coldspot_status == "Attenuated coldspot",
      na.rm = TRUE
    ),
    persistent = sum(
      coldspot_status == "Persistent coldspot",
      na.rm = TRUE
    ),
    gained = sum(
      coldspot_status == "Newly identified coldspot",
      na.rm = TRUE
    ),
    baseline_denominator = sum(
      cluster_before == "Coldspot",
      na.rm = TRUE
    ),
    gained_denominator = sum(
      cluster_before != "Coldspot",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

coldspot_summary <- bind_rows(
  coldspot_brazil,
  coldspot_states
) |>
  transmute(
    region_name,
    state_name,
    Coldspot_resolved = format_n_over_N(
      resolved,
      baseline_denominator
    ),
    Coldspot_persistent = format_n_over_N(
      persistent,
      baseline_denominator
    ),
    Coldspot_gained = format_n_over_N(
      gained,
      gained_denominator
    )
  )

table_s16_data <- hotspot_summary |>
  left_join(
    coldspot_summary,
    by = c(
      "region_name",
      "state_name"
    )
  ) |>
  mutate(
    region_name = factor(
      region_name,
      levels = c(
        "Brazil",
        REGION_ORDER
      )
    ),
    state_order = if_else(
      state_name == "Brazil",
      0L,
      1L
    )
  ) |>
  arrange(
    region_name,
    state_order,
    state_name
  ) |>
  select(
    -state_order
  )

writexl::write_xlsx(
  table_s16_data |>
    mutate(
      region_name = as.character(region_name)
    ),
  file.path(
    RESULTS_TABLE_DIR,
    "Table_S16_hotspot_coldspot_changes.xlsx"
  )
)

table_s16 <- table_s16_data |>
  gt(
    groupname_col = "region_name"
  ) |>
  cols_label(
    state_name = "Region/State",
    Hotspot_resolved = "Attenuated\nn/N (%)",
    Hotspot_persistent = "Persistent\nn/N (%)",
    Hotspot_gained = "Newly\nn/N (%)",
    Coldspot_resolved = "Attenuated\nn/N (%)",
    Coldspot_persistent = "Persistent\nn/N (%)",
    Coldspot_gained = "Newly\nn/N (%)"
  ) |>
  tab_spanner(
    label = "Hotspots",
    columns = starts_with("Hotspot_")
  ) |>
  tab_spanner(
    label = "Coldspots",
    columns = starts_with("Coldspot_")
  ) |>
  cols_align(
    align = "center",
    columns = -state_name
  ) |>
  tab_header(
    title = md(
      "**Table S16.** Changes in hotspot and coldspot classification after covariate adjustment, by region and federative units, Brazil, 2010–2024."
    )
  ) |>
  tab_source_note(
    source_note = md(
      paste(
        "*Note:* Percentages for persistent classifications and classifications attenuated",
        "after adjustment were calculated using the number of hotspots or coldspots in the",
        "unadjusted model as the denominator. Percentages for classifications newly identified",
        "after adjustment were calculated using municipalities not classified as hotspots or",
        "coldspots in the unadjusted model as the denominator. Persistent refers to",
        "classifications retained after adjustment."
      )
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) |>
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(3)
  )

gt::gtsave(
  table_s16,
  filename = file.path(
    RESULTS_TABLE_DIR,
    "Table_S16_hotspot_coldspot_changes.docx"
  )
)

# Save the enriched comparative database for the next analysis ----

readr::write_csv(
  comp,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_classification_stability.csv"
  )
)

saveRDS(
  comp,
  file.path(
    RESULTS_DATA_DIR,
    "municipal_classification_stability.rds"
  )
)

# Reproducibility information ----

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    RESULTS_DATA_DIR,
    "session_info_cluster_trend_stability.txt"
  )
)
