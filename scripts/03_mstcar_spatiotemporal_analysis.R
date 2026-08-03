# Article: Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024:
# a Bayesian spatiotemporal modelling study
# Repository: mortalitytb
# Script: MSTCAR spatiotemporal smoothing, diagnostics and estimate suppression
# Author: José Mário Nunes da Silva, PhD
# Last updated: 2026-08-03

# Packages ----

library(RSTr)
library(coda)
library(tidyverse)
library(sf)
library(spdep)
library(geobr)
library(matrixStats)
library(flextable)
library(officer)
library(ggspatial)
library(patchwork)

# Output directories ----

MODEL_DIR <- "models/mstcar"
OUTPUT_DIR <- "results/mstcar"
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
OBJECT_DIR <- "objects/mstcar"

for (directory in c(
  MODEL_DIR,
  OUTPUT_DIR,
  TABLE_DIR,
  FIGURE_DIR,
  OBJECT_DIR
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# Prepare MSTCAR input data ----

DATA_FILE <- "data/processed/tb_mortality.csv"

STUDY_YEARS <- 2010:2024
AGE_GROUPS <- c("0-19", "20-59", "60+")
PERIODS <- c("2010-2014", "2015-2019", "2020-2024")
RATE_MULTIPLIER <- 1e5

mortality_data <- readr::read_csv(DATA_FILE, show_col_types = FALSE)

mortality_data <- mortality_data |>
  dplyr::transmute(
    code_muni = stringr::str_pad(
      as.character(code_muni),
      6,
      pad = "0"
    ),
    year = suppressWarnings(
      as.integer(year)
    ),
    age_group = as.character(age_group),
    deaths = suppressWarnings(
      as.numeric(deaths)
    ),
    population = suppressWarnings(
      as.numeric(population)
    )
  )

mortality_data <- mortality_data |>
  dplyr::filter(
    year %in% STUDY_YEARS
  ) |>
  dplyr::mutate(
    period = dplyr::case_when(
      year <= 2014 ~ "2010-2014",
      year <= 2019 ~ "2015-2019",
      year <= 2024 ~ "2020-2024",
      TRUE ~ NA_character_
    ),
    age_group_3 = dplyr::case_when(
      age_group %in% c(
        "0-4",
        "5-9",
        "10-14",
        "15-19"
      ) ~
        "0-19",

      age_group %in% c(
        "20-24",
        "25-29",
        "30-34",
        "35-39",
        "40-44",
        "45-49",
        "50-54",
        "55-59"
      ) ~
        "20-59",

      age_group %in% c(
        "60-64",
        "65-69",
        "70-74",
        "75-79",
        "80+"
      ) ~
        "60+",

      TRUE ~
        NA_character_
    )
  )

unclassified_age_groups <- mortality_data |>
  dplyr::filter(is.na(age_group_3)) |>
  dplyr::distinct(age_group) |>
  dplyr::pull(age_group)

if (length(unclassified_age_groups) > 0L) {
  stop(
    "Unclassified age groups: ",
    paste(unclassified_age_groups, collapse = ", ")
  )
}

source_totals <- mortality_data |>
  dplyr::summarise(
    deaths_source = sum(deaths),
    person_years_source = sum(population)
  )

mstcar_long <- mortality_data |>
  dplyr::mutate(
    age_group_3 = factor(age_group_3, levels = AGE_GROUPS),
    period = factor(period, levels = PERIODS)
  ) |>
  dplyr::group_by(code_muni, age_group_3, period) |>
  dplyr::summarise(
    events = sum(deaths),
    population = sum(population),
    .groups = "drop"
  )

aggregated_totals <- mstcar_long |>
  dplyr::summarise(
    deaths_aggregated = sum(events),
    person_years_aggregated = sum(population)
  )

aggregation_check <- dplyr::bind_cols(source_totals, aggregated_totals) |>
  dplyr::mutate(
    deaths_difference = deaths_aggregated - deaths_source,
    person_years_difference = person_years_aggregated - person_years_source,
    totals_preserved = abs(deaths_difference) < 1e-8 &
      abs(person_years_difference) < 1e-4
  )

readr::write_csv(
  aggregation_check,
  file.path(OUTPUT_DIR, "period_aggregation_check.csv")
)

municipality_codes <- sort(unique(mstcar_long$code_muni))

mstcar_long <- mstcar_long |>
  dplyr::mutate(
    code_muni = factor(code_muni, levels = municipality_codes),
    age_group_3 = factor(age_group_3, levels = AGE_GROUPS),
    period = factor(period, levels = PERIODS)
  ) |>
  tidyr::complete(
    code_muni,
    age_group_3,
    period,
    fill = list(events = 0, population = 0)
  ) |>
  dplyr::arrange(code_muni, age_group_3, period)

N_MUNICIPALITIES <- length(municipality_codes)
N_AGE_GROUPS <- length(AGE_GROUPS)
N_PERIODS <- length(PERIODS)

Y <- array(0, dim = c(N_MUNICIPALITIES, N_AGE_GROUPS, N_PERIODS))
N <- array(0, dim = c(N_MUNICIPALITIES, N_AGE_GROUPS, N_PERIODS))

array_index <- cbind(
  as.integer(mstcar_long$code_muni),
  as.integer(mstcar_long$age_group_3),
  as.integer(mstcar_long$period)
)

Y[array_index] <- mstcar_long$events
N[array_index] <- mstcar_long$population

dimnames(Y) <- list(
  region = municipality_codes,
  group = AGE_GROUPS,
  time = PERIODS
)
dimnames(N) <- dimnames(Y)

mstcar_data <- list(Y = Y, n = N)

saveRDS(
  mstcar_data,
  file.path(OBJECT_DIR, "mstcar_input_3_age_groups_3_periods.rds")
)

readr::write_csv(
  mstcar_long,
  file.path(OUTPUT_DIR, "mstcar_input_3_age_groups_3_periods.csv")
)

period_descriptive_summary <- mstcar_long |>
  dplyr::group_by(period, age_group_3) |>
  dplyr::summarise(
    municipalities = dplyr::n_distinct(code_muni),
    events = sum(events),
    person_years = sum(population),
    crude_average_annual_rate_per_100k = dplyr::if_else(
      person_years > 0,
      RATE_MULTIPLIER * events / person_years,
      NA_real_
    ),
    .groups = "drop"
  )

readr::write_csv(
  period_descriptive_summary,
  file.path(OUTPUT_DIR, "period_descriptive_summary.csv")
)

# Build or load the symmetric six-nearest-neighbour graph ----

K_NEIGHBOURS <- 6L
MUNICIPAL_BOUNDARY_YEAR <- 2022L

ADJACENCY_FILE <- file.path(
  OBJECT_DIR,
  sprintf(
    "adj_knn%d_municipalities_%d.rds",
    K_NEIGHBOURS,
    MUNICIPAL_BOUNDARY_YEAR
  )
)

if (file.exists(ADJACENCY_FILE)) {
  adjacency <- readRDS(ADJACENCY_FILE)

  if (length(adjacency) != N_MUNICIPALITIES) {
    stop("The stored adjacency does not match the number of municipalities.")
  }

  if (is.null(names(adjacency))) {
    names(adjacency) <- municipality_codes
  }

  if (!identical(names(adjacency), municipality_codes)) {
    stop(
      "The municipality order in the stored adjacency differs from the data. ",
      "Delete the adjacency file and rerun the script."
    )
  }
} else {
  municipalities_sf <- geobr::read_municipality(
    code_muni = "all",
    year = MUNICIPAL_BOUNDARY_YEAR,
    showProgress = FALSE
  ) |>
    dplyr::mutate(
      code_muni = substr(sprintf("%07.0f", as.numeric(code_muni)), 1, 6)
    ) |>
    dplyr::filter(code_muni %in% municipality_codes) |>
    dplyr::arrange(match(code_muni, municipality_codes))

  if (nrow(municipalities_sf) != N_MUNICIPALITIES ||
      !identical(as.character(municipalities_sf$code_muni), municipality_codes)) {
    stop("Municipal boundaries could not be matched and ordered to the analytical data.")
  }

  municipality_points <- municipalities_sf |>
    sf::st_transform(5880) |>
    sf::st_point_on_surface()

  coordinates <- sf::st_coordinates(municipality_points)
  adjacency <- spdep::knearneigh(coordinates, k = K_NEIGHBOURS) |>
    spdep::knn2nb() |>
    spdep::make.sym.nb()

  if (min(spdep::card(adjacency)) < 1L) {
    stop("At least one municipality has no neighbour in the adjacency graph.")
  }

  names(adjacency) <- municipality_codes
  saveRDS(adjacency, ADJACENCY_FILE)
}

# Fit or load the three MSTCAR chains ----

MODEL_PREFIX <- "tb_mstcar_2010_2024_3groups_3periods_50k"

# Set to TRUE only when fitting new independent chains.
RUN_MULTIPLE_CHAINS <- FALSE

N_CHAINS <- 3L
SEEDS <- c(123L, 2026L, 4517L)

ITERATIONS <- 50000L
BURN_IN <- 25000L
ITERATIONS_PER_BATCH <- 100L
MODEL_CREDIBLE_INTERVAL <- 0.95

if (N_CHAINS < 2L) {
  stop(
    "At least two independent chains are required for convergence diagnostics."
  )
}

if (length(SEEDS) != N_CHAINS || anyDuplicated(SEEDS)) {
  stop("Provide one distinct seed for each MCMC chain.")
}

if (
  ITERATIONS %% ITERATIONS_PER_BATCH != 0L ||
    BURN_IN %% ITERATIONS_PER_BATCH != 0L
) {
  stop(
    "ITERATIONS and BURN_IN must be divisible by ",
    "ITERATIONS_PER_BATCH."
  )
}

if (BURN_IN >= ITERATIONS) {
  stop("BURN_IN must be lower than ITERATIONS.")
}

EXPECTED_BATCHES <- ITERATIONS %/% ITERATIONS_PER_BATCH
BURN_BATCHES <- BURN_IN %/% ITERATIONS_PER_BATCH

CHAIN_NAMES <- sprintf(
  "%s_chain_%02d",
  MODEL_PREFIX,
  seq_len(N_CHAINS)
)

CHAIN_DIRS <- file.path(
  MODEL_DIR,
  CHAIN_NAMES
)

models <- vector("list", N_CHAINS)

for (chain_id in seq_len(N_CHAINS)) {
  chain_name <- CHAIN_NAMES[chain_id]
  chain_path <- CHAIN_DIRS[chain_id]

  if (RUN_MULTIPLE_CHAINS) {
    if (dir.exists(chain_path)) {
      stop(
        "Model directory already exists: ", chain_path,
        ". Set RUN_MULTIPLE_CHAINS to FALSE or change MODEL_PREFIX."
      )
    }

    message("Fitting ", chain_name, "...")
    models[[chain_id]] <- RSTr::mstcar(
      name = chain_name,
      data = mstcar_data,
      adjacency = adjacency,
      dir = MODEL_DIR,
      seed = SEEDS[chain_id],
      method = "poisson",
      iterations = ITERATIONS,
      burn = BURN_IN,
      perc_ci = MODEL_CREDIBLE_INTERVAL,
      show_plots = FALSE,
      verbose = TRUE
    )
  } else {
    if (!dir.exists(chain_path)) {
      stop("Model directory not found: ", chain_path)
    }

    message("Loading ", chain_name, "...")
    models[[chain_id]] <- RSTr::load_model(
      name = chain_name,
      dir = MODEL_DIR
    )
  }
}

names(models) <- CHAIN_NAMES
reference_dimensions <- dim(models[[1L]]$data$Y)
reference_dimnames <- dimnames(models[[1L]]$data$Y)

for (chain_id in seq_len(N_CHAINS)) {
  if (!identical(dim(models[[chain_id]]$data$Y), reference_dimensions) ||
      !identical(dimnames(models[[chain_id]]$data$Y), reference_dimnames)) {
    stop("Data dimensions or dimension names differ across chains.")
  }

  same_events <- isTRUE(all.equal(
    as.numeric(models[[chain_id]]$data$Y),
    as.numeric(mstcar_data$Y),
    tolerance = 0,
    check.attributes = FALSE
  ))

  same_population <- isTRUE(all.equal(
    as.numeric(models[[chain_id]]$data$n),
    as.numeric(mstcar_data$n),
    tolerance = 1e-8,
    check.attributes = FALSE
  ))

  if (!same_events || !same_population) {
    stop(
      "The data stored in chain ", chain_id,
      " do not match the data prepared by this script."
    )
  }
}

region_names <- reference_dimnames[[1L]]
group_names <- reference_dimnames[[2L]]
time_names <- reference_dimnames[[3L]]

if (is.null(region_names)) region_names <- as.character(seq_len(reference_dimensions[1L]))
if (is.null(group_names)) group_names <- as.character(seq_len(reference_dimensions[2L]))
if (is.null(time_names)) time_names <- as.character(seq_len(reference_dimensions[3L]))

# Essential functions for batched RSTr output ----

USE_PARAMETER_CACHE <- TRUE

SCALAR_CACHE_SIGNATURE <- paste0(
  MODEL_PREFIX,
  "_",
  N_CHAINS,
  "chains_",
  ITERATIONS,
  "iter_",
  BURN_IN,
  "burn"
)

SCALAR_CACHE_DIR <- file.path(
  OUTPUT_DIR,
  "cache",
  "selected_parameters",
  SCALAR_CACHE_SIGNATURE
)

dir.create(
  SCALAR_CACHE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

find_parameter_batch_files <- function(chain_dir, parameter, post_burn_only = FALSE) {
  parameter_dir <- file.path(chain_dir, parameter)

  if (!dir.exists(parameter_dir)) {
    stop("Parameter directory not found: ", parameter_dir)
  }

  pattern <- paste0("^", parameter, "_out_([0-9]+)\\.Rds$")
  files <- list.files(parameter_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0L) {
    stop("No batches were found for parameter ", parameter, ".")
  }

  file_table <- tibble::tibble(
    file = files,
    batch = as.integer(sub(pattern, "\\1", basename(files)))
  ) |>
    dplyr::filter(!is.na(batch), batch <= EXPECTED_BATCHES) |>
    dplyr::arrange(batch)

  duplicated_batches <- unique(file_table$batch[duplicated(file_table$batch)])
  missing_batches <- setdiff(seq_len(EXPECTED_BATCHES), file_table$batch)

  if (length(duplicated_batches) > 0L) {
    stop("Duplicated batches for ", parameter, ": ", paste(duplicated_batches, collapse = ", "))
  }

  if (length(missing_batches) > 0L) {
    stop(
      "Missing batches for ", parameter, " in ", basename(chain_dir), ": ",
      paste(missing_batches, collapse = ", ")
    )
  }

  if (post_burn_only) {
    file_table <- file_table |>
      dplyr::filter(batch > BURN_BATCHES)
  }

  file_table
}

axis_labels_for_parameter <- function(parameter, model, object_dimensions) {
  labels_from_model <- model$params$dimnames
  if (is.null(labels_from_model)) labels_from_model <- dimnames(model$data$Y)

  labels <- switch(
    parameter,
    lambda = list(labels_from_model[[1L]], labels_from_model[[2L]], labels_from_model[[3L]]),
    beta = list(seq_len(object_dimensions[1L]), labels_from_model[[2L]], labels_from_model[[3L]]),
    tau2 = list(labels_from_model[[2L]]),
    Ag = list(labels_from_model[[2L]], labels_from_model[[2L]]),
    G = list(labels_from_model[[2L]], labels_from_model[[2L]], labels_from_model[[3L]]),
    lapply(object_dimensions[-length(object_dimensions)], seq_len)
  )

  purrr::map2(labels, object_dimensions[-length(object_dimensions)], function(x, n) {
    if (is.null(x) || length(x) != n) as.character(seq_len(n)) else as.character(x)
  })
}

make_selection_grid <- function(parameter, model, first_object, mode = "all") {
  object_dimensions <- dim(first_object)

  if (is.null(object_dimensions) || length(object_dimensions) < 2L) {
    stop("Unexpected object structure for parameter ", parameter, ".")
  }

  non_draw_dimensions <- object_dimensions[-length(object_dimensions)]
  labels <- axis_labels_for_parameter(parameter, model, object_dimensions)

  if (mode == "all") {
    selection <- expand.grid(
      lapply(non_draw_dimensions, seq_len),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    names(selection) <- paste0("d", seq_along(non_draw_dimensions))
  } else if (parameter == "G") {
    group_positions <- unique(c(1L, ceiling(non_draw_dimensions[1L] / 2), non_draw_dimensions[1L]))
    time_positions <- unique(c(1L, ceiling(non_draw_dimensions[3L] / 2), non_draw_dimensions[3L]))

    selection <- tidyr::crossing(d1 = group_positions, d3 = time_positions) |>
      dplyr::mutate(d2 = d1) |>
      dplyr::select(d1, d2, d3)

    if (non_draw_dimensions[1L] > 1L) {
      selection <- dplyr::bind_rows(
        selection,
        tibble::tibble(
          d1 = 1L,
          d2 = non_draw_dimensions[2L],
          d3 = ceiling(non_draw_dimensions[3L] / 2)
        )
      )
    }
  } else if (parameter == "Ag") {
    group_positions <- unique(c(1L, ceiling(non_draw_dimensions[1L] / 2), non_draw_dimensions[1L]))
    selection <- tibble::tibble(d1 = group_positions, d2 = group_positions)

    if (non_draw_dimensions[1L] > 1L) {
      selection <- dplyr::bind_rows(
        selection,
        tibble::tibble(d1 = 1L, d2 = non_draw_dimensions[2L])
      )
    }
  } else {
    stop("Selected mode is available only for G and Ag.")
  }

  selection <- dplyr::distinct(selection)

  axis_names <- switch(
    parameter,
    lambda = c("region", "group", "time"),
    beta = c("island", "group", "time"),
    tau2 = "group",
    Ag = c("group1", "group2"),
    G = c("group1", "group2", "time"),
    paste0("dim", seq_along(non_draw_dimensions))
  )

  selection$series <- purrr::map_chr(seq_len(nrow(selection)), function(i) {
    values <- purrr::map_chr(seq_along(non_draw_dimensions), function(k) {
      labels[[k]][selection[[paste0("d", k)]][i]]
    })

    paste0(
      parameter, "[",
      paste0(axis_names, "=", values, collapse = "; "),
      "]"
    )
  })

  selection
}

extract_selected_from_object <- function(object, selection) {
  object_dimensions <- dim(object)
  draw_dimension <- length(object_dimensions)
  n_draws <- object_dimensions[draw_dimension]
  n_non_draw_dimensions <- draw_dimension - 1L

  output <- vapply(seq_len(nrow(selection)), function(i) {
    indices <- as.list(as.integer(unlist(
      selection[i, paste0("d", seq_len(n_non_draw_dimensions)), drop = FALSE],
      use.names = FALSE
    )))

    as.numeric(do.call(
      `[`,
      c(list(object), indices, list(seq_len(n_draws)), list(drop = TRUE))
    ))
  }, numeric(n_draws))

  if (is.null(dim(output))) output <- matrix(output, ncol = 1L)
  colnames(output) <- selection$series
  output
}

load_selected_parameter_chain <- function(
    chain_id,
    parameter,
    selection,
    scale = 1,
    use_cache = USE_PARAMETER_CACHE
) {
  selection_key <- sum(utf8ToInt(paste(selection$series, collapse = "|"))) %% 100000000L
  cache_file <- file.path(
    SCALAR_CACHE_DIR,
    sprintf("%s_chain_%02d_full_%08d.rds", parameter, chain_id, selection_key)
  )

  if (use_cache && file.exists(cache_file)) {
    cached_matrix <- readRDS(cache_file)
    cached_batches <- attr(cached_matrix, "batch_numbers")
    cached_draws_per_batch <- attr(cached_matrix, "draws_per_batch")

    cache_is_valid <- is.matrix(cached_matrix) &&
      !is.null(cached_batches) &&
      !is.null(cached_draws_per_batch) &&
      identical(as.integer(cached_batches), seq_len(EXPECTED_BATCHES)) &&
      nrow(cached_matrix) == EXPECTED_BATCHES * as.integer(cached_draws_per_batch) &&
      all(selection$series %in% colnames(cached_matrix))

    if (cache_is_valid) return(cached_matrix)
    message("Ignoring incompatible cache: ", basename(cache_file))
  }

  file_table <- find_parameter_batch_files(
    CHAIN_DIRS[chain_id],
    parameter,
    post_burn_only = FALSE
  )

  batch_matrices <- vector("list", nrow(file_table))

  for (batch_id in seq_len(nrow(file_table))) {
    batch_object <- readRDS(file_table$file[batch_id])
    batch_matrices[[batch_id]] <- extract_selected_from_object(batch_object, selection)

    if (batch_id %% 50L == 0L || batch_id == nrow(file_table)) {
      message(
        "  ", parameter, ", chain ", chain_id, ": ",
        batch_id, "/", nrow(file_table), " batches"
      )
    }
  }

  draws_per_batch <- vapply(batch_matrices, nrow, integer(1))
  if (length(unique(draws_per_batch)) != 1L) {
    stop("Inconsistent number of draws per batch for parameter ", parameter, ".")
  }

  parameter_matrix <- do.call(rbind, batch_matrices) * scale
  usable_columns <- vapply(seq_len(ncol(parameter_matrix)), function(column_id) {
    values <- parameter_matrix[, column_id]
    all(is.finite(values)) && length(unique(values)) > 1L && stats::sd(values) > 0
  }, logical(1))

  parameter_matrix <- parameter_matrix[, usable_columns, drop = FALSE]
  if (ncol(parameter_matrix) == 0L) {
    stop("No usable scalar series remained for parameter ", parameter, ".")
  }

  attr(parameter_matrix, "draws_per_batch") <- unique(draws_per_batch)
  attr(parameter_matrix, "batch_numbers") <- file_table$batch
  saveRDS(parameter_matrix, cache_file)
  parameter_matrix
}

subset_post_burn <- function(parameter_matrix) {
  draws_per_batch <- attr(parameter_matrix, "draws_per_batch")
  batch_numbers <- attr(parameter_matrix, "batch_numbers")

  if (is.null(draws_per_batch) || is.null(batch_numbers)) {
    stop("Batch attributes required to apply burn-in are missing.")
  }

  keep_rows <- rep(batch_numbers > BURN_BATCHES, each = draws_per_batch)
  output <- parameter_matrix[keep_rows, , drop = FALSE]
  attr(output, "draws_per_batch") <- draws_per_batch
  attr(output, "retained_batches") <- sum(batch_numbers > BURN_BATCHES)
  output
}

make_mcmc_list <- function(matrices, parameter) {
  common_columns <- purrr::reduce(purrr::map(matrices, colnames), intersect)
  if (length(common_columns) == 0L) {
    stop("No common scalar series across chains for parameter ", parameter, ".")
  }

  matrices <- purrr::map(matrices, ~ .x[, common_columns, drop = FALSE])
  common_n <- min(purrr::map_int(matrices, nrow))
  matrices <- purrr::map(matrices, ~ tail(.x, common_n))

  do.call(
    coda::mcmc.list,
    lapply(matrices, function(x) coda::mcmc(x, start = 1L, thin = 1L))
  )
}

# Select representative mortality-rate series ----

RATE_GROUP_POSITIONS <- c(
  "first",
  "middle",
  "last"
)

RATE_TIME_POSITIONS <- c(
  "first",
  "middle",
  "last"
)

LOW_INFORMATION_QUANTILE <- 0.10

Y_model <- models[[1L]]$data$Y
N_model <- models[[1L]]$data$n
position_index <- function(positions, n) {
  unname(c(first = 1L, middle = ceiling(n / 2), last = n)[positions])
}

group_indices <- unique(position_index(RATE_GROUP_POSITIONS, dim(Y_model)[2L]))
time_indices <- unique(position_index(RATE_TIME_POSITIONS, dim(Y_model)[3L]))
selected_rates <- list()
selection_counter <- 1L

for (group_id in group_indices) {
  for (time_id in time_indices) {
    population_values <- N_model[, group_id, time_id]
    event_values <- Y_model[, group_id, time_id]
    valid_regions <- which(is.finite(population_values) & population_values > 0)

    if (length(valid_regions) == 0L) next

    high_region <- valid_regions[which.max(population_values[valid_regions])]
    low_target <- as.numeric(stats::quantile(
      population_values[valid_regions],
      probs = LOW_INFORMATION_QUANTILE,
      na.rm = TRUE
    ))

    zero_event_regions <- valid_regions[event_values[valid_regions] == 0]
    low_pool <- if (length(zero_event_regions) > 0L) zero_event_regions else valid_regions
    low_region <- low_pool[which.min(abs(population_values[low_pool] - low_target))]

    for (selection in list(
      list(region_id = low_region, information = "Low information"),
      list(region_id = high_region, information = "High information")
    )) {
      selected_rates[[selection_counter]] <- tibble::tibble(
        d1 = selection$region_id,
        d2 = group_id,
        d3 = time_id,
        region = region_names[selection$region_id],
        group = group_names[group_id],
        time = time_names[time_id],
        information = selection$information,
        events = as.numeric(event_values[selection$region_id]),
        population = as.numeric(population_values[selection$region_id])
      )
      selection_counter <- selection_counter + 1L
    }
  }
}

lambda_selection <- dplyr::bind_rows(selected_rates) |>
  dplyr::distinct(d1, d2, d3, .keep_all = TRUE) |>
  dplyr::mutate(
    series = paste0(
      "lambda[region=", region,
      "; group=", group,
      "; time=", time,
      "; ", information, "]"
    )
  )

if (nrow(lambda_selection) < 2L) {
  stop("Fewer than two representative lambda series were selected.")
}

readr::write_csv(
  lambda_selection |>
    dplyr::select(region, group, time, information, events, population, series),
  file.path(OUTPUT_DIR, "selected_mortality_rates.csv")
)

# Load selected scalar series ----

first_batch_object <- function(parameter) {
  first_file <- find_parameter_batch_files(
    CHAIN_DIRS[1L], parameter, post_burn_only = FALSE
  )$file[1L]
  readRDS(first_file)
}

beta_selection <- make_selection_grid(
  "beta", models[[1L]], first_batch_object("beta"), mode = "all"
)
tau2_selection <- make_selection_grid(
  "tau2", models[[1L]], first_batch_object("tau2"), mode = "all"
)
G_selection <- make_selection_grid(
  "G", models[[1L]], first_batch_object("G"), mode = "selected"
)
Ag_selection <- make_selection_grid(
  "Ag", models[[1L]], first_batch_object("Ag"), mode = "selected"
)

selections <- list(
  lambda = lambda_selection |>
    dplyr::select(d1, d2, d3, series),
  beta = beta_selection,
  tau2 = tau2_selection,
  G = G_selection,
  Ag = Ag_selection
)

parameter_scales <- c(lambda = RATE_MULTIPLIER, beta = 1, tau2 = 1, G = 1, Ag = 1)

loaded_full <- purrr::map(names(selections), function(parameter) {
  message("Loading selected ", parameter, " series...")
  purrr::map(seq_len(N_CHAINS), function(chain_id) {
    load_selected_parameter_chain(
      chain_id = chain_id,
      parameter = parameter,
      selection = selections[[parameter]],
      scale = parameter_scales[[parameter]],
      use_cache = USE_PARAMETER_CACHE
    )
  })
})
names(loaded_full) <- names(selections)

loaded_post <- purrr::map(loaded_full, ~ purrr::map(.x, subset_post_burn))
mcmc_full <- purrr::imap(loaded_full, make_mcmc_list)
mcmc_post <- purrr::imap(loaded_post, make_mcmc_list)

# Diagnostic reference values ----

RHAT_THRESHOLD <- 1.01
ESS_THRESHOLD <- 400

# Table S2: model and MCMC configuration ----

draws_per_batch <- attr(loaded_full$lambda[[1L]], "draws_per_batch")
full_saved_draws_per_chain <- nrow(loaded_full$lambda[[1L]])
retained_saved_draws_per_chain <- nrow(loaded_post$lambda[[1L]])
total_retained_draws <- retained_saved_draws_per_chain * N_CHAINS

model_configuration <- tibble::tibble(
  Item = c(
    "Likelihood",
    "Age groups",
    "Periods",
    "Spatial structure",
    "Number of chains",
    "Iterations per chain",
    "Burn-in per chain",
    "Stored draws per chain before burn-in exclusion",
    "Retained stored draws per chain",
    "Total retained draws",
    "Posterior interval",
    "Diagnostic reference values",
    "Age standardisation"
  ),
  Specification = c(
    "Poisson",
    "0–19, 20–59 and ≥60 years",
    "2010–2014, 2015–2019 and 2020–2024",
    if (K_NEIGHBOURS == 6L) {
      "Symmetric six-nearest-neighbour graph"
    } else {
      paste0("Symmetric ", K_NEIGHBOURS, "-nearest-neighbour graph")
    },
    format(N_CHAINS, big.mark = ","),
    format(ITERATIONS, big.mark = ","),
    format(BURN_IN, big.mark = ","),
    format(full_saved_draws_per_chain, big.mark = ","),
    format(retained_saved_draws_per_chain, big.mark = ","),
    format(total_retained_draws, big.mark = ","),
    paste0(100 * MODEL_CREDIBLE_INTERVAL, "% credible interval"),
    paste0("R̂ ≤ ", RHAT_THRESHOLD, "; ESS ≥ ", ESS_THRESHOLD),
    "WHO World Standard Population"
  )
)

readr::write_csv(
  model_configuration,
  file.path(TABLE_DIR, "Table_S2_MSTCAR_configuration.csv")
)

model_configuration_ft <- flextable::flextable(model_configuration) |>
  flextable::theme_booktabs() |>
  flextable::bold(part = "header") |>
  flextable::set_caption(caption = "Table S2. MSTCAR model and MCMC configuration.") |>
  flextable::add_footer_lines(values = "Abbreviation: WHO, World Health Organisation.") |>
  flextable::autofit()

flextable::save_as_docx(
  model_configuration_ft,
  path = file.path(TABLE_DIR, "Table_S2_MSTCAR_configuration.docx")
)

# Convergence and sampling-efficiency diagnostics ----

safe_rhat_table <- function(mcmc_list) {
  gelman_result <- tryCatch(
    coda::gelman.diag(
      mcmc_list,
      confidence = 0.95,
      autoburnin = FALSE,
      multivariate = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(gelman_result)) {
    psrf <- as.data.frame(gelman_result$psrf)
    return(tibble::tibble(
      scalar_parameter = rownames(psrf),
      rhat_point = psrf[, 1L],
      rhat_upper_95 = psrf[, 2L]
    ))
  }

  scalar_names <- colnames(as.matrix(mcmc_list[[1L]]))
  purrr::map_dfr(scalar_names, function(scalar_name) {
    scalar_result <- tryCatch(
      coda::gelman.diag(
        mcmc_list[, scalar_name, drop = FALSE],
        confidence = 0.95,
        autoburnin = FALSE,
        multivariate = FALSE
      ),
      error = function(e) NULL
    )

    tibble::tibble(
      scalar_parameter = scalar_name,
      rhat_point = if (is.null(scalar_result)) NA_real_ else scalar_result$psrf[1L, 1L],
      rhat_upper_95 = if (is.null(scalar_result)) NA_real_ else scalar_result$psrf[1L, 2L]
    )
  })
}

diagnose_parameter <- function(mcmc_list, parameter) {
  chains <- unclass(mcmc_list)
  combined_ess <- coda::effectiveSize(mcmc_list)

  per_chain_ess <- purrr::map2_dfr(chains, seq_along(chains), function(chain, chain_id) {
    ess <- coda::effectiveSize(chain)
    tibble::tibble(
      scalar_parameter = names(ess),
      chain_id = chain_id,
      ess_chain = as.numeric(ess)
    )
  })

  ess_table <- tibble::tibble(
    scalar_parameter = names(combined_ess),
    ess_combined = as.numeric(combined_ess)
  ) |>
    dplyr::left_join(
      per_chain_ess |>
        dplyr::group_by(scalar_parameter) |>
        dplyr::summarise(
          ess_min_chain = min(ess_chain, na.rm = TRUE),
          ess_median_chain = median(ess_chain, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "scalar_parameter"
    )

  combined_matrix <- do.call(rbind, lapply(chains, as.matrix))
  posterior_means <- colMeans(combined_matrix)

  ess_table |>
    dplyr::left_join(safe_rhat_table(mcmc_list), by = "scalar_parameter") |>
    dplyr::mutate(
      parameter_group = parameter,
      posterior_mean = posterior_means[scalar_parameter],
      rhat_flag = dplyr::case_when(
        is.na(rhat_point) ~ "Not available",
        rhat_point <= RHAT_THRESHOLD ~ "Adequate",
        TRUE ~ "Review"
      ),
      ess_flag = dplyr::case_when(
        is.na(ess_combined) ~ "Not available",
        ess_combined >= ESS_THRESHOLD ~ "Adequate",
        TRUE ~ "Low"
      )
    ) |>
    dplyr::select(
      parameter_group, scalar_parameter, posterior_mean,
      ess_combined, ess_min_chain, ess_median_chain,
      rhat_point, rhat_upper_95, rhat_flag, ess_flag
    )
}

scalar_diagnostics <- purrr::imap_dfr(mcmc_post, diagnose_parameter)

readr::write_csv(
  scalar_diagnostics,
  file.path(OUTPUT_DIR, "selected_scalar_diagnostics.csv")
)

safe_statistic <- function(x, statistic) {
  if (all(is.na(x))) NA_real_ else statistic(x, na.rm = TRUE)
}

compact_diagnostic_summary <- scalar_diagnostics |>
  dplyr::group_by(parameter_group) |>
  dplyr::summarise(
    n_scalar_series = dplyr::n(),
    min_ess_combined = safe_statistic(ess_combined, min),
    median_ess_combined = safe_statistic(ess_combined, median),
    min_ess_single_chain = safe_statistic(ess_min_chain, min),
    median_rhat = safe_statistic(rhat_point, median),
    max_rhat = safe_statistic(rhat_point, max),
    max_rhat_upper_95 = safe_statistic(rhat_upper_95, max),
    rhat_le_threshold_pct = 100 * safe_statistic(
      rhat_point <= RHAT_THRESHOLD,
      mean
    ),
    ess_ge_threshold_pct = 100 * safe_statistic(
      ess_combined >= ESS_THRESHOLD,
      mean
    ),
    .groups = "drop"
  )

readr::write_csv(
  compact_diagnostic_summary,
  file.path(OUTPUT_DIR, "compact_convergence_summary.csv")
)

# Table S3: convergence diagnostics ----

table_s3 <- compact_diagnostic_summary |>
  dplyr::mutate(
    Parameter = dplyr::recode(
      parameter_group,
      Ag = "Ag",
      G = "G",
      beta = "Beta",
      lambda = "Lambda",
      tau2 = "Tau²"
    ),
    `Number of series` = as.character(n_scalar_series),
    `Minimum combined ESS` = formatC(min_ess_combined, format = "f", digits = 2),
    `Median combined ESS` = formatC(median_ess_combined, format = "f", digits = 2),
    `Minimum ESS (per chain)` = formatC(min_ess_single_chain, format = "f", digits = 2),
    `Median R̂` = formatC(median_rhat, format = "f", digits = 3),
    `Maximum R̂` = formatC(max_rhat, format = "f", digits = 3),
    `Max R̂ (upper 95%)` = formatC(max_rhat_upper_95, format = "f", digits = 3)
  ) |>
  dplyr::select(
    Parameter,
    `Number of series`,
    `Minimum combined ESS`,
    `Median combined ESS`,
    `Minimum ESS (per chain)`,
    `Median R̂`,
    `Maximum R̂`,
    `Max R̂ (upper 95%)`
  ) |>
  tidyr::pivot_longer(
    -Parameter,
    names_to = "Item",
    values_to = "Value"
  ) |>
  tidyr::pivot_wider(names_from = Parameter, values_from = Value) |>
  dplyr::select(Item, Ag, G, Beta, Lambda, `Tau²`)

readr::write_csv(table_s3, file.path(TABLE_DIR, "Table_S3_MSTCAR_diagnostics.csv"))

table_s3_ft <- flextable::flextable(table_s3) |>
  flextable::theme_booktabs() |>
  flextable::bold(part = "header") |>
  flextable::align(j = 2:6, align = "center", part = "all") |>
  flextable::set_caption(
    caption = "Table S3. Convergence and sampling-efficiency diagnostics for selected MSTCAR parameters."
  ) |>
  flextable::add_footer_lines(
    values = "Abbreviations: ESS, effective sample size; R̂, potential scale reduction factor."
  ) |>
  flextable::autofit()

flextable::save_as_docx(
  table_s3_ft,
  path = file.path(TABLE_DIR, "Table_S3_MSTCAR_diagnostics.docx")
)

# Representative traceplots ----

TRACE_FIGURE_WIDTH_IN <- 13
TRACE_FIGURE_HEIGHT_IN <- 10
TRACE_PNG_DPI <- 300

available_lambda_selection <- lambda_selection |>
  dplyr::filter(
    series %in% colnames(
      as.matrix(
        mcmc_full$lambda[[1L]]
      )
    )
  )

if (nrow(available_lambda_selection) < 2L) {
  stop(
    "Fewer than two representative mortality-rate series remained ",
    "available for the traceplot figure."
  )
}

trace_selection <- c(
  lambda_low = available_lambda_selection$series[
    which.min(available_lambda_selection$population)
  ],
  lambda_high = available_lambda_selection$series[
    which.max(available_lambda_selection$population)
  ],
  beta_first = colnames(as.matrix(mcmc_full$beta[[1L]]))[1L],
  beta_last = tail(colnames(as.matrix(mcmc_full$beta[[1L]])), 1L),
  tau2_first = colnames(as.matrix(mcmc_full$tau2[[1L]]))[1L],
  tau2_last = tail(colnames(as.matrix(mcmc_full$tau2[[1L]])), 1L),
  G_selected = colnames(as.matrix(mcmc_full$G[[1L]]))[1L],
  Ag_selected = colnames(as.matrix(mcmc_full$Ag[[1L]]))[1L]
)

readr::write_csv(
  tibble::tibble(panel = LETTERS[1:8], scalar_parameter = unname(trace_selection)),
  file.path(OUTPUT_DIR, "representative_traceplot_panels.csv")
)

chain_colours <- grDevices::hcl.colors(N_CHAINS, palette = "Dark 3")

plot_one_trace <- function(
    mcmc_list,
    scalar_name,
    panel_title,
    y_label,
    show_legend = FALSE
) {
  matrices <- lapply(unclass(mcmc_list), as.matrix)
  n_draws <- nrow(matrices[[1L]])
  x <- seq_len(n_draws)
  y_range <- range(unlist(lapply(matrices, function(x) x[, scalar_name])), finite = TRUE)

  graphics::plot(
    x,
    matrices[[1L]][, scalar_name],
    type = "l",
    col = chain_colours[1L],
    lwd = 0.7,
    xlab = "Saved posterior draw",
    ylab = y_label,
    main = panel_title,
    ylim = y_range,
    las = 1
  )

  if (length(matrices) > 1L) {
    for (chain_id in 2:length(matrices)) {
      graphics::lines(
        x,
        matrices[[chain_id]][, scalar_name],
        col = chain_colours[chain_id],
        lwd = 0.7
      )
    }
  }

  burn_saved_draw <- BURN_BATCHES * attr(loaded_full$lambda[[1L]], "draws_per_batch")
  graphics::abline(v = burn_saved_draw, lty = 2, lwd = 1)

  if (show_legend) {
    graphics::legend(
      "topright",
      legend = c(paste("Chain", seq_len(N_CHAINS)), "Burn-in cutoff"),
      col = c(chain_colours, "black"),
      lty = c(rep(1L, N_CHAINS), 2L),
      lwd = 1,
      bty = "n",
      cex = 0.72
    )
  }
}

make_trace_figure <- function(device = c("pdf", "png")) {
  device <- match.arg(device)
  output_file <- file.path(
    FIGURE_DIR,
    if (device == "pdf") {
      "Fig_S3_MSTCAR_traceplots.pdf"
    } else {
      "Fig_S3_MSTCAR_traceplots_300dpi.png"
    }
  )

  if (device == "pdf") {
    grDevices::pdf(
      output_file,
      width = TRACE_FIGURE_WIDTH_IN,
      height = TRACE_FIGURE_HEIGHT_IN,
      onefile = TRUE
    )
  } else {
    grDevices::png(
      output_file,
      width = TRACE_FIGURE_WIDTH_IN,
      height = TRACE_FIGURE_HEIGHT_IN,
      units = "in",
      res = TRACE_PNG_DPI
    )
  }

  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_parameters)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(
    mfrow = c(4, 2),
    mar = c(3.5, 4.2, 3.0, 1.0),
    oma = c(2.0, 0.5, 1.0, 0.5),
    mgp = c(2.2, 0.7, 0),
    cex.axis = 0.78,
    cex.lab = 0.85,
    cex.main = 0.82
  )

  panel_parameters <- list(
    list(mcmc_full$lambda, trace_selection[["lambda_low"]], "A", "Rate per 100,000", TRUE),
    list(mcmc_full$lambda, trace_selection[["lambda_high"]], "B", "Rate per 100,000", FALSE),
    list(mcmc_full$beta, trace_selection[["beta_first"]], "C", expression(beta), FALSE),
    list(mcmc_full$beta, trace_selection[["beta_last"]], "D", expression(beta), FALSE),
    list(mcmc_full$tau2, trace_selection[["tau2_first"]], "E", expression(tau^2), FALSE),
    list(mcmc_full$tau2, trace_selection[["tau2_last"]], "F", expression(tau^2), FALSE),
    list(mcmc_full$G, trace_selection[["G_selected"]], "G", "G", FALSE),
    list(mcmc_full$Ag, trace_selection[["Ag_selected"]], "H", expression(A[G]), FALSE)
  )

  for (panel in panel_parameters) {
    plot_one_trace(
      mcmc_list = panel[[1L]],
      scalar_name = panel[[2L]],
      panel_title = paste0(
        panel[[3L]], ". ",
        stringr::str_trunc(panel[[2L]], width = 68L, side = "right")
      ),
      y_label = panel[[4L]],
      show_legend = panel[[5L]]
    )
  }

  graphics::mtext(
    paste0("Representative traceplots from ", N_CHAINS, " independent MSTCAR chains"),
    side = 1,
    outer = TRUE,
    line = 0.5,
    cex = 0.95
  )

  invisible(output_file)
}

make_trace_figure("pdf")
make_trace_figure("png")

low_rate_row <- available_lambda_selection |>
  dplyr::slice_min(population, n = 1, with_ties = FALSE)
high_rate_row <- available_lambda_selection |>
  dplyr::slice_max(population, n = 1, with_ties = FALSE)

figure_caption <- paste0(
  "Fig S3. Representative traceplots for the MSTCAR model. Traceplots are shown for ",
  "selected municipality-, age-group-, and period-specific mortality rates (λ), mean ",
  "parameters (β), residual variance parameters (τ²), period-specific cross-age covariance ",
  "components (G), and the common covariance hyperparameter (Ag). The representative ",
  "mortality-rate panels show municipalities with low and high information: code ",
  low_rate_row$region, " and code ", high_rate_row$region, ", respectively. The three ",
  "colours represent independent MCMC chains. Each chain was run for ",
  format(ITERATIONS, big.mark = ","), " iterations, with the first ",
  format(BURN_IN, big.mark = ","), " iterations discarded as burn-in. The vertical dashed ",
  "line marks the end of the burn-in period."
)

writeLines(figure_caption, file.path(FIGURE_DIR, "Fig_S3_caption.txt"))

# Agreement between chain-specific posterior medians ----

lambda_post_matrices <- lapply(unclass(mcmc_post$lambda), as.matrix)
lambda_series <- colnames(lambda_post_matrices[[1L]])

lambda_chain_summaries <- purrr::map2_dfr(
  lambda_post_matrices,
  seq_along(lambda_post_matrices),
  function(chain_matrix, chain_id) {
    tibble::tibble(
      scalar_parameter = lambda_series,
      chain_id = chain_id,
      posterior_median = matrixStats::colMedians(chain_matrix),
      posterior_mean = colMeans(chain_matrix)
    )
  }
)

readr::write_csv(
  lambda_chain_summaries,
  file.path(OUTPUT_DIR, "lambda_chain_specific_summaries.csv")
)

lambda_median_wide <- lambda_chain_summaries |>
  dplyr::select(scalar_parameter, chain_id, posterior_median) |>
  dplyr::mutate(chain = sprintf("chain_%02d", chain_id)) |>
  dplyr::select(-chain_id) |>
  tidyr::pivot_wider(names_from = chain, values_from = posterior_median)

chain_columns <- grep("^chain_", names(lambda_median_wide), value = TRUE)

lambda_between_chain_stability <- lambda_median_wide |>
  dplyr::rowwise() |>
  dplyr::mutate(
    median_across_chains = median(c_across(dplyr::all_of(chain_columns)), na.rm = TRUE),
    min_across_chains = min(c_across(dplyr::all_of(chain_columns)), na.rm = TRUE),
    max_across_chains = max(c_across(dplyr::all_of(chain_columns)), na.rm = TRUE),
    relative_range_pct = 100 *
      (max_across_chains - min_across_chains) /
      pmax(abs(median_across_chains), .Machine$double.eps)
  ) |>
  dplyr::ungroup()

readr::write_csv(
  lambda_between_chain_stability,
  file.path(OUTPUT_DIR, "lambda_between_chain_stability.csv")
)

lambda_pairwise_correlations <- combn(chain_columns, 2L, simplify = FALSE) |>
  purrr::map_dfr(function(pair) {
    x <- lambda_median_wide[[pair[1L]]]
    y <- lambda_median_wide[[pair[2L]]]

    tibble::tibble(
      chain_1 = pair[1L],
      chain_2 = pair[2L],
      pearson = stats::cor(x, y, method = "pearson", use = "complete.obs"),
      spearman = stats::cor(x, y, method = "spearman", use = "complete.obs"),
      median_relative_difference_pct = median(
        100 * abs(x - y) /
          pmax(abs((x + y) / 2), .Machine$double.eps),
        na.rm = TRUE
      )
    )
  })

readr::write_csv(
  lambda_pairwise_correlations,
  file.path(OUTPUT_DIR, "lambda_pairwise_chain_comparisons.csv")
)

# Table S4: agreement between chain-specific rates ----

table_s4 <- lambda_pairwise_correlations |>
  dplyr::transmute(
    Comparison = paste0(
      "Chains ",
      as.integer(stringr::str_remove(chain_1, "chain_")),
      "–",
      as.integer(stringr::str_remove(chain_2, "chain_"))
    ),
    Pearson = formatC(pearson, format = "f", digits = 4),
    Spearman = formatC(spearman, format = "f", digits = 4),
    `Median relative difference` = paste0(
      formatC(median_relative_difference_pct, format = "f", digits = 2),
      "%"
    )
  )

readr::write_csv(
  table_s4,
  file.path(TABLE_DIR, "Table_S4_chain_specific_rate_agreement.csv")
)

table_s4_ft <- flextable::flextable(table_s4) |>
  flextable::theme_booktabs() |>
  flextable::bold(part = "header") |>
  flextable::align(j = 2:4, align = "center", part = "all") |>
  flextable::set_caption(
    caption = "Table S4. Agreement between chain-specific posterior median mortality rates."
  ) |>
  flextable::autofit()

flextable::save_as_docx(
  table_s4_ft,
  path = file.path(TABLE_DIR, "Table_S4_chain_specific_rate_agreement.docx")
)

# Suppression rules for 95% and 90% intervals ----

SUPPRESSION_CREDIBLE_INTERVALS <- c(
  0.95,
  0.90
)

POPULATION_SUPPRESSION_THRESHOLD <- 1000

CI_95_ALPHA <- (
  1 - SUPPRESSION_CREDIBLE_INTERVALS[1L]
) / 2

CI_90_ALPHA <- (
  1 - SUPPRESSION_CREDIBLE_INTERVALS[2L]
) / 2

SUPPRESSION_QUANTILE_PROBABILITIES <- c(
  CI_95_ALPHA,
  CI_90_ALPHA,
  0.50,
  1 - CI_90_ALPHA,
  1 - CI_95_ALPHA
)

apply_suppression_rules <- function(rate_data) {
  rate_data |>
    dplyr::mutate(
      interval_width_95 = ci_upper_95 - ci_lower_95,
      interval_width_90 = ci_upper_90 - ci_lower_90,

      relative_precision_95 = dplyr::if_else(
        is.finite(interval_width_95) &
          interval_width_95 > 0,
        posterior_median / interval_width_95,
        NA_real_
      ),

      relative_precision_90 = dplyr::if_else(
        is.finite(interval_width_90) &
          interval_width_90 > 0,
        posterior_median / interval_width_90,
        NA_real_
      ),

      suppressed_low_precision_95 =
        is.na(relative_precision_95) |
        relative_precision_95 < 1,

      suppressed_low_precision_90 =
        is.na(relative_precision_90) |
        relative_precision_90 < 1,

      suppressed_low_population =
        average_annual_population <
        POPULATION_SUPPRESSION_THRESHOLD,

      suppressed_95 =
        suppressed_low_precision_95 |
        suppressed_low_population,

      suppressed_90 =
        suppressed_low_precision_90 |
        suppressed_low_population,

      posterior_median_suppressed_95 =
        dplyr::if_else(
          suppressed_95,
          NA_real_,
          posterior_median
        ),

      posterior_median_suppressed_90 =
        dplyr::if_else(
          suppressed_90,
          NA_real_,
          posterior_median
        ),

      medians_suppressed_95 =
        posterior_median_suppressed_95,

      medians_suppressed_90 =
        posterior_median_suppressed_90,

      suppression_reason_95 = dplyr::case_when(
        suppressed_low_precision_95 &
          suppressed_low_population ~
          paste0(
            "Relative precision based on the 95% interval < 1 ",
            "and average annual population < ",
            POPULATION_SUPPRESSION_THRESHOLD
          ),

        suppressed_low_precision_95 ~
          "Relative precision based on the 95% interval < 1",

        suppressed_low_population ~
          paste0(
            "Average annual population < ",
            POPULATION_SUPPRESSION_THRESHOLD
          ),

        TRUE ~
          "Not suppressed"
      ),

      suppression_reason_90 = dplyr::case_when(
        suppressed_low_precision_90 &
          suppressed_low_population ~
          paste0(
            "Relative precision based on the 90% interval < 1 ",
            "and average annual population < ",
            POPULATION_SUPPRESSION_THRESHOLD
          ),

        suppressed_low_precision_90 ~
          "Relative precision based on the 90% interval < 1",

        suppressed_low_population ~
          paste0(
            "Average annual population < ",
            POPULATION_SUPPRESSION_THRESHOLD
          ),

        TRUE ~
          "Not suppressed"
      )
    )
}

# Combine posterior draws and calculate period-specific age-standardised rates ----

WHO_STANDARD_17 <- c(
  8.86,
  8.69,
  8.60,
  8.47,
  8.22,
  7.93,
  7.61,
  7.15,
  6.59,
  6.04,
  5.37,
  4.55,
  3.72,
  2.96,
  2.21,
  1.52,
  0.91 + 0.44 + 0.15 + 0.04 + 0.005
)

WHO_STANDARD_3 <- c(
  "0-19" = sum(WHO_STANDARD_17[1:4]),
  "20-59" = sum(WHO_STANDARD_17[5:12]),
  "60+" = sum(WHO_STANDARD_17[13:17])
)

WHO_STANDARD_3 <- WHO_STANDARD_3 /
  sum(WHO_STANDARD_3)

YEARS_PER_PERIOD <- 5L
REGION_CHUNK_SIZE <- 500L
USE_RATE_CACHE <- TRUE

# Use a compact cache identifier to avoid Windows path-length limitations.
RATE_CACHE_ID <- paste0(
  "c", N_CHAINS,
  "_i", ITERATIONS,
  "_b", BURN_IN,
  "_p", POPULATION_SUPPRESSION_THRESHOLD,
  "_w", paste(round(WHO_STANDARD_3 * 1e6), collapse = "-")
)

RATE_CACHE_DIR <- file.path(
  OUTPUT_DIR,
  "cache",
  RATE_CACHE_ID
)

PERIOD_CHUNK_DIR <- file.path(
  RATE_CACHE_DIR,
  "period"
)

OVERALL_CHUNK_DIR <- file.path(
  RATE_CACHE_DIR,
  "overall"
)

create_cache_directory <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = TRUE
    )
  }

  if (!dir.exists(directory)) {
    stop(
      "The cache directory could not be created: ",
      normalizePath(directory, winslash = "/", mustWork = FALSE),
      "\nCheck write permission and the total Windows path length."
    )
  }

  invisible(directory)
}

create_cache_directory(RATE_CACHE_DIR)
create_cache_directory(PERIOD_CHUNK_DIR)
create_cache_directory(OVERALL_CHUNK_DIR)

standardise_lambda_period_batch <- function(object, region_indices, weights) {
  object_dimensions <- dim(object)

  if (length(object_dimensions) != 4L ||
      object_dimensions[2L] != length(weights) ||
      object_dimensions[3L] != N_PERIODS) {
    stop(
      "Unexpected lambda batch dimensions: ",
      paste(object_dimensions, collapse = " × ")
    )
  }

  selected_object <- object[
    region_indices,
    seq_len(object_dimensions[2L]),
    seq_len(object_dimensions[3L]),
    seq_len(object_dimensions[4L]),
    drop = FALSE
  ]

  reordered <- aperm(selected_object, c(1L, 3L, 4L, 2L))
  group_matrix <- matrix(reordered, ncol = length(weights))
  standardised_vector <- drop(group_matrix %*% weights)
  standardised_array <- array(
    standardised_vector,
    dim = c(length(region_indices), N_PERIODS, object_dimensions[4L])
  )

  matrix(
    aperm(standardised_array, c(3L, 1L, 2L)),
    nrow = object_dimensions[4L],
    ncol = length(region_indices) * N_PERIODS
  )
}

process_period_chunk <- function(region_indices, chunk_id) {
  # Recreate the directory when this section is run independently.
  create_cache_directory(PERIOD_CHUNK_DIR)

  chunk_cache <- file.path(
    PERIOD_CHUNK_DIR,
    sprintf("period_%03d.rds", chunk_id)
  )

  if (USE_RATE_CACHE && file.exists(chunk_cache)) return(readRDS(chunk_cache))

  pooled_by_chain <- vector("list", N_CHAINS)

  for (chain_id in seq_len(N_CHAINS)) {
    file_table <- find_parameter_batch_files(
      CHAIN_DIRS[chain_id], "lambda", post_burn_only = TRUE
    )
    chain_batches <- vector("list", nrow(file_table))

    for (batch_id in seq_len(nrow(file_table))) {
      chain_batches[[batch_id]] <- standardise_lambda_period_batch(
        readRDS(file_table$file[batch_id]),
        region_indices,
        WHO_STANDARD_3
      )

      if (batch_id %% 50L == 0L || batch_id == nrow(file_table)) {
        message(
          "  period chunk ", chunk_id, ", chain ", chain_id, ": ",
          batch_id, "/", nrow(file_table), " lambda batches"
        )
      }
    }

    pooled_by_chain[[chain_id]] <- do.call(rbind, chain_batches)
  }

  chain_draw_counts <- vapply(pooled_by_chain, nrow, integer(1))
  chain_column_counts <- vapply(pooled_by_chain, ncol, integer(1))

  if (length(unique(chain_draw_counts)) != 1L ||
      length(unique(chain_column_counts)) != 1L) {
    stop("Chains differ in their number of retained draws or rate cells.")
  }

  pooled_draws <- do.call(rbind, pooled_by_chain) * RATE_MULTIPLIER
  posterior_quantiles <- matrixStats::colQuantiles(
    pooled_draws,
    probs = SUPPRESSION_QUANTILE_PROBABILITIES,
    useNames = FALSE,
    drop = FALSE
  )

  cell_grid <- expand.grid(
    code_muni = region_names[region_indices],
    period = time_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  events_matrix <- apply(
    Y[region_indices, , , drop = FALSE],
    c(1L, 3L),
    sum
  )
  population_matrix <- apply(
    N[region_indices, , , drop = FALSE],
    c(1L, 3L),
    sum
  )

  expected_cells <- length(region_indices) * N_PERIODS
  if (nrow(posterior_quantiles) != expected_cells ||
      ncol(posterior_quantiles) != 5L ||
      nrow(cell_grid) != expected_cells) {
    stop("Posterior quantile dimensions do not match the municipality-period grid.")
  }

  result <- tibble::tibble(
    code_muni = cell_grid$code_muni,
    period = cell_grid$period,
    group = "ASMR_WHO",
    posterior_median = posterior_quantiles[, 3L],
    ci_lower_95 = posterior_quantiles[, 1L],
    ci_upper_95 = posterior_quantiles[, 5L],
    ci_lower_90 = posterior_quantiles[, 2L],
    ci_upper_90 = posterior_quantiles[, 4L],
    events = as.vector(events_matrix),
    population_person_years = as.vector(population_matrix),
    average_annual_population = as.vector(population_matrix) / YEARS_PER_PERIOD,
    pooled_posterior_draws = nrow(pooled_draws)
  ) |>
    apply_suppression_rules()

  saveRDS(result, chunk_cache)
  result
}

region_chunks <- split(
  seq_len(N_MUNICIPALITIES),
  ceiling(seq_len(N_MUNICIPALITIES) / REGION_CHUNK_SIZE)
)

message("Combining posterior draws for period-specific age-standardised rates...")

combined_period_rates <- purrr::map2_dfr(
  region_chunks,
  seq_along(region_chunks),
  process_period_chunk
) |>
  dplyr::mutate(period = factor(period, levels = PERIODS)) |>
  dplyr::arrange(period, code_muni) |>
  dplyr::mutate(period = as.character(period))

readr::write_csv(
  combined_period_rates,
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_by_period.csv")
)

readr::write_csv(
  combined_period_rates |>
    dplyr::filter(!suppressed_95),
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_by_period_reliable_95.csv")
)

readr::write_csv(
  combined_period_rates |>
    dplyr::filter(!suppressed_90),
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_by_period_reliable_90.csv")
)

suppression_summary_by_period <- combined_period_rates |>
  dplyr::group_by(period) |>
  dplyr::summarise(
    total_estimates = dplyr::n(),
    reliable_estimates_95 = sum(!suppressed_95),
    suppressed_estimates_95 = sum(suppressed_95),
    reliable_pct_95 = 100 * mean(!suppressed_95),
    suppressed_low_precision_95 = sum(suppressed_low_precision_95),
    reliable_estimates_90 = sum(!suppressed_90),
    suppressed_estimates_90 = sum(suppressed_90),
    reliable_pct_90 = 100 * mean(!suppressed_90),
    suppressed_low_precision_90 = sum(suppressed_low_precision_90),
    suppressed_low_population = sum(suppressed_low_population),
    .groups = "drop"
  )

readr::write_csv(
  suppression_summary_by_period,
  file.path(OUTPUT_DIR, "suppression_summary_by_period_95_and_90.csv")
)

# Municipal age-standardised smoothed mortality rates by period ----

MAP_MUNICIPAL_BOUNDARY_YEAR <- 2022L
MAP_STATE_BOUNDARY_YEAR <- 2020L
MAP_SUPPRESSION_INTERVAL <- 0.90

PANEL_A_FILE <- file.path(
  "objects",
  "spatial",
  "map_rates_A.RData"
)

if (!file.exists(PANEL_A_FILE)) {
  stop(
    "Panel A map file not found: ",
    PANEL_A_FILE,
    "\nRun the descriptive spatial-analysis script before this section."
  )
}

panel_a_environment <- new.env(
  parent = globalenv()
)

loaded_panel_a_objects <- load(
  PANEL_A_FILE,
  envir = panel_a_environment
)

if (!"map_rates_A" %in% loaded_panel_a_objects) {
  stop(
    "The file ",
    PANEL_A_FILE,
    " does not contain an object named `map_rates_A`."
  )
}

map_rates_A <- panel_a_environment$map_rates_A

municipal_boundaries <- geobr::read_municipality(
  code_muni = "all",
  year = MAP_MUNICIPAL_BOUNDARY_YEAR,
  showProgress = FALSE
) |>
  dplyr::mutate(
    code_muni = stringr::str_sub(
      as.character(code_muni),
      1,
      6
    )
  )

state_boundaries <- geobr::read_state(
  code_state = "all",
  year = MAP_STATE_BOUNDARY_YEAR,
  showProgress = FALSE
)

municipality_count <- dplyr::n_distinct(
  municipal_boundaries$code_muni
)

if (
  municipality_count != N_MUNICIPALITIES ||
    !setequal(
      municipal_boundaries$code_muni,
      municipality_codes
    )
) {
  stop(
    "The municipal geometries do not match the municipalities ",
    "in the MSTCAR analytical dataset."
  )
}

if (MAP_SUPPRESSION_INTERVAL == 0.90) {
  map_rate_data <- combined_period_rates |>
    dplyr::select(
      code_muni,
      period,
      rate = posterior_median_suppressed_90
    )
} else if (MAP_SUPPRESSION_INTERVAL == 0.95) {
  map_rate_data <- combined_period_rates |>
    dplyr::select(
      code_muni,
      period,
      rate = posterior_median_suppressed_95
    )
} else {
  stop(
    "MAP_SUPPRESSION_INTERVAL must be 0.90 or 0.95."
  )
}

municipality_period_grid <- tidyr::crossing(
  code_muni = municipal_boundaries$code_muni,
  period = PERIODS
)

map_data <- municipal_boundaries |>
  dplyr::left_join(
    municipality_period_grid,
    by = "code_muni"
  ) |>
  dplyr::left_join(
    map_rate_data,
    by = c(
      "code_muni",
      "period"
    )
  )

rate_levels <- c(
  "Suppressed",
  "< 1.2",
  "1.2–1.9",
  "1.9–2.7",
  "2.7–4.2",
  "> 4.2"
)

rate_colours <- c(
  "Suppressed" = "#F1F1F1",
  "< 1.2" = "#FEF0D9",
  "1.2–1.9" = "#FDCC8A",
  "1.9–2.7" = "#FC8D59",
  "2.7–4.2" = "#E34A33",
  "> 4.2" = "#B30000"
)

map_data <- map_data |>
  dplyr::mutate(
    rate_category = dplyr::case_when(
      is.na(rate) ~ "Suppressed",
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

rate_legend_counts <- map_data |>
  sf::st_drop_geometry() |>
  dplyr::count(
    period,
    rate_category,
    .drop = FALSE,
    name = "n"
  ) |>
  dplyr::mutate(
    legend_label = paste0(
      rate_category,
      " [",
      scales::comma(
        n,
        big.mark = ",",
        decimal.mark = "."
      ),
      "]"
    )
  )

map_theme <- ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 18
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 18,
      hjust = 0
    ),
    legend.title = ggplot2::element_text(
      size = 12
    ),
    legend.text = ggplot2::element_text(
      size = 12
    ),
    legend.key.size = grid::unit(
      0.7,
      "cm"
    )
  )

make_period_rate_map <- function(
    period_value,
    panel_title = NULL
) {
  period_labels <- rate_legend_counts |>
    dplyr::filter(
      period == period_value
    ) |>
    dplyr::select(
      rate_category,
      legend_label
    ) |>
    tibble::deframe()

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = map_data |>
        dplyr::filter(
          period == period_value
        ),
      ggplot2::aes(
        fill = rate_category
      ),
      colour = NA
    ) +
    ggplot2::geom_sf(
      data = state_boundaries,
      ggplot2::aes(
        colour = "Federative units"
      ),
      fill = NA,
      linewidth = 0.8
    ) +
    ggplot2::scale_fill_manual(
      name = paste0(
        "Tuberculosis deaths per\n",
        "100,000 inhabitants\n",
        "[",
        scales::comma(
          municipality_count,
          big.mark = ","
        ),
        " municipalities]"
      ),
      values = rate_colours,
      breaks = rate_levels,
      labels = period_labels[rate_levels],
      drop = FALSE,
      guide = ggplot2::guide_legend(
        override.aes = list(
          colour = "black",
          linewidth = 0.5
        )
      )
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c(
        "Federative units" = "black"
      )
    ) +
    ggplot2::labs(
      title = panel_title,
      subtitle = stringr::str_replace_all(
        period_value,
        "-",
        "–"
      )
    ) +
    map_theme
}

map_2010_2014 <- make_period_rate_map(
  "2010-2014",
  panel_title = paste(
    "B) Age-standardised, spatiotemporally",
    "smoothed mortality rates"
  )
)

map_2015_2019 <- make_period_rate_map(
  "2015-2019"
)

map_2020_2024 <- make_period_rate_map(
  "2020-2024"
) +
  ggspatial::annotation_scale(
    location = "br",
    text_cex = 1,
    height = grid::unit(0.3, "cm"),
    pad_y = grid::unit(0.9, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    style = ggspatial::north_arrow_fancy_orienteering(),
    location = "bl",
    width = grid::unit(1.5, "cm"),
    height = grid::unit(2, "cm"),
    pad_x = grid::unit(1.5, "cm"),
    pad_y = grid::unit(1.5, "cm")
  )

map_rates_B <- (
  map_2010_2014 +
    map_2015_2019 +
    map_2020_2024
) +
  patchwork::plot_layout(
    ncol = 3
  )

final_map <- (
  map_rates_A /
    map_rates_B
) 

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Fig_2_asmr_rates.tiff"
  ),
  plot = final_map,
  device = "tiff",
  width = 26,
  height = 12,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Fig_2_asmr_rates.pdf"
  ),
  plot = final_map,
  device = "pdf",
  width = 26,
  height = 12,
  dpi = 300,
  bg = "white"
)

writeLines(
  paste(
    "Fig 2. Spatiotemporal distribution of age-standardised tuberculosis",
    "mortality rates before (A) and after spatiotemporal smoothing (B)",
    "in Brazil, 2010–2024."
  ),
  file.path(
    ARTICLE_FIGURE_DIR,
    "Fig_2_caption.txt"
  )
)

# Calculate age-standardised rates for the complete 2010–2024 period ----

standardise_lambda_overall_batch <- function(
    object,
    region_indices,
    population_array,
    age_weights
) {
  object_dimensions <- dim(object)

  if (length(object_dimensions) != 4L ||
      object_dimensions[2L] != length(age_weights) ||
      object_dimensions[3L] != N_PERIODS) {
    stop(
      "Unexpected lambda batch dimensions: ",
      paste(object_dimensions, collapse = " × ")
    )
  }

  lambda_subset <- object[region_indices, , , , drop = FALSE]
  population_subset <- population_array[region_indices, , , drop = FALSE]
  total_population_by_age <- apply(population_subset, c(1L, 2L), sum)

  standardised_draws <- matrix(
    0,
    nrow = object_dimensions[4L],
    ncol = length(region_indices)
  )

  for (local_region_id in seq_along(region_indices)) {
    for (group_id in seq_len(object_dimensions[2L])) {
      denominator <- total_population_by_age[local_region_id, group_id]
      if (!is.finite(denominator) || denominator <= 0) next

      time_weights <- population_subset[local_region_id, group_id, ] / denominator
      lambda_time_draws <- lambda_subset[local_region_id, group_id, , , drop = TRUE]

      if (is.null(dim(lambda_time_draws))) {
        lambda_time_draws <- matrix(lambda_time_draws, nrow = object_dimensions[3L])
      }

      weighted_age_specific_rate <- colSums(
        sweep(lambda_time_draws, MARGIN = 1L, STATS = time_weights, FUN = "*")
      )

      standardised_draws[, local_region_id] <-
        standardised_draws[, local_region_id] +
        age_weights[group_id] * weighted_age_specific_rate
    }
  }

  standardised_draws
}

process_overall_chunk <- function(region_indices, chunk_id) {
  # Recreate the directory when this section is run independently.
  create_cache_directory(OVERALL_CHUNK_DIR)

  chunk_cache <- file.path(
    OVERALL_CHUNK_DIR,
    sprintf("overall_%03d.rds", chunk_id)
  )

  if (USE_RATE_CACHE && file.exists(chunk_cache)) return(readRDS(chunk_cache))

  pooled_by_chain <- vector("list", N_CHAINS)

  for (chain_id in seq_len(N_CHAINS)) {
    file_table <- find_parameter_batch_files(
      CHAIN_DIRS[chain_id], "lambda", post_burn_only = TRUE
    )
    chain_batches <- vector("list", nrow(file_table))

    for (batch_id in seq_len(nrow(file_table))) {
      chain_batches[[batch_id]] <- standardise_lambda_overall_batch(
        readRDS(file_table$file[batch_id]),
        region_indices,
        N,
        WHO_STANDARD_3
      )

      if (batch_id %% 50L == 0L || batch_id == nrow(file_table)) {
        message(
          "  overall chunk ", chunk_id, ", chain ", chain_id, ": ",
          batch_id, "/", nrow(file_table), " lambda batches"
        )
      }
    }

    pooled_by_chain[[chain_id]] <- do.call(rbind, chain_batches)
  }

  chain_draw_counts <- vapply(pooled_by_chain, nrow, integer(1))
  if (length(unique(chain_draw_counts)) != 1L) {
    stop("Chains differ in their number of retained draws.")
  }

  pooled_draws <- do.call(rbind, pooled_by_chain) * RATE_MULTIPLIER
  posterior_quantiles <- matrixStats::colQuantiles(
    pooled_draws,
    probs = SUPPRESSION_QUANTILE_PROBABILITIES,
    useNames = FALSE,
    drop = FALSE
  )

  total_events <- apply(Y[region_indices, , , drop = FALSE], 1L, sum)
  total_population <- apply(N[region_indices, , , drop = FALSE], 1L, sum)

  result <- tibble::tibble(
    code_muni = region_names[region_indices],
    period = "2010-2024",
    group = "ASMR_WHO",
    posterior_median = posterior_quantiles[, 3L],
    ci_lower_95 = posterior_quantiles[, 1L],
    ci_upper_95 = posterior_quantiles[, 5L],
    ci_lower_90 = posterior_quantiles[, 2L],
    ci_upper_90 = posterior_quantiles[, 4L],
    events = as.numeric(total_events),
    population_person_years = as.numeric(total_population),
    average_annual_population = as.numeric(total_population) / length(STUDY_YEARS),
    pooled_posterior_draws = nrow(pooled_draws)
  ) |>
    apply_suppression_rules()

  saveRDS(result, chunk_cache)
  result
}

message("Combining posterior draws for age-standardised rates in 2010–2024...")

overall_rates <- purrr::map2_dfr(
  region_chunks,
  seq_along(region_chunks),
  process_overall_chunk
) |>
  dplyr::arrange(code_muni)

readr::write_csv(
  overall_rates,
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_2010_2024.csv")
)

readr::write_csv(
  overall_rates |>
    dplyr::filter(!suppressed_95),
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_2010_2024_reliable_95.csv")
)

readr::write_csv(
  overall_rates |>
    dplyr::filter(!suppressed_90),
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_2010_2024_reliable_90.csv")
)

overall_suppression_summary <- overall_rates |>
  dplyr::summarise(
    total_estimates = dplyr::n(),
    reliable_estimates_95 = sum(!suppressed_95),
    suppressed_estimates_95 = sum(suppressed_95),
    reliable_pct_95 = 100 * mean(!suppressed_95),
    suppressed_low_precision_95 = sum(suppressed_low_precision_95),
    reliable_estimates_90 = sum(!suppressed_90),
    suppressed_estimates_90 = sum(suppressed_90),
    reliable_pct_90 = 100 * mean(!suppressed_90),
    suppressed_low_precision_90 = sum(suppressed_low_precision_90),
    suppressed_low_population = sum(suppressed_low_population)
  )

readr::write_csv(
  overall_suppression_summary,
  file.path(OUTPUT_DIR, "suppression_summary_2010_2024_95_and_90.csv")
)

# Final checks and reproducibility files ----

configuration_details <- tibble::tibble(
  metric = c(
    "Number of chains",
    "Iterations per chain",
    "Burn-in per chain",
    "Expected batches per chain",
    "Burn-in batches per chain",
    "Saved draws per batch",
    "Stored draws per chain before burn-in exclusion",
    "Retained stored draws per chain",
    "Total retained draws",
    "Age groups",
    "Periods",
    "Population suppression threshold",
    "Precision intervals used for suppression"
  ),
  value = c(
    N_CHAINS,
    ITERATIONS,
    BURN_IN,
    EXPECTED_BATCHES,
    BURN_BATCHES,
    draws_per_batch,
    full_saved_draws_per_chain,
    retained_saved_draws_per_chain,
    total_retained_draws,
    paste(AGE_GROUPS, collapse = "; "),
    paste(PERIODS, collapse = "; "),
    POPULATION_SUPPRESSION_THRESHOLD,
    paste0(100 * SUPPRESSION_CREDIBLE_INTERVALS, "%", collapse = "; ")
  )
)

readr::write_csv(
  configuration_details,
  file.path(OUTPUT_DIR, "mcmc_configuration_details.csv")
)

max_selected_rhat <- max(scalar_diagnostics$rhat_point, na.rm = TRUE)
min_selected_ess <- min(scalar_diagnostics$ess_combined, na.rm = TRUE)

if (is.finite(max_selected_rhat) && max_selected_rhat > RHAT_THRESHOLD) {
  warning(
    "At least one selected scalar parameter has R-hat > ", RHAT_THRESHOLD,
    ". Combining chains does not correct non-convergence."
  )
}

if (is.finite(min_selected_ess) && min_selected_ess < ESS_THRESHOLD) {
  warning("At least one selected scalar parameter has ESS < ", ESS_THRESHOLD, ".")
}

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "session_info.txt")
)

print(model_configuration)
print(table_s3)
print(table_s4)
print(suppression_summary_by_period)
print(overall_suppression_summary)

message(
  "MSTCAR processing completed. Main outputs:\n",
  "  - Table S2: ", file.path(TABLE_DIR, "Table_S2_MSTCAR_configuration.docx"), "\n",
  "  - Table S3: ", file.path(TABLE_DIR, "Table_S3_MSTCAR_diagnostics.docx"), "\n",
  "  - Table S4: ", file.path(TABLE_DIR, "Table_S4_chain_specific_rate_agreement.docx"), "\n",
  "  - Fig S3 traceplots: ", file.path(FIGURE_DIR, "Fig_S3_MSTCAR_traceplots.pdf"), "\n",
  "  - Period estimates with 95% and 90% suppression: ",
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_by_period.csv"), "\n",
  "  - Overall estimates with 95% and 90% suppression: ",
  file.path(OUTPUT_DIR, "combined_age_standardised_rates_2010_2024.csv"), "\n",
  "  - Fig 2: ", file.path(ARTICLE_FIGURE_DIR, "Fig_2.tiff")
)

