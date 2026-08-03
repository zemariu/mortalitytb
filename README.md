# mortalitytb

[![R](https://img.shields.io/badge/R-analysis-276DC3?logo=r&logoColor=white)](https://www.r-project.org/) [![Bayesian modelling](https://img.shields.io/badge/Bayesian-INLA%20%7C%20MSTCAR-7B2CBF)](#software-and-packages) [![Study period](https://img.shields.io/badge/study%20period-2010--2024-1F7A8C)](#study-overview) [![Geographic unit](https://img.shields.io/badge/unit-5%2C570%20municipalities-2E8B57)](#study-overview) [![Manuscript](https://img.shields.io/badge/manuscript-under%20review-orange)](#citation) [![Licence](https://img.shields.io/badge/licence-MIT-green)](LICENSE)

## Overview

This repository contains the analysis code supporting the manuscript:

> **Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024: a Bayesian spatiotemporal modelling study**

The study aimed to characterise the temporal and spatial distribution of tuberculosis (TB) mortality across Brazilian municipalities between 2010 and 2024, examine contextual factors associated with mortality, and develop a Bayesian municipal prioritisation framework integrating risk magnitude, temporal dynamics, spatial patterns, posterior uncertainty, and robustness to support public health decision-making.

The repository provides the complete analytical workflow, from preparation of mortality and population data to descriptive epidemiology, spatial analysis, Bayesian spatiotemporal modelling, classification-stability assessment, and construction of the final municipal TB mortality priority index.

> **Manuscript status:** the manuscript is currently under review at *Scientific Reports*. The citation and repository documentation will be updated after completion of the peer-review process.

------------------------------------------------------------------------

## Interactive dashboard

An interactive web dashboard was developed to facilitate the exploration and visualisation of the main study findings at the municipal level.

The dashboard allows users to examine:

- final municipal priority classifications;
- adjusted priority scores;
- estimated relative risks and posterior probabilities;
- municipality-specific temporal trends;
- hotspot persistence and emergence;
- spatial–temporal regime stability;
- municipality-level estimates and 95% credible intervals.

The dashboard includes geographic filters by region and state, interactive municipal maps, evidence plots, searchable municipal results, and options to download the filtered data.

The dashboard is publicly available at: [Open the Municipal Tuberculosis Mortality Priority Dashboard](https://cenariostb.unb.br/mortalitytb/)

The dashboard is intended to support the interpretation and dissemination of the study findings. It should not replace local epidemiological assessment, data-quality evaluation, or public-health decision-making.

------------------------------------------------------------------------

## Study overview

- **Design:** national ecological study.

- **Geographic unit:** 5,570 Brazilian municipalities.

- **Study period:** 2010–2024.

- **Outcome:** deaths with TB as the underlying cause.

- **Analytical unit:** municipality-year for the principal Bayesian models.

- **Spatial framework:** municipal adjacency and small-area spatial smoothing.

- **Main modelling framework:** Bayesian spatiotemporal count models fitted using INLA.

- **Additional smoothing analysis:** multivariate spatiotemporal conditional autoregressive modelling across age groups and periods.

- **Primary output:** municipal priority classification integrating central risk, temporal dynamics, spatial persistence, posterior evidence, uncertainty, and robustness.

## Repository structure

The recommended GitHub repository structure is:

``` text
mortalitytb/
├── README.md
├── LICENSE
├── CITATION.cff
├──.gitignore
│
│── scripts/
│   ├── 01_descriptive_analysis.R
│   ├── 02_spatial_analysis.R
│   ├── 03_mstcar_spatiotemporal_analysis.R
│   ├── 04_bayesian_spatiotemporal_analysis.R
│   ├── 05_cluster_trend_stability.R
│   ├── 06_tb_priority_index.R
│
├── data/
│   ├── SIM_TB.csv
│   ├── pop_br.csv    # Not included; download separately
│   ├── pop_race.csv
│   ├── municipal_indicators.csv
│   └── processed/
│
├── objects/
│   ├── spatial/
│   ├── mstcar/
│   └── bayesian_models/
│
├── models/
│   └── mstcar/
│
└── results/
    ├── descriptive/
    │   ├── figures/
    │   └── tables/
    ├── spatial/
    │   ├── figures/
    │   └── tables/
    ├── mstcar/
    │   ├── figures/
    │   └── tables/
    ├── spatiotemporal/
    │   ├── outputs/
    │   ├── figures/
    │   └── tables/
    └── priority_index/
        ├── outputs/
        ├── figures/
        └── tables/
```

The scripts create most output directories automatically. Raw data, fitted models, intermediate objects, and large output files do not need to be version-controlled unless their public redistribution is permitted and scientifically justified.

------------------------------------------------------------------------

## Analytical workflow

Run the scripts in numerical order:

``` text
01_descriptive_analysis.R
        ↓
02_spatial_analysis.R
        ↓
03_mstcar_spatiotemporal_analysis.R
        ↓
04_bayesian_spatiotemporal_analysis.R
        ↓
05_cluster_trend_stability.R
        ↓
06_tb_priority_index.R
```
The scripts are not fully independent. Later stages use processed data or model outputs generated by earlier stages.

------------------------------------------------------------------------

## Script overview

### `01_descriptive_analysis.R`

Prepares mortality and population data and produces the descriptive epidemiological indicators used in subsequent analyses.

Main tasks include:

- importing TB mortality records;
- standardising six-digit municipality codes;
- recoding age, sex, race/ethnicity, marital status, education, place of death, and anatomical TB classification;
- defining five-year study periods;
- combining mortality and population denominators;
- preparing age-, sex-, race-, municipality-, and year-specific analytical datasets;
- producing descriptive tables and temporal figures;
- exporting the principal processed dataset used by later scripts.

Principal inputs:

``` text
data/SIM_TB.csv
data/pop_br.csv
data/pop_race.csv
```

Principal outputs:

``` text
data/processed/tb_mortality.csv
data/processed/tb_mortality_race.csv
data/processed/population_by_race.csv
data/processed/joinpoint_mortality_rates.csv
results/descriptive/
```

### `02_spatial_analysis.R`

Calculates municipal mortality rates, evaluates global spatial autocorrelation, and produces descriptive maps.

Main tasks include:

- calculating age-standardised municipal mortality rates by five-year period;
- obtaining municipal and state boundaries;
- constructing spatial neighbourhood objects;
- estimating global Moran's I;
- producing maps of Brazil and municipal mortality patterns.

Principal input:

``` text
data/processed/tb_mortality.csv
```

Principal outputs:

``` text
data/processed/municipal_asmr_by_period.csv
objects/
results/spatial/
```

### `03_mstcar_spatiotemporal_analysis.R`

Performs multivariate spatiotemporal smoothing, convergence diagnostics, between-chain agreement assessment, and estimate-suppression procedures.

Main tasks include:

- aggregating deaths and population into three age groups and three five-year periods;
- constructing a symmetric six-nearest-neighbour graph;
- fitting or loading three independent MSTCAR chains;
- assessing convergence and effective sample size;
- comparing chain-specific posterior mortality-rate estimates;
- generating age-standardised smoothed mortality rates;
- applying reliability and suppression criteria;
- exporting diagnostics, tables, figures, and municipal estimates.

Principal input:

``` text
data/processed/tb_mortality.csv
```

Principal outputs:

``` text
models/mstcar/
objects/mstcar/
results/mstcar/
```

Because this stage can be computationally intensive, fitted chains and cached parameter summaries should normally be stored locally rather than committed to GitHub.

### `04_bayesian_spatiotemporal_analysis.R`

Fits and compares Bayesian spatiotemporal models and estimates adjusted municipal risks and temporal trends.

Main tasks include:

- calculating observed and expected deaths;
- describing zero-death patterns;
- constructing a connected municipal adjacency graph;
- fitting and comparing Poisson, negative-binomial, and zero-inflated model specifications;
- fitting BYM2 spatial effects, temporal effects, municipality-specific trends, and space–time interactions;
- importing and evaluating contextual covariates;
- performing univariable and multivariable analyses;
- generating adjusted and unadjusted municipal classifications;
- exporting model-comparison statistics, fixed effects, hyperparameters, diagnostics, maps, and analytical datasets.

Principal inputs:

``` text
data/processed/tb_mortality.csv
data/municipal_indicators.csv
```

Principal outputs:

``` text
objects/spatial/
objects/bayesian_models/
results/spatiotemporal/
```

### `05_cluster_trend_stability.R`

Compares municipal classifications before and after covariate adjustment and evaluates spatial-cluster, temporal-trend, and complete-regime stability.

Main tasks include:

- aligning adjusted and unadjusted classifications;
- constructing transition matrices;
- calculating national, regional, state-level, and regime-specific persistence measures;
- evaluating hotspot emergence and resolution;
- calculating entropy-based stability measures;
- producing stability tables and figures;
- generating the enriched municipal stability dataset used by the priority-index script.

Principal inputs:

``` text
results/spatiotemporal/outputs/spatial_temporal_classification_before_adjustment.csv
results/spatiotemporal/outputs/spatial_temporal_classification_after_adjustment.csv
```

Principal outputs:

``` text
results/spatiotemporal/outputs/municipal_classification_stability.csv
results/spatiotemporal/outputs/municipal_classification_stability.rds
results/spatiotemporal/tables/
results/spatiotemporal/figures/
```

### `06_tb_priority_index.R`

Constructs the municipal TB mortality priority index, incorporates posterior uncertainty and robustness, and generates final priority classifications.

Main tasks include:

- assigning component scores for spatial classification, hotspot persistence or emergence, relative risk, exceedance probability, temporal trend, temporal evidence, regime stability, and state context;
- rescaling the central priority score to 0–100;
- constructing an explicit uncertainty score;
- applying an uncertainty penalty;
- simulating score and rank uncertainty;
- assessing sensitivity across alternative scenarios;
- assigning final priority and operational categories;
- generating the final map, municipal decision tables, and supplementary public-use dataset.

Principal inputs:

``` text
results/spatiotemporal/outputs/municipal_classification_stability.rds
results/mstcar/combined_age_standardised_rates_2010_2024.csv
```

Principal outputs:

``` text
results/priority_index/outputs/municipal_priority_index.csv
results/priority_index/outputs/municipal_priority_index.rds
results/priority_index/outputs/Supplementary_Table_S1.csv
results/priority_index/outputs/Supplementary_Table_S1.xlsx
results/priority_index/tables/
results/priority_index/figures/
```

------------------------------------------------------------------------

## Input data

The analysis expects the following source files.

| File | Purpose |
|------------------------------------|------------------------------------|
| `data/SIM_TB.csv` | De-identified TB mortality records used to construct the mortality numerator. |
| `data/pop_br.csv` | Population denominators by municipality, year, age group, and sex.  |
| `data/pop_race.csv` | Population denominators or population composition by race/ethnicity. |
| `data/municipal_indicators.csv` | Municipality-level contextual covariates used in the adjusted Bayesian models. |

Municipality identifiers are handled as **six-digit codes** in the analytical datasets. Geometry files may contain seven-digit IBGE codes; the scripts retain the first six digits when matching spatial and analytical data.

### Population denominator data

The population denominator file used in the analyses:

```text
data/pop_br.csv
```
is not included in this repository because its size exceeds GitHub's 100 MB file-size limit.

The population estimates can be obtained from the Brazilian Ministry of Health population data portal:

[Population data — Brazilian Ministry of Health](https://www.gov.br/saude/pt-br/composicao/seidigi/demas/dados-populacionais).

After downloading and preparing the population data, save the resulting file as:

```
data/pop_br.csv
```

The file must contain the municipality-, year-, sex-, and age-specific population denominators required by `01_descriptive_analysis.R.` Variable names and formats must match those referenced in the script.

The absence of `pop_br.csv` from the GitHub repository does not affect access to the analysis code, but the file is required to reproduce the mortality-rate calculations and subsequent analyses.


------------------------------------------------------------------------

## De-identified public-use data

This repository is designed for reproducible analysis of de-identified, public-use, or aggregated data.

The analytical workflow does not require names, national identification numbers, addresses, or other direct personal identifiers. Nevertheless:

- raw source data should only be redistributed when permitted by the original data provider;
- individual-level mortality extracts should not be committed to a public GitHub repository unless their public redistribution is explicitly authorised;
- public repository releases should preferentially contain aggregated municipality-level analytical data;
- small-cell disclosure and data-provider requirements should be reviewed before publication;
- users remain responsible for complying with the terms of use of the Brazilian Mortality Information System (SIM), population datasets (IBGE), contextual indicators, and spatial data providers.

The final supplementary municipal dataset generated by `06_tb_priority_index.R` is intended to contain municipality-level analytical results rather than identifiable individual records.

The `data/README.md` file should describe how authorised users can obtain the original source data and how the local files must be named, rather than redistributing restricted raw files.

------------------------------------------------------------------------

## Variable codebook

The following table summarises the principal variables used across the workflow. It is not a substitute for the complete codebook included with the final supplementary dataset.

### Core mortality and population variables

| Variable | Description | Typical level or unit |
|------------------------|------------------------|------------------------|
| `code_muni` | Six-digit municipality identifier used across analytical datasets. | Municipality |
| `year` | Calendar year. | 2010–2024 |
| `period` | Five-year analytical period. | 2010–2014, 2015–2019, 2020–2024 |
| `deaths` | Number of TB deaths. | Count |
| `population` | Population or person-years at risk. | Count |
| `sex` | Sex category used in standardisation and descriptive analyses. | Male, Female |
| `age_group` | Five-year age group. | 0–4 to 80+ |
| `age_group_3` | Broad age group used in the MSTCAR analysis. | 0–19, 20–59, 60+ |
| `race` | Race/ethnicity category used in descriptive analyses. | White, Black, Brown, Yellow, Indigenous |
| `region` | Brazilian macroregion. | North, Northeast, Southeast, South, Central-West |

### Bayesian modelling variables

| Variable | Description | Interpretation |
|------------------------|------------------------|------------------------|
| `O` | Observed municipal TB deaths. | Count |
| `E` | Expected municipal TB deaths based on national age- and sex-specific rates. | Positive continuous value |
| `SMR` | Standardised mortality ratio, calculated as `O / E`. | Values above 1 indicate excess mortality |
| `RR_mean` or `RR_after` | Posterior mean adjusted relative risk. | Values above 1 indicate risk above the reference |
| `RR_lower_95` | Lower bound of the 95% credible interval for relative risk. | Posterior uncertainty |
| `RR_upper_95` | Upper bound of the 95% credible interval for relative risk. | Posterior uncertainty |
| `Pr_RR_gt1` or `Pr_RR_after` | Posterior probability that relative risk exceeds 1. | 0–1 |
| `slope_mean` or `slope_after` | Posterior mean municipality-specific temporal coefficient. | Positive values indicate increasing risk |
| `slope_lower_95` | Lower bound of the 95% credible interval for the temporal coefficient. | Posterior uncertainty |
| `slope_upper_95` | Upper bound of the 95% credible interval for the temporal coefficient. | Posterior uncertainty |
| `Pr_slope_gt0` or `Pr_slope_after` | Posterior probability that the temporal coefficient is positive. | 0–1 |

### Classification and stability variables

| Variable | Description |
|------------------------------------|------------------------------------|
| `cluster_before` | Spatial classification before covariate adjustment. |
| `cluster_after` | Spatial classification after covariate adjustment. |
| `trend_before` | Temporal classification before covariate adjustment. |
| `trend_after` | Temporal classification after covariate adjustment. |
| `regime_before` | Combined spatial–temporal classification before adjustment. |
| `regime_after` | Combined spatial–temporal classification after adjustment. |
| `hotspot_status` | Persistence, emergence, resolution, or absence of hotspot status. |
| `cluster_persistent` | Indicator that the spatial classification persisted after adjustment. |
| `trend_persistent` | Indicator that the temporal classification persisted after adjustment. |
| `strict_persistent` | Indicator that the complete spatial–temporal regime persisted. |
| `regime_entropy_normalized` | Normalised entropy measure representing regime uncertainty or instability. |
| `state_strict_persistence` | State-level proportion of municipalities retaining the complete regime. |

### Priority-index variables

| Variable | Description |
|------------------------------------|------------------------------------|
| `central_raw_priority_score` | Sum of the central priority components before rescaling. |
| `central_priority_score` | Central priority score rescaled to 0–100. |
| `uncertainty_score` | Composite score representing posterior and structural uncertainty. |
| `uncertainty_penalty` | Penalty subtracted from the central score. |
| `adjusted_priority_score` | Priority score after accounting for uncertainty. |
| `median_rank` | Median simulated municipal rank. |
| `rank_p025` | 2.5th percentile of the simulated rank distribution. |
| `rank_p975` | 97.5th percentile of the simulated rank distribution. |
| `Pr_score_very_high` | Probability that a municipality meets the very-high-score criterion. |
| `final_priority_class` | Final municipal priority category. |
| `operational_priority` | Public-health action category linked to the final priority class. |

The complete codebook is provided in:

``` text
results/priority_index/outputs/Supplementary_Table_S1.xlsx
```

------------------------------------------------------------------------

## Priority framework

The central priority score combines evidence from:

1.  current spatial classification;
2.  hotspot persistence or emergence;
3.  adjusted relative-risk magnitude;
4.  posterior probability that relative risk exceeds 1;
5.  municipality-specific temporal trend;
6.  posterior probability that the temporal trend is positive;
7.  stability or entropy of the spatial–temporal regime;
8.  persistence within the corresponding federative unit.

The central score is rescaled to 0–100. An explicit uncertainty score is then constructed from credible-interval widths, posterior probabilities, and structural stability measures. The final adjusted score subtracts an uncertainty penalty from the central score.

The final outputs distinguish evidence magnitude from evidence certainty and include robustness and sensitivity analyses rather than relying solely on a deterministic rank.

------------------------------------------------------------------------

## Software and packages 

The analyses were written in R. The exact software versions used for each stage are recorded by the `sessionInfo()` files generated by the scripts.

### Main R packages

| Purpose | Packages |
|------------------------------------|------------------------------------|
| Data management | `tidyverse`, `dplyr`, `tidyr`, `tibble`, `purrr`, `readr`, `stringr` |
| Tables and reporting | `gtsummary`, `gt`, `flextable`, `officer`, `writexl` |
| Visualisation | `ggplot2`, `patchwork`, `ggrepel`, `scales` |
| Spatial data | `sf`, `spdep`, `geobr`, `ggspatial`, `rnaturalearth`, `rnaturalearthdata` |
| Bayesian modelling | `INLA`, `RSTr`, `coda`, `matrixStats` |

### Installing CRAN packages

``` r
install.packages(c(
  "tidyverse",
  "gtsummary",
  "flextable",
  "officer",
  "patchwork",
  "sf",
  "spdep",
  "geobr",
  "ggspatial",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggrepel",
  "gt",
  "writexl",
  "coda",
  "matrixStats"
))
```

`INLA` and `RSTr` may require installation outside the standard CRAN workflow. Install them using the instructions provided by their official maintainers.

For reproducibility, use a project-level dependency manager such as `renv`:

``` r
install.packages("renv")
renv::init()
renv::snapshot()
```

A generated `renv.lock` file can then be committed to the repository.

------------------------------------------------------------------------

## Computational requirements

The descriptive and conventional spatial analyses can generally be run on a standard desktop computer.

The following stages may require substantially more time, memory, and storage:

- fitting multiple Bayesian model specifications with INLA;
- fitting three independent MSTCAR chains;
- reading and summarising batched posterior parameter files;
- posterior simulation of score and rank uncertainty;
- generating high-resolution national municipal maps.

Recommended practice:

- store fitted models and intermediate posterior objects locally;
- avoid refitting models when valid cached objects are available;
- use fixed random seeds where specified;
- preserve the generated `sessionInfo()` files;
- document any changes to iterations, burn-in, priors, model families, or sensitivity scenarios.

------------------------------------------------------------------------

## Reproducing the analysis

1.  Clone the repository:

``` bash
git clone https://github.com/<USERNAME>/mortalitytb.git
cd mortalitytb
```

2.  Create the required directory structure:

``` r
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("objects", recursive = TRUE, showWarnings = FALSE)
dir.create("models", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
```

3. Obtain the source datasets from their original providers and place the locally prepared files in `data/`. The `pop_br.csv` file is not included because it exceeds GitHub's file-size limit and must be obtained from the [Brazilian Ministry of Health population data portal](https://www.gov.br/saude/pt-br/composicao/seidigi/demas/dados-populacionais).

4.  Install the required R packages.

5.  Check file paths and analysis-control parameters near the beginning of each script.

6.  Run the scripts in numerical order.

``` r
source("01_descriptive_analysis.R")
source("02_spatial_analysis.R")
source("03_mstcar_spatiotemporal_analysis.R")
source("04_bayesian_spatiotemporal_analysis.R")
source("05_cluster_trend_stability.R")
source("06_tb_priority_index.R")
```

For computationally intensive stages, running scripts interactively or through a batch scheduler is preferable to sourcing the complete pipeline in a single R session.

------------------------------------------------------------------------

## Reproducibility safeguards

The scripts include several checks intended to protect the analytical workflow, including:

- municipality-code harmonisation;
- validation of required columns;
- detection of duplicated municipalities;
- comparison of municipality sets across data sources;
- preservation checks after temporal and age-group aggregation;
- validation of connected spatial graphs;
- validation of model files and chain dimensions;
- checks for missing MCMC batches;
- convergence and effective-sample-size diagnostics;
- comparison of chain-specific posterior mortality-rate estimates;
- explicit output of model configuration and session information;
- sensitivity analyses for the priority framework.

Users who modify the workflow should retain or extend these checks.

------------------------------------------------------------------------

## Citation 

The manuscript is currently under review:

> **da Silva JMN, et al.** *Identifying priority areas for tuberculosis mortality in Brazil, 2010–2024: a Bayesian spatiotemporal modelling study.* Manuscript under review at *Scientific Reports*.

The citation will be updated after completion of peer review and publication.

Users of the code or derived municipal outputs should cite both the repository and the associated article when the final citation becomes available.

------------------------------------------------------------------------

## Contact

For questions, reproducibility issues, or technical problems, please contact:

**José Mário Nunes da Silva, PhD**\
Email: [zemariu\@hotmail.com](mailto:zemariu@hotmail.com)

When reporting a problem, please include:

- the script name and relevant line number;
- the complete R error message;
- your R version and operating system;
- package versions or the generated `sessionInfo()` file;
- the steps needed to reproduce the issue.

------------------------------------------------------------------------

## Licence

The source code in this repository may be distributed under the [MIT License](LICENSE).

The MIT License applies to the original repository code and documentation. It does not override the terms of use, attribution requirements, or redistribution restrictions of the original mortality, population, contextual, or spatial datasets.

------------------------------------------------------------------------

## Disclaimer

The priority framework is intended to support population-level public-health planning and technical review. It should not be interpreted as a substitute for local epidemiological investigation, programme knowledge, data-quality assessment, or resource-allocation deliberation.

Municipal rankings and categories depend on the study period, analytical specifications, input data, posterior uncertainty, and decision rules documented in the manuscript and code.
