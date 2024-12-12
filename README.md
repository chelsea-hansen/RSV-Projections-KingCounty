README
================

#### Scenario projections of RSV hospitalizations averted due to new immunizations in King County, Washington, October 2023 to May 2025

This work was done in collaboration with Public Health - Seattle & King
County as part of a CSTE/CDC - supported initiative, “Development of
forecast, analytic, and visualization tools to improve outbreak response
and support public health decision making.”

Please see our manuscript and documentation for our R Package
[R.Scenario.Vax](https://chelsea-hansen.github.io/R.Scenario.Vax/) For
more information regarding the model structure.

We used RSV hospitalization data for the greater Seattle area from two
sources:

1.  For age groups ≥10 years, we used RSV-diagnosed hospital admissions
    among King County residents or at King County hospitals from
    [Washington State’s syndromic surveillance
    platform](https://doh.wa.gov/public-health-provider-resources/healthcare-professions-and-facilities/data-exchange/syndromic-surveillance-rhino).
2.  For children \<10 years, we used RSV hospitalizations with a
    positive RSV test (during or preceding admission) seen at [Seattle
    Children’s Hospital](https://www.seattlechildrens.org/) for the same
    period.

**Note:** Raw data are not publicly available. Aggregated datasets have
been provided here to reproduce our study results. Researchers wishing
to access Washington’s syndromic surveillance data should contact
[RHINO@doh.wa.gov](RHINO@doh.wa.gov). Researchers wishing to access
aggregate data from Seattle Children’s Hospital should contact
[amanda.adler@seattlechildrens.org](amanda.adler@seattlechildrens.org)
for information on data sharing agreements.

## Data

#### Calibration datasets and parameter values

- `time_series_public.rds`: Time series of weekly RSV hospitalization
  rate per 100,000 population in King County, WA. These data have been
  adjusted to account for changes in RSV testing during the COVID-19
  pandemic, and a 3-week moving average has been applied. For more
  information, see the eMethods from the corresponding manuscript.
- `age_distribution_public.rds`: Age distribution of RSV
  hospitalizations in King County, WA.
- `weekly_immunizations_public.rds`: Estimated weekly and cumulative
  immunization coverage for:
  - Infants born to vaccinated mothers.
  - Infants receiving birth doses of nirsevimab.
  - Infants receiving catch-up doses of nirsevimab.
  - Adults ≥60 years receiving RSV vaccines.
- `fixed_parameters.rds`: Fixed parameter values for running the model.
  These values are taken from previous modeling studies and publicly
  available census data:
  - `[[1]]`: Fixed parameter values (see eTable 2 of the manuscript for
    full list).
  - `[[2]]`: Starting values for model compartments, based on publicly
    available age-stratified census data (vector form).
  - `[[3]]`: Starting values for model compartments, based on publicly
    available age-stratified census data (matrix form).
- `fitted_paramaters_part1.rds`: Fitted parameter values for the
  baseline transmission rate, timing and amplitude of seasonal forcing,
  and age-specific reporting rates for RSV hospitalizations, fit to data
  from 2017 to 2020. These values are fitted using maximum likelihood
  estimation (`R/MLE.R`).
- `fitted_parameters_part2.rds`: Fitted parameters for contact
  reductions during the COVID-19 period (see eMethods).
- `fitted_parameters_100.rds`: 100 parameter sets based on Latin
  Hypercube Sampling from the upper and lower bounds of the fitted
  parameter values (after adding noise). See eMethods.

## R

#### R Scripts for fitting model parameters and running scenario projections.

- `MLE.R`: R function for fitting parameter values using Maximum
  Likelihood Estimation.
- `MSIRS_immunization_dynamics.R`: Differential equations describing the
  transmission dynamics of RSV.
- `projection_function.R`: Function to produce scenario projections.
- `Projections_2023_24.R`: Application of `projection_function.R` to
  scenarios for the 2023–2024 RSV season.
- `Projections_2024_25.R`: Application of `projection_function.R` to
  scenarios for the 2024–2025 RSV season.

## Results

#### 100 Trajectories of RSV Hospitalizations for each Scenario

- `results_23-24_100replicates`: Scenario projections for the 2023–2024
  RSV season.
- `results_24-25_100replicates`: Scenario projections for the 2024–2025
  RSV season.

#### Contact

For questions regarding this repository or the manuscript please contact
[chelsea.hansen@nih.gov](chelsea.hansen@nih.gov)
