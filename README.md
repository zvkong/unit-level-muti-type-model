The code uses Illinois ACS PUMS data stored in `source data/IL23.csv`.

The responses are:
- Gaussian response: scaled log income
- Binary response: poverty indicator

The spatial effect is represented using basis functions built from the positive eigenvectors of the PUMA adjacency matrix.

The analysis is done at the PUMA level, with poststratification cells defined by:
- `PUMA`
- `SEX`
- `BACH`

## Requirements

The code is written in R.

All required packages are handled by:

```r
source("packages.r")

This script installs any missing packages and then loads them.

Before running

Set the working directory to the root of this repository.

Make sure the following file is present:

source data/IL23.csv
Recommended workflow

Run the scripts from the repository root in the following order.

1. Load packages
source("packages.r")
2. Load functions
source("functions.r")
3. Run the region-level analysis
source("region basis.r")

This script:

reads the Illinois PUMS data,
constructs poststratification cells,
builds the spatial basis,
fits the univariate Gaussian model,
fits the univariate Bernoulli model,
fits the joint multi-type model,
saves the fitted objects to data/region_basis_tau_in_gauss_results.RData.
4. Summarize the region-level results
source("region basis results.r")

This script:

loads data/region_basis_tau_in_gauss_results.RData,
computes region-level summaries,
compares posterior variances across models,
writes summary tables,
saves figures to figs/.

Typical outputs include:

basis_region_summary.csv
basis_region_variance_summary.csv
data/region_result_complete_basis.RData
figures in figs/
5. Run the empirical simulation study
source("empirical basis.r")

This script:

uses the Illinois data to define the area-level truth,
builds the same spatial basis,
repeatedly samples survey-like datasets,
fits direct, univariate, and multi-type estimators for each replicate,
stores the simulation output in data/empirical_basis_tau_in_gauss_results.RData.

Current default settings in the script are:

n_sim = 100
sample_size = 1000
nburn = 1000
nsim = 1000
nthin = 1
6. Summarize the empirical simulation results
source("empirical basis results.r")

This script:

loads the saved simulation results,
computes MSE, interval score, and coverage summaries,
produces both replicate-level and area-level comparison tables,
saves figures to figs/.

Typical outputs include:

basis_summary_table_area.csv
basis_summary_table_rep.csv
basis_absolute_table.csv
basis_paired_table.csv
basis_better_prop_table.csv
basis_paired_replicate_metrics.csv
basis_paired_area_metrics.csv
basis_better_prop_area_table.csv
figures in figs/
7. Optional: direct-estimator benchmark
source("HT_result.r")

This script computes direct-estimator benchmark metrics for the same simulation design.

Output folders

The main scripts create these folders if they do not already exist:

data/
figs/

Typical outputs are:

saved .RData objects in data/
.csv summary tables in the repository root
.png figures in figs/
Notes
region basis.r will use IL_poststrat_cells.csv if it is already available.
Otherwise, it constructs the poststratification cells directly from the sample data.
The basis-function specification keeps the eigenvectors associated with positive eigenvalues of the adjacency matrix.
The same modeling functions in functions.r are used by both the region analysis and the empirical simulation.
Minimal run example

If you want the main analysis pipeline only, run:

source("packages.r")
source("functions.r")
source("region basis.r")
source("region basis results.r")
source("empirical basis.r")
source("empirical basis results.r")
Reproducibility

For a clean run:

Clone the repository.
Set the working directory to the repository root.
Make sure source data/IL23.csv is available.
Run the scripts in the order listed above.
Contact

This repository contains the code used for the unit-level multi-type modeling project.