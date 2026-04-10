# Reproducibility code for the paper

This repository contains the R code used to reproduce the simulation study, scenario figures, and applied example reported in the paper. The project is organized into two main components:

1. `simulation/`: code for the simulation study and the comparisons with competing methods.
2. `applied_example/`: code for the applied data analysis.


## Repository structure

```text
.
├── applied_example/
|   ├── figures/
│   ├── fitting.R
│   ├── init.R
│   └── results.R
├── simulation/
│   ├── data/
|   ├── figures/
│   ├── functions.R
│   ├── init.R
│   ├── inspection.R
│   ├── run_cluster.R
│   ├── run_single.R
│   └── scenarios.R
└── README.md
```

## Requirements

The code is written in R and relies on the following packages:

- `tidyverse`
- `kDGLM`
- `qs`
- `Rfast`
- `cubature`
- `stringr`
- `tictoc`
- `readr`
- `clubSandwich`
- `extrafont`
- `lme4`
- `geepack`
- `latex2exp`
- `parallel`
- `plotly`
- `knitr`

Install the packages before running the scripts. For example:

```r
install.packages(c(
  "tidyverse", "qs", "Rfast", "cubature", "stringr", "tictoc",
  "readr", "clubSandwich", "extrafont", "lme4", "geepack",
  "latex2exp", "plotly", "knitr"
))
```

Some scripts call `loadfonts(device = "win")`. This reflects the original Windows setup used for figure generation. On non-Windows systems, you may need to adapt or remove the font-loading lines.
If you want to use the same font workflow, run `font_import()` once before the first execution.
This is purely aesthetic and completely optional.

## Data availability

The applied example uses an external dataset that is **not included** in this repository. The script `applied_example/init.R` expects the file:

```text
restricted_data/dataset2/data_long_stream.csv
```

Because this dataset is restricted, the applied example cannot be reproduced from a clean clone of the repository unless the user has access to the original data and places it at the expected location.

The simulation study does not require external data files; it generates synthetic datasets internally.

## Reproducing the simulation study

The simulation code defines the study settings, generates synthetic data under multiple scenarios, fits the proposed method and several competitors, and summarizes the results.

### Step 1: set up the simulation environment

`simulation/init.R` loads dependencies and defines the main simulation parameters, including:

- study lengths in `T.list = c(10, 5, 3, 2, 1)`
- sample sizes in `N.list = c(300, 250, 200, 150, 100, 50)`
- bootstrap size `N.bootstrap = 2000`
- scenario-specific total effects and true PTE values

### Step 2: run simulation replications

For a single batch of replications, run:

```r
source("simulation/run_single.R")
```

This script sources `simulation/init.R`, runs a batch of replications through `run_set()`, and writes the results to a timestamped CSV file.

For parallel execution, run:

```r
source("simulation/run_cluster.R")
```

This optional script distributes batches across multiple cores using the `parallel` package. In the original paper run, 20 parallel runs were used.

### Step 3: generate scenario illustrations

To recreate the scenario plots, run:

```r
source("simulation/scenarios.R")
```

This script generates the figure showing the treatment-effect dynamics for Study 1 and the additional scenario plots used for Study 2.

### Step 4: summarize simulation results

After the simulation output files have been created, run:

```r
source("simulation/inspection.R")
```

This script reads the simulation result files, computes summary metrics such as bias, empirical coverage, confidence-interval width, and rejection rates, and saves the corresponding figures.

## Reproducing the applied example

The applied example is organized as a three-step workflow.

### Step 1: prepare the data and objects

Run:

```r
source("applied_example/init.R")
```

This script:

- loads the dataset
- reshapes the data into long format
- standardizes the surrogate and outcome variables
- defines the design matrices and analysis objects used later in the model-fitting step

### Step 2: fit the models

Run:

```r
source("applied_example/fitting.R")
```

This script fits the marginal and conditional state-space models for multiple maximum surrogate lags, performs the bootstrap calculations, and saves the workspace image for later use.

### Step 3: produce the applied-example results

Run:

```r
source("applied_example/results.R")
```

This script loads the saved workspace from the fitting step, computes the treatment-effect summaries, derives LPTE and CPTE curves, performs the time-homogeneity test, and exports the applied-example figures.

## Main outputs

Depending on which scripts are run, the repository produces:

- simulation output CSV files from the batch runs
- scenario figures from `simulation/scenarios.R`
- simulation summary figures from `simulation/inspection.R`
- a saved workspace from `applied_example/fitting.R`
- applied-example figures from `applied_example/results.R`

## Minimal reproduction order

For convenience, the recommended execution order is:

### Simulation study

```r
source("simulation/init.R")
source("simulation/run_single.R")   # or simulation/run_cluster.R
source("simulation/scenarios.R")
source("simulation/inspection.R")
```

### Applied example

```r
source("applied_example/init.R")
source("applied_example/fitting.R")
source("applied_example/results.R")
```

## Citation

If you use this code, please cite the corresponding paper.
