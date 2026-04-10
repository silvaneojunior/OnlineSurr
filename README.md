# OnlineSurr <img src="hex_badge.png" align="right" height="220" alt="OnlineSurr hex logo" />

`OnlineSurr` is an R package for surrogate evaluation when both the primary outcome and the surrogate marker are measured longitudinally. The package implements a state-space approach based on the methodology developed in *A Causal Framework for Evaluating Jointly Longitudinal Outcomes and Surrogate Markers: A State-Space Approach*.

The current implementation fits two Gaussian state-space models:

- a **marginal model** for the longitudinal outcome as a function of treatment and time;
- a **conditional model** that additionally adjusts for a user-specified surrogate structure.

From these two fitted models, the package produces time-specific treatment-effect estimates, bootstrap-based uncertainty summaries, estimates of the local and cumulative proportion of treatment effect explained (LPTE and CPTE), and a test for temporal homogeneity of the PTE.

## Installation

The package is not yet on CRAN. Install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("silvaneojunior/OnlineSurr")
```

Depending on your local setup, you may also need to install the package dependencies first.

## Dependencies

The current source imports functionality from:

- `kDGLM`
- `dplyr`
- `tidyr`
- `rlang`
- `ggplot2`
- `Rfast`
- `latex2exp`

## Basic usage

```r
library(OnlineSurr)

fit <- fit.surr(
  formula   = y ~ 1,
  id        = id,
  surrogate = ~ s1 + s2,
  treat     = trt,
  data      = dat,
  time      = time,
  N.boots   = 500
)

summary(fit, t = fit$T, cumulative = TRUE)
plot(fit, type = "LPTE")
plot(fit, type = "CPTE")
plot(fit, type = "Delta")
time_homo_test(fit)
```

Go here to view a tutorial for this package: [OnlineSurr vignette](https://htmlpreview.github.io/?https://github.com/silvaneojunior/OnlineSurr/vignettes/onlinesurr-vignette.html)

## Notes and limitations

At its current stage, the package is focused on the core methodology and assumes that the user supplies an appropriate surrogate specification. In practice, users should pay careful attention to:

- the choice of surrogate history or lag structure;
- whether treatment effects are well defined at each time point;
- study design assumptions required for causal interpretation.

The current implementation is best viewed as a research package accompanying the methodological paper.

## Citation

If you use this repository, please cite the associated paper:

```text
Silvaneo V. dos Santos Jr. and Layla Parast.
A Causal Framework for Evaluating Jointly Longitudinal Outcomes and Surrogate Markers: A State-Space Approach.
```

You may also want to cite the package repository directly once a public versioned release is available.
