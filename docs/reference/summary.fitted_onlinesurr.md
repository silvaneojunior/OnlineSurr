# Summarize a `"fitted_onlinesurr"` object

Prints a human-readable report for an object of class
`"fitted_onlinesurr"` returned by `fit.surr`. The report includes
marginal and conditional treatment-effect estimates at a selected time
point (or cumulatively up to that time), an estimate of the LPTE/CPTE,
and a time-homogeneity test of the LPTE.

## Usage

``` r
# S3 method for class 'fitted_onlinesurr'
summary(object, t = object$T, cumulative = T, signif.level = 0.05, ...)
```

## Arguments

- object:

  A `"fitted_onlinesurr"` object.

- t:

  Integer time index at which to evaluate treatment effects and the PTE.
  If `cumulative = TRUE`, effects are aggregated over times `1:t`. If
  `cumulative = FALSE`, effects are evaluated at time `t` only.

- cumulative:

  Logical; if `TRUE` (default), the report uses cumulative (up to time
  `t`) marginal and conditional treatment effects. If `FALSE`, the
  report uses the effects at time `t` only.

- signif.level:

  Numeric in \\(0,1)\\ giving the significance level for the
  time-homogeneity test that is reported (e.g., via `time_homo_test`).

- ...:

  Additional arguments passed to downstream summary/print utilities (if
  any).

## Value

No return value. Called for its side effect of printing a summary
report.

## Details

The `"fitted_onlinesurr"` object stores point estimates and bootstrap
samples for marginal and surrogate-adjusted (conditional) models in
`object$Marginal` and `object$Conditional`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit.surr(y ~ 1,
  id = id, surrogate = ~ s1 + s2, treat = trt,
  data = dat, time = time, N.boots = 2000
)

# Cumulative up to time 5
summary(fit, t = 5, cumulative = TRUE, signif.level = 0.05)

# Time-specific at time 5
summary(fit, t = 5, cumulative = FALSE)
} # }
```
