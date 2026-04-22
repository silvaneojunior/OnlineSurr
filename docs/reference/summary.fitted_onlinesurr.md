# Summarize a `"fitted_onlinesurr"` object

Prints a human-readable report for an object of class
`"fitted_onlinesurr"` returned by `fit.surr`. The report includes
marginal and conditional treatment-effect estimates at a selected time
point (or cumulatively up to that time), an estimate of the LPTE/CPTE,
and a time-homogeneity test of the LPTE.

## Usage

``` r
# S3 method for class 'fitted_onlinesurr'
summary(object, t = object$T, cumulative = TRUE, signif.level = 0.05, ...)
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
fit <- fit.surr(y ~ 1,
  id = id,
  surrogate = ~s,
  treat = trt,
  data = sim_onlinesurr, # This dataset is included in the OnlineSurr package
  time = time,
  verbose = 0,
  N.boots = 500 # Generally, this value would be too small.
  # Remember to increase it for your dataset.
)

# Cumulative up to time 5
summary(fit, t = 5, cumulative = TRUE, signif.level = 0.05)
#> Fitted Online Surrogate
#> 
#> Cummulated effects at time 5:
#>         Estimate Std. Error t value   Pr(>|t|)   
#> Delta    6.65540  0.04904   135.72003 0.0000e+00 ***
#> Delta.R  2.27557  0.35247     6.45608 1.0745e-10 ***
#> CPTE     0.65809  0.05284       -         -       
#> 
#> Time homogeneity test: 
#> 
#> Test stat.   Crit. value   p-value     
#>    1.11084       2.43794    0.61612    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Time-specific at time 5
summary(fit, t = 5, cumulative = FALSE)
#> Fitted Online Surrogate
#> 
#> Local effects at time 5:
#>         Estimate Std. Error t value  Pr(>|t|)   
#> Delta    2.01395  0.02172   92.71419 0.0000e+00 ***
#> Delta.R  0.62847  0.11065    5.67994 1.3474e-08 ***
#> LPTE     0.68794  0.05494      -         -       
#> 
#> Time homogeneity test: 
#> 
#> Test stat.   Crit. value   p-value     
#>    1.11084       2.43875    0.61710    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
