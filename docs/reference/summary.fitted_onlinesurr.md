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
#> Delta    6.65706  0.04906   135.70256 0.0000e+00 ***
#> Delta.R  2.27495  0.34700     6.55610 5.5232e-11 ***
#> CPTE     0.65827  0.05204       -         -       
#> 
#> Time homogeneity test: 
#> 
#> Test stat.   Crit. value   p-value     
#>    1.14145       2.44355    0.59894    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Time-specific at time 5
summary(fit, t = 5, cumulative = FALSE)
#> Fitted Online Surrogate
#> 
#> Local effects at time 5:
#>         Estimate Std. Error t value  Pr(>|t|)   
#> Delta    2.01452  0.02173   92.71791 0.0000e+00 ***
#> Delta.R  0.62921  0.10903    5.77078 7.8907e-09 ***
#> LPTE     0.68766  0.05412      -         -       
#> 
#> Time homogeneity test: 
#> 
#> Test stat.   Crit. value   p-value     
#>    1.14145       2.44569    0.59910    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
