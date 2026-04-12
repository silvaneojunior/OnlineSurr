# Test time-homogeneity of the PTE

Tests the null hypothesis that the LPTE is constant over time. The test
is based on the difference between the conditional and marginal
treatment-effect trajectories implied by a fitted `"fitted_onlinesurr"`
object, standardized by an estimated covariance, and uses a max-type
statistic to control the family wise error across time points.

## Usage

``` r
time_homo_test(model, signif.level = 0.05, N.boots = 50000)
```

## Arguments

- model:

  A fitted object of class `"fitted_onlinesurr"`, typically returned by
  `fit.surr`. Must contain `$T`, `$n.fixed`, and the elements
  `$Marginal` and `$Conditional` with `point` and `smp` components.

- signif.level:

  Numeric in (0,1) giving the test significance level used to form the
  critical value from the bootstrap distribution. Default is `0.05`.

- N.boots:

  Integer number of Monte Carlo draws used to approximate the null
  distribution of the max standardized deviation statistic and to
  compute the p-value. Default is `50000`.

## Value

A named list with:

- `T`: the observed test statistic (maximum absolute standardized
  deviation).

- `T.crit`: the 1-signif.level critical value.

- `p.value`: the Monte Carlo p-value `mean(T_null > T_obs)`.

## Details

Notes:

- The function assumes the first `T` time-specific treatment-effect
  parameters are stored contiguously at the beginning of
  `model$Marginal$point` and `model$Conditional$point` (and similarly
  for `smp`). It uses the index `1:(n.fixed)` as implemented in the
  code: `1:(T + n.fixed - T)`.

- `N.boots` here is a Monte Carlo size for the null simulation (distinct
  from the bootstrap size used when fitting `model`).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit.surr(y ~ 1,
  id = id, surrogate = ~ s1 + s2, treat = trt,
  data = dat, time = time, N.boots = 2000
)

time_homo_test(fit, signif.level = 0.05, N.boots = 50000)
} # }
```
