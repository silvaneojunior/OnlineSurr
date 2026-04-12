# Plot time-varying PTE measures and treatment effects from a `"fitted_onlinesurr"` object

Produces a `ggplot2` figure showing, over time, either the Local PTE
(LPTE), the Cumulative PTE (CPTE), or the marginal and residual
treatment effects \\\Delta(t)\\ and \\\Delta_R(t)\\ (labeled \\\Delta\\
and \\\Delta_R\\ in the plot). Point estimates are taken from
`object$Marginal$point` and `object$Conditional$point`, with uncertainty
bands computed from the stored bootstrap draws.

## Usage

``` r
# S3 method for class 'fitted_onlinesurr'
plot(x, type = "LPTE", conf.level = 0.95, one.sided = TRUE, ...)
```

## Arguments

- x:

  A `"fitted_onlinesurr"` object, typically returned by `fit.surr`. It
  must contain `$T`, `$n.fixed`, and the components `$Marginal` and
  `$Conditional`, each with `point` and `smp`.

- type:

  Character string specifying what to plot. One of `"LPTE"`, `"CPTE"`,
  or `"Delta"` (case-insensitive). `"Delta"` plots both \\\Delta(t)\\
  and \\\Delta_R(t)\\ with separate colors.

- conf.level:

  Numeric in \\(0,1)\\ giving the confidence level for the plotted
  intervals. Default is `0.95`.

- one.sided:

  Logical; if `TRUE` (default), uses `signif.level = (1-conf.level)/2`
  when taking quantiles, so each tail excludes `1-conf.level` (i.e., a
  wider interval than the usual two-sided `conf.level` interval). This
  is convenient when visually assessing one-sided surrogate validation
  criteria. If `FALSE`, uses the standard two-sided construction
  `signif.level = 1-conf.level`.

- ...:

  Additional arguments (currently unused) included for S3 method
  compatibility.

## Value

A `ggplot` object.

## Details

The function extracts time-indexed treatment-effect estimates
\\\Delta(t)\\ (marginal) and \\\Delta_R(t)\\ (residual/conditional) from
the fitted object, along with bootstrap draws for each. It then
constructs:

- **LPTE:** \\\mathrm{LPTE}(t) = 1 - \Delta_R(t)/\Delta(t)\\.

- **CPTE:** \\\mathrm{CPTE}(t) = 1 - \sum\_{u\le
  t}\Delta_R(u)/\sum\_{u\le t}\Delta(u)\\.

- **Delta:** plots \\\Delta(t)\\ and \\\Delta_R(t)\\ directly.

Point estimates are plotted as points; intervals are empirical quantile
intervals computed from the bootstrap sample matrices stored in
`object`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit.surr(y ~ 1,
  id = id, surrogate = ~ s1 + s2, treat = trt,
  data = dat, time = time, N.boots = 2000
)

plot(fit, type = "LPTE")
plot(fit, type = "CPTE", conf.level = 0.90, one.sided = FALSE)
plot(fit, type = "Delta")
} # }
```
