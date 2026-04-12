# Compute a spline configuration

Derive the configuration needed to evaluate
[`s()`](https://silvaneojunior.github.io/OnlineSurr/reference/s.md) for
a given input vector. This helper resolves boundary limits, chooses or
validates knot locations, augments the knot vector with repeated
boundary knots, and returns the arguments in a form that can be used to
rebuild or update spline calls.

## Usage

``` r
get.config.s(
  x,
  P = 3,
  K = min(7, max(3, floor(log2(length(unique(x)))))),
  limits = c(NA, NA),
  knots = "eq"
)
```

## Arguments

- x:

  A numeric vector of predictor values.

- P:

  A non-negative integer giving the spline degree. `P = 3` corresponds
  to a cubic B-spline basis.

- K:

  An integer giving the number of basis functions implied by the spline
  specification.

- limits:

  A numeric vector of length 2 giving the lower and upper boundary
  limits for the spline basis. Missing values are replaced by `min(x)`
  and `max(x)`.

- knots:

  Either a numeric vector of knot locations, or one of `"eq"` or
  `"quantile"`. If `"eq"`, knots are placed uniformly between
  `limits[1]` and `limits[2]`. If `"quantile"`, knots are placed at
  equally spaced empirical quantiles of `x`.

## Value

A named list containing:

- x:

  The unevaluated expression supplied as `x`, returned via
  `substitute(x)`.

- knots:

  The full augmented knot vector, including repeated boundary knots.

- limits:

  The resolved lower and upper boundary limits.

- P:

  The spline degree.

- K:

  The number of basis functions.

## Details

This function is intended for programmatic use, for example when
rewriting model formulas that contain calls to `s(...)`.
