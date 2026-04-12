# Construct a B-spline basis matrix

Build a B-spline basis for a numeric vector using a Cox-de Boor style
recursion. By default, the function constructs a cubic spline basis
(`P = 3`) and chooses the number of basis functions from the number of
unique values in `x`.

## Usage

``` r
s(
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

  An integer giving the number of basis functions to return. The default
  increases slowly with the number of unique values in `x`.

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

A numeric matrix with one row per element of `x` and one column per
spline basis function.

## Details

Boundary limits are taken from `x` unless supplied explicitly. Knot
locations may be given directly as a numeric vector, or generated either
at equally spaced locations (`"eq"`) or at empirical quantiles
(`"quantile"`).

The returned basis has `length(x)` rows and `k` columns.

When `knots` is generated internally, the function first creates
`K - P + 1` knot locations and then augments them with repeated boundary
knots so the recursion can be evaluated.

## Examples

``` r
x <- seq(0, 1, length.out = 10)

# Default cubic basis
B <- s(x)
dim(B)
#> [1] 10  3

# Equally spaced knots with custom basis size
B2 <- s(x, K = 5, knots = "eq")

# Quantile-based knots
B3 <- s(x, knots = "quantile")
```
