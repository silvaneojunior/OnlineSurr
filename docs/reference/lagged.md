# Compute lagged values of a vector

Returns a lagged version of `x` by shifting its values forward by `k`
positions and padding the first `k` entries with zeros.

## Usage

``` r
lagged(x, k = 1)
```

## Arguments

- x:

  A vector to be lagged.

- k:

  A non-negative integer giving the lag order.

## Value

A vector of the same length as `x`, with the first `k` values set to `0`
and the remaining values taken from `x` shifted by `k` positions.

## Details

This function is intended for use in model formulas when delayed effects
of a predictor should be included explicitly.

## Examples

``` r
x <- 1:5

lagged(x)
#> [1] 0 1 2 3 4
lagged(x, k = 2)
#> [1] 0 0 1 2 3
```
