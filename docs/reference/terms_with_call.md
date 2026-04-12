# Extract formula terms containing a function call

Returns the term labels from a formula whose parsed expressions contain
a call to `fun`.

## Usage

``` r
terms_with_call(formula, fun = "s")
```

## Arguments

- formula:

  A model formula.

- fun:

  A character string giving the function name to look for.

## Value

A character vector of term labels containing a call to `fun`.
