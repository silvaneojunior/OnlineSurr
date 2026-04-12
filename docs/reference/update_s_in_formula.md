# Update spline calls in a formula

Rewrites calls to `fun` in a formula using `config_fun`, evaluating
arguments in an environment built from `data`.

## Usage

``` r
update_s_in_formula(formula, data, config_fun = get.config.s, fun = "s")
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame or list providing variables used in the formula.

- config_fun:

  A function returning the replacement argument list for each matched
  call.

- fun:

  A character string giving the function name to rewrite.

## Value

A formula with updated calls to `fun`.
