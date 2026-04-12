# Rewrite function calls in an expression

Recursively traverses `expr` and replaces calls to `fun` using the
configuration returned by `config_fun`.

## Usage

``` r
rewrite_calls(
  expr,
  fun = "s",
  config_fun = get.config.s,
  eval_env = parent.frame()
)
```

## Arguments

- expr:

  An R language object.

- fun:

  A character string giving the function name to rewrite.

- config_fun:

  A function returning the replacement argument list for each matched
  call.

- eval_env:

  An environment used to evaluate arguments passed to `config_fun`.

## Value

A modified language object with matching calls rewritten.
