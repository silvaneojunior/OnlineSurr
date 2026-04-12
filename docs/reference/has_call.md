# Test whether an expression contains a function call

Recursively checks whether `expr` contains a call to `fun`.

## Usage

``` r
has_call(expr, fun = "s")
```

## Arguments

- expr:

  An R language object.

- fun:

  A character string giving the function name to look for.

## Value

`TRUE` if `expr` contains a call to `fun`, otherwise `FALSE`.
