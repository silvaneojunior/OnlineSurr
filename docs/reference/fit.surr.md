# Fit marginal and conditional state-space models for longitudinal surrogate evaluation

Fits two Gaussian state-space models (Dynamic Linear Models) to jointly
longitudinal outcome data: (i) a marginal model for the outcome
trajectory given treatment and time, and (ii) a conditional model that
additionally adjusts for a user-specified surrogate structure. The
function returns per-time treatment-effect estimates from both models
and subject-level bootstrap draws obtained via subject-level resampling.

## Usage

``` r
fit.surr(
  formula,
  id,
  surrogate,
  treat,
  data = NULL,
  time = NULL,
  N.boots = 2000,
  verbose = 1
)
```

## Arguments

- formula:

  An object of class `formula` describing the fixed-effects mean
  structure for the primary outcome. The left-hand side must be the
  outcome variable. Internally, the right-hand side is augmented to
  include treatment-by-time fixed effects.

- id:

  A variable (unquoted) identifying subjects. Each subject must have at
  most one measurement per `time` value.

- surrogate:

  A formula describing the surrogate structure to be included in the
  conditional model. May be provided either as a `formula` (e.g.,
  `~ s1 + s2`) or as a string that can be coerced to a formula.

- treat:

  A variable (unquoted) indicating treatment assignment. Must encode
  exactly two treatment levels after coercion to a factor.

- data:

  A `data.frame` containing all variables referenced in `formula`, `id`,
  `treat`, `surrogate`, and (optionally) `time`.

- time:

  Optional variable (unquoted) giving the measurement time index. Must
  be numeric and equally spaced across observed time points. If `NULL`,
  an equally spaced within-subject index is created in the current row
  order (with a warning).

- N.boots:

  Integer number of subject-level bootstrap replicates. Each replicate
  resamples subjects with replacement and recombines subject-specific
  sufficient quantities to form bootstrap draws of the fixed effects.

- verbose:

  Logical scalar indicating whether to print progress information during
  model fitting. If `TRUE`, progress updates are shown; if `FALSE`, no
  progress output is produced.

## Value

An object of class `"fitted_onlinesurr"`: a named list with elements
`$Marginal` and `$Conditional`. Each of these contains:

- `point`: the point estimate vector of the fixed effects (excluding
  subject-specific random-walk states) at the final time point.

- `smp`: a matrix of bootstrap draws for those fixed effects, with one
  column per bootstrap replicate.

The object also includes:

- `T`: number of unique time points.

- `N`: number of subjects.

- `n.fixed`: number of fixed-effect coefficients implied by `formula`
  for a single subject prior to stacking across subjects.

## Details

The implementation follows a two-model decomposition used for estimating
longitudinal treatment effects and surrogate-adjusted (residual)
treatment effects in a state-space framework.

**Data requirements.** The data must have at most one row per
subject-time pair; time must be numeric and equally spaced (or omitted,
in which case an index is created). Treatment and subject identifiers
are coerced to factors with sorted levels.

**Model structure.** The marginal model includes treatment-by-time fixed
effects and a subject-specific random-walk component to capture
within-subject correlation. The conditional model adds the
user-specified surrogate structure to the design, and checks that
treatment is not a linear combination of the surrogate design (rank
check).

**Bootstrap.** Subjects are resampled with replacement. Subject-specific
filtered quantities are computed once and recombined in each bootstrap
iteration to reduce computational cost, consistent with a subject-level
nonparametric bootstrap strategy for replicated time series.

## Examples

``` r
if (FALSE) { # \dontrun{
# data columns: y (outcome), id (subject id), trt (0/1 or two-level factor),
# time (numeric equally spaced), s1 and s2 (surrogates)

fit <- fit.surr(
  formula    = y ~ 1, # baseline fixed effects; function adds trt*time terms
  id         = id,
  surrogate  = ~ s1 + s2,
  treat      = trt,
  data       = dat,
  time       = time,
  N.boots    = 500
)

# Access point estimates and bootstrap samples
fit$Marginal$point
fit$Conditional$smp[, 1:10]
} # }
```
