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
  verbose = 1,
  D.local = 0.8
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

- D.local:

  Numeric, a number between 0 and 1 indicating the discount factor to be
  used for the random effect block. This factor controls how smooth the
  random effect evolve over time. A discount factor of 1 means that the
  random effects do not change over time, so that each individual has
  its own local level, but that level is the same for all times. A
  discount factor of 0 is not acceptable (the kDGLM package will replace
  it by 1), but values closer to 0 imply in a more flexible dynamic. See
  West and Harrison (1997) or the appendix in dos Santos Jr. and
  Parast (2026) for instructions on how to specify the discount factor.

## Value

An object of class `"fitted_onlinesurr"`: a named list with elements
`$Marginal` and `$Conditional`. Each of these contains:

- `point`: the point estimate vector of the treatment effect at each
  time point.

- `smp`: a matrix of bootstrap draws for the treatment effect at each
  time point, with one column per bootstrap replicate. The draws are
  generated from the joint distribution of the full vector, thereby
  accounting for the dependence among different time points. The samples
  from the marginal (total effect) and conditional (residual effect)
  models are paired, so that the i-th samples from both models are drawn
  jointly from the distribution of the estimators.

The object also includes:

- `T`: number of unique time points.

- `N`: number of subjects.

- `n.fixed`: number of fixed-effect coefficients implied by `formula`
  for a single subject prior to stacking across subjects.

## Details

The implementation follows a two-model decomposition used for estimating
longitudinal treatment effects and surrogate-adjusted (residual)
treatment effects in a state-space framework.

See dos Santos Jr. and Parast (2026) for details on the methodology.

See West and Harrison (1997) for best practices on model specification
in the state-space model setting.

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

## References

Silvaneo V. dos Santos Jr., Layla Parast (2026). “A Causal Framework for
Evaluating Jointly Longitudinal Outcomes and Surrogate Markers: A
State-Space Approach.” 2604.12882, <https://arxiv.org/abs/2604.12882>.  
  
Mike West, Jeff Harrison (1997). *Bayesian Forecasting and Dynamic
Models (Springer Series in Statistics)*. Springer-Verlag. ISBN
0387947256.

## Examples

``` r
fit <- fit.surr(y ~ 1,
  id = id,
  surrogate = ~s,
  treat = trt,
  data = sim_onlinesurr, # This dataset is included in the OnlineSurr package
  time = time,
  verbose = 0,
  N.boots = 500 # Generally, this value would be too small.
  # Remember to increase it for your dataset.
)
summary(fit)
#> Fitted Online Surrogate
#> 
#> Cummulated effects at time 6:
#>         Estimate Std. Error t value   Pr(>|t|)   
#> Delta    8.99606  0.04861   185.06592 0.0000e+00 ***
#> Delta.R  2.99490  0.47800     6.26547 3.7169e-10 ***
#> CPTE     0.66709  0.05306       -         -       
#> 
#> Time homogeneity test: 
#> 
#> Test stat.   Crit. value   p-value     
#>    1.11091       2.43248    0.62236    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
