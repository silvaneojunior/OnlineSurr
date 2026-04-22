# Simulated longitudinal surrogate dataset for \`OnlineSurr\`

A simulated long-format dataset illustrating the input structure
expected by \[fit.surr()\] for surrogate evaluation with jointly
longitudinal outcomes and surrogate markers.

## Usage

``` r
sim_onlinesurr
```

## Format

A data frame with 600 rows and 5 variables:

- id:

  Integer subject identifier.

- trt:

  Binary treatment indicator: \`0\` for control and \`1\` for treated.

- time:

  Numeric measurement time index taking values \`1\` through \`6\`.

- s:

  Continuous longitudinal surrogate marker.

- y:

  Continuous longitudinal primary outcome.

## Source

Simulated data generated within the package; not based on an external
study.

## Details

The dataset contains 100 subjects observed at 6 equally spaced time
points. Treatment assignment is binary and constant within subject. The
surrogate marker \`s\` varies over time and is affected by treatment.
The primary outcome \`y\` depends on treatment, time, and the surrogate
marker.

Rows are ordered by subject identifier and time.

This dataset was generated for package examples and testing. It
represents a balanced longitudinal design with one observation per
subject-time pair. Measurement times are equally spaced, which is a
requirement for use with \[fit.surr()\].

In the data-generating mechanism, the surrogate marker is affected by
time and treatment, and the outcome depends on time, treatment, and the
surrogate.
