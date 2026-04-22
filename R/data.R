#' Simulated longitudinal surrogate dataset for `OnlineSurr`
#'
#' A simulated long-format dataset illustrating the input structure expected by
#' [fit.surr()] for surrogate evaluation with jointly longitudinal outcomes and
#' surrogate markers.
#'
#' The dataset contains 100 subjects observed at 6 equally spaced time points.
#' Treatment assignment is binary and constant within subject. The surrogate
#' marker `s` varies over time and is affected by treatment. The primary outcome
#' `y` depends on treatment, time, and the surrogate marker.
#'
#' Rows are ordered by subject identifier and time.
#'
#' @format A data frame with 600 rows and 5 variables:
#' \describe{
#'   \item{id}{Integer subject identifier.}
#'   \item{trt}{Binary treatment indicator: `0` for control and `1` for treated.}
#'   \item{time}{Numeric measurement time index taking values `1` through `6`.}
#'   \item{s}{Continuous longitudinal surrogate marker.}
#'   \item{y}{Continuous longitudinal primary outcome.}
#' }
#'
#' @details
#' This dataset was generated for package examples and testing. It represents a balanced longitudinal design with one observation per subject-time pair.
#' Measurement times are equally spaced, which is a requirement for use with [fit.surr()].
#'
#' In the data-generating mechanism, the surrogate marker is affected by time and treatment, and the outcome depends on time, treatment, and the surrogate.
#'
#' @source Simulated data generated within the package; not based on an external study.
#'
#' @keywords datasets
"sim_onlinesurr"
