#' print.fitted_onlinesurr
#'
#' This method is wrapper for the summary.fitted_onlinesurr method.
#'
#' @param x A fitted_onlinesurr object.
#' @param ... Arguments passed to summary.fitted_onlinesurr
#'
#' @return No return value, called to print a summary of the fitted kDGLM model.
#'
#' @export
#' @keywords internal
print.fitted_onlinesurr <- function(x, ...) {
  summary.fitted_onlinesurr(x, ...)
}

#' Summarize a \code{"fitted_onlinesurr"} object
#'
#' Prints a human-readable report for an object of class \code{"fitted_onlinesurr"} returned by \code{fit.surr}. The report includes marginal and conditional treatment-effect estimates at a selected time point (or cumulatively up to that time), an estimate of the LPTE/CPTE, and a time-homogeneity test of the LPTE.
#'
#' @param object A \code{"fitted_onlinesurr"} object.
#' @param t Integer time index at which to evaluate treatment effects and the PTE. If \code{cummulative = TRUE}, effects are aggregated over times \code{1:t}. If \code{cummulative = FALSE}, effects are evaluated at time \code{t} only.
#' @param cummulative Logical; if \code{TRUE} (default), the report uses cumulative (up to time \code{t}) marginal and conditional treatment effects. If \code{FALSE}, the report uses the effects at time \code{t} only.
#' @param signif.level Numeric in \eqn{(0,1)} giving the significance level for the time-homogeneity test that is reported (e.g., via \code{time_homo_test}).
#' @param ... Additional arguments passed to downstream summary/print utilities (if any).
#'
#' @return No return value. Called for its side effect of printing a summary report.
#'
#' @details
#' The \code{"fitted_onlinesurr"} object stores point estimates and bootstrap samples for marginal and surrogate-adjusted (conditional) models in \code{object$Marginal} and \code{object$Conditional}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit.surr(y ~ 1,
#'   id = id, surrogate = ~ s1 + s2, treat = trt,
#'   data = dat, time = time, N.boots = 2000
#' )
#'
#' # Cumulative up to time 5
#' summary(fit, t = 5, cummulative = TRUE, signif.level = 0.05)
#'
#' # Time-specific at time 5
#' summary(fit, t = 5, cummulative = FALSE)
#' }
summary.fitted_onlinesurr <- function(object, t = object$T, cummulative = T, signif.level = 0.05, ...) {
  T <- object$T
  n <- object$n.fixed

  delta.est <- object$Marginal$point[1:T + n - T]
  delta.R.est <- object$Conditional$point[1:T + n - T]

  if (cummulative) {
    delta.est <- cumsum(delta.est)
    delta.R.est <- cumsum(delta.R.est)
  }

  format.num <- function(x) {
    ifelse(abs(x) < 0.001,
      format(x, digits = 5, justify = "l", scientific = TRUE),
      format(round(x, 5), width = 8, justify = "l", scientific = FALSE)
    )
  }

  names <- c("Delta", "Delta.R", ifelse(cummulative, "CPTE", "LPTE"))
  vals <- c(delta.est[t], delta.R.est[t], 1 - delta.R.est[t] / delta.est[t]) |>
    format.num() |>
    as.character()
  len.names <- sapply(1:3, function(x) {
    max(nchar(names[x]), nchar(vals[x]))
  })

  homo.test <- time_homo_test(object, signif.level)
  names.test <- c("Test stat.", "Crit. value", "p-value")
  vals.test <- as.numeric(homo.test) |>
    format.num() |>
    as.character()
  len.names.test <- sapply(seq_along(names.test), function(x) {
    max(nchar(names.test[x]), nchar(vals.test[x]))
  })

  cat(paste0(
    "Fitted Online Surrogate\n\n",
    ifelse(cummulative, "Cummulated", "Local"), " effects at time ", t, ":\n",
    paste(sapply(seq_along(names), function(x) {
      format(names[x], width = len.names[x], justify = "l")
    }), collapse = "   "), "\n",
    paste(sapply(seq_along(names), function(x) {
      format(vals[x], width = len.names[x], justify = "r")
    }), collapse = "   "), "\n",
    "---\n",
    "Time homogeneity test: \n",
    paste(sapply(seq_along(names.test), function(x) {
      format(names.test[x], width = len.names.test[x], justify = "l")
    }), collapse = "   "), "\n",
    paste(sapply(seq_along(names.test), function(x) {
      format(vals.test[x], width = len.names.test[x], justify = "r")
    }), collapse = "   "), "\n",
    "Signif. level: ", signif.level, "\n",
    "---\n"
  ))
}

#' Plot time-varying PTE measures and treatment effects from a \code{"fitted_onlinesurr"} object
#'
#' Produces a \code{ggplot2} figure showing, over time, either the Local PTE (LPTE), the Cumulative PTE (CPTE), or the marginal and residual treatment effects \eqn{\Delta(t)} and \eqn{\Delta_R(t)} (labeled \eqn{\Delta} and \eqn{\Delta_R} in the plot). Point estimates are taken from \code{object$Marginal$point} and \code{object$Conditional$point}, with uncertainty bands computed from the stored bootstrap draws.
#'
#' @param x A \code{"fitted_onlinesurr"} object, typically returned by \code{fit.surr}. It must contain \code{$T}, \code{$n.fixed}, and the components \code{$Marginal} and \code{$Conditional}, each with \code{point} and \code{smp}.
#' @param type Character string specifying what to plot. One of \code{"LPTE"}, \code{"CPTE"}, or \code{"Delta"} (case-insensitive). \code{"Delta"} plots both \eqn{\Delta(t)} and \eqn{\Delta_R(t)} with separate colors.
#' @param conf.level Numeric in \eqn{(0,1)} giving the confidence level for the plotted intervals. Default is \code{0.95}.
#' @param one.sided Logical; if \code{TRUE} (default), uses \code{signif.level = (1-conf.level)/2} when taking quantiles, so each tail excludes \code{1-conf.level} (i.e., a wider interval than the usual two-sided \code{conf.level} interval). This is convenient when visually assessing one-sided surrogate validation criteria. If \code{FALSE}, uses the standard two-sided construction \code{signif.level = 1-conf.level}.
#' @param ... Additional arguments (currently unused) included for S3 method compatibility.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The function extracts time-indexed treatment-effect estimates \eqn{\Delta(t)} (marginal) and \eqn{\Delta_R(t)} (residual/conditional) from the fitted object, along with bootstrap draws for each. It then constructs:
#' \itemize{
#'   \item \strong{LPTE:} \eqn{\mathrm{LPTE}(t) = 1 - \Delta_R(t)/\Delta(t)}.
#'   \item \strong{CPTE:} \eqn{\mathrm{CPTE}(t) = 1 - \sum_{u\le t}\Delta_R(u)/\sum_{u\le t}\Delta(u)}.
#'   \item \strong{Delta:} plots \eqn{\Delta(t)} and \eqn{\Delta_R(t)} directly.
#' }
#' Point estimates are plotted as points; intervals are empirical quantile intervals computed from the bootstrap sample matrices stored in \code{object}.
#'
#' @import ggplot2
#' @import latex2exp
#' @importFrom Rfast colCumSums
#' @importFrom kDGLM rowQuantile
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit.surr(y ~ 1,
#'   id = id, surrogate = ~ s1 + s2, treat = trt,
#'   data = dat, time = time, N.boots = 2000
#' )
#'
#' plot(fit, type = "LPTE")
#' plot(fit, type = "CPTE", conf.level = 0.90, one.sided = FALSE)
#' plot(fit, type = "Delta")
#' }
plot.fitted_onlinesurr <- function(x, type = "LPTE", conf.level = 0.95, one.sided = TRUE, ...) {
  signif.level <- 1 - conf.level
  if (one.sided) {
    signif.level <- signif.level / 2
  }

  type <- tolower(type)
  T <- x$T
  n <- x$n.fixed

  delta.est <- x$Marginal$point[1:T + n - T]
  delta.R.est <- x$Conditional$point[1:T + n - T]

  delta.smp <- x$Marginal$smp[1:T + n - T, ]
  delta.R.smp <- x$Conditional$smp[1:T + n - T, ]

  if (type == "lpte") {
    pte <- 1 - delta.R.est / delta.est
    pte.smp <- 1 - delta.R.smp / delta.smp

    plot.data <- data.frame(
      Time = 1:T,
      point = pte,
      icl = rowQuantile(pte.smp, signif.level),
      icu = rowQuantile(pte.smp, 1 - signif.level)
    )
  } else if (type == "cpte") {
    pte <- 1 - cumsum(delta.R.est) / cumsum(delta.est)
    pte.smp <- 1 - colCumSums(delta.R.smp) / colCumSums(delta.smp)

    plot.data <- data.frame(
      Time = 1:T,
      point = pte,
      icl = rowQuantile(pte.smp, signif.level),
      icu = rowQuantile(pte.smp, 1 - signif.level)
    )
  } else if (type == "delta") {
    plot.data <- data.frame(
      Time = rep(1:T, 2),
      point = c(delta.est, delta.R.est),
      icl = c(rowQuantile(delta.smp, signif.level), rowQuantile(delta.R.smp, signif.level)),
      icu = c(rowQuantile(delta.smp, 1 - signif.level), rowQuantile(delta.R.smp, 1 - signif.level)),
      color = rep(c("$Delta$", "$Delta_R$"), each = T)
    )
  }

  p <- ggplot(
    plot.data,
    aes(x = .data$Time, color = if (type == "delta") {
      .data$color
    } else {
      NULL
    })
  ) +
    geom_point(aes(y = .data$point)) +
    geom_errorbar(aes(ymin = .data$icl, ymax = .data$icu)) +
    theme_bw()
  if (type == "delta") {
    p <- p +
      scale_y_continuous("Treatment effect") +
      scale_color_hue("", labels = TeX)
  } else {
    p <- p + scale_y_continuous(toupper(type))
  }
  p
}
