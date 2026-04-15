#' rmvnorm
#'
#' Obtains a sample from a multivariate normal distribution.
#'
#' @param n integer: The sample size.
#' @param mu numeric: The mean vector
#' @param Sigma matrix: The Covariance matrix.
#'
#' @importFrom Rfast matrnorm transpose
#' @importFrom stats runif
#'
#' @keywords internal
rmvnorm <- function(n, mu, Sigma,
                    norm.x = matrnorm(k, n, seed = round(runif(1) * 1e15))) {
  k <- length(mu)
  chol.Sigma <- var_decomp(Sigma)
  transpose(chol.Sigma) %*% norm.x + c(mu)
}

#' var_decomp
#'
#' This function receives a covariance matrix S and creates a matrix Q, so that t(Q) \%*\% Q = S.
#'
#' @param S A covariance matrix
#'
#' @keywords internal
var_decomp <- function(S) {
  n <- dim(S)[1]
  chol.decomp <- suppressWarnings({
    chol(S, pivot = TRUE)
  })
  pivot <- attr(chol.decomp, "pivot")
  oo <- order(pivot)
  chol.decomp <- chol.decomp[, oo, drop = FALSE]
  return(chol.decomp)
}

#' formula.to.structure
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which glm is called.
#' @param label An optional character naming the linear predictor.
#'
#' @importFrom stats update.formula model.matrix
#' @importFrom rlang is_formula
#'
#' @keywords internal
formula.to.structure <- function(formula, data, label = "mu") {
  terms <- attr(terms(formula), "term.labels")
  intercept.add <- attr(terms(formula), "intercept") & !any(grepl("pol(", terms, fixed = TRUE))
  intercept.flag <- attr(terms(formula), "intercept")
  terms <- attr(terms(formula), "term.labels")
  args <- list()
  terms.mat <- c()

  if (length(terms) >= 1) {
    for (term in terms) {
      if (check_is_dlm_block(term)) {
        args <- append(args, list(eval(parse(text = term), envir = data)))
      } else {
        terms.mat <- append(terms.mat, term)
      }
    }
  }
  if (length(terms.mat) > 0 | intercept.add) {
    mat.formula <- update.formula(formula, as.formula(paste0("~", paste(c(ifelse(intercept.flag, 1, 0), terms.mat), collapse = "+"))))

    X <- model.matrix(mat.formula, data = data)
    if (!intercept.add & intercept.flag) {
      X <- X[, -1, drop = FALSE]
    }
    for (j in 1:dim(X)[2]) {
      args <- append(args, list(reg(X[, j], name = colnames(X)[j], D = 1, R1 = 900)))
    }
  }
  block_rename(do.call(block_superpos, args), label)
}


#' ginv
#'
#' This function receives a covariance matrix S and calculates the generalized inverse of S.
#'
#' @param S A covariance matrix
#'
#' @keywords internal
ginv <- function(S) {
  if (length(S) == 1) {
    return(1 / S)
  } else {
    S <- as.matrix(S)
    n <- dim(S)[1]
    chol.decomp <- suppressWarnings({
      chol(S, pivot = TRUE)
    })
    rank <- attr(chol.decomp, "rank")
    pivot <- attr(chol.decomp, "pivot")
    oo <- order(pivot)
    inv <- matrix(0, n, n)
    if (rank > 0) {
      inv[1:rank, 1:rank] <- chol2inv(chol.decomp[1:rank, 1:rank])
    }
    inv <- inv[oo, oo, drop = FALSE]
    return(inv)
  }
}


#' Compute lagged values of a vector
#'
#' Returns a lagged version of \code{x} by shifting its values forward by \code{k} positions and padding the first \code{k} entries with zeros.
#'
#' @param x A vector to be lagged.
#' @param k A non-negative integer giving the lag order.
#'
#' @returns A vector of the same length as \code{x}, with the first \code{k} values set to \code{0} and the remaining values taken from \code{x} shifted by \code{k} positions.
#'
#' @details
#' This function is intended for use in model formulas when delayed effects of a predictor should be included explicitly.
#'
#' @examples
#' x <- 1:5
#'
#' lagged(x)
#' lagged(x, k = 2)
#'
#' @export
lagged <- function(x, k = 1) {
  c(rep(0, k), x[1:(length(x) - k)])
}

#' Test whether an expression contains a function call
#'
#' Recursively checks whether \code{expr} contains a call to \code{fun}.
#'
#' @param expr An R language object.
#' @param fun A character string giving the function name to look for.
#'
#' @returns \code{TRUE} if \code{expr} contains a call to \code{fun}, otherwise \code{FALSE}.
#' @keywords internal
has_call <- function(expr, fun = "s") {
  if (is.call(expr)) {
    if (identical(expr[[1]], as.name(fun))) {
      return(TRUE)
    }
    return(any(vapply(as.list(expr)[-1], has_call, logical(1), fun = fun)))
  }
  FALSE
}

#' Extract formula terms containing a function call
#'
#' Returns the term labels from a formula whose parsed expressions contain a call to \code{fun}.
#'
#' @param formula A model formula.
#' @param fun A character string giving the function name to look for.
#'
#' @returns A character vector of term labels containing a call to \code{fun}.
#'
#' @importFrom stats terms
#' @keywords internal
terms_with_call <- function(formula, fun = "s") {
  labs <- attr(terms(formula), "term.labels")
  exprs <- lapply(labs, str2lang)
  labs[vapply(exprs, has_call, logical(1), fun = fun)]
}

#' Rewrite function calls in an expression
#'
#' Recursively traverses \code{expr} and replaces calls to \code{fun} using the configuration returned by \code{config_fun}.
#'
#' @param expr An R language object.
#' @param fun A character string giving the function name to rewrite.
#' @param config_fun A function returning the replacement argument list for each matched call.
#' @param eval_env An environment used to evaluate arguments passed to \code{config_fun}.
#'
#' @returns A modified language object with matching calls rewritten.
#' @keywords internal
rewrite_calls <- function(expr, fun = "s", config_fun = get.config.s, eval_env = parent.frame()) {
  if (is.call(expr)) {
    # First rewrite children
    parts <- lapply(as.list(expr), rewrite_calls, fun = fun, eval_env = eval_env, config_fun = config_fun)
    expr2 <- as.call(parts)

    # Then, if this node is fun(...), replace its arguments
    if (identical(expr2[[1]], as.name(fun))) {
      old_args <- as.list(expr2)[-1]
      names(old_args) <- names(as.list(expr2))[-1]

      new_args <- do.call(config_fun, old_args, envir = eval_env)

      if (!is.list(new_args)) {
        stop("config_fun must return a list of arguments for the rebuilt call")
      }

      return(as.call(c(list(as.name(fun)), new_args)))
    }

    return(expr2)
  }

  if (is.pairlist(expr)) {
    out <- lapply(as.list(expr), rewrite_calls, fun = fun, eval_env = eval_env, config_fun = config_fun)
    return(as.pairlist(out))
  }

  expr
}

#' Update spline calls in a formula
#'
#' Rewrites calls to \code{fun} in a formula using \code{config_fun}, evaluating arguments in an environment built from \code{data}.
#'
#' @param formula A model formula.
#' @param data A data frame or list providing variables used in the formula.
#' @param config_fun A function returning the replacement argument list for each matched call.
#' @param fun A character string giving the function name to rewrite.
#'
#' @returns A formula with updated calls to \code{fun}.
#' @keywords internal
update_s_in_formula <- function(formula, data, config_fun = get.config.s, fun = "s") {
  env <- environment(formula)
  data_env <- list2env(as.list(data), parent = env)
  out <- rewrite_calls(formula,
    fun = fun, config_fun = config_fun,
    eval_env = data_env
  )
  as.formula(out, env = env)
}

#' Construct a B-spline basis matrix
#'
#' Build a B-spline basis for a numeric vector using a Cox-de Boor style recursion. By default, the function constructs a cubic spline basis (\code{P = 3}) and chooses the number of basis functions from the number of unique values in \code{x}.
#'
#' Boundary limits are taken from \code{x} unless supplied explicitly. Knot locations may be given directly as a numeric vector, or generated either at equally spaced locations (\code{"eq"}) or at empirical quantiles (\code{"quantile"}).
#'
#' @param x A numeric vector of predictor values.
#' @param P A non-negative integer giving the spline degree. \code{P = 3} corresponds to a cubic B-spline basis.
#' @param K An integer giving the number of basis functions to return. The default increases slowly with the number of unique values in \code{x}.
#' @param limits A numeric vector of length 2 giving the lower and upper boundary limits for the spline basis. Missing values are replaced by \code{min(x)} and \code{max(x)}.
#' @param knots Either a numeric vector of knot locations, or one of \code{"eq"} or \code{"quantile"}. If \code{"eq"}, knots are placed uniformly between \code{limits[1]} and \code{limits[2]}. If \code{"quantile"}, knots are placed at equally spaced empirical quantiles of \code{x}.
#'
#' @details
#' The returned basis has \code{length(x)} rows and \code{k} columns.
#'
#' When \code{knots} is generated internally, the function first creates \code{K - P + 1} knot locations and then augments them with repeated boundary knots so the recursion can be evaluated.
#'
#' @return A numeric matrix with one row per element of \code{x} and one column per spline basis function.
#'
#' @examples
#' x <- seq(0, 1, length.out = 10)
#'
#' # Default cubic basis
#' B <- s(x)
#' dim(B)
#'
#' # Equally spaced knots with custom basis size
#' B2 <- s(x, K = 5, knots = "eq")
#'
#' # Quantile-based knots
#' B3 <- s(x, knots = "quantile")
#'
#' @export
s <- function(x, P = 3, K = min(7, max(3, floor(log2(length(unique(x)))))), limits = c(NA, NA), knots = "eq") {
  config <- get.config.s(x, P, K, limits, knots)
  P <- config$P
  K <- config$K
  limits <- config$limits
  knots <- config$knots

  knots <- c(rep(limits[1], P), knots, rep(limits[2], P))
  N <- length(x)
  spline <- array(x, c(N, K + P))
  for (k in 1:(K + P)) {
    spline[, k] <- (x >= knots[k]) & ((x < knots[k + 1]) | ((k == (K + P)) & (x == knots[k + 1])))
  }
  for (p in 1:P) {
    for (k in 1:(K + P - p)) {
      spline[, k] <-
        spline[, k] * (if (knots[k + p] == knots[k]) {
          0
        } else {
          (x - knots[k]) / (knots[k + p] - knots[k])
        }) +
        spline[, k + 1] * (if (knots[k + p + 1] == knots[k + 1]) {
          0
        } else {
          (knots[k + p + 1] - x) / (knots[k + p + 1] - knots[k + 1])
        })
    }
    spline <- spline[, -(K + P - p + 1), drop = FALSE]
  }
  spline
}

#' Compute a spline configuration
#'
#' Derive the configuration needed to evaluate \code{s()} for a given input vector. This helper resolves boundary limits, chooses or validates knot locations, augments the knot vector with repeated boundary knots, and returns the arguments in a form that can be used to rebuild or update spline calls.
#'
#' This function is intended for programmatic use, for example when rewriting model formulas that contain calls to \code{s(...)}.
#'
#' @param x A numeric vector of predictor values.
#' @param P A non-negative integer giving the spline degree. \code{P = 3} corresponds to a cubic B-spline basis.
#' @param K An integer giving the number of basis functions implied by the spline specification.
#' @param limits A numeric vector of length 2 giving the lower and upper boundary limits for the spline basis. Missing values are replaced by \code{min(x)} and \code{max(x)}.
#' @param knots Either a numeric vector of knot locations, or one of \code{"eq"} or \code{"quantile"}. If \code{"eq"}, knots are placed uniformly between \code{limits[1]} and \code{limits[2]}. If \code{"quantile"}, knots are placed at equally spaced empirical quantiles of \code{x}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{x}{The unevaluated expression supplied as \code{x}, returned via \code{substitute(x)}.}
#'   \item{knots}{The full augmented knot vector, including repeated boundary knots.}
#'   \item{limits}{The resolved lower and upper boundary limits.}
#'   \item{P}{The spline degree.}
#'   \item{K}{The number of basis functions.}
#' }
#'
#' @keywords internal
get.config.s <- function(x, P = 3, K = min(7, max(3, floor(log2(length(unique(x)))))), limits = c(NA, NA), knots = "eq") {
  tol <- sqrt(.Machine$double.eps) * max(1, abs(x))
  limits <- if.na(limits, c(min(x) - tol, max(x) + tol))

  if (is.numeric(knots)) {
    if (length(knots) != K - P + 1) {
      stop(paste0("Incorrect number of knots. Expected ", K - P + 1, ", got ", length(knots), "."))
    } else if ((max(knots) < max(x)) | (min(knots) > min(x))) {
      stop(paste0("The range of the knots does not cover the data. knots range from ", min(knots), " to ", max(knots), ", data ranges from ", min(x) - tol, " to ", max(x) + tol, "."))
    }
  } else {
    if (knots == "eq") {
      knots <- seq.int(limits[1], limits[2], length.out = K - P + 1)
      knots[1] <- limits[1]
      knots[K - P + 1] <- limits[2]
    } else if (knots == "quantile") {
      if (!all(limits == c(min(x) - tol, max(x) + tol))) {
        warning("Knots defined by quantiles. Ignoring the limits argument.")
        limits <- c(min(x), max(x))
      }
      knots <- seq.int(0, 1, length.out = K - P + 1)
      knots <- quantile(x, knots, na.rm = TRUE)
    } else {
      stop("Invalid knots argument.")
    }
  }
  return(list(x = substitute(x), knots = knots, limits = limits, P = P, K = K))
}

#' Check if a dlm block has the treatment as covariate
#'
#' @param formula A formula describing the model passed to the fit.surr function.
#' @param name.G The name of the surrogate
#'
#' @keywords internal
check_has_G <- function(formula, name.G, inside_dlm = FALSE) {
  terms <- as.character(attr(terms(formula), "variables"))[-1]
  for (term in terms) {
    dlm_flag <- check_is_dlm_block(term)
    if (inside_dlm | dlm_flag) {
      if (dlm_flag) {
        call.term <- as.call(str2lang(term))
      } else {
        call.term <- list(X = term)
      }
      for (arg in names(call.term)) {
        if (arg == "X") {
          if (is_formula(call.term[[arg]])) {
            check_has_G(as.formula(call.term[[arg]]), name.G, inside_dlm = TRUE)
          } else if (call.term[[arg]] == name.G) {
            stop("Treatment cannot be used inside DLM blocks. If you need time dependent treatment effect or iterations, use a iteraction with the times index.")
          }
        }
      }
    }
  }
}

#' Check if a formula term is a dlm block
#'
#' @param term A string.
#'
#' @keywords internal
check_is_dlm_block <- function(term) {
  any(sapply(c("har(", "ffs(", "reg(", "pol(", "AR(", "TF(", "noise("), function(x) {
    grepl(x, term, fixed = TRUE)
  }))
}
