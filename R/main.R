#' Fit marginal and conditional state-space models for longitudinal surrogate evaluation
#'
#' Fits two Gaussian state-space models (Dynamic Linear Models) to jointly longitudinal outcome data: (i) a marginal model for the outcome trajectory given treatment and time, and (ii) a conditional model that additionally adjusts for a user-specified surrogate structure. The function returns per-time treatment-effect estimates from both models and subject-level bootstrap draws obtained via subject-level resampling.
#'
#' The implementation follows a two-model decomposition used for estimating longitudinal treatment effects and surrogate-adjusted (residual) treatment effects in a state-space framework.
#'
#' See \insertCite{santos2026causalframeworkevaluatingjointly}{OnlineSurr} for details on the methodology.
#'
#' See \insertCite{West-DLM}{OnlineSurr} for best practices on model specification in the state-space model setting.
#'
#' @param formula An object of class \code{formula} describing the fixed-effects mean structure for the primary outcome. The left-hand side must be the outcome variable. Internally, the right-hand side is augmented to include treatment-by-time fixed effects.
#' @param id A variable (unquoted) identifying subjects. Each subject must have at most one measurement per \code{time} value.
#' @param surrogate A formula describing the surrogate structure to be included in the conditional model. May be provided either as a \code{formula} (e.g., \code{~ s1 + s2}) or as a string that can be coerced to a formula.
#' @param treat A variable (unquoted) indicating treatment assignment. Must encode exactly two treatment levels after coercion to a factor.
#' @param data A \code{data.frame} containing all variables referenced in \code{formula}, \code{id}, \code{treat}, \code{surrogate}, and (optionally) \code{time}.
#' @param time Optional variable (unquoted) giving the measurement time index. Must be numeric and equally spaced across observed time points. If \code{NULL}, an equally spaced within-subject index is created in the current row order (with a warning).
#' @param N.boots Integer number of subject-level bootstrap replicates. Each replicate resamples subjects with replacement and recombines subject-specific sufficient quantities to form bootstrap draws of the fixed effects.
#' @param verbose Logical scalar indicating whether to print progress information during model fitting. If \code{TRUE}, progress updates are shown; if \code{FALSE}, no progress output is produced.
#' @param D.local Numeric, a number between 0 and 1 indicating the discount factor to be used for the random effect block. This factor controls how smooth the random effect evolve over time. A discount factor of 1 means that the random effects do not change over time, so that each individual has its own local level, but that level is the same for all times. A discount factor of 0 is not acceptable (the kDGLM package will replace it by 1), but values closer to 0 imply in a more flexible dynamic. See \insertCite{West-DLM}{OnlineSurr} or the appendix in \insertCite{santos2026causalframeworkevaluatingjointly}{OnlineSurr} for instructions on how to specify the discount factor.
#'
#' @return An object of class \code{"fitted_onlinesurr"}: a named list with elements \code{$Marginal} and \code{$Conditional}. Each of these contains:
#'   \itemize{
#'     \item \code{point}: the point estimate vector of the fixed effects (excluding subject-specific random-walk states) at the final time point.
#'     \item \code{smp}: a matrix of bootstrap draws for those fixed effects, with one column per bootstrap replicate.
#'   }
#'   The object also includes:
#'   \itemize{
#'     \item \code{T}: number of unique time points.
#'     \item \code{N}: number of subjects.
#'     \item \code{n.fixed}: number of fixed-effect coefficients implied by \code{formula} for a single subject prior to stacking across subjects.
#'   }
#'
#' @details
#' \strong{Data requirements.} The data must have at most one row per subject-time pair; time must be numeric and equally spaced (or omitted, in which case an index is created). Treatment and subject identifiers are coerced to factors with sorted levels.
#'
#' \strong{Model structure.} The marginal model includes treatment-by-time fixed effects and a subject-specific random-walk component to capture within-subject correlation. The conditional model adds the user-specified surrogate structure to the design, and checks that treatment is not a linear combination of the surrogate design (rank check).
#'
#' \strong{Bootstrap.} Subjects are resampled with replacement. Subject-specific filtered quantities are computed once and recombined in each bootstrap iteration to reduce computational cost, consistent with a subject-level nonparametric bootstrap strategy for replicated time series.
#'
#' @examples
#' \dontrun{
#' # data columns: y (outcome), id (subject id), trt (0/1 or two-level factor),
#' # time (numeric equally spaced), s1 and s2 (surrogates)
#'
#' fit <- fit.surr(
#'   formula    = y ~ 1, # baseline fixed effects; function adds trt*time terms
#'   id         = id,
#'   surrogate  = ~ s1 + s2,
#'   treat      = trt,
#'   data       = dat,
#'   time       = time,
#'   N.boots    = 500
#' )
#'
#' # Access point estimates and bootstrap samples
#' fit$Marginal$point
#' fit$Conditional$smp[, 1:10]
#' }
#'
#' @import kDGLM
#' @import dplyr
#' @import tidyr
#' @import rlang
#' @importFrom stats update.formula model.matrix as.formula model.frame na.pass var terms
#' @export
fit.surr <- function(formula, id, surrogate, treat, data = NULL, time = NULL, N.boots = 2000, verbose = 1, D.local = 0.8) {
  family <- Normal

  # Possible errors:
  # variables are not what they're supposed to be


  Y <- model.frame(update.formula(formula, ~1), data = data, na.action = na.pass)
  name.Y <- names(Y)[1]
  name.id <- deparse(substitute(id))
  name.G <- deparse(substitute(treat))
  if (name.Y == "") {
    stop("Argument 'Y' is missing.")
  }
  if (name.id == "") {
    stop("Argument 'id' is missing.")
  } else {
    if (is.null(data[[name.id]])) {
      stop("Invalid 'id' argument.")
    } else {
      if (any(is.na(data[[name.id]]))) {
        stop("Cannot have NAs as id.")
      }
    }
  }
  if (name.G == "") {
    stop("Argument 'treat' is missing.")
  } else if (is.null(data[[name.G]])) {
    stop("Invalid 'treat' argument.")
  } else if (any(is.na(data[[name.G]]))) {
    stop("Cannot have NAs as treatment.")
  }

  data <- data %>%
    mutate(
      !!name.id := factor(!!sym(name.id), levels = sort(unique(!!sym(name.id)))),
      !!name.G := factor(!!sym(name.G), levels = sort(unique(!!sym(name.G))))
    )

  if (is.null(time)) {
    data <- data %>%
      arrange(name.id) %>%
      group_by(!!sym(name.id)) %>%
      mutate(Time = seq_along(!!sym(name.id)))
    name.T <- "Time"
    warning("Time index was not specified. We assume that the data is equally spaced AND already ordered.")
  } else {
    name.T <- deparse(substitute(time))
  }

  counts <- data %>%
    mutate(n = 1) %>%
    group_by(!!sym(name.id), !!sym(name.T)) %>%
    summarize(n = sum(n), .groups = "keep")
  if (any(counts$n > 1)) {
    stop("Multiple measurements for the same individual-time combination. Have you passed the right variables?")
  }

  if (is.numeric(data[[name.T]])) {
    time.var <- data %>%
      group_by(!!sym(name.id)) %>%
      mutate(diff = c(0, diff(!!sym(name.T)))) %>%
      filter(diff > 0)
    if (var(time.var$diff) > 0) {
      stop("Data is not equally spaced. If this is due to missing measurements, you can fill the missing times as NA.")
    }
    data <- data %>%
      mutate(!!name.T := factor(!!sym(name.T), levels = sort(unique(!!sym(name.T))))) %>%
      arrange(name.id, name.T)
  } else {
    stop("Time index is not numeric.")
  }


  counts.G <- table(data[[name.G]])
  n.treat <- length(counts.G[counts.G > 0])
  if (n.treat != 2) {
    if (n.treat > 2) {
      stop(paste0("Incorrect number of treatments. Expected 2, but ", n.treat, " were found."))
    } else {
      if (n.treat == 1) {
        stop(paste0("Incorrect number of treatments. Expected 2, but only ", n.treat, " was found."))
      } else {
        stop(paste0("Incorrect number of treatments. Expected 2, but none were found."))
      }
    }
  }


  Y.mat <- (data %>% select(name.Y, name.id, name.T) %>%
    pivot_wider(names_from = name.id, values_from = name.Y))[, -1]
  treat.list <- (data %>% group_by(!!sym(name.id)) %>%
    summarize(treat = (!!sym(name.G))[1]))$treat
  control.name <- (unique(treat.list)[1])
  treat.name <- (unique(treat.list)[2])
  treat.list <- (treat.list == treat.name)

  T <- dim(Y.mat)[1]
  N <- dim(Y.mat)[2]
  pat.names <- names(Y.mat)

  random.base <- polynomial_block(mu = 1, order = 1, D = D.local, R1 = 900)
  random.effects <- random.base * N

  pred.names <- random.effects$pred.names

  check_has_G(formula, name.G)

  formula <- update_s_in_formula(formula, data = data)
  formula <- update.formula(formula, eval(paste0(". ~ . + ", name.G, "*", name.T, "  - ", name.G, " -1")))
  fixed.effects <- list()
  fixed.effects[[1]] <- formula.to.structure(formula, data %>% filter(!!sym(name.id) == pat.names[1]), label = pred.names[1])
  structure.base <- random.effects + fixed.effects[[1]]


  if (!is_formula(surrogate)) {
    surrogate <- as.formula(paste0("~", surrogate))
  }
  check_has_G(surrogate, name.G)
  surrogate <- update_s_in_formula(surrogate, data = data)
  surrogate <- update.formula(surrogate, eval(paste0("~ . - 1")))
  surrogate.effects <- list()

  surrogate.effects[[1]] <- formula.to.structure(surrogate, data %>% filter(!!sym(name.id) == pat.names[1]), pred.names[1])


  surrogate.base <- surrogate.effects[[1]]
  surrogate.base$FF <- array(surrogate.base$FF, c(dim(surrogate.base$FF)[c(1, 3)], N)) %>% aperm(c(1, 3, 2))

  colnames(surrogate.base$FF) <- pred.names
  surrogate.base$FF.labs <- matrix(surrogate.base$FF.labs, dim(surrogate.base$FF.labs)[1], N)
  surrogate.base$k <- N
  surrogate.base$pred.names <- pred.names

  n.ref <- fixed.effects[[1]]$n
  for (i in 2:N) {
    fixed.effects[[i]] <- formula.to.structure(formula, data %>% filter(!!sym(name.id) == pat.names[i]), label = pred.names[i])
    surrogate.effects[[i]] <- formula.to.structure(surrogate, data %>% filter(!!sym(name.id) == pat.names[i]), pred.names[i])

    structure.base$FF[-(1:N), i, ] <- fixed.effects[[i]]$FF
    surrogate.base$FF[, i, ] <- surrogate.effects[[i]]$FF
  }

  test.F <- surrogate.base$FF
  dim(test.F) <- c(dim(test.F)[1], dim(test.F)[3] * dim(test.F)[2])
  test.F <- cbind(t(test.F), 1)
  test.G <- (data %>% filter(!!sym(name.T) == (!!sym(name.T))[1]))[[name.G]] %>% as.numeric()

  tol_rank <- 1e-12

  rA <- qr(test.F, tol = tol_rank)$rank
  rAv <- qr(cbind(test.F, test.G), tol = tol_rank)$rank

  if (rA == rAv) {
    stop("The treatment is a linear combination of the surrogate.")
  }

  models <- list()
  K <- 1
  boots.idx.mat <- matrix(sample.int(N, N * N.boots, replace = TRUE), N, N.boots)
  for (marg in c(TRUE, FALSE)) {
    args <- list()
    args[[name.Y]] <- Normal(mu = pred.names, V = diag(N), data = Y.mat)
    args$smooth <- FALSE
    args$safe.mode <- FALSE
    if (marg) {
      structure <- structure.base
    } else {
      structure <- structure.base + surrogate.base
    }
    args$structure <- structure
    model <- do.call(fit_model, args)
    W.est <- model$W
    R1.est.control <- var(model$mt[(1:N)[!treat.list], T])
    R1.est.treat <- var(model$mt[(1:N)[treat.list], T])


    mt <- matrix(NA, structure$n - N + 1, N)
    St <- array(NA, c(structure$n - N + 1, structure$n - N + 1, N))

    start <- Sys.time()
    for (i in 1:N) {
      args <- list()
      args[[name.Y]] <- Normal(mu = pred.names[i], V = 1, data = Y.mat[, i])
      args$smooth <- FALSE
      args$safe.mode <- FALSE

      structure.cur <- (random.base %>% block_rename(pred.names[i])) + fixed.effects[[i]]


      if (!marg) {
        structure.cur <- structure.cur + surrogate.effects[[i]]
      }
      idx <- -((1:N)[-i])
      structure.cur$D[, , ] <- 0
      structure.cur$H[, , ] <- W.est[idx, idx, ]
      structure.cur$R1[1, 1] <- if (treat.list[i]) {
        R1.est.treat
      } else {
        R1.est.control
      }
      structure.cur$R1[-1, -1] <- structure.cur$R1[-1, -1] * N
      args$structure <- structure.cur
      model.cur <- do.call(fit_model, args)

      # max(abs(model$W[,,T]-diag(diag(model$W[,,T]))))

      St[, , i] <- S.cur <- model.cur$Ct[, , T] %>% ginv()
      mt[, i] <- S.cur %*% model.cur$mt[, T]

      if (verbose) {
        spent <- difftime(Sys.time(), start, units = "mins")
        ETA <- (N - i) * spent / i
        cat(paste0(format(paste0(ifelse(marg, "Marginal", "Conditional"), " model"), width = 17, justify = "r"), " - pre-processing: ", format(round(100 * i / N, 2), width = 6, nsmall = 2), "% - ETA: ", round(ETA, 2), "mins                   \r"))
      }
    }

    mt.smp <- matrix(NA, structure$n - N, N.boots)
    rownames(mt.smp) <- row.names(model.cur$mt)[-1]
    start <- Sys.time()

    boots.idx <- 1:N

    mt.remix <- mt[, boots.idx, drop = FALSE]
    St.remix <- St[, , boots.idx, drop = FALSE]

    vals.diag <- ifelse(St.remix[1, 1, ] == 0, 1, St.remix[1, 1, ])
    A_11 <- rowSums(St.remix[-1, -1, ], dims = 2)
    A_12 <- St.remix[1, -1, ]
    S_11 <- A_11 - tcrossprod(sweep(A_12, 2, vals.diag, FUN = "/"), A_12)

    a_1 <- rowSums(mt.remix[-1, ])
    a_2 <- mt.remix[1, ] / vals.diag

    point.est <- solve(S_11, a_1 - A_12 %*% a_2)

    for (i in 1:N.boots) {
      boots.idx <- boots.idx.mat[, i]
      # boots.idx=1:N

      mt.remix <- mt[, boots.idx, drop = FALSE]
      St.remix <- St[, , boots.idx, drop = FALSE]

      vals.diag <- ifelse(St.remix[1, 1, ] == 0, 1, St.remix[1, 1, ])
      A_11 <- rowSums(St.remix[-1, -1, ], dims = 2)
      A_12 <- St.remix[1, -1, ]
      S_11 <- A_11 - tcrossprod(sweep(A_12, 2, vals.diag, FUN = "/"), A_12)

      a_1 <- rowSums(mt.remix[-1, ])
      a_2 <- mt.remix[1, ] / vals.diag

      mt.smp[, i] <- solve(S_11, a_1 - A_12 %*% a_2)

      if (verbose) {
        spent <- difftime(Sys.time(), start, units = "mins")
        ETA <- (N.boots - i) * spent / i
        cat(paste0(format(paste0(ifelse(marg, "Marginal", "Conditional"), " model"), width = 17, justify = "r"), " - bootstraping:   ", format(round(100 * i / N.boots, 2), width = 6, nsmall = 2), "% - ETA: ", round(ETA, 2), "mins                   \r"))
      }
    }

    if (marg) {
      treat.effect <- rep(NA, T)
      treat.effect.smp <- matrix(NA, N.boots, T)
      theta <- point.est
      theta.smp <- mt.smp
      vals.control <- model$FF[-(1:N), !treat.list, , drop = FALSE]
      vals.treat <- model$FF[-(1:N), treat.list, , drop = FALSE]
      for (t in 1:T) {
        treat.effect[t] <- mean(t(vals.treat[, , t]) %*% theta) -
          mean(t(vals.control[, , t]) %*% theta)
        treat.effect.smp[, t] <- colMeans(t(vals.treat[, , t]) %*% theta.smp) -
          colMeans(t(vals.control[, , t]) %*% theta.smp)
      }
    } else {
      N.control <- 2 * sum(!treat.list)
      data.control <- rbind(
        data %>%
          filter(!!sym(name.G) == control.name),
        data %>%
          filter(!!sym(name.G) == control.name) %>%
          mutate(!!sym(name.G) := treat.name)
      )
      char.len <- floor(log10(N.control)) + 1
      pred.control <- paste0("mu.", formatC(1:N.control, width = char.len, flag = "0"))
      data.control[[name.id]] <- rep(paste0("mu.", formatC(1:N.control, width = char.len, flag = "0")), each = T)
      data.control[[name.id]] <- factor(data.control[[name.id]], levels = unique(data.control[[name.id]]))
      pat.names.control <- levels(data.control[[name.id]])
      treat.list.control <- (data.control %>% group_by(!!sym(name.id)) %>%
        summarize(treat = (!!sym(name.G))[1]))$treat
      treat.list.control <- (treat.list.control == treat.name)

      fixed.effects.control <- list()
      fixed.effects.control[[1]] <- formula.to.structure(formula, data.control %>% filter(!!sym(name.id) == pat.names.control[1]), label = pred.control[1])

      structure.base.control <- fixed.effects.control[[1]]
      structure.base.control$FF <- array(structure.base.control$FF, c(dim(structure.base.control$FF)[c(1, 3)], N.control)) %>% aperm(c(1, 3, 2))


      colnames(structure.base.control$FF) <- pred.control
      structure.base.control$FF.labs <- matrix(structure.base.control$FF.labs, dim(structure.base.control$FF.labs)[1], N.control)
      structure.base.control$k <- N.control
      structure.base.control$pred.names <- pred.control

      surrogate.effects.control <- list()
      surrogate.effects.control[[1]] <- formula.to.structure(surrogate, data.control %>% filter(!!sym(name.id) == pat.names.control[1]), pred.control[1])

      surrogate.base.control <- surrogate.effects.control[[1]]
      surrogate.base.control$FF <- array(surrogate.base.control$FF, c(dim(surrogate.base.control$FF)[c(1, 3)], N.control)) %>% aperm(c(1, 3, 2))

      colnames(surrogate.base.control$FF) <- pred.control
      surrogate.base.control$FF.labs <- matrix(surrogate.base.control$FF.labs, dim(surrogate.base.control$FF.labs)[1], N.control)
      surrogate.base.control$k <- N.control
      surrogate.base.control$pred.names <- pred.control

      for (i in 2:N.control) {
        fixed.effects.control[[i]] <- formula.to.structure(formula, data.control %>% filter(!!sym(name.id) == pat.names.control[i]), label = pred.control[i])
        surrogate.effects.control[[i]] <- formula.to.structure(surrogate, data.control %>% filter(!!sym(name.id) == pat.names.control[i]), label = pred.control[i])

        structure.base.control$FF[, i, ] <- fixed.effects.control[[i]]$FF
        surrogate.base.control$FF[, i, ] <- surrogate.effects.control[[i]]$FF
      }
      structure.control <- structure.base.control + surrogate.base.control

      treat.effect <- rep(NA, T)
      treat.effect.smp <- matrix(NA, N.boots, T)
      theta <- point.est
      theta.smp <- mt.smp
      vals.control <- structure.control$FF[, !treat.list.control, , drop = FALSE]
      vals.treat <- structure.control$FF[, treat.list.control, , drop = FALSE]

      for (t in 1:T) {
        treat.effect[t] <- mean(t(vals.treat[, , t]) %*% theta) -
          mean(t(vals.control[, , t]) %*% theta)
        treat.effect.smp[, t] <- colMeans(t(vals.treat[, , t]) %*% theta.smp) -
          colMeans(t(vals.control[, , t]) %*% theta.smp)
      }
    }


    # models[[ifelse(marg, "Marginal", "Conditional")]] <- list(point = model$mt[-(1:N), T], smp = mt.smp)
    models[[ifelse(marg, "Marginal", "Conditional")]] <- list(point = treat.effect, smp = t(treat.effect.smp))
  }

  models$marg

  models <- append(models, list(
    "T" = T,
    "N" = N,
    "n.fixed" = n.ref
  ))
  class(models) <- "fitted_onlinesurr"
  return(models)
}


#' Test time-homogeneity of the PTE
#'
#' Tests the null hypothesis that the LPTE is constant over time. The test is based on the difference between the conditional and marginal treatment-effect trajectories implied by a fitted \code{"fitted_onlinesurr"} object, standardized by an estimated covariance, and uses a max-type statistic to control the family wise error across time points.
#'
#' See \insertCite{santos2026causalframeworkevaluatingjointly}{OnlineSurr} for the theoretical details about this test.
#'
#' @param model A fitted object of class \code{"fitted_onlinesurr"}, typically returned by \code{fit.surr}. Must contain \code{$T}, \code{$n.fixed}, and the elements \code{$Marginal} and \code{$Conditional} with \code{point} and \code{smp} components.
#' @param signif.level Numeric in (0,1) giving the test significance level used to form the critical value from the bootstrap distribution. Default is \code{0.05}.
#'
#' @param N.boots Integer number of Monte Carlo draws used to approximate the null distribution of the max standardized deviation statistic and to compute the p-value. Default is \code{50000}.
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{T}: the observed test statistic (maximum absolute standardized deviation).
#'     \item \code{T.crit}: the 1-signif.level critical value.
#'     \item \code{p.value}: the Monte Carlo p-value \code{mean(T_null > T_obs)}.
#'   }
#'
#' @details
#'
#' Notes:
#' \itemize{
#'   \item The function assumes the first \code{T} time-specific treatment-effect parameters are stored contiguously at the beginning of \code{model$Marginal$point} and \code{model$Conditional$point} (and similarly for \code{smp}). It uses the index \code{1:(n.fixed)} as implemented in the code: \code{1:(T + n.fixed - T)}.
#'   \item \code{N.boots} here is a Monte Carlo size for the null simulation (distinct from the bootstrap size used when fitting \code{model}).
#' }
#'
#' @importFrom Rfast colMaxs
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit.surr(y ~ 1,
#'   id = id, surrogate = ~ s1 + s2, treat = trt,
#'   data = dat, time = time, N.boots = 2000
#' )
#'
#' time_homo_test(fit, signif.level = 0.05, N.boots = 50000)
#' }
time_homo_test <- function(model, signif.level = 0.05, N.boots = 50000) {
  T <- model$T
  n <- model$n.fixed

  delta.est <- model$Marginal$point
  delta.R.est <- model$Conditional$point

  delta.smp <- t(model$Marginal$smp)
  delta.R.smp <- t(model$Conditional$smp)

  pte <- 1 - sum(delta.R.est) / sum(delta.est)
  pte.smp <- 1 - colSums(delta.R.smp) / colSums(delta.smp)

  err <- delta.R.smp - (1 - pte.smp) * delta.smp

  th <- delta.R.est - (1 - pte) * delta.est
  var.th <- var(err)

  # _______________________________ MSD __________________________________
  A <- rmvnorm(N.boots, rep(0, T), var.th)
  sd.th <- sqrt(diag(var.th))
  sigma <- matrix(sd.th, T, N.boots)
  test.smp <- colMaxs(abs(A / sigma), value = TRUE)
  crit.val <- quantile(test.smp, 1 - signif.level)
  test.obs <- max(abs(th / sd.th))
  return(list(T = test.obs, T.crit = crit.val, p.value = mean(test.smp > test.obs)))
}
