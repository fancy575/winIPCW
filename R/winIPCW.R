#' Restricted win statistics with conditional tie weighting
#'
#' Estimate restricted win statistics for prioritized time-to-event composite
#' endpoints under right censoring. The function computes the verified-tie IPCW
#' estimator and/or the proposed conditional tie-weighting estimator.
#'
#' @param formula A right-censored survival formula, for example
#'   `survival::Surv(time, status) ~ Z1 + Z2 + Z3`.
#' @param data A long-format data frame. It must contain one follow-up row
#'   per subject with `endpoint == 0`, plus event rows for prioritized
#'   endpoints. The follow-up row encodes the censoring/follow-up time, and
#'   positive endpoint rows encode component-specific event times.
#' @param id,trt,endpoint Column names for subject id, treatment arm, and
#'   endpoint type. Treatment must be coded `1` for treatment and `0` for
#'   control.
#' @param tau Numeric vector of restriction times.
#' @param method Estimators to compute: `"ipcw"` for the verified-tie IPCW
#'   estimator and `"m-ipcw"` for the conditional tie-weighting estimator.
#' @param estimand Summaries to report: `"NB"` for net benefit, `"WR"` for
#'   win ratio, and/or `"WO"` for win odds.
#' @param priority Endpoint priority order, most important first. If `NULL`,
#'   all positive endpoint values are used in increasing order.
#' @param censor Censoring working model: `"cox"` or `"km"`.
#' @param margin Event-margin working model. Only `"cox"` is currently
#'   supported.
#' @param copula Working copula family for the conditional tie estimator.
#'   Available options are `"gumbel"`, `"clayton"`, `"frank"`, and
#'   `"plackett"`.
#' @param eps Lower bound used in inverse censoring weights.
#' @param conf.level Confidence level for Wald intervals.
#'
#' @return An object of class `"winIPCW"` with components:
#' \describe{
#'   \item{table}{A data frame containing estimates, standard errors, Wald
#'   confidence intervals, and the estimated treatment/control win
#'   probabilities.}
#'   \item{fits}{A list of fitted nuisance models.}
#'   \item{spec}{A list recording the analysis specification.}
#' }
#'
#' @examples
#' utils::data("windat", package = "winIPCW")
#'
#' ## Use a small subset so the example runs quickly during R CMD check.
#' follow <- windat[windat$event_type == 0, ]
#' id_t <- follow$id[follow$A == 1]
#' id_c <- follow$id[follow$A == 0]
#' id_keep <- c(
#'   id_t[seq_len(min(80, length(id_t)))],
#'   id_c[seq_len(min(80, length(id_c)))]
#' )
#' windat_ex <- windat[windat$id %in% id_keep, ]
#'
#' fit <- winIPCW(
#'   survival::Surv(time, status) ~ Z1 + Z2 + Z3,
#'   data = windat_ex,
#'   id = "id",
#'   trt = "A",
#'   endpoint = "event_type",
#'   tau = 36,
#'   method = c("ipcw", "m-ipcw"),
#'   estimand = "NB",
#'   priority = c(1, 2),
#'   censor = "cox",
#'   margin = "cox",
#'   copula = "gumbel"
#' )
#'
#' results(fit)
#'
#' @export
winIPCW <- function(formula,
                  data,
                  id,
                  trt,
                  endpoint,
                  tau,
                  method = c("ipcw", "m-ipcw"),
                  estimand = c("NB", "WR", "WO"),
                  priority = NULL,
                  censor = c("cox", "km"),
                  margin = "cox",
                  copula = c("gumbel", "clayton", "frank", "plackett"),
                  eps = 0.01,
                  conf.level = 0.95) {
  if (is.null(method)) {
    method <- c("ipcw", "m-ipcw")
  }
  method <- tolower(as.character(method))
  method[method == "ct"] <- "m-ipcw"
  method <- .arg_match(method, c("ipcw", "m-ipcw"), several = TRUE, name = "method")
  estimand <- .arg_match(estimand, c("NB", "WR", "WO"), several = TRUE, name = "estimand")
  censor <- .arg_match(censor, c("cox", "km"), name = "censor")
  copula <- .std_copula(copula)

  eps <- as.numeric(eps)
  if (length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a positive finite scalar.", call. = FALSE)
  }
  conf.level <- as.numeric(conf.level)
  if (length(conf.level) != 1L || !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be a finite scalar between 0 and 1.", call. = FALSE)
  }
  if (!identical(margin, "cox")) {
    stop("Only margin = 'cox' is currently supported. KM event margins require new event-nuisance influence-function and variance terms.", call. = FALSE)
  }
  tau <- sort(unique(as.numeric(tau)))
  if (!length(tau) || any(!is.finite(tau)) || any(tau <= 0)) {
    stop("`tau` must be a positive numeric vector.", call. = FALSE)
  }

  rows <- list()
  fits <- list()
  for (tt in tau) {
    pack <- .pack_winIPCW(formula, data, id = id, trt = trt, endpoint = endpoint,
                        tau = tt, priority = priority)
    cens_fit <- .fit_censor(pack, censor)
    margins <- NULL

    if ("ipcw" %in% method) {
      est <- .pair_ipcw(pack, cens_fit, eps)
      rows[[length(rows) + 1L]] <- .rows_from_est(est, tt, "ipcw", "none", estimand, conf.level)
    }

    if ("m-ipcw" %in% method) {
      margins <- .fit_margins(pack)
      for (cf in copula) {
        theta <- .fit_copula_theta(pack, margins, cf)
        est <- .pair_m_ipcw(pack, cens_fit, margins, cf, theta, eps)
        rows[[length(rows) + 1L]] <- .rows_from_est(est, tt, "m-ipcw", cf, estimand, conf.level)
        fits[[paste(tt, "m-ipcw", cf, sep = "_")]] <- list(censor = cens_fit, margins = margins, theta = theta)
      }
    }
    fits[[paste(tt, "ipcw", sep = "_")]] <- list(censor = cens_fit)
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  out <- list(
    call = match.call(),
    table = tab,
    fits = fits,
    spec = list(
      tau = tau,
      method = method,
      estimand = estimand,
      priority = priority,
      censor = censor,
      margin = margin,
      copula = copula,
      eps = eps,
      conf.level = conf.level
    )
  )
  class(out) <- "winIPCW"
  out
}

.rows_from_est <- function(est, tau, method, copula, estimand, conf.level) {
  do.call(rbind, lapply(estimand, function(e) {
    estimate <- est[[e]]
    se <- switch(e,
      NB = est$se_NB,
      WR = est$se_logWR,
      WO = est$se_logWO
    )
    ci <- .ci_one(estimate, se, e, conf.level)
    data.frame(
      tau = tau,
      method = method,
      copula = copula,
      estimand = e,
      estimate = as.numeric(estimate)[1L],
      se = as.numeric(se)[1L],
      lower = as.numeric(ci[1L]),
      upper = as.numeric(ci[2L]),
      pi_t = as.numeric(est$pi_t)[1L],
      pi_c = as.numeric(est$pi_c)[1L],
      stringsAsFactors = FALSE
    )
  }))
}

#' Extract result tables from a winIPCW fit
#'
#' @description
#' `results()` returns the analysis table stored in a fitted [`winIPCW()`]
#' object. The table contains the requested restricted-time win summaries,
#' standard errors, confidence intervals, and the estimated treatment and
#' control win probabilities. `tidy.winIPCW()` is provided as a convenience
#' method with the same return value.
#'
#' @param x A fitted object of class `"winIPCW"`.
#' @param ... Additional arguments, currently unused.
#'
#' @return A data frame with one row per restriction time, estimator, copula
#'   specification, and requested estimand.
#'
#' @examples
#' utils::data("windat", package = "winIPCW")
#'
#' follow <- windat[windat$event_type == 0, ]
#' id_t <- follow$id[follow$A == 1]
#' id_c <- follow$id[follow$A == 0]
#' id_keep <- c(
#'   id_t[seq_len(min(40, length(id_t)))],
#'   id_c[seq_len(min(40, length(id_c)))]
#' )
#' windat_ex <- windat[windat$id %in% id_keep, ]
#'
#' fit <- winIPCW(
#'   survival::Surv(time, status) ~ Z1 + Z2 + Z3,
#'   data = windat_ex,
#'   id = "id",
#'   trt = "A",
#'   endpoint = "event_type",
#'   tau = 36,
#'   method = "ipcw",
#'   estimand = "NB",
#'   priority = c(1, 2),
#'   censor = "cox"
#' )
#' results(fit)
#'
#' @export
results <- function(x, ...) {
  UseMethod("results")
}

#' @rdname results
#' @export
results.winIPCW <- function(x, ...) {
  x$table
}

#' @rdname results
#' @export
tidy.winIPCW <- function(x, ...) {
  x$table
}
