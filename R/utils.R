.arg_match <- function(x, choices, several = FALSE, name = deparse(substitute(x))) {
  if (missing(x) || is.null(x)) {
    return(choices[1L])
  }
  x <- as.character(x)
  ok <- x %in% choices
  if (!all(ok)) {
    stop(sprintf(
      "`%s` must be one of %s.",
      name,
      paste(sprintf("'%s'", choices), collapse = ", ")
    ), call. = FALSE)
  }
  if (!several && length(x) != 1L) {
    stop(sprintf("`%s` must have length 1.", name), call. = FALSE)
  }
  unique(x)
}

.as_name <- function(x, name = deparse(substitute(x))) {
  if (is.character(x) && length(x) == 1L) {
    return(x)
  }
  stop(sprintf("`%s` must be a single column name.", name), call. = FALSE)
}

.clip01 <- function(x, eps = 1e-10) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

.wr <- function(pi_t, pi_c) {
  ifelse(is.finite(pi_c) & pi_c > 0, pi_t / pi_c, NA_real_)
}

.wo <- function(pi_t, pi_c) {
  nb <- pi_t - pi_c
  ifelse(is.finite(nb) & abs(nb) < 1, (1 + nb) / (1 - nb), NA_real_)
}

.ci_one <- function(est, se, estimand, level) {
  ## `se` is often extracted from a named vector, e.g. se["NB"].
  ## If names are not removed, R creates names such as "lower.NB" and
  ## "upper.NB", so later extraction by ci["lower"] returns NA.
  est <- as.numeric(est)[1L]
  se <- as.numeric(se)[1L]
  level <- as.numeric(level)[1L]
  estimand <- as.character(estimand)[1L]

  if (!is.finite(est) || !is.finite(se) || !is.finite(level) ||
      level <= 0 || level >= 1 || se < 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z <- stats::qnorm(1 - (1 - level) / 2)

  if (estimand %in% c("WR", "WO")) {
    if (est <= 0) {
      return(c(lower = NA_real_, upper = NA_real_))
    }
    if (se == 0) {
      return(c(lower = est, upper = est))
    }
    return(c(
      lower = exp(log(est) - z * se),
      upper = exp(log(est) + z * se)
    ))
  }

  c(lower = est - z * se, upper = est + z * se)
}

.pretty_method <- function(x) {
  out <- x
  out[out == "ipcw"] <- "IPCW"
  out[out == "m-ipcw"] <- "M-IPCW"
  out
}

.pretty_copula <- function(x) {
  out <- x
  out[is.na(out) | out == "none"] <- "-"
  labs <- c(indep = "Independence", gumbel = "Gumbel", clayton = "Clayton",
            frank = "Frank", plackett = "Plackett")
  hit <- out %in% names(labs)
  out[hit] <- unname(labs[out[hit]])
  out
}
