.parse_formula <- function(formula, data, trt) {
  trt <- .as_name(trt, "trt")
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv") || NCOL(y) != 2L) {
    stop("`formula` must have a right-censored response, e.g. Surv(time, status) ~ covariates.", call. = FALSE)
  }
  tt <- attr(stats::terms(formula), "term.labels")
  covars <- setdiff(tt, trt)
  list(time = as.numeric(y[, 1L]), status = as.integer(y[, 2L]), covars = covars)
}

.design_matrix <- function(data, covars) {
  if (length(covars) == 0L) {
    out <- matrix(0, nrow(data), 0L)
    colnames(out) <- character()
    return(out)
  }
  mm <- stats::model.matrix(stats::reformulate(covars), data = data)
  if (ncol(mm) && colnames(mm)[1L] == "(Intercept)") {
    mm <- mm[, -1L, drop = FALSE]
  }
  storage.mode(mm) <- "double"
  mm
}

.pack_winIPCW <- function(formula, data, id, trt, endpoint, tau, priority = NULL) {
  id <- .as_name(id, "id")
  trt <- .as_name(trt, "trt")
  endpoint <- .as_name(endpoint, "endpoint")
  parsed <- .parse_formula(formula, data, trt)

  needed <- c(id, trt, endpoint)
  miss <- setdiff(needed, names(data))
  if (length(miss)) {
    stop(sprintf("`data` is missing required column(s): %s.", paste(miss, collapse = ", ")), call. = FALSE)
  }

  dat <- data
  dat$.winIPCW_time <- parsed$time
  dat$.winIPCW_status <- parsed$status
  dat$.winIPCW_id <- dat[[id]]
  dat$.winIPCW_trt <- as.integer(dat[[trt]])
  dat$.winIPCW_endpoint <- as.integer(dat[[endpoint]])

  if (!all(dat$.winIPCW_trt %in% c(0L, 1L))) {
    stop("`trt` must be coded 0 for control and 1 for treatment.", call. = FALSE)
  }
  if (is.null(priority)) {
    priority <- sort(unique(dat$.winIPCW_endpoint[dat$.winIPCW_endpoint > 0L]))
  }
  priority <- as.integer(priority)
  if (length(priority) < 1L) {
    stop("At least one positive endpoint level is required.", call. = FALSE)
  }

  follow <- dat[dat$.winIPCW_endpoint == 0L, , drop = FALSE]
  if (!nrow(follow)) {
    stop("`data` must contain one follow-up/censoring row per subject with endpoint == 0.", call. = FALSE)
  }
  follow <- follow[order(follow$.winIPCW_id, follow$.winIPCW_time), , drop = FALSE]
  follow <- follow[!duplicated(follow$.winIPCW_id), , drop = FALSE]

  arm_counts <- table(factor(follow$.winIPCW_trt, levels = c(0L, 1L)))
  if (any(arm_counts < 2L)) {
    stop(
      "At least two control subjects and two treatment subjects are required for standard errors and confidence intervals.",
      call. = FALSE
    )
  }

  X <- .design_matrix(follow, parsed$covars)
  ids <- follow$.winIPCW_id
  n <- length(ids)
  K <- length(priority)
  Y <- matrix(rep(pmin(follow$.winIPCW_time, tau), K), nrow = n)
  D <- matrix(0L, nrow = n, ncol = K)

  for (k in seq_along(priority)) {
    sub <- dat[dat$.winIPCW_endpoint == priority[k] & dat$.winIPCW_status == 1L,
               c(".winIPCW_id", ".winIPCW_time"), drop = FALSE]
    if (!nrow(sub)) next
    sub <- sub[order(sub$.winIPCW_id, sub$.winIPCW_time), , drop = FALSE]
    sub <- sub[!duplicated(sub$.winIPCW_id), , drop = FALSE]
    idx <- match(sub$.winIPCW_id, ids)
    ok <- !is.na(idx)
    idx <- idx[ok]
    evtime <- sub$.winIPCW_time[ok]
    Y[idx, k] <- pmin(evtime, follow$.winIPCW_time[idx], tau)
    D[idx, k] <- as.integer(evtime <= pmin(follow$.winIPCW_time[idx], tau))
  }

  colnames(Y) <- paste0("Y", priority)
  colnames(D) <- paste0("d", priority)
  list(
    id = ids,
    trt = as.integer(follow$.winIPCW_trt),
    tau = as.numeric(tau),
    priority = priority,
    covars = parsed$covars,
    X = X,
    follow_time = pmin(as.numeric(follow$.winIPCW_time), tau),
    cens_status = as.integer(as.numeric(follow$.winIPCW_time) < tau & follow$.winIPCW_status == 1L),
    Y = Y,
    D = D,
    data = follow
  )
}
