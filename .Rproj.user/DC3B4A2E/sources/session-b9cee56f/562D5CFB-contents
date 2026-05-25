.fit_surv_arm <- function(time, status, X, model) {
  time <- as.numeric(time)
  status <- as.integer(status)
  X <- as.matrix(X)
  n <- length(time)
  if (model == "km" || ncol(X) == 0L) {
    sf <- survival::survfit(survival::Surv(time, status) ~ 1)
    return(list(type = "km", fit = sf, beta = numeric(), X = X))
  }
  df <- data.frame(.time = time, .status = status, X, check.names = FALSE)
  names(df)[-(1:2)] <- paste0("x", seq_len(ncol(X)))
  frm <- stats::as.formula(paste("survival::Surv(.time, .status) ~", paste(names(df)[-(1:2)], collapse = " + ")))
  fit <- survival::coxph(frm, data = df, ties = "breslow", x = TRUE, y = TRUE, model = FALSE)
  list(type = "cox", fit = fit, beta = stats::coef(fit), X = X)
}

.surv_predict <- function(obj, time, Xnew) {
  time <- as.numeric(time)
  Xnew <- as.matrix(Xnew)
  if (obj$type == "km") {
    ss <- summary(obj$fit, times = time, extend = TRUE)$surv
    return(pmax(ss, 1e-8))
  }
  bh <- survival::basehaz(obj$fit, centered = FALSE)
  idx <- findInterval(time, bh$time)
  H0 <- ifelse(idx == 0L, 0, bh$hazard[idx])
  beta <- stats::coef(obj$fit)
  eta <- if (length(beta)) as.vector(Xnew %*% beta) else rep(0, nrow(Xnew))
  pmax(exp(-H0 * exp(eta)), 1e-8)
}

.fit_censor <- function(pack, model) {
  out <- vector("list", 2L)
  names(out) <- c("0", "1")
  for (a in c(0L, 1L)) {
    idx <- which(pack$trt == a)
    out[[as.character(a)]] <- .fit_surv_arm(
      time = pack$follow_time[idx],
      status = pack$cens_status[idx],
      X = pack$X[idx, , drop = FALSE],
      model = model
    )
  }
  out
}

.fit_margins <- function(pack) {
  out <- vector("list", 2L)
  names(out) <- c("0", "1")
  for (a in c(0L, 1L)) {
    idx <- which(pack$trt == a)
    arm <- vector("list", ncol(pack$Y))
    for (k in seq_len(ncol(pack$Y))) {
      arm[[k]] <- .fit_surv_arm(
        time = pack$Y[idx, k],
        status = pack$D[idx, k],
        X = pack$X[idx, , drop = FALSE],
        model = "cox"
      )
    }
    out[[as.character(a)]] <- arm
  }
  out
}

.fit_copula_theta <- function(pack, margins, family) {
  if (family == "indep" || ncol(pack$Y) != 2L) {
    return(c("0" = 1, "1" = 1))
  }
  theta_grid <- switch(family,
    gumbel = seq(1.05, 8, length.out = 40),
    clayton = seq(0.05, 8, length.out = 40),
    frank = seq(0.05, 20, length.out = 40),
    plackett = seq(1.05, 30, length.out = 40)
  )
  ans <- c("0" = NA_real_, "1" = NA_real_)
  for (a in c(0L, 1L)) {
    idx <- which(pack$trt == a)
    if (!length(idx)) {
      ans[as.character(a)] <- .copula_theta_default(family)
      next
    }
    X <- pack$X[idx, , drop = FALSE]
    u <- .surv_predict(margins[[as.character(a)]][[1L]], pack$Y[idx, 1L], X)
    v <- .surv_predict(margins[[as.character(a)]][[2L]], pack$Y[idx, 2L], X)
    ll <- vapply(theta_grid, function(th) {
      d <- pmax(.copula_density(u, v, th, family), 1e-8)
      sum(log(d[pack$D[idx, 1L] == 1L & pack$D[idx, 2L] == 1L]), na.rm = TRUE)
    }, numeric(1L))
    ans[as.character(a)] <- theta_grid[which.max(ll)]
  }
  ans
}
