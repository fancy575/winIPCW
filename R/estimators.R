.pair_ipcw <- function(pack, censor, eps) {
  trt <- which(pack$trt == 1L)
  ctl <- which(pack$trt == 0L)
  K <- ncol(pack$Y)
  pt <- pc <- numeric(K)
  pair <- matrix(0, nrow = length(trt), ncol = length(ctl))
  pair_loss <- pair

  for (ii in seq_along(trt)) {
    i <- trt[ii]
    Xi <- pack$X[i, , drop = FALSE]
    for (jj in seq_along(ctl)) {
      j <- ctl[jj]
      Xj <- pack$X[j, , drop = FALSE]
      for (q in seq_len(K)) {
        gate <- TRUE
        if (q > 1L) {
          gate <- all(pack$Y[i, seq_len(q - 1L)] >= pack$tau - eps &
                        pack$Y[j, seq_len(q - 1L)] >= pack$tau - eps)
        }
        if (!gate) next
        if (pack$Y[i, q] > pack$Y[j, q] && pack$D[j, q] == 1L) {
          G1 <- .surv_predict(censor[["1"]], pack$Y[j, q], Xi)
          G0 <- .surv_predict(censor[["0"]], pack$Y[j, q], Xj)
          val <- 1 / pmax(G1 * G0, eps)
          pt[q] <- pt[q] + val
          pair[ii, jj] <- val
          break
        }
        if (pack$Y[j, q] > pack$Y[i, q] && pack$D[i, q] == 1L) {
          G1 <- .surv_predict(censor[["1"]], pack$Y[i, q], Xi)
          G0 <- .surv_predict(censor[["0"]], pack$Y[i, q], Xj)
          val <- 1 / pmax(G1 * G0, eps)
          pc[q] <- pc[q] + val
          pair_loss[ii, jj] <- val
          break
        }
      }
    }
  }
  denom <- length(trt) * length(ctl)
  if (denom == 0L) {
    stop("Both treatment arms must contain at least one subject.", call. = FALSE)
  }
  .finish_estimate(pt / denom, pc / denom, pair, pair_loss, pack$priority)
}

.ct_ratio_gt <- function(margins, family, theta, u, t, X) {
  S1tau <- .surv_predict(margins[[1L]], rep(max(u, t), length(t)), X)
  S1u <- .surv_predict(margins[[1L]], u, X)
  S2t <- .surv_predict(margins[[2L]], t, X)
  .copula_C(S1tau, S2t, theta, family) / pmax(.copula_C(S1u, S2t, theta, family), 1e-8)
}

.ct_ratio_eq <- function(margins, family, theta, u, t, X) {
  S1tau <- .surv_predict(margins[[1L]], rep(max(u, t), length(t)), X)
  S1u <- .surv_predict(margins[[1L]], u, X)
  S2t <- .surv_predict(margins[[2L]], t, X)
  .copula_C2(S1tau, S2t, theta, family) / pmax(.copula_C2(S1u, S2t, theta, family), 1e-8)
}

.pair_m_ipcw <- function(pack, censor, margins, copula, theta, eps) {
  if (ncol(pack$Y) != 2L) {
    stop("The current M-IPCW implementation supports two prioritized survival endpoints.", call. = FALSE)
  }
  trt <- which(pack$trt == 1L)
  ctl <- which(pack$trt == 0L)
  pt <- pc <- numeric(2L)
  pair <- matrix(0, nrow = length(trt), ncol = length(ctl))
  pair_loss <- pair

  for (ii in seq_along(trt)) {
    i <- trt[ii]
    Xi <- pack$X[i, , drop = FALSE]
    for (jj in seq_along(ctl)) {
      j <- ctl[jj]
      Xj <- pack$X[j, , drop = FALSE]

      if (pack$Y[i, 1L] > pack$Y[j, 1L] && pack$D[j, 1L] == 1L) {
        G1 <- .surv_predict(censor[["1"]], pack$Y[j, 1L], Xi)
        G0 <- .surv_predict(censor[["0"]], pack$Y[j, 1L], Xj)
        val <- 1 / pmax(G1 * G0, eps)
        pt[1L] <- pt[1L] + val
        pair[ii, jj] <- val
        next
      }
      if (pack$Y[j, 1L] > pack$Y[i, 1L] && pack$D[i, 1L] == 1L) {
        G1 <- .surv_predict(censor[["1"]], pack$Y[i, 1L], Xi)
        G0 <- .surv_predict(censor[["0"]], pack$Y[i, 1L], Xj)
        val <- 1 / pmax(G1 * G0, eps)
        pc[1L] <- pc[1L] + val
        pair_loss[ii, jj] <- val
        next
      }

      no_death_i <- pack$D[i, 1L] == 0L
      no_death_j <- pack$D[j, 1L] == 0L
      if (!no_death_i || !no_death_j) next

      if (pack$Y[i, 2L] > pack$Y[j, 2L] && pack$D[j, 2L] == 1L) {
        t <- pack$Y[j, 2L]
        ri <- .ct_ratio_gt(margins[["1"]], copula, theta[["1"]], pack$Y[i, 1L], t, Xi)
        rj <- .ct_ratio_eq(margins[["0"]], copula, theta[["0"]], pack$Y[j, 1L], t, Xj)
        G1 <- .surv_predict(censor[["1"]], t, Xi)
        G0 <- .surv_predict(censor[["0"]], t, Xj)
        val <- .clip01(ri, 0) * .clip01(rj, 0) / pmax(G1 * G0, eps)
        pt[2L] <- pt[2L] + val
        pair[ii, jj] <- val
        next
      }
      if (pack$Y[j, 2L] > pack$Y[i, 2L] && pack$D[i, 2L] == 1L) {
        t <- pack$Y[i, 2L]
        ri <- .ct_ratio_eq(margins[["1"]], copula, theta[["1"]], pack$Y[i, 1L], t, Xi)
        rj <- .ct_ratio_gt(margins[["0"]], copula, theta[["0"]], pack$Y[j, 1L], t, Xj)
        G1 <- .surv_predict(censor[["1"]], t, Xi)
        G0 <- .surv_predict(censor[["0"]], t, Xj)
        val <- .clip01(ri, 0) * .clip01(rj, 0) / pmax(G1 * G0, eps)
        pc[2L] <- pc[2L] + val
        pair_loss[ii, jj] <- val
      }
    }
  }
  denom <- length(trt) * length(ctl)
  if (denom == 0L) {
    stop("Both treatment arms must contain at least one subject.", call. = FALSE)
  }
  .finish_estimate(pt / denom, pc / denom, pair, pair_loss, pack$priority)
}

.finish_estimate <- function(ptk, pck, pair_win, pair_loss, priority) {
  pi_t <- sum(ptk)
  pi_c <- sum(pck)
  nb <- pi_t - pi_c
  out <- list(
    pi_t = pi_t,
    pi_c = pi_c,
    NB = nb,
    WR = .wr(pi_t, pi_c),
    WO = .wo(pi_t, pi_c),
    pi_tk = stats::setNames(ptk, paste0("pi_t", priority)),
    pi_ck = stats::setNames(pck, paste0("pi_c", priority))
  )
  se <- .projection_se(pair_win, pair_loss, out)
  out$se_NB <- se["NB"]
  out$se_logWR <- se["logWR"]
  out$se_logWO <- se["logWO"]
  out
}

.projection_se <- function(pair_win, pair_loss, est) {
  n1 <- nrow(pair_win)
  n0 <- ncol(pair_win)
  if (n1 < 2L || n0 < 2L) {
    return(c(NB = NA_real_, logWR = NA_real_, logWO = NA_real_))
  }
  h <- pair_win - pair_loss
  psi1 <- rowMeans(h) - est$NB
  psi0 <- colMeans(h) - est$NB
  v_nb <- stats::var(psi1) / n1 + stats::var(psi0) / n0
  se_nb <- sqrt(pmax(v_nb, 0))

  wt <- pair_win
  wl <- pair_loss
  p1t <- rowMeans(wt) - est$pi_t
  p0t <- colMeans(wt) - est$pi_t
  p1c <- rowMeans(wl) - est$pi_c
  p0c <- colMeans(wl) - est$pi_c
  v_t <- stats::var(p1t) / n1 + stats::var(p0t) / n0
  v_c <- stats::var(p1c) / n1 + stats::var(p0c) / n0
  c_tc <- stats::cov(p1t, p1c) / n1 + stats::cov(p0t, p0c) / n0

  se_lwr <- if (est$pi_t > 0 && est$pi_c > 0) {
    sqrt(pmax(v_t / est$pi_t^2 + v_c / est$pi_c^2 - 2 * c_tc / (est$pi_t * est$pi_c), 0))
  } else NA_real_
  g <- 2 / pmax(1 - est$NB^2, 1e-8)
  se_lwo <- abs(g) * se_nb
  c(NB = se_nb, logWR = se_lwr, logWO = se_lwo)
}
