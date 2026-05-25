.copulas <- c("indep", "gumbel", "clayton", "frank", "plackett")

.std_copula <- function(x) {
  x <- tolower(as.character(x))
  x[x %in% c("independence", "independent")] <- "indep"
  .arg_match(x, .copulas, several = TRUE, name = "copula")
}

.copula_theta_default <- function(family) {
  switch(family,
    indep = 1,
    gumbel = 2,
    clayton = 2,
    frank = 5.736283,
    plackett = 5
  )
}

.copula_C <- function(u, v, theta, family) {
  u <- .clip01(u)
  v <- .clip01(v)
  if (family == "indep") {
    return(u * v)
  }
  if (family == "gumbel") {
    if (abs(theta - 1) < 1e-8) return(u * v)
    A <- (-log(u))^theta + (-log(v))^theta
    return(exp(-A^(1 / theta)))
  }
  if (family == "clayton") {
    if (theta < 1e-8) return(u * v)
    B <- pmax(u^(-theta) + v^(-theta) - 1, 1e-12)
    return(B^(-1 / theta))
  }
  if (family == "frank") {
    if (abs(theta) < 1e-8) return(u * v)
    eu <- exp(-theta * u)
    ev <- exp(-theta * v)
    e <- exp(-theta)
    B <- 1 + (eu - 1) * (ev - 1) / (e - 1)
    return(-log(pmax(B, 1e-12)) / theta)
  }
  if (family == "plackett") {
    if (abs(theta - 1) < 1e-8) return(u * v)
    B <- 1 + (theta - 1) * (u + v)
    D <- pmax(B^2 - 4 * theta * (theta - 1) * u * v, 1e-12)
    return((B - sqrt(D)) / (2 * (theta - 1)))
  }
}

.copula_C2 <- function(u, v, theta, family) {
  u <- .clip01(u)
  v <- .clip01(v)
  if (family == "indep") {
    return(u)
  }
  if (family == "gumbel") {
    if (abs(theta - 1) < 1e-8) return(u)
    A <- (-log(u))^theta + (-log(v))^theta
    C <- exp(-A^(1 / theta))
    return(C * A^(1 / theta - 1) * (-log(v))^(theta - 1) / v)
  }
  if (family == "clayton") {
    if (theta < 1e-8) return(u)
    B <- pmax(u^(-theta) + v^(-theta) - 1, 1e-12)
    return(v^(-theta - 1) * B^(-1 / theta - 1))
  }
  if (family == "frank") {
    if (abs(theta) < 1e-8) return(u)
    eu <- exp(-theta * u)
    ev <- exp(-theta * v)
    e <- exp(-theta)
    B <- 1 + (eu - 1) * (ev - 1) / (e - 1)
    return(ev * (eu - 1) / ((e - 1) * pmax(B, 1e-12)))
  }
  if (family == "plackett") {
    if (abs(theta - 1) < 1e-8) return(u)
    B <- 1 + (theta - 1) * (u + v)
    D <- pmax(B^2 - 4 * theta * (theta - 1) * u * v, 1e-12)
    return(0.5 * (1 - (B - 2 * theta * u) / sqrt(D)))
  }
}

.copula_density <- function(u, v, theta, family) {
  h <- 1e-5
  lo <- pmax(v - h, 1e-8)
  hi <- pmin(v + h, 1 - 1e-8)
  (.copula_C2(u, hi, theta, family) - .copula_C2(u, lo, theta, family)) /
    pmax(hi - lo, 1e-8)
}
