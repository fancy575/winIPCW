#' Example data for restricted win statistics
#'
#' A simulated two-component prioritized survival data set for illustrating
#' restricted win statistics with right censoring. The data are stored in long
#' format. Each subject has one follow-up/censoring row with `event_type = 0`;
#' subjects with observed component events before the restriction time have
#' additional event rows with `event_type = 1` or `event_type = 2`.
#'
#' @name windat
#' @docType data
#' @usage data(windat)
#' @format A data frame in long format with columns:
#' \describe{
#'   \item{\code{id}}{Subject identifier.}
#'   \item{\code{A}}{Treatment arm indicator; \code{1} for treatment and
#'     \code{0} for control.}
#'   \item{\code{Z1}}{Binary baseline covariate.}
#'   \item{\code{Z2}}{Continuous baseline covariate.}
#'   \item{\code{Z3}}{Binary baseline covariate.}
#'   \item{\code{event_type}}{Endpoint type. \code{0} is the follow-up or
#'     censoring row, \code{1} is the first-priority endpoint, and \code{2} is
#'     the second-priority endpoint.}
#'   \item{\code{priority}}{Endpoint priority order; \code{0} for the
#'     follow-up/censoring row.}
#'   \item{\code{label}}{Text label for the endpoint type.}
#'   \item{\code{time}}{Observed row-specific time used in
#'     \code{survival::Surv(time, status)}.}
#'   \item{\code{status}}{Row-specific event indicator. For event rows,
#'     \code{1} indicates that the endpoint event was observed. For the
#'     follow-up row, \code{1} indicates censoring before the restriction time.}
#'   \item{\code{tau}}{Restriction time used to construct the observed data.}
#' }
#'
#' @details
#' The data were simulated from Weibull proportional hazards margins linked by
#' a Gumbel--Hougaard dependence structure, with covariate-dependent censoring.
#' The object is intended as a reproducible example for the package API.
#'
#' @source Simulated for the package examples.
#'
#' @examples
#' data(windat)
#' head(windat)
#'
#' @keywords datasets
"windat"
