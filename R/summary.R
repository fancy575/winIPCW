#' @export
print.winIPCW <- function(x, digits = 3, ...) {
  cat("Restricted win statistics with conditional ties\n\n")
  cat("Methods:  ", paste(.pretty_method(x$spec$method), collapse = ", "), "\n", sep = "")
  cat("Censor:   ", x$spec$censor, "\n", sep = "")
  cat("Margins:  ", x$spec$margin, " only\n", sep = "")
  if ("m-ipcw" %in% x$spec$method) {
    cat("Copulas:  ", paste(.pretty_copula(x$spec$copula), collapse = ", "), "\n", sep = "")
  }
  cat("Tau:      ", paste(x$spec$tau, collapse = ", "), "\n\n", sep = "")
  print(utils::head(.format_winIPCW_table(x$table, digits), 12L), row.names = FALSE)
  invisible(x)
}

#' @export
summary.winIPCW <- function(object, ...) {
  out <- list(
    call = object$call,
    spec = object$spec,
    table = object$table
  )
  class(out) <- "summary.winIPCW"
  out
}

#' @export
print.summary.winIPCW <- function(x, digits = 3, ...) {
  cat("Summary of restricted win statistics\n\n")
  print(.format_winIPCW_table(x$table, digits), row.names = FALSE)
  invisible(x)
}

.format_winIPCW_table <- function(tab, digits) {
  cols <- c("tau", "method", "copula", "estimand", "estimate", "se", "lower", "upper")
  missing_cols <- setdiff(cols, names(tab))
  if (length(missing_cols)) {
    stop(
      sprintf("Result table is missing required column(s): %s.", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- tab[, cols, drop = FALSE]
  out$method <- .pretty_method(out$method)
  out$copula <- .pretty_copula(out$copula)

  num <- vapply(out, is.numeric, logical(1L))
  out[num] <- lapply(out[num], function(z) round(z, digits = digits))

  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      z <- out[[nm]]
      formatted <- format(z, trim = TRUE, scientific = FALSE)
      formatted[is.na(z)] <- "not estimable"
      formatted[is.infinite(z) & z > 0] <- "Inf"
      formatted[is.infinite(z) & z < 0] <- "-Inf"
      out[[nm]] <- formatted
    } else {
      z <- as.character(out[[nm]])
      z[is.na(z)] <- ""
      out[[nm]] <- z
    }
  }
  out
}
