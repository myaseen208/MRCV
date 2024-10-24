#' @name    print.summary.genloglin
#' @aliases print.summary.genloglin
#' @title Control Printed Display of MRCV Summary Modeling Objects
#' @description
#'   A method function that controls the display of \code{summary.genloglin}
#' @param x An object of class \code{"genloglin"}, produced by the
#'   \code{\link{genloglin}} function.
#' @param digits Minimum number of digits; see \code{\link{print.default}} for additional explanation.
#' @param symbolic.cor A logical value indicating whether correlations should be printed in symbolic form; see the \code{\link{summary}} method for class \code{\link{glm}} for additional explanation.
#' @param signif.stars A logical value indicating whether significance stars should be printed; see the \code{\link{summary}} method for class \code{\link{glm}} for additional explanation.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @export
#' @rdname print
print.summary.genloglin <- function(x, digits = max(3, getOption("digits") - 3),
                                    symbolic.cor = x$symbolic.cor, signif.stars =
                                      getOption("show.signif.stars"), ...) {
  op <- options()
  on.exit(options(op))
  options(scipen = 5)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
    sep = ""
  )
  cat("Deviance Residuals: \n")
  if (x$df.residual > 5) {
    x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
    names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
  }
  xx <- zapsmall(x$deviance.resid, digits + 1)
  print.default(xx, digits = digits, na.print = "", print.gap = 2)
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    df <- if ("df" %in% names(x)) {
      x[["df"]]
    } else {
      NULL
    }
    if (!is.null(df) && (nsingular <- df[3L] - df[1L])) {
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
        sep = ""
      )
    } else {
      cat("\nCoefficients:\n")
    }
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs,
      digits = digits, signif.stars = signif.stars,
      na.print = "NA", ...
    )
  }
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
    format(x$dispersion), ")\n\n", apply(cbind(paste(format(c(
      "Null",
      "    Residual"
    ), justify = "right"), "deviance:"), format(unlist(x[c(
      "null.deviance",
      "deviance"
    )]), digits = max(5, digits + 1))), 1L, paste, collapse = " "),
    sep = ""
  )
  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }
  cat("\n", "        Number of Fisher Scoring iterations: ", x$iter, "\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      } else {
        correl <- format(round(correl, 2), nsmall = 2, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}
