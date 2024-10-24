#' @name    print.genloglin
#' @aliases print.genloglin
#' @title Control Printed Display of MRCV Regression Modeling Objects
#' @description
#' Method functions that control the printed display of MRCV regression
#' modeling objects.
#' 
#' @details
#' The \code{print.genloglin} function is based on the \code{\link{print}}
#' method for class \code{"glm"}.
#'
#' The \code{print.summary.genloglin} function is based on the
#' \code{\link{print}} method for class \code{"summary.glm"}.
#'
#' @aliases print.genloglin print.summary.genloglin print.anova.genloglin
#' print.predict.genloglin
#' @param x An object of class \code{"genloglin"} produced by the
#' \code{\link{genloglin}} function.
#' @param digits Minimum number of digits; see \code{\link{print.default}} for
#' additional explanation.
#' @param ...  Additional arguments passed to or from other methods.
#' @return A print out of selected results.
#' @export
#' @rdname print
print.genloglin <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  op <- options()
  on.exit(options(op))
  options(scipen = 5)
  x <- x$mod.fit
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
    sep = ""
  )
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts)) {
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    }
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
      print.gap = 2,
      quote = FALSE
    )
  } else {
    cat("No coefficients\n\n")
  }
  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }
  cat(
    "Null Deviance:\t   ", format(signif(x$null.deviance, digits)),
    "\nResidual Deviance:", format(signif(x$deviance, digits)), "\n"
  )
  invisible(x)
}
