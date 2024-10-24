# summary.genloglin()
# A method function that summarizes results given by genloglin()



#' Summarize Two or Three MRCV Model Fit Information
#'
#' The \code{summary.genloglin} function summarizes model fit information
#' provided by the \code{\link{genloglin}} function.
#'
#' The \code{summary.genloglin} function is based on the \code{\link{summary}}
#' method for class \code{"glm"} with a few modifications.  The
#' \code{coefficients} object contains Rao-Scott second-order adjusted standard
#' errors, z-values, and p-values.  The \code{cov.unscaled} object contains the
#' Rao-Scott second-order adjusted covariance matrix of the estimated
#' coefficients.
#'
#' The deviance information printed by \code{summary.genloglin} should not be
#' used to conduct traditional model comparison tests.  The
#' \code{\link{anova.genloglin}} function offers adjusted tests.
#'
#' @param object An object of class \code{'genloglin'} produced by the
#' \code{\link{genloglin}} function.
#' @param \dots Additional arguments passed to or from other methods.
#' @return The \code{summary.genloglin} function returns the same list returned
#' by the \code{\link{summary}} method for class \code{"glm"} with the
#' exception of AIC.
#' @examples
#'
#' ## For examples see help(genloglin).
#'
#' @export
#' @rdname summary
summary.genloglin <- function(object, ...) {
  op <- options()
  on.exit(options(op))
  options(scipen = 5)
  sum.fit <- object$sum.fit
  class(sum.fit) <- "summary.genloglin"
  sum.fit
}
