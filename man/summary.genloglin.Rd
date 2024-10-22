\name{summary.genloglin}
\alias{summary.genloglin}
\title{
Summarize Two or Three MRCV Model Fit Information
}
\description{
The \code{summary.genloglin} function summarizes model fit information provided by the \code{\link{genloglin}} function.
}
\usage{
\method{summary}{genloglin}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{'genloglin'} produced by the \code{\link{genloglin}} function.
}
  \item{\dots}{
Additional arguments passed to or from other methods.
}
}
\details{
The \code{summary.genloglin} function is based on the \code{\link{summary}} method for class \code{"glm"} with a few modifications.  The \code{coefficients} object contains Rao-Scott second-order adjusted standard errors, z-values, and p-values.  The \code{cov.unscaled} object contains the Rao-Scott second-order adjusted covariance matrix of the estimated coefficients.

The deviance information printed by \code{summary.genloglin} should not be used to conduct traditional model comparison tests.  The \code{\link{anova.genloglin}} function offers adjusted tests.
}
\value{
The \code{summary.genloglin} function returns the same list returned by the \code{\link{summary}} method for class \code{"glm"} with the exception of AIC.
}
\examples{
## For examples see help(genloglin).
}