\name{residuals.genloglin}
\alias{residuals.genloglin}
\title{
Calculate Standardized Pearson Residuals for MRCV Data
}
\description{
The \code{residuals.genloglin} method function calculates standardized Pearson residuals for the model specified in the \code{\link{genloglin}} function. It offers an asymptotic approximation and a bootstrap approximation for estimating the variance of the residuals.
}
\usage{
\method{residuals}{genloglin}(object, ...)
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
The bootstrap results are only available when \code{boot = TRUE} in the call to the \code{\link{genloglin}} function.

The \code{residuals.genloglin} function uses \code{\link[tables:tabular]{tables:tabular()}} to display the results for the two MRCV case.

See Bilder and Loughin (2007) for additional details about calculating the residuals.
}
\value{
--- A list containing at least \code{std.pearson.res.asymp.var}.  For the two MRCV case, the object is a 2Ix2J table of class \code{'tabular'} containing the standardized Pearson residuals based on the estimated asymptotic variance.  For the three MRCV case, the object is a data frame containing the 2Ix2Jx2K residuals.

--- For \code{boot = TRUE} in the call to the \code{\link{genloglin}} function, the list additionally includes:
\itemize{
\item{\code{B.use}: The number of bootstrap resamples used.}
\item{\code{B.discard}: The number of bootstrap resamples discarded due to having at least one item with all positive or negative responses.}
\item{\code{std.pearson.res.boot.var}: For the two MRCV case, a 2Ix2J table of class \code{'tabular'} containing the standardized Pearson residuals based on the bootstrap variance.  For the three MRCV case, a data frame containing the 2Ix2Jx2K residuals.}
}
}
\references{
Bilder, C. and Loughin, T. (2007)  Modeling association between two or more categorical variables that allow for multiple category choices.  \emph{Communications in Statistics--Theory and Methods}, \bold{36}, 433--451.
}
\examples{
## For examples see help(genloglin).
}