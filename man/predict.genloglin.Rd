\name{predict.genloglin}
\alias{predict.genloglin}
\title{
Calculate Observed and Model-Predicted Odds Ratios for MRCV Data
}
\description{
The \code{predict.genloglin} method function calculates observed and model-predicted odds ratios and their confidence intervals using results from \code{\link{genloglin}}.  It offers an asymptotic normal approximation for estimating the confidence intervals for the observed and model-predicted odds ratios, and a bootstrap approach for estimating the confidence intervals for the model-predicted odds ratios.
}
\usage{
\method{predict}{genloglin}(object, alpha = 0.05, pair = "WY", print.status = TRUE, ...)
}
\arguments{
  \item{object}{
An object of class \code{'genloglin'} produced by the \code{\link{genloglin}} function.
}
  \item{alpha}{
The desired alpha level.  The \code{predict.genloglin} function provides two-sided (1-\code{alpha})x100\% confidence intervals.
}
  \item{pair}{
For the case of three MRCVs, a character string specifying the pair of items for which odds ratios will be calculated: \code{"WY"} indicates odds ratios should be calculated for each (Wi, Yj) pair conditional on the response for each Zk, \code{"WZ"} indicates conditional odds ratios should be calculated for each (Wi, Zk) pair, and \code{"YZ"} indicates conditional odds ratios should be calculated for each (Yj, Zk) pair.
}
  \item{print.status}{
A logical value indicating whether bootstrap progress updates should be provided.
}
  \item{\dots}{
Additional arguments passed to or from other methods.
}
}
\details{
Wald confidence intervals are estimated for both model-based (see Appendix A of Bilder and Loughin, 2007) and observed (see Agresti, 2013, p. 70) odds ratios.

A bootstrap method is also available which provides bias-corrected accelerated (BCa) confidence intervals for the model-predicted odds ratios.  See Efron (1987) for more information about BCa intervals.  The \code{predict.genloglin} function uses a jackknife approximation for estimating the empirical influence values.

The bootstrap confidence intervals are available only when \code{boot = TRUE} in the original call to the \code{\link{genloglin}} function.
}
\value{
--- A list containing at least \code{original.arg}, \code{OR.obs}, and \code{OR.model.asymp}. 

\code{original.arg} is a list containing the following objects:
\itemize{
\item{\code{data}: The original data frame supplied to the \code{data} argument.}
\item{\code{I}: The original value supplied to the \code{I} argument.}
\item{\code{J}: The original value supplied to the \code{J} argument.}
\item{\code{K}: The original value supplied to the \code{K} argument.}
\item{\code{nvars}: The number of MRCVs.}
\item{\code{alpha}: The original value supplied to the \code{alpha} argument.}
}

\code{OR.obs} is a numeric matrix. For the two MRCV case, the matrix contains IxJ rows corresponding to the IxJ possible pairs (Wi, Yj) and 3 columns, where column 1 corresponds to the observed odds ratio for (Wi, Yj) and columns 2 and 3 correspond to the estimated lower and upper confidence bounds, respectively.  For the three MRCV case, the matrix contains 2xIxJxK rows corresponding to all possible combinations of \code{pair} conditional on the response for each item of the 3rd MRCV, and 3 columns as described for the 2 MRCV case.

\code{OR.model.asymp} is a numeric matrix similar to \code{OR.obs} but where column 1 corresponds to the model-predicted odds ratios and columns 2 and 3 correspond to the estimated lower and upper confidence bounds, respectively, using an asymptotic normal approximation.

--- For \code{boot = TRUE} in the call to the \code{\link{genloglin}} function, the primary list additionally includes \code{boot.results}, a list containing the following objects:
\itemize{
\item{\code{B.use}: The number of bootstrap resamples used.}
\item{\code{B.discard}: The number of bootstrap resamples discarded due to having at least one item with all positive or negative responses.}
\item{\code{OR.model.BCa}: A numeric matrix similar to \code{OR.obs} but where column 1 corresponds to the model-predicted odds ratios and columns 2 and 3 correspond to the estimated lower and upper confidence bounds, respectively, of the BCa intervals.}
}
}
\references{
Agresti, A. (2013)  \emph{Categorical data analysis (3rd ed.)}.  Hoboken, New Jersey: John Wiley & Sons.

Bilder, C. and Loughin, T. (2007)  Modeling association between two or more categorical variables that allow for multiple category choices.  \emph{Communications in Statistics--Theory and Methods}, \bold{36}, 433--451.

Efron, B. (1987)  Better bootstrap confidence intervals.  \emph{Journal of the American Statistical Association}, \bold{82}, 171--185.
}
\examples{
## For examples see help(genloglin).
}