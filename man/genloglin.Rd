\name{genloglin}
\alias{genloglin}
\title{
Model the Association Among Two or Three MRCVs
}
\description{
The \code{genloglin} function uses a generalized loglinear modeling approach to estimate the association among two or three MRCVs.  Standard errors are adjusted using a second-order Rao-Scott approach.
}
\usage{
genloglin(data, I, J, K = NULL, model, add.constant = 0.5, boot = TRUE, 
    B = 1999, B.max = B, print.status = TRUE)
}
\arguments{
  \item{data}{
A data frame containing the raw data where rows correspond to the individual item response vectors, and columns correspond to the binary items, W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
}
  \item{I}{
The number of items corresponding to row variable W.
}
  \item{J}{
The number of items corresponding to column variable Y.
}
  \item{K}{
The number of items corresponding to strata variable Z.
}
  \item{model}{
For the two MRCV case, a character string specifying one of the following models: \code{"spmi"} (the complete independence model), \code{"homogeneous"} (the homogeneous association model), \code{"w.main"} (the w-main effects model), \code{"y.main"} (the y-main effects model), \code{"wy.main"} (the w-main and y-main effects model), or \code{"saturated"}.  Alternatively, a user-supplied formula can be specified, where the formula is limited to the generic variables \code{W}, \code{Y}, \code{wi}, \code{yj}, \code{count}, \code{W1},\ldots, \code{WI}, and \code{Y1},\ldots, \code{YJ}.  For the three MRCV case, only user-supplied formulas are accepted. In addition to the generic variables defined for two MRCVs, the formula may include the generic variables \code{Z}, \code{zk}, and \code{Z1},\ldots, \code{ZK}.
}
  \item{add.constant}{
A positive constant to be added to all zero marginal cell counts.
}
  \item{boot}{
A logical value indicating whether bootstrap resamples should be taken.
}
  \item{B}{
The desired number of bootstrap resamples.
}
  \item{B.max}{
The maximum number of bootstrap resamples.  Resamples for which at least one item has all positive or negative responses are thrown out; \code{genloglin} uses the first \code{B} valid resamples or all valid resamples if that number is less than \code{B}.
}
  \item{print.status}{
A logical value indicating whether progress updates should be provided.  When \code{print.status = TRUE}, the status of the IPF algorithm is printed after every 5 iterations.  Upon completion of the IPF algorithm, a progress bar appears that documents progress of the bootstrap.
}
}
\details{
The \code{genloglin} function first converts the raw data into a form that can be used for estimation.  For the two MRCV case, the reformatted data frame contains 2Ix2J rows and 5 columns generically named \code{W}, \code{Y}, \code{wi}, \code{yj}, and \code{count}.  For the three MRCV case, the reformatted data frame contains 2Ix2Jx2K rows and 7 columns generically named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, \code{zk}, and \code{count}.  Then, the model of interest is estimated by calling the \code{\link{glm}} function where the \code{family} argument is specified as \code{poisson}.  For all predictor variables, the first level is the reference group (i.e., 1 is the reference for variables \code{W}, \code{Y}, and \code{Z}, and 0 is the reference for variables \code{wi}, \code{yj}, and \code{zj}).

The \code{boot} argument must equal \code{TRUE} in order to obtain bootstrap results with the \code{genloglin} method functions.
}
\value{
--- \code{genloglin} returns an object of class \code{'genloglin'}.  The object is a list containing at least the following objects: \code{original.arg}, \code{mod.fit}, \code{sum.fit}, and \code{rs.results}.

\code{original.arg} is a list containing the following objects:
\itemize{
\item{\code{data}: The original data frame supplied to the \code{data} argument.}
\item{\code{I}: The original value supplied to the \code{I} argument.}
\item{\code{J}: The original value supplied to the \code{J} argument.}
\item{\code{K}: The original value supplied to the \code{K} argument.}
\item{\code{nvars}: The number of MRCVs.}
\item{\code{model}: The original value supplied to the \code{model} argument.}
\item{\code{add.constant}: The original value supplied to the \code{add.constant} argument.}
\item{\code{boot}: The original value supplied to the \code{boot} argument.}
}

\code{mod.fit} is a list containing the same objects returned by \code{\link{glm}} with a few modifications as described in \code{\link{summary.genloglin}}.

\code{sum.fit} is a list containing the same objects returned by the \code{\link{summary}} method for class \code{"glm"} with a few modifications as described in \code{\link{summary.genloglin}}.

\code{rs.results} is a list containing the following objects (see Appendix A of Bilder and Loughin, 2007, for more detail):
\itemize{
\item{\code{cov.mu}: The covariance matrix for the estimated cell counts.}
\item{\code{E}: The covariance matrix for the residuals.}
\item{\code{gamma}: Eigenvalues used in computing second-order Rao-Scott adjusted statistics.}
}

--- For \code{boot = TRUE}, the primary list additionally includes \code{boot.results}, a list containing the following objects:
\itemize{
\item{\code{B.use}: The number of bootstrap resamples used.}
\item{\code{B.discard}: The number of bootstrap resamples discarded due to having at least one item with all positive or negative responses.}
\item{\code{model.data.star}: For the two MRCV case, a numeric matrix containing 2Ix2J rows and \code{B.use}+4 columns, where the first 4 columns correspond to the model variables \code{W}, \code{Y}, \code{wi}, and \code{yj}, and the last \code{B.use} columns correspond to the observed counts for each resample.  For the three MRCV case, a numeric matrix containing 2Ix2Jx2K rows and \code{B.use+6} columns, where the first 6 columns correspond to the model variables \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, and \code{zk}, and the last \code{B.use} columns correspond to the observed counts for each resample.}
\item{\code{mod.fit.star}: For the two MRCV case, a numeric matrix containing \code{B.use} rows and 2Ix2J +1 columns, where the first 2Ix2J columns correspond to the model-predicted counts for each resample, and the last column corresponds to the residual deviance for each resample.  For the three MRCV case, a numeric matrix containing \code{B.use} rows and 2Ix2Jx2K+1 columns, where the first 2Ix2Jx2K columns correspond to the model-predicted counts for each resample, and the last column corresponds to the residual deviance for each resample.}
\item{\code{chisq.star}: A numeric vector of length \code{B.use} containing the Pearson statistics (comparing \code{model} to the saturated model) calculated for each resample.}
\item{\code{lrt.star}: A numeric vector of length \code{B.use} containing the LRT statistics calculated for each resample.}
\item{\code{residual.star}: A numeric matrix with 2Ix2J rows (or 2Ix2Jx2K rows for the three MRCV case) and \code{B.use} columns containing the residuals calculated for each resample.}
}
}
\references{
Bilder, C. and Loughin, T. (2007)  Modeling association between two or more categorical variables that allow for multiple category choices.  \emph{Communications in Statistics--Theory and Methods}, \bold{36}, 433--451.
}
\seealso{
The \code{genloglin} methods \code{\link{summary.genloglin}}, \code{\link{residuals.genloglin}}, \code{\link{anova.genloglin}}, and \code{\link{predict.genloglin}}, and the corresponding generic functions \code{\link{summary}}, \code{\link{residuals}}, \code{\link{anova}}, and \code{\link{predict}}.

The \code{\link{glm}} function for fitting generalized linear models.

The \code{\link{MI.test}} function for testing for MMI (one MRCV case) or SPMI (two MRCV case).
}
\examples{
# Estimate the y-main effects model for 2 MRCVs
mod.fit <- genloglin(data = farmer2, I = 3, J = 4, model = "y.main", boot = FALSE)
# Summarize model fit information
summary(mod.fit)
# Examine standardized Pearson residuals
residuals(mod.fit)
# Compare the y-main effects model to the saturated model
anova(mod.fit, model.HA = "saturated", type = "rs2")
# Obtain observed and model-predicted odds ratios
predict(mod.fit)

# Estimate a model that is not one of the named models
# Note that this was the final model chosen by Bilder and Loughin (2007)
mod.fit.other <- genloglin(data = farmer2, I = 3, J = 4, model = count ~ -1 + W:Y + 
    wi\%in\%W:Y + yj\%in\%W:Y + wi:yj + wi:yj\%in\%Y + wi:yj\%in\%W3:Y1, boot = 
    FALSE)
# Compare this model to the y-main effects model
anova(mod.fit, model.HA = count ~ -1 + W:Y + wi\%in\%W:Y + yj\%in\%W:Y + wi:yj + 
    wi:yj\%in\%Y + wi:yj\%in\%W3:Y1, type = "rs2", gof = TRUE)

# To obtain bootstrap results from the method functions the genloglin() boot 
# argument must be specified as TRUE (the default)
# A small B is used for demonstration purposes; normally, a larger B should be used
\donttest{
mod.fit <- genloglin(data = farmer2, I = 3, J = 4, model = "y.main", boot = TRUE,
    B = 99)
residuals(mod.fit)
anova(mod.fit, model.HA = "saturated", type = "all")
predict(mod.fit)
}

# Estimate a model for 3 MRCVs
\donttest{
mod.fit.three <- genloglin(data = farmer3, I = 3, J = 4, K = 5, model = count ~ 
    -1 + W:Y:Z + wi\%in\%W:Y:Z + yj\%in\%W:Y:Z + zk\%in\%W:Y:Z + wi:yj + 
    wi:yj\%in\%Y + wi:yj\%in\%W + wi:yj\%in\%Y:W + yj:zk + yj:zk\%in\%Z + 
    yj:zk\%in\%Y + yj:zk\%in\%Z:Y, boot = TRUE, B = 99)
residuals(mod.fit.three)
anova(mod.fit.three, model.HA = "saturated", type = "all")
predict(mod.fit.three, pair = "WY")}
}