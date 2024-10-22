\name{anova.genloglin}
\alias{anova.genloglin}
\title{
Perform MRCV Model Comparison Tests
}
\description{
The \code{anova.genloglin} method function offers second-order Rao-Scott and bootstrap adjusted model comparison and goodness-of-fit (Pearson and LRT) statistics appropriate for evaluating models estimated by the \code{\link{genloglin}} function.
}
\usage{
\method{anova}{genloglin}(object, model.HA = "saturated", type = "all", gof = TRUE, 
    print.status = TRUE, ...)
}
\arguments{
  \item{object}{
An object of class \code{'genloglin'} produced by the \code{\link{genloglin}} function.
}
  \item{model.HA}{
For the two MRCV case, a character string specifying one of the following models to be compared to the null model (where the null model should be nested within the alternative model): \code{"homogeneous"} (the homogeneous association model), \code{"w.main"} (the w-main effects model), \code{"y.main"} (the y-main effects model), \code{"wy.main"} (the w- and y-main effects model), or \code{"saturated"}.  Alternatively, a user-supplied formula can be specified.  For the three MRCV case, only \code{"saturated"} or user-supplied formulas are accepted. 
}
  \item{type}{
A character string specifying one of the following approaches for performing adjusted model comparison tests: \code{"boot"} specifies a bootstrapping procedure; \code{"rs2"} specifies a Rao-Scott second-order adjustment; \code{"all"} specifies both approaches.
}
  \item{gof}{
A logical value indicating whether goodness-of-fit statistics should be calculated in addition to model comparison statistics.  For \code{model.HA = "saturated"}, model comparison statistics and goodness-of-fit statistics are identical, so only one set of statistics is presented.
}
  \item{print.status}{
A logical value indicating whether bootstrap progress updates should be provided.
}
  \item{\dots}{
Additional arguments passed to or from other methods.
}
}
\details{
The Rao-Scott approach applies a second-order adjustment to the model comparison statistic and its sampling distribution.  Formulas are provided in Appendix A of Bilder and Loughin (2007).

The bootstrap approach empirically estimates the sampling distribution of the model comparison statistic.  Gange's (1995) method for generating correlated binary data is used for taking resamples under the null hypothesis.  Bootstrap results are available only when \code{boot = TRUE} in the call to the \code{\link{genloglin}} function.
}
\value{
--- A list containing at least the following objects: \code{original.arg} and \code{test.statistics}.

\code{original.arg} is a list containing the following objects:
\itemize{
\item{\code{model}: The original model specified in the call to the \code{\link{genloglin}} function.}
\item{\code{model.HA}: The alternative model specified for the \code{model.HA} argument.}
\item{\code{gof}: The original value supplied to the \code{gof} argument.}
}

\code{test.statistics} is a list containing at least the following objects:
\itemize{
\item{\code{Pearson.chisq}: The Pearson model comparison statistic calculated using the observed data.}
\item{\code{lrt}: The LRT model comparison statistic calculated using the observed data.}
}
If \code{gof = TRUE}, \code{test.statistics} additionally contains
\itemize{
\item{\code{Pearson.chisq.gof}: The Pearson goodness-of-fit statistic calculated using the observed data.}
\item{\code{lrt.gof}: The LRT goodness-of-fit statistic calculated using the observed data.}
}

--- For \code{type = "boot"}, the primary list additionally includes \code{boot.results}, a list containing at least the following objects:
\itemize{
\item{\code{B.use}: The number of bootstrap resamples used.}
\item{\code{B.discard}: The number of bootstrap resamples discarded due to having at least one item with all positive or negative responses.}
\item{\code{p.chisq.boot}: The bootstrap p-value for the Pearson model comparison test.}
\item{\code{p.lrt.boot}: The bootstrap p-value for the LRT model comparison test.}
}
If \code{gof = TRUE}, \code{boot.results} additionally contains
\itemize{
\item{\code{p.chisq.gof.boot}: The bootstrap p-value for the Pearson goodness-of-fit test.}
\item{\code{p.lrt.gof.boot}: The bootstrap p-value for the LRT goodness-of-fit test.}
}

--- For \code{type = "rs2"}, the primary list additionally includes \code{rs.results}, a list that includes at least \code{Pearson.chisq.rs} and \code{lrt.rs}.

\code{Pearson.chisq.rs} is a list containing the following objects:
\itemize{
\item{\code{Pearson.chisq.rs}: The Rao-Scott second-order adjusted Pearson model comparison statistic.}
\item{\code{df}: The Rao-Scott second-order adjusted degrees of freedom for the model comparison statistic.}
\item{\code{p.value}: The p-value for the Rao-Scott second-order adjusted Pearson model comparison test.}
}

\code{lrt.rs} is a list containing the following objects:
\itemize{
\item{\code{lrt.rs}: The Rao-Scott second-order adjusted LRT model comparison statistic.}
\item{\code{df}: Same as \code{df} given above.}
\item{\code{p.value}: The p-value for the Rao-Scott second-order adjusted LRT model comparison test.}
}

If \code{gof = TRUE}, \code{rs.results} additionally includes \code{Pearson.chisq.gof.rs} and \code{lrt.gof.rs}.

\code{Pearson.chisq.gof.rs} is a list containing the following objects:
\itemize{
\item{\code{Pearson.chisq.gof.rs}: The Rao-Scott second-order adjusted Pearson goodness-of-fit statistic.}
\item{\code{df}: Same as \code{df} given above.}
\item{\code{p.value}: The p-value for the Rao-Scott second-order adjusted Pearson goodness-of-fit test.}
}

\code{lrt.gof.rs} is a list containing the following objects:
\itemize{
\item{\code{lrt.gof.rs}: The Rao-Scott second-order adjusted LRT goodness-of-fit statistic.}
\item{\code{df}: Same as \code{df} given above.}
\item{\code{p.value}: The p-value for the Rao-Scott second-order adjusted LRT goodness-of-fit test.}
}

--- For \code{type = "all"}, the original list includes the \code{boot.results} and \code{rs.results} output.
}
\references{
Bilder, C. and Loughin, T. (2007)  Modeling association between two or more categorical variables that allow for multiple category choices.  \emph{Communications in Statistics--Theory and Methods}, \bold{36}, 433--451.

Gange, S. (1995)  Generating multivariate categorical variates using the iterative proportional fitting algorithm.  \emph{The American Statistician}, \bold{49}, 134--138.
}
\examples{
## For examples see help(genloglin).
}