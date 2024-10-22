\name{MI.test}
\alias{MI.test}
\alias{MI.stat}
\title{
Test for Marginal Independence
}
\description{
The \code{MI.test} function offers three approaches for testing multiple marginal independence (MMI) between one SRCV and one MRCV, or simultaneous pairwise marginal independence (SPMI) between two MRCVs.
}
\usage{
MI.test(data, I, J, type = "all", B = 1999, B.max = B, summary.data = 
    FALSE, add.constant = 0.5, plot.hist = FALSE, print.status = TRUE)
    
MI.stat(data, I, J, summary.data = FALSE, add.constant = 0.5)
}
\arguments{
  \item{data}{
For \code{summary.data = FALSE}: a data frame containing the raw data where rows correspond to the individual item response vectors, and columns correspond to the items, W1, \ldots, WI and Y1, \ldots, YJ (in this order).

For \code{summary.data = TRUE}: a data frame containing 4 columns generically named \code{W}, \code{Y}, \code{yj}, and \code{count} (one MRCV case), or 5 columns named \code{W}, \code{Y}, \code{wi}, \code{yj}, and \code{count} (two MRCV case).
}
  \item{I}{
The number of items corresponding to row variable W.  I = 1 for the one MRCV case.
}
  \item{J}{
The number of items corresponding to column variable Y.
}
  \item{type}{
A character string specifying one of the following approaches for testing for MI: \code{"boot"} specifies a nonparametric bootstrap procedure; \code{"rs2"} specifies a Rao-Scott second-order adjustment; \code{"bon"} specifies a Bonferroni adjustment; \code{"all"} specifies all three approaches.
}
  \item{B}{
The desired number of bootstrap resamples.
}
  \item{B.max}{
The maximum number of bootstrap resamples.  A resample is thrown out if at least one of the J (one MRCV case) or IxJ (two MRCV case) contingency tables does not have the correct dimension; \code{MI.test} uses the first \code{B} valid resamples or all valid resamples if that number is less than \code{B}.
}
  \item{summary.data}{
A logical value indicating whether \code{data} is a summary file containing the item response data instead of the raw data.  Only \code{type = "bon"} is available for \code{summary.data = TRUE}.
}
  \item{add.constant}{
A positive constant to be added to all zero marginal cell counts.
}
  \item{plot.hist}{
A logical value indicating whether plots of the emprical bootstrap sampling distributions should be provided.
}
  \item{print.status}{
A logical value indicating whether bootstrap progress updates should be provided.
}
}
\details{
The \code{MI.test} function calls \code{MI.stat} to calculate a modified Pearson statistic (see Bilder, Loughin, and Nettleton (2000) and Bilder and Loughin (2004)), and then performs the testing of MMI or SPMI.  Three sets of testing methods are implemented:
\itemize{
\item{The nonparametric bootstrap resamples under the null hypothesis by independently sampling the W and Y vectors with replacement from the observed data.  Fixed row counts (i.e., fixed counts for each level of the SRCV) are maintained for the one MRCV case.  A modified Pearson statistic is calculated for each resample.  In addition, bootstrap p-value combination methods are available to take advantage of the modified Pearson statistic's decomposition into J (one MRCV case) or IxJ (two MRCV case) individual Pearson statistics.  The minimum or the product of p-values is the combination for each resample.} {}
\item{The Rao-Scott approach applies a second-order adjustment to the modified Pearson statistic and its sampling distribution.  Formulas are provided in Appendix A of Bilder, Loughin, and Nettleton (2000) and Bilder and Loughin (2004).  Note that this test can be conservative at times.} {}
\item{The Bonferroni adjustment multiplies each p-value (using a standard chi-square approximation) from the individual Pearson statistics by J (one MRCV case) or IxJ (two MRCV case).  If a resulting p-value is greater than 1 after the multiplication, it is set to a value of 1.  The overall p-value for the test then is the minimum of these adjusted p-values.  Note that the Bonferroni adjustment tends to produce an overly conservative test when the number of individual Pearson statistics is large.} {}
}

Agresti and Liu (1999) discuss a marginal logit model approach that uses generalized estimation equations (GEE) to test for MMI.  As shown in the example given below, this approach can be performed via functions available from the \bold{geepack} package.  However, Bilder, Loughin, and Nettleton (2000) caution that the Wald test produced by this approach does not hold the correct size, particularly when the sample size is not large and marginal probabilities are small.
}
\value{
--- \code{MI.test} returns a list containing at least \code{general}, a list containing the following objects:
\itemize{
\item{\code{data}: The original data frame supplied to the \code{data} argument.}
\item{\code{I}: The original value supplied to the \code{I} argument.}
\item{\code{J}: The original value supplied to the \code{J} argument.}
\item{\code{summary.data}: The original value supplied to the \code{summary.data} argument.}
\item{\code{X.sq.S}: The modified Pearson statistic; NA if at least one of the J (one MRCV case) or IxJ (two MRCV case) contingency tables does not have the correct dimension.}
\item{\code{X.sq.S.ij}: A matrix containing the individual Pearson statistics.}
}

--- For \code{type = "boot"}, the primary list additionally includes \code{boot}, a list containing the following objects:
\itemize{
\item{\code{B.use}: The number of bootstrap resamples used.}
\item{\code{B.discard}: The number of bootstrap resamples discarded due to having at least one contingency table with incorrect dimension.}
\item{\code{p.value.boot}: The bootstrap p-value for the test of MMI or SPMI.}
\item{\code{p.combo.min.boot}: The bootstrap p-value for the minimum p-value combination method.}
\item{\code{p.combo.prod.boot}: The bootstrap p-value for the product p-value combination method.}
\item{\code{X.sq.S.star}: A numeric vector containing the modified Pearson statistics calculated for each resample.}
\item{\code{X.sq.S.ij.star}: A matrix containing the individual Pearson statistics calculated for each resample.}
\item{\code{p.combo.min.star}: A numeric vector containing the minimum p-value calculated for each resample.}
\item{\code{p.combo.prod.star}: A numeric vector containing the product p-value calculated for each resample.}
}

--- For \code{type = "rs2"}, the primary list additionally includes \code{rs2}, a list containing the following objects:
\itemize{
\item{\code{X.sq.S.rs2}: The Rao-Scott second-order adjusted Pearson statistic.}
\item{\code{df.rs2}: The degrees of freedom for testing the second-order Rao-Scott adjusted Pearson statistic.}
\item{\code{p.value.rs2}: The p-value based on the Rao-Scott second-order adjustment.}
}

--- For \code{type = "bon"}, the primary list additionally includes \code{bon}, a list containing the following objects:
\itemize{
\item{\code{p.value.bon}: The Bonferroni adjusted p-value for the test of MMI or SPMI.}
\item{\code{X.sq.S.ij.p.bon}: A matrix containing the Bonferroni adjusted p-values for the individual Pearson statistics.}
}

--- For \code{type = "all"}, the list includes all of the above objects.

--- \code{MI.stat} returns a list containing the following objects:
\itemize{
\item{\code{X.sq.S}: Defined above.}
\item{\code{X.sq.S.ij}: Defined above.}
\item{\code{valid.margins}: The number of contingency tables with correct dimension.}
}
}
\references{
Agresti, A. and Liu, I.-M. (1999)  Modeling a categorical variable allowing arbitrarily many category choices.  \emph{Biometrics}, \bold{55}, 936--943.

Bilder, C. and Loughin, T. (2004)  Testing for marginal independence between two categorical variables with multiple responses.  \emph{Biometrics}, \bold{36}, 433--451.

Bilder, C., Loughin, T., and Nettleton, D. (2000)  Multiple marginal independence testing for pick any/c variables.  \emph{Communications in Statistics--Theory and Methods}, \bold{29}, 1285--1316.
}
\seealso{
The \code{\link{genloglin}} function offers a generalized loglinear modeling approach for testing the relationship among two or three MRCVs.
}
\examples{
# Test for MMI using the second-order Rao-Scott adjustment
test.mmi.rs2 <- MI.test(data = farmer1, I = 1, J = 5, type = "rs2")
test.mmi.rs2

# Test for MMI using all three approaches
# A small B is used for demonstration purposes; normally, a larger B should be used
test.mmi.all <- MI.test(data = farmer1, I = 1, J = 5, type = "all", B = 99,
    plot.hist = TRUE)
test.mmi.all

# Use MI.test() with summary data
# Convert raw data file to summary file for this example 
farmer1.irdframe <- item.response.table(data = farmer1, I = 1, J = 5, create.dataframe = 
    TRUE)
# Test for MMI using the Bonferroni adjustment
test.mmi.bon <- MI.test(data = farmer1.irdframe, I = 1, J = 5, type = "bon", 
    summary.data = TRUE)
test.mmi.bon

# Test for SPMI using the second-order Rao-Scott adjustment
test.spmi.rs2 <- MI.test(data = farmer2, I = 3, J = 4, type = "rs2")
test.spmi.rs2

# Test for MMI using the marginal logit model approach
library(geepack)
n<-nrow(farmer1)
farmer1.id<-cbind(case=1:n, farmer1)
# Reshape raw data into long format as required by geeglm() function
# Assumes 3:ncol(farmer1.id) corresponds to MRCV items
farmer1.gee<-reshape(data = farmer1.id, 
                 varying = names(farmer1.id)[3:ncol(farmer1.id)], 
                 v.names = "response", timevar = "item", idvar = "case", 
                 direction = "long") 
row.names(farmer1.gee)<-NULL
farmer1.gee[,2:3]<-lapply(farmer1.gee[,2:3], factor)
# Data frame must be ordered by case
farmer1.gee<-farmer1.gee[order(farmer1.gee$case),]
head(farmer1.gee)
tail(farmer1.gee)
mod.fit.H0<-geeglm(formula = response ~ item, family = binomial(link = logit), 
                  data = farmer1.gee, na.action = na.omit, id = case, 
                  corstr = "unstructured")
mod.fit.HA<-geeglm(formula = response ~ Ed*item, family = binomial(link = logit), 
                   data = farmer1.gee, na.action = na.omit, id = case, 
                   corstr = "unstructured")
# Compute Wald test
anova(mod.fit.HA, mod.fit.H0)
}