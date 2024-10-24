#' @name MRCV-package
#' @aliases MRCV-package MRCV
#' @title Methods for Analyzing Multiple Response Categorical Variables
#' @author Natalie Koziol, Chris Bilder
#' Maintainer: Chris Bilder <bilder@@unl.edu>
#' @description
#' Provides functions for analyzing the association between one single response
#' categorical variable (SRCV) and one multiple response categorical variable
#' (MRCV), or between two or three MRCVs.  A modified Pearson chi-square
#' statistic can be used to test for marginal independence for the one or two
#' MRCV case, or a more general loglinear modeling approach can be used to
#' examine various other structures of association for the two or three MRCV
#' case.  Bootstrap- and asymptotic-based standardized residuals and
#' model-predicted odds ratios are available, in addition to other descriptive
#' information.
#'
#' \tabular{ll}{ Package: \tab MRCV\cr Version: \tab 0.4-0\cr Date: \tab
#' 2024-10-18\cr Depends: \tab R (>= 4.4.0)\cr Imports: \tab tables\cr
#' Suggests: \tab geepack\cr LazyData: \tab TRUE\cr License: \tab GPL (>= 3)\cr
#' }
#'
#' \strong{Notation:}\cr For the two or three MRCV case, define row variable,
#' W, column variable, Y, and strata variable, Z, as MRCVs with binary items
#' (i.e., categories) Wi for i = 1, \ldots, I, Yj for j = 1, \ldots, J, and Zk
#' for k = 1, \ldots, K, respectively.  Define a marginal count as the number
#' of subjects who responded (Wi = a, Yj = b, Zk = c) for a, b, and c belonging
#' to the set \{0, 1\}.  For the one MRCV case, let W be an SRCV such that I = 1
#' and 'a' corresponds to one of r levels of W, and let Y be the MRCV as
#' previously defined.
#'
#' \strong{Format of Data Frame:}\cr Many of the functions require a data frame
#' containing the raw data structured such that the \emph{n} rows correspond to
#' the individual item response vectors, and the columns correspond to the
#' items, W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
#' Some of the functions use a summary version of the raw data frame (converted
#' automatically without need for user action) formatted to have rx2J rows and
#' 4 columns generically named \code{W}, \code{Y}, \code{yj}, and \code{count}
#' (one MRCV case), 2Ix2J rows and 5 columns named \code{W}, \code{Y},
#' \code{wi}, \code{yj}, and \code{count} (two MRCV case), or 2Ix2Jx2K rows and
#' 7 columns named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj},
#' \code{zk}, and \code{count} (three MRCV case). The column named \code{count}
#' contains the marginal counts defined above.
#'
#' \strong{Descriptive Functions:}\cr Users can call the
#' \code{\link{item.response.table}} function to obtain a cross-tabulation of
#' the positive and negative responses for each combination of items, or the
#' \code{\link{marginal.table}} function to obtain a cross-tabulation of only
#' the positive responses.
#'
#' \strong{Functions to Test for Marginal Independence:}\cr Methods proposed by
#' Agresti and Liu (1999), Bilder and Loughin (2004), Bilder, Loughin, and
#' Nettleton (2000), and Thomas and Decady (2004) are implemented using the
#' \code{\link{MI.test}} function.  This function calculates a modified Pearson
#' chi-square statistic that can be used to test for multiple marginal
#' independence (MMI; one MRCV case) or simultaneous pairwise marginal
#' independence (SPMI; two MRCV case).  MMI is a test of whether the SRCV, W,
#' is marginally independent of each Yj, where the modified statistic is the
#' sum of the J Pearson statistics used to test for independence of each (W,
#' Yj) pair.  SPMI is a test of whether each Wi is pairwise independent of each
#' Yj, where the modified statistic is the sum of the IxJ Pearson statistics
#' used to test for independence of each (Wi, Yj) pair.  The asymptotic
#' distribution of the modified statistics is a linear combination of
#' independent chi-square(1) random variables, so traditional methods for
#' analyzing the association between categorical variables W and Y are
#' inappropriate.  The \code{\link{MI.test}} function offers three sets of
#' testing methods: a nonparametric bootstrap approach, a Rao-Scott
#' second-order adjustment, and a Bonferroni adjustment, that can be used in
#' conjunction with the modified statistic to construct an appropriate test for
#' independence.
#'
#' \strong{Functions for Performing Regression Modeling:}\cr Regression
#' modeling methods described by Bilder and Loughin (2007) are implemented
#' using \code{\link{genloglin}} and methods \code{\link{summary.genloglin}},
#' \code{\link{residuals.genloglin}}, \code{\link{anova.genloglin}}, and
#' \code{\link{predict.genloglin}}.  The \code{\link{genloglin}} function
#' provides parameter estimates and Rao-Scott adjusted standard errors for
#' models involving two or three MRCVs.  The \code{\link{anova.genloglin}}
#' function offers second-order Rao-Scott and bootstrap adjusted model
#' comparison and goodness-of-fit (Pearson and LRT) statistics.  The
#' \code{\link{residuals.genloglin}} and \code{\link{predict.genloglin}}
#' functions provide bootstrap- and asymptotic-based standardized Pearson
#' residuals and model-based odds ratios, respectively.
#'
#' \strong{General Notes:}\cr Rao-Scott adjustments may not be feasible when
#' the total number of MRCV items is large.  In this case, an error message
#' will be returned describing a memory allocation issue.
#'
#' @references Agresti, A. and Liu, I.-M. (1999) Modeling a categorical
#' variable allowing arbitrarily many category choices.  \emph{Biometrics},
#' \bold{55}, 936--943.
#'
#' Bilder, C. and Loughin, T. (2004) Testing for marginal independence between
#' two categorical variables with multiple responses.  \emph{Biometrics},
#' \bold{36}, 433--451.
#'
#' Bilder, C. and Loughin, T. (2007) Modeling association between two or more
#' categorical variables that allow for multiple category choices.
#' \emph{Communications in Statistics--Theory and Methods}, \bold{36},
#' 433--451.
#'
#' Bilder, C., Loughin, T., and Nettleton, D. (2000) Multiple marginal
#' independence testing for pick any/c variables.  \emph{Communications in
#' Statistics--Theory and Methods}, \bold{29}, 1285--1316.
#'
#' Koziol, N. and Bilder, C. (2014) MRCV: A package for analyzing categorical
#' variables with multiple response options.  \emph{The R Journal}, \bold{6},
#' 144--150.
#'
#' Thomas, D. and Decady, Y. (2004) Testing for association using multiple
#' response survey data: Approximate procedures based on the Rao-Scott
#' approach.  \emph{International Journal of Testing}, \bold{4}, 43--59.
#' 
#' @import tables
#' @importFrom graphics abline hist layout
#' @importFrom stats chisq.test coef glm model.matrix naprint pchisq pnorm poisson printCoefmat qnorm quantile rmultinom sd symnum
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
NULL

.onLoad <- function(libname, pkgname) {
  assign("MRCV_globals", new.env(), envir = parent.env(environment()))
}
