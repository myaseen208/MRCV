\name{print.genloglin}
\alias{print.genloglin}
\alias{print.summary.genloglin}
\alias{print.anova.genloglin}
\alias{print.predict.genloglin}
\title{
Control Printed Display of MRCV Regression Modeling Objects
}
\description{
Method functions that control the printed display of MRCV regression modeling objects.
}
\usage{
\method{print}{genloglin}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{summary.genloglin}(x, digits = max(3, getOption("digits") - 3), symbolic.cor = 
    x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...)
\method{print}{anova.genloglin}(x, ...)
\method{print}{predict.genloglin}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"genloglin"} produced by the \code{\link{genloglin}} function.
}
  \item{digits}{
Minimum number of digits; see \code{\link{print.default}} for additional explanation.
}
  \item{symbolic.cor}{
A logical value indicating whether correlations should be printed in symbolic form; see the \code{\link{summary}} method for class \code{"glm"} for additional explanation.
}
  \item{signif.stars}{
A logical value indicating whether significance stars should be printed; see the \code{\link{summary}} method for class \code{"glm"} for additional explanation.
}
  \item{...}{
Additional arguments passed to or from other methods.
}
}
\details{
The \code{print.genloglin} function is based on the \code{\link{print}} method for class \code{"glm"}.

The \code{print.summary.genloglin} function is based on the \code{\link{print}} method for class \code{"summary.glm"}.
}
\value{
A print out of selected results.
}