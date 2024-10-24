#' @name residuals.genloglin
#' @aliases residuals.genloglin
#' @title Calculate Standardized Pearson Residuals for MRCV Data
#' @description
#' A method function that calculates standardized residuals
#' @details
#' The \code{residuals.genloglin} method function calculates standardized
#' Pearson residuals for the model specified in the \code{\link{genloglin}}
#' function. It offers an asymptotic approximation and a bootstrap
#' approximation for estimating the variance of the residuals.
#'
#' The bootstrap results are only available when \code{boot = TRUE} in the call
#' to the \code{\link{genloglin}} function.
#'
#' The \code{residuals.genloglin} function uses
#' \code{\link[tables:tabular]{tables:tabular()}} to display the results for
#' the two MRCV case.
#'
#' See Bilder and Loughin (2007) for additional details about calculating the
#' residuals.
#'
#' @param object An object of class \code{'genloglin'} produced by the
#' \code{\link{genloglin}} function.
#' @param \dots Additional arguments passed to or from other methods.
#' @return --- A list containing at least \code{std.pearson.res.asymp.var}.
#' For the two MRCV case, the object is a 2Ix2J table of class \code{'tabular'}
#' containing the standardized Pearson residuals based on the estimated
#' asymptotic variance.  For the three MRCV case, the object is a data frame
#' containing the 2Ix2Jx2K residuals.
#'
#' --- For \code{boot = TRUE} in the call to the \code{\link{genloglin}}
#' function, the list additionally includes: \itemize{ \item\code{B.use}: The
#' number of bootstrap resamples used. \item\code{B.discard}: The number of
#' bootstrap resamples discarded due to having at least one item with all
#' positive or negative responses. \item\code{std.pearson.res.boot.var}: For
#' the two MRCV case, a 2Ix2J table of class \code{'tabular'} containing the
#' standardized Pearson residuals based on the bootstrap variance.  For the
#' three MRCV case, a data frame containing the 2Ix2Jx2K residuals. }
#' @references Bilder, C. and Loughin, T. (2007) Modeling association between
#' two or more categorical variables that allow for multiple category choices.
#' \emph{Communications in Statistics--Theory and Methods}, \bold{36},
#' 433--451.
#' @examples
#'
#' ## For examples see help(genloglin).
#'
#' @export
#' @rdname residuals
residuals.genloglin <- function(object, ...) {
  data <- object$original.arg$data
  boot <- object$original.arg$boot
  model <- object$original.arg$model
  nvars <- object$original.arg$nvars
  I <- object$original.arg$I
  J <- object$original.arg$J
  K <- object$original.arg$K
  model.data <- object$mod.fit$data[, 1:(2 * nvars + 1)]
  mu.hat <- object$mod.fit$fitted.values
  mu.hat.sat <- model.data[, (2 * nvars + 1)]
  E <- object$rs.results$E

  # Calculate standardized pearson residuals using the estimated asymptotic variance
  resid.num <- mu.hat.sat - mu.hat
  for (i in 1:length(resid.num)) {
    if (abs(resid.num[i]) < .000000001) {
      resid.num[i] <- 0
    }
  }
  std.pearson.res.asymp.var <- resid.num / sqrt(diag(E))
  std.pearson.res.asymp.var <- data.frame(model.data[, 1:(2 * nvars)],
    res = round(std.pearson.res.asymp.var, 2)
  )
  if (nvars == 2) {
    res.table.asymp <- tabular(Heading() * res * Heading() * (mean) * W * Heading()
      * Factor(wi, wi, c(1, 0)) ~ Heading() * Y * Heading()
      * Factor(yj, yj, c(1, 0)), data = std.pearson.res.asymp.var)
    output <- list(std.pearson.res.asymp.var = res.table.asymp)
  }
  if (nvars == 3) {
    output <- list(std.pearson.res.asymp.var = std.pearson.res.asymp.var)
  }

  if (boot) {
    B.use <- object$boot.results$B.use
    B.discard <- object$boot.results$B.discard
    residual.star <- object$boot.results$residual.star
    resid.denom <- apply(X = residual.star, MARGIN = 2, FUN = sd)
    # Calculate standardized pearson residuals using the bootstrap variance
    std.pearson.res.boot.var <- resid.num / resid.denom
    std.pearson.res.boot.var <- data.frame(model.data[, 1:(2 * nvars)],
      res = round(std.pearson.res.boot.var, 2)
    )
    if (nvars == 2) {
      res.table.boot <- tabular(Heading() * res * Heading() * (mean) * W * Heading()
        * Factor(wi, wi, c(1, 0)) ~ Heading() * Y * Heading()
        * Factor(yj, yj, c(1, 0)), data = std.pearson.res.boot.var)
      output <- list(
        std.pearson.res.asymp.var = res.table.asymp, B.use = B.use,
        B.discard = B.discard, std.pearson.res.boot.var = res.table.boot
      )
    }
    if (nvars == 3) {
      output <- list(
        std.pearson.res.asymp.var = std.pearson.res.asymp.var,
        B.use = B.use, B.discard = B.discard, std.pearson.res.boot.var =
          std.pearson.res.boot.var
      )
    }
  }
  output
}
