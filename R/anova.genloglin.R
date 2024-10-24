#' @name    anova.genloglin
#' @aliases anova.genloglin
#' @title   Perform MRCV Model Comparison Tests
#' @description The \code{anova.genloglin} method function offers second-order Rao-Scott and
#' bootstrap adjusted model comparison and goodness-of-fit (Pearson and LRT)
#' statistics appropriate for evaluating models estimated by the
#' \code{\link{genloglin}} function.
#'
#' @details
#' The Rao-Scott approach applies a second-order adjustment to the model
#' comparison statistic and its sampling distribution.  Formulas are provided
#' in Appendix A of Bilder and Loughin (2007).
#'
#' The bootstrap approach empirically estimates the sampling distribution of
#' the model comparison statistic.  Gange's (1995) method for generating
#' correlated binary data is used for taking resamples under the null
#' hypothesis.  Bootstrap results are available only when \code{boot = TRUE} in
#' the call to the \code{\link{genloglin}} function.
#'
#' @param object An object of class \code{'genloglin'} produced by the
#' \code{\link{genloglin}} function.
#' @param model.HA For the two MRCV case, a character string specifying one of
#' the following models to be compared to the null model (where the null model
#' should be nested within the alternative model): \code{"homogeneous"} (the
#' homogeneous association model), \code{"w.main"} (the w-main effects model),
#' \code{"y.main"} (the y-main effects model), \code{"wy.main"} (the w- and
#' y-main effects model), or \code{"saturated"}.  Alternatively, a
#' user-supplied formula can be specified.  For the three MRCV case, only
#' \code{"saturated"} or user-supplied formulas are accepted.
#' @param type A character string specifying one of the following approaches
#' for performing adjusted model comparison tests: \code{"boot"} specifies a
#' bootstrapping procedure; \code{"rs2"} specifies a Rao-Scott second-order
#' adjustment; \code{"all"} specifies both approaches.
#' @param gof A logical value indicating whether goodness-of-fit statistics
#' should be calculated in addition to model comparison statistics.  For
#' \code{model.HA = "saturated"}, model comparison statistics and
#' goodness-of-fit statistics are identical, so only one set of statistics is
#' presented.
#' @param print.status A logical value indicating whether bootstrap progress
#' updates should be provided.
#' @param \dots Additional arguments passed to or from other methods.
#'
#' @return --- A list containing at least the following objects:
#' \code{original.arg} and \code{test.statistics}.
#'
#' \code{original.arg} is a list containing the following objects: \itemize{
#' \item\code{model}: The original model specified in the call to the
#' \code{\link{genloglin}} function. \item\code{model.HA}: The alternative
#' model specified for the \code{model.HA} argument. \item\code{gof}: The
#' original value supplied to the \code{gof} argument. }
#'
#' \code{test.statistics} is a list containing at least the following objects:
#' \itemize{ \item\code{Pearson.chisq}: The Pearson model comparison statistic
#' calculated using the observed data. \item\code{lrt}: The LRT model
#' comparison statistic calculated using the observed data. } If \code{gof =
#' TRUE}, \code{test.statistics} additionally contains \itemize{
#' \item\code{Pearson.chisq.gof}: The Pearson goodness-of-fit statistic
#' calculated using the observed data. \item\code{lrt.gof}: The LRT
#' goodness-of-fit statistic calculated using the observed data. }
#'
#' --- For \code{type = "boot"}, the primary list additionally includes
#' \code{boot.results}, a list containing at least the following objects:
#' \itemize{ \item\code{B.use}: The number of bootstrap resamples used.
#' \item\code{B.discard}: The number of bootstrap resamples discarded due to
#' having at least one item with all positive or negative responses.
#' \item\code{p.chisq.boot}: The bootstrap p-value for the Pearson model
#' comparison test. \item\code{p.lrt.boot}: The bootstrap p-value for the LRT
#' model comparison test. } If \code{gof = TRUE}, \code{boot.results}
#' additionally contains \itemize{ \item\code{p.chisq.gof.boot}: The bootstrap
#' p-value for the Pearson goodness-of-fit test. \item\code{p.lrt.gof.boot}:
#' The bootstrap p-value for the LRT goodness-of-fit test. }
#'
#' --- For \code{type = "rs2"}, the primary list additionally includes
#' \code{rs.results}, a list that includes at least \code{Pearson.chisq.rs} and
#' \code{lrt.rs}.
#'
#' \code{Pearson.chisq.rs} is a list containing the following objects:
#' \itemize{ \item\code{Pearson.chisq.rs}: The Rao-Scott second-order adjusted
#' Pearson model comparison statistic. \item\code{df}: The Rao-Scott
#' second-order adjusted degrees of freedom for the model comparison statistic.
#' \item\code{p.value}: The p-value for the Rao-Scott second-order adjusted
#' Pearson model comparison test. }
#'
#' \code{lrt.rs} is a list containing the following objects: \itemize{
#' \item\code{lrt.rs}: The Rao-Scott second-order adjusted LRT model comparison
#' statistic. \item\code{df}: Same as \code{df} given above.
#' \item\code{p.value}: The p-value for the Rao-Scott second-order adjusted LRT
#' model comparison test. }
#'
#' If \code{gof = TRUE}, \code{rs.results} additionally includes
#' \code{Pearson.chisq.gof.rs} and \code{lrt.gof.rs}.
#'
#' \code{Pearson.chisq.gof.rs} is a list containing the following objects:
#' \itemize{ \item\code{Pearson.chisq.gof.rs}: The Rao-Scott second-order
#' adjusted Pearson goodness-of-fit statistic. \item\code{df}: Same as
#' \code{df} given above. \item\code{p.value}: The p-value for the Rao-Scott
#' second-order adjusted Pearson goodness-of-fit test. }
#'
#' \code{lrt.gof.rs} is a list containing the following objects: \itemize{
#' \item\code{lrt.gof.rs}: The Rao-Scott second-order adjusted LRT
#' goodness-of-fit statistic. \item\code{df}: Same as \code{df} given above.
#' \item\code{p.value}: The p-value for the Rao-Scott second-order adjusted LRT
#' goodness-of-fit test. }
#'
#' --- For \code{type = "all"}, the original list includes the
#' \code{boot.results} and \code{rs.results} output.
#' @references Bilder, C. and Loughin, T. (2007) Modeling association between
#' two or more categorical variables that allow for multiple category choices.
#' \emph{Communications in Statistics--Theory and Methods}, \bold{36},
#' 433--451.
#'
#' Gange, S. (1995) Generating multivariate categorical variates using the
#' iterative proportional fitting algorithm.  \emph{The American Statistician},
#' \bold{49}, 134--138.
#'
#' @examples
#'
#' ## For examples see help(genloglin).
#'
#' @export
#' @rdname anova
anova.genloglin <- function(object, model.HA = "saturated", type = "all", gof = TRUE,
                            print.status = TRUE, ...) {
  op <- options()
  on.exit(options(op))
  options(warn = 1)
  if (!is.logical(print.status) & (print.status != 1) & (print.status != 0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status <- TRUE
  }
  if (!is.logical(gof) & (gof != 1) & (gof != 0)) {
    warning("The \"gof\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    gof <- TRUE
  }



  print.status <- as.logical(print.status)
  gof <- as.logical(gof)
  if (length(type) > 1) {
    warning("The \"type\" argument requires an object of length 1. \n  Only the first element is used.")
    type <- type[1]
  }
  if (!(type %in% c("all", "boot", "rs2"))) {
    warning("The \"type\" argument can only take on values of \"boot\", \"rs2\", and \"all\". \n  The input value has been changed to the default value of \"all\".")
    type <- "all"
  }
  boot <- object$original.arg$boot
  if (!boot) {
    if (type != "rs2") {
      warning("You must specify the boot option in genloglin() in order to obtain bootstrap results. \n  The \"type\" argument has been changed to \"rs2\".")
      type <- "rs2"
    }
  }
  data <- object$original.arg$data
  nvars <- object$original.arg$nvars
  I <- object$original.arg$I
  J <- object$original.arg$J
  K <- object$original.arg$K
  if (!is.null(K)) {
    if (!inherits(model.HA, "formula") & (model.HA != "saturated")) {
      warning("For the 3 MRCV case, only \"saturated\" or user-supplied formulas are accepted by the \"model.HA\" argument. \n  The input value has been changed to the default value of \"saturated\".")
      model.HA <- "saturated"
    }
  }
  if (!inherits(model.HA, "formula") & (model.HA != "homogeneous") & (model.HA != "w.main") & (model.HA != "y.main") & (model.HA != "wy.main") & (model.HA != "saturated")) {
    warning("The \"model.HA\" argument requires a formula or one of 5 recognized character strings. \n  See help(anova.genloglin) for details. \n  The input value has been changed to the default value of \"saturated\".")
    model.HA <- "saturated"
  }
  model <- object$original.arg$model
  model.data <- object$mod.fit$data[, 1:(2 * nvars + 1)]
  mu.hat <- object$mod.fit$fitted.values
  mu.hat.HA <- model.data[, (2 * nvars + 1)]
  deviance <- object$mod.fit$deviance
  gamma <- object$rs.results$gamma
  MRCV_globals$print.status <- print.status
  # Compute observed test statistics
  chisq.obs <- sum(((mu.hat.HA - mu.hat)^2) / mu.hat)
  lrt.obs <- deviance
  chisq.gof.obs <- chisq.obs
  lrt.gof.obs <- lrt.obs
  # If alternative model is not saturated model then need to estimate model
  if (model.HA != "saturated") {
    # if (class(model.HA)=="formula") {
    if (inherits(model.HA, "formula")) {
      for (i in 1:I) {
        parm <- paste("W", i, sep = "")
        if (length(agrep(parm, model.HA, max.distance = 0)) > 0) {
          model.data <- data.frame(
            model.data,
            as.numeric((model.data[, 1] == names(data)[i]))
          )
          colnames(model.data)[ncol(model.data)] <- parm
        }
      }
      for (j in 1:J) {
        parm <- paste("Y", j, sep = "")
        if (length(agrep(parm, model.HA, max.distance = 0)) > 0) {
          model.data <- data.frame(
            model.data,
            as.numeric((model.data[, 2] == names(data)[(I + j)]))
          )
          colnames(model.data)[ncol(model.data)] <- parm
        }
      }
      if (nvars == 3) {
        for (k in 1:K) {
          parm <- paste("Z", k, sep = "")
          if (length(agrep(parm, model.HA, max.distance = 0)) > 0) {
            model.data <- data.frame(
              model.data,
              as.numeric((model.data[, 3] == names(data)[(I + J + k)]))
            )
            colnames(model.data)[ncol(model.data)] <- parm
          }
        }
      }
    }
    mod.fit.HA <- genloglin.fit(data = model.data, model = model.HA, nvars = nvars)
    model.data <- model.data[, 1:(2 * nvars + 1)]
    mu.hat.HA <- mod.fit.HA$fitted.values
    deviance.HA <- mod.fit.HA$deviance
    chisq.obs <- sum(((mu.hat.HA - mu.hat)^2) / mu.hat)
    lrt.obs <- deviance - deviance.HA
  }
  original.arg <- list(model = model, model.HA = model.HA, gof = gof)
  test.statistics <- list(Pearson.chisq = chisq.obs, lrt = lrt.obs)
  # If gof=TRUE and alternative model is not saturated then provide gof statistics
  if (gof & (model.HA != "saturated")) {
    test.statistics <- list(
      Pearson.chisq = chisq.obs, lrt = lrt.obs,
      Pearson.chisq.gof = chisq.gof.obs,
      lrt.gof = lrt.gof.obs
    )
  }

  if (any(type == "rs2" | type == "all")) {
    sum.gamma <- sum(gamma)
    sum.gamma.sq <- sum(gamma^2)
    chisq.rs <- (sum.gamma * chisq.obs) / sum.gamma.sq
    lrt.rs <- (sum.gamma * lrt.obs) / sum.gamma.sq
    chisq.gof.rs <- (sum.gamma * chisq.gof.obs) / sum.gamma.sq
    lrt.gof.rs <- (sum.gamma * lrt.gof.obs) / sum.gamma.sq
    df.rs <- (sum.gamma^2) / sum.gamma.sq
    p.chisq.rs <- 1 - pchisq(q = chisq.rs, df = df.rs)
    p.lrt.rs <- 1 - pchisq(q = lrt.rs, df = df.rs)
    p.chisq.gof.rs <- 1 - pchisq(q = chisq.gof.rs, df = df.rs)
    p.lrt.gof.rs <- 1 - pchisq(q = lrt.gof.rs, df = df.rs)
    rs.results <- list(
      Pearson.chisq.rs = list(
        Pearson.chisq.rs = chisq.rs,
        df = df.rs, p.value = p.chisq.rs
      ), lrt.rs =
        list(lrt.rs = lrt.rs, df = df.rs, p.value = p.lrt.rs)
    )
    if (gof & (model.HA != "saturated")) {
      rs.results <- list(
        Pearson.chisq.rs = list(
          Pearson.chisq.rs = chisq.rs,
          df = df.rs, p.value = p.chisq.rs
        ), lrt.rs = list(
          lrt.rs =
            lrt.rs, df = df.rs, p.value = p.lrt.rs
        ),
        Pearson.chisq.gof.rs = list(
          Pearson.chisq.gof.rs =
            chisq.gof.rs, df = df.rs, p.value = p.chisq.gof.rs
        ),
        lrt.gof.rs = list(
          lrt.gof.rs = lrt.gof.rs, df = df.rs,
          p.value = p.lrt.gof.rs
        )
      )
    }
    output <- list(
      original.arg = original.arg, test.statistics = test.statistics,
      rs.results = rs.results
    )
  }

  if (any(type == "boot" | type == "all")) {
    nrows <- (2^nvars) * I * J * max(1, K)
    B.use <- object$boot.results$B.use
    B.discard <- object$boot.results$B.discard
    model.data.star <- object$boot.results$model.data.star
    mu.hat.star <- object$boot.results$mod.fit.star[, 1:nrows]
    deviance.star <- object$boot.results$mod.fit.star[, (nrows + 1)]
    chisq.star <- object$boot.results$chisq.star
    lrt.star <- object$boot.results$lrt.star
    chisq.gof.star <- chisq.star
    lrt.gof.star <- lrt.star
    # If alternative model is not saturated then need to estimate t*'s for HA
    if (model.HA != "saturated") {
      # if (class(model.HA)=="formula") {
      if (inherits(model.HA, "formula")) {
        for (i in 1:I) {
          parm <- paste("W", i, sep = "")
          if (length(agrep(parm, model.HA, max.distance = 0)) > 0) {
            model.data.star <- data.frame(
              model.data.star,
              as.numeric((model.data.star[, 1] == names(data)[i]))
            )
            colnames(model.data.star)[ncol(model.data.star)] <- parm
          }
        }
        for (j in 1:J) {
          parm <- paste("Y", j, sep = "")
          if (length(agrep(parm, model.HA, max.distance = 0)) > 0) {
            model.data.star <- data.frame(
              model.data.star,
              as.numeric((model.data.star[, 2] == names(data)[(I + j)]))
            )
            colnames(model.data.star)[ncol(model.data.star)] <- parm
          }
        }
        if (nvars == 3) {
          for (k in 1:K) {
            parm <- paste("Z", k, sep = "")
            if (length(agrep(parm, model.HA, max.distance = 0)) > 0) {
              model.data.star <- data.frame(
                model.data.star,
                as.numeric((model.data.star[, 3] == names(data)[(I + J + k)]))
              )
              colnames(model.data.star)[ncol(model.data.star)] <- parm
            }
          }
        }
      }
      var.cols <- c(1:(2 * nvars))
      if (ncol(model.data.star) > (B.use + (2 * nvars))) {
        var.cols <- c(1:(2 * nvars), (B.use + (2 * nvars + 1)):ncol(model.data.star))
      }
      # Create progress bar for bootstrapping
      if (print.status) {
        cat("Bootstrap Progress:", "\n")
        MRCV_globals$pb <- txtProgressBar(min = 0, max = B.use, style = 3)
        MRCV_globals$pb.counter <- 0
      }
      mod.fit.star.HA <- matrix(NA, nrow = B.use, ncol = (nrows + 1))
      mod.fit.star.HA <- t(apply(
        X =
          as.matrix(model.data.star[, (2 * nvars + 1):(B.use + (2 * nvars))]),
        MARGIN = 2, FUN = genloglin.fit, model = model.HA,
        nvars = nvars, limit.output = TRUE, model.vars =
          model.data.star[, var.cols]
      ))
      if (print.status) {
        setTxtProgressBar(MRCV_globals$pb, B.use)
        close(MRCV_globals$pb)
      }
      model.data.star <- model.data.star[, 1:(B.use + 2 * nvars)]
      mu.hat.star.HA <- mod.fit.star.HA[, (1:nrows)]
      deviance.star.HA <- mod.fit.star.HA[, ncol(mod.fit.star.HA)]
      chisq.star <- rowSums(((mu.hat.star.HA - mu.hat.star)^2) / mu.hat.star)
      lrt.star <- deviance.star - deviance.star.HA
    }
    p.chisq.boot <- (1 / B.use) * sum(chisq.star >= chisq.obs)
    p.lrt.boot <- (1 / B.use) * sum(lrt.star >= lrt.obs)
    p.chisq.gof.boot <- (1 / B.use) * sum(chisq.gof.star >= chisq.gof.obs)
    p.lrt.gof.boot <- (1 / B.use) * sum(lrt.gof.star >= lrt.gof.obs)
    boot.results <- list(
      B.use = B.use, B.discard = B.discard, p.chisq.boot =
        p.chisq.boot, p.lrt.boot = p.lrt.boot
    )
    if (gof & (model.HA != "saturated")) {
      boot.results <- list(
        B.use = B.use, B.discard = B.discard, p.chisq.boot =
          p.chisq.boot, p.lrt.boot = p.lrt.boot, p.chisq.gof.boot =
          p.chisq.gof.boot, p.lrt.gof.boot = p.lrt.gof.boot
      )
    }
    output <- list(
      original.arg = original.arg, test.statistics = test.statistics,
      boot.results = boot.results
    )
  }
  if (type == "all") {
    output <- list(
      original.arg = original.arg, test.statistics = test.statistics,
      rs.results = rs.results, boot.results = boot.results
    )
  }
  if (model == "saturated") {
    test.statistics <- list(Pearson.chisq = 0, lrt = 0)
    rs.results <- list(
      Pearson.chisq.rs = list(p.value = 1),
      lrt.rs = list(p.value = 1)
    )
    boot.results <- list(
      Pearson.chisq.boot = list(p.value = 1), lrt.boot =
        list(p.value = 1)
    )
    output <- list(
      original.arg = original.arg, test.statistics = test.statistics,
      rs.results = rs.results, boot.results = boot.results
    )
  }
  class(output) <- "anova.genloglin"
  output
}
