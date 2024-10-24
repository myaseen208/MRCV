#' @name    predict.genloglin
#' @aliases predict.genloglin
#' @title   Calculate Observed and Model-Predicted Odds Ratios for MRCV Data
#' @description The \code{predict.genloglin} method calculates observed and model-predicted odds ratios
#' and their confidence intervals using results from \code{\link{genloglin}}. It uses an
#' asymptotic normal approximation for the confidence intervals for both the observed and
#' model-predicted odds ratios, and provides a bootstrap approach for estimating the
#' confidence intervals of the model-predicted odds ratios.
#'  
#' @details
#' Wald confidence intervals are estimated for both model-based (see Appendix A of Bilder
#' and Loughin, 2007) and observed (see Agresti, 2013, p. 70) odds ratios.
#'
#' A bootstrap method is also available, providing bias-corrected accelerated (BCa)
#' confidence intervals for the model-predicted odds ratios. For more information on BCa
#' intervals, see Efron (1987). The \code{predict.genloglin} function uses a jackknife
#' approximation to estimate empirical influence values.
#'
#' Note: The bootstrap confidence intervals are only available when \code{boot = TRUE}
#' in the original call to the \code{\link{genloglin}} function.
#'
#' @param object An object of class \code{"genloglin"}, produced by the
#'   \code{\link{genloglin}} function.
#' @param alpha A numeric value specifying the desired alpha level for the confidence
#'   intervals. The function provides two-sided (1-\code{alpha})x100\% confidence
#'   intervals.
#' @param pair For the case of three MRCVs, a character string specifying the pair of
#'   items for which odds ratios will be calculated: \code{"WY"} for (Wi, Yj) conditional
#'   on Zk, \code{"WZ"} for (Wi, Zk) conditional on Yj, and \code{"YZ"} for (Yj, Zk)
#'   conditional on Wi.
#' @param print.status A logical value indicating whether to provide bootstrap progress
#'   updates.
#' @param \dots Additional arguments passed to or from other methods.
#'
#' @return A list containing at least \code{original.arg}, \code{OR.obs}, and
#'   \code{OR.model.asymp}:
#'   \itemize{
#'     \item \code{original.arg}: A list containing the original arguments supplied to
#'           \code{genloglin}.
#'     \item \code{OR.obs}: A numeric matrix containing observed odds ratios and their
#'           confidence intervals.
#'     \item \code{OR.model.asymp}: A numeric matrix containing model-predicted odds
#'           ratios and their confidence intervals using asymptotic normal approximation.
#'   }
#'
#' For \code{boot = TRUE} in the \code{\link{genloglin}} function, the returned list
#' additionally includes:
#'   \itemize{
#'     \item \code{boot.results}: A list containing bootstrap-related results such as
#'           \code{B.use}, \code{B.discard}, and \code{OR.model.BCa}.
#'   }
#'
#' @references
#' Agresti, A. (2013). \emph{Categorical Data Analysis} (3rd ed.). Hoboken, NJ: John Wiley & Sons.
#'
#' Bilder, C., & Loughin, T. (2007). Modeling association between two or more categorical variables
#' that allow for multiple category choices. \emph{Communications in Statistics--Theory and Methods},
#' \bold{36}, 433-451.
#'
#' Efron, B. (1987). Better bootstrap confidence intervals. \emph{Journal of the American Statistical
#' Association}, \bold{82}, 171-185.
#'
#' @examples
#' ## For examples see help(genloglin).
#'
#' @export
#' @rdname predict
#'
predict.genloglin <- function(object, alpha = .05, pair = "WY", print.status = TRUE, ...) {
  op <- options()
  on.exit(options(op))
  options(warn = 1)
  if (!is.numeric(alpha)) {
    warning("The \"alpha\" argument only accepts numeric values. \n  The input value has been changed to the default value of .05.")
    alpha <- .05
  }
  if ((alpha <= 0) | (alpha >= 1)) {
    warning("The \"alpha\" argument must be between 0 and 1. \n  The input value has been changed to the default value of .05.")
    alpha <- .05
  }
  if (length(pair) > 1) {
    warning("The \"pair\" argument requires an object of length 1. \n  Only the first element is used.")
    pair <- pair[1]
  }
  if (!(pair %in% c("WY", "WZ", "YZ", "YW", "ZW", "ZY", "wy", "wz", "yz", "yw", "zw", "zy"))) {
    warning("The \"pair\" argument can only take on values of \"WY\", \"WZ\", and \"YZ\". \n  The input value has been changed to the default value of \"WY\".")
    pair <- "WY"
  }

  # if ((class(print.status)!="logical")&(print.status!=1)&(print.status!=0)) {
  #  warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
  #  print.status<-TRUE
  # }
  if (!is.logical(print.status) & (print.status != 1) & (print.status != 0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status <- TRUE
  }



  print.status <- as.logical(print.status)
  boot <- object$original.arg$boot
  data <- object$original.arg$data
  n <- nrow(data)
  nvars <- object$original.arg$nvars
  I <- object$original.arg$I
  J <- object$original.arg$J
  K <- object$original.arg$K
  add.constant <- object$original.arg$add.constant
  model <- object$original.arg$model
  model.data <- object$mod.fit$data[, 1:(2 * nvars + 1)]
  mu.hat <- object$mod.fit$fitted.values
  cov.mu <- object$rs.results$cov.mu
  nrows <- I * J * max(1, K)
  MRCV_globals$print.status <- print.status
  if (nvars == 3) {
    if (any(pair == "WZ" | pair == "ZW" | pair == "wz" | pair == "zw")) {
      cov.mu <- cbind(cov.mu, model.data[, 1:6])
      cov.mu <- cov.mu[order(-cov.mu$yj, -cov.mu$wi, -cov.mu$zk, cov.mu$Y, cov.mu$W, cov.mu$Z), ]
      cov.mu <- cbind(t(cov.mu[, 1:(ncol(cov.mu) - 6)]), model.data[, 1:6])
      cov.mu <- cov.mu[order(-cov.mu$yj, -cov.mu$wi, -cov.mu$zk, cov.mu$Y, cov.mu$W, cov.mu$Z), ]
      cov.mu <- t(cov.mu[, 1:(ncol(cov.mu) - 6)])
      mu.hat <- cbind(as.data.frame(mu.hat), model.data[, 1:6])
      mu.hat <- mu.hat[order(-mu.hat$yj, -mu.hat$wi, -mu.hat$zk, mu.hat$Y, mu.hat$W, mu.hat$Z), ]
      mu.hat <- as.matrix(mu.hat[, 1])
      model.data <- model.data[order(
        -model.data$yj, -model.data$wi, -model.data$zk,
        model.data$Y, model.data$W, model.data$Z
      ), ]
    }
    if (any(pair == "WY" | pair == "YW" | pair == "wy" | pair == "yw")) {
      cov.mu <- cbind(cov.mu, model.data[, 1:6])
      cov.mu <- cov.mu[order(-cov.mu$zk, -cov.mu$wi, -cov.mu$yj, cov.mu$Z, cov.mu$W, cov.mu$Y), ]
      cov.mu <- cbind(t(cov.mu[, 1:(ncol(cov.mu) - 6)]), model.data[, 1:6])
      cov.mu <- cov.mu[order(-cov.mu$zk, -cov.mu$wi, -cov.mu$yj, cov.mu$Z, cov.mu$W, cov.mu$Y), ]
      cov.mu <- t(cov.mu[, 1:(ncol(cov.mu) - 6)])
      mu.hat <- cbind(as.data.frame(mu.hat), model.data[, 1:6])
      mu.hat <- mu.hat[order(-mu.hat$zk, -mu.hat$wi, -mu.hat$yj, mu.hat$Z, mu.hat$W, mu.hat$Y), ]
      mu.hat <- as.matrix(mu.hat[, 1])
      model.data <- model.data[order(
        -model.data$zk, -model.data$wi, -model.data$yj,
        model.data$Z, model.data$W, model.data$Y
      ), ]
    }
  }
  # Get q (used for computing asymptotic variance of model estimated ORs)
  if (nvars == 2) {
    q <- matrix(data = NA, nrow = nrows, ncol = ((2^nvars) * nrows))
    for (i in 1:nrows) {
      delta <- as.matrix(c(
        rep(0, (i - 1)), (1 / mu.hat[i]), rep(0, (nrows - 1)),
        (-1 / mu.hat[(nrows + i)]), rep(0, (nrows - 1)),
        (-1 / mu.hat[(2 * nrows + i)]), rep(0, (nrows - 1)),
        (1 / mu.hat[(3 * nrows + i)]), rep(0, (nrows - i))
      ))
      q[i, ] <- delta
    }
    logOR.obs <- matrix(data = NA, nrow = nrows, ncol = 3)
    var.logOR.obs <- matrix(data = NA, nrow = nrows, ncol = 1)
    logOR.est <- matrix(data = NA, nrow = nrows, ncol = 3)
  }
  if (nvars == 3) {
    q <- matrix(data = NA, nrow = 2 * nrows, ncol = ((2^nvars) * nrows))
    for (i in 1:nrows) {
      delta <- as.matrix(c(
        rep(0, (i - 1)), (1 / mu.hat[i]), rep(0, (nrows - 1)),
        (-1 / mu.hat[(nrows + i)]), rep(0, (nrows - 1)),
        (-1 / mu.hat[(2 * nrows + i)]), rep(0, (nrows - 1)),
        (1 / mu.hat[(3 * nrows + i)]), rep(0, (nrows - i)),
        rep(0, (4 * nrows))
      ))
      q[i, ] <- delta
    }
    for (i in 1:nrows) {
      delta <- as.matrix(c(
        rep(0, (4 * nrows)), rep(0, (i - 1)), (1 / mu.hat[(i + 4 * nrows)]), rep(0, (nrows - 1)),
        (-1 / mu.hat[(i + 5 * nrows)]), rep(0, (nrows - 1)),
        (-1 / mu.hat[(i + 6 * nrows)]), rep(0, (nrows - 1)),
        (1 / mu.hat[(i + 7 * nrows)]), rep(0, (nrows - i))
      ))
      q[(nrows + i), ] <- delta
    }
    logOR.obs <- matrix(data = NA, nrow = 2 * nrows, ncol = 3)
    var.logOR.obs <- matrix(data = NA, nrow = 2 * nrows, ncol = 1)
    logOR.est <- matrix(data = NA, nrow = 2 * nrows, ncol = 3)
  }
  # Asymptotic covariance matrix for model estimated ORs
  var.logOR.asymp <- q %*% tcrossprod(cov.mu, q)
  var.logOR.asymp <- as.matrix(diag(var.logOR.asymp))
  for (i in 1:nrows) {
    # Observed ORs and corresponding CIs
    logOR.obs[i, 1] <- (log(model.data[i, (2 * nvars + 1)])
    - log(model.data[(nrows + i), (2 * nvars + 1)])
      - log(model.data[(2 * nrows + i), (2 * nvars + 1)])
      + log(model.data[(3 * nrows + i), (2 * nvars + 1)]))
    var.logOR.obs[i, 1] <- (1 / model.data[i, (2 * nvars + 1)]
      + 1 / model.data[(nrows + i), (2 * nvars + 1)]
      + 1 / model.data[(2 * nrows + i), (2 * nvars + 1)]
      + 1 / model.data[(3 * nrows + i), (2 * nvars + 1)])
    logOR.obs[i, 2] <- (logOR.obs[i, 1] - qnorm(1 - alpha / 2) * sqrt(var.logOR.obs[i, 1]))
    logOR.obs[i, 3] <- (logOR.obs[i, 1] + qnorm(1 - alpha / 2) * sqrt(var.logOR.obs[i, 1]))
    # Model estimated ORs and corresponding CIs based on asymptotic variance
    logOR.est[i, 1] <- (log(mu.hat[i]) - log(mu.hat[(nrows + i)])
      - log(mu.hat[(2 * nrows + i)]) + log(mu.hat[(3 * nrows + i)]))
    # For SPMI, prevent warnings/NAs by setting lower and upper CI bounds to 1
    if (model == "spmi") {
      logOR.est[i, 2] <- log(1)
      logOR.est[i, 3] <- log(1)
    }
    if (model != "spmi") {
      logOR.est[i, 2] <- (logOR.est[i, 1] - qnorm(1 - alpha / 2) * sqrt(var.logOR.asymp[i, 1]))
      logOR.est[i, 3] <- (logOR.est[i, 1] + qnorm(1 - alpha / 2) * sqrt(var.logOR.asymp[i, 1]))
    }
  }
  if (nvars == 3) {
    for (i in (4 * nrows + 1):(5 * nrows)) {
      # Observed ORs and corresponding CIs
      logOR.obs[(i - 3 * nrows), 1] <- (log(model.data[i, (2 * nvars + 1)])
      - log(model.data[(nrows + i), (2 * nvars + 1)])
        - log(model.data[(2 * nrows + i), (2 * nvars + 1)])
        + log(model.data[(3 * nrows + i), (2 * nvars + 1)]))
      var.logOR.obs[(i - 3 * nrows), 1] <- (1 / model.data[i, (2 * nvars + 1)]
        + 1 / model.data[(nrows + i), (2 * nvars + 1)]
        + 1 / model.data[(2 * nrows + i), (2 * nvars + 1)]
        + 1 / model.data[(3 * nrows + i), (2 * nvars + 1)])
      logOR.obs[(i - 3 * nrows), 2] <- (logOR.obs[(i - 3 * nrows), 1]
      - qnorm(1 - alpha / 2) * sqrt(var.logOR.obs[(i - 3 * nrows), 1]))
      logOR.obs[(i - 3 * nrows), 3] <- (logOR.obs[(i - 3 * nrows), 1]
      + qnorm(1 - alpha / 2) * sqrt(var.logOR.obs[(i - 3 * nrows), 1]))
      # Model estimated ORs and corresponding CIs based on asymptotic variance
      logOR.est[(i - 3 * nrows), 1] <- ((log(mu.hat[i]) - log(mu.hat[(nrows + i)])
        - log(mu.hat[(2 * nrows + i)]) + log(mu.hat[(3 * nrows + i)])))
      logOR.est[(i - 3 * nrows), 2] <- (logOR.est[(i - 3 * nrows), 1]
      - qnorm(1 - alpha / 2) * sqrt(var.logOR.asymp[(i - 3 * nrows), 1]))
      logOR.est[(i - 3 * nrows), 3] <- (logOR.est[(i - 3 * nrows), 1]
      + qnorm(1 - alpha / 2) * sqrt(var.logOR.asymp[(i - 3 * nrows), 1]))
    }
    for (i in 1:(2 * nrows)) {
      if (abs(exp(logOR.est[i, 1]) - 1) < .00000001) {
        logOR.est[i, 2] <- log(1)
        logOR.est[i, 3] <- log(1)
      }
    }
  }
  OR.obs <- exp(logOR.obs)
  OR.model.asymp <- exp(logOR.est)
  counter <- 0
  if (nvars == 2) {
    rowlabels <- matrix(data = NA, nrow = 1, ncol = nrows)
    for (i in 1:I) {
      for (j in (I + 1):(I + J)) {
        counter <- counter + 1
        cell <- paste(colnames(data)[i], colnames(data)[j], sep = "")
        rowlabels[, counter] <- cell
      }
    }
  }
  if (nvars == 3) {
    rowlabels <- matrix(data = NA, nrow = 1, ncol = 2 * nrows)
    if (any(pair == "YZ" | pair == "ZY" | pair == "yz" | pair == "zy")) {
      for (h in 1:0) {
        for (i in 1:I) {
          for (j in (I + 1):(I + J)) {
            for (k in (I + J + 1):(I + J + K)) {
              counter <- counter + 1
              cell <- paste(colnames(data)[i], "=", h, ",", colnames(data)[j], colnames(data)[k], sep = "")
              rowlabels[, counter] <- cell
            }
          }
        }
      }
    }
    if (any(pair == "WZ" | pair == "ZW" | pair == "wz" | pair == "zw")) {
      for (h in 1:0) {
        for (j in (I + 1):(I + J)) {
          for (i in 1:I) {
            for (k in (I + J + 1):(I + J + K)) {
              counter <- counter + 1
              cell <- paste(colnames(data)[j], "=", h, ",", colnames(data)[i], colnames(data)[k], sep = "")
              rowlabels[, counter] <- cell
            }
          }
        }
      }
    }
    if (any(pair == "WY" | pair == "YW" | pair == "wy" | pair == "yw")) {
      for (h in 1:0) {
        for (k in (I + J + 1):(I + J + K)) {
          for (i in 1:I) {
            for (j in (I + 1):(I + J)) {
              counter <- counter + 1
              cell <- paste(colnames(data)[k], "=", h, ",", colnames(data)[i], colnames(data)[j], sep = "")
              rowlabels[, counter] <- cell
            }
          }
        }
      }
    }
  }
  colnames(OR.obs) <- c("OR", "lower.bound", "upper.bound")
  rownames(OR.obs) <- rowlabels
  colnames(OR.model.asymp) <- c("OR", "lower.bound", "upper.bound")
  rownames(OR.model.asymp) <- rowlabels
  original.arg <- list(data = data, I = I, J = J, K = K, nvars = nvars, alpha = alpha)
  output <- list(
    original.arg = original.arg, OR.obs = OR.obs,
    OR.model.asymp = OR.model.asymp
  )

  if (boot) {
    B.use <- object$boot.results$B.use
    B.discard <- object$boot.results$B.discard
    if (model == "spmi") {
      bca.ci.lower <- matrix(data = 1, nrow = 1, ncol = nrows)
      bca.ci.upper <- matrix(data = 1, nrow = 1, ncol = nrows)
    }
    if (model != "spmi") {
      mu.hat.star <- object$boot.results$mod.fit.star[, (1:((2^nvars) * nrows))]
      # Perform jackknife calculations
      data.n_1 <- merge(data, 1:n, by = NULL, all = TRUE)
      data.n_1 <- data.n_1[-seq(from = 1, to = n * n, by = (n + 1)), ]
      # Create progress bar for bootstrapping
      if (print.status) {
        cat("Bootstrap Progress:", "\n")
        MRCV_globals$pb <- txtProgressBar(min = 0, max = (1.7 * n + n), style = 3)
        MRCV_globals$pb.counter <- 0
      }
      if (nvars == 2) {
        model.data.unsorted.n_1 <- by(
          data = data.n_1[, 1:(I + J)],
          INDICES = data.n_1[, (I + J + 1)],
          FUN = data.format, I = I, J = J, K = K,
          nvars = nvars, add.constant =
            add.constant, predict.func = TRUE
        )
      }
      if (nvars == 3) {
        model.data.unsorted.n_1 <- by(
          data = data.n_1[, 1:(I + J + K)],
          INDICES = data.n_1[, (I + J + K + 1)],
          FUN = data.format, I = I, J = J, K = K,
          nvars = nvars, add.constant =
            add.constant, predict.func = TRUE
        )
      }
      model.data.unsorted.n_1 <- as.data.frame(do.call(rbind, model.data.unsorted.n_1))
      model.data.unsorted.n_1 <- data.frame(model.data.unsorted.n_1, index = rep(c(1:n),
        each = ((2^nvars) * nrows)
      ))
      # if (class(model)=="formula") {
      if (inherits(model, "formula")) {
        for (i in 1:I) {
          parm <- paste("W", i, sep = "")
          if (length(agrep(parm, model, max.distance = 0)) > 0) {
            model.data.unsorted.n_1 <- data.frame(
              model.data.unsorted.n_1,
              as.numeric((model.data.unsorted.n_1[, 1] == names(data)[i]))
            )
            colnames(model.data.unsorted.n_1)[ncol(model.data.unsorted.n_1)] <- parm
          }
        }
        for (j in 1:J) {
          parm <- paste("Y", j, sep = "")
          if (length(agrep(parm, model, max.distance = 0)) > 0) {
            model.data.unsorted.n_1 <- data.frame(
              model.data.unsorted.n_1,
              as.numeric((model.data.unsorted.n_1[, 2] == names(data)[(I + j)]))
            )
            colnames(model.data.unsorted.n_1)[ncol(model.data.unsorted.n_1)] <- parm
          }
        }
        if (nvars == 3) {
          for (k in 1:K) {
            parm <- paste("Z", k, sep = "")
            if (length(agrep(parm, model, max.distance = 0)) > 0) {
              model.data.unsorted.n_1 <- data.frame(
                model.data.unsorted.n_1,
                as.numeric((model.data.unsorted.n_1[, 3] == names(data)[(I + J + k)]))
              )
              colnames(model.data.unsorted.n_1)[ncol(model.data.unsorted.n_1)] <- parm
            }
          }
        }
      }
      var.cols <- c(1:(2 * nvars + 1))
      if (ncol(model.data.unsorted.n_1) > (2 * nvars + 2)) {
        var.cols <- c(1:(2 * nvars + 1), (2 * nvars + 3):ncol(model.data.unsorted.n_1))
      }
      mu.hat.n_1 <- by(
        data = model.data.unsorted.n_1[, var.cols],
        INDICES = model.data.unsorted.n_1[, (2 * nvars + 2)],
        FUN = genloglin.fit, model = model, nvars = nvars,
        limit.output = TRUE
      )
      if (print.status) {
        setTxtProgressBar(MRCV_globals$pb, (1.7 * n + n))
        close(MRCV_globals$pb)
      }
      model.data.unsorted.n_1 <- model.data.unsorted.n_1[, 1:(2 * nvars + 2)]
      mu.hat.n_1 <- t(as.data.frame(matrix(do.call(rbind, mu.hat.n_1),
        nrow = n,
        ncol = ((2^nvars) * nrows)
      )))
      mu.hat.n_1 <- data.frame(
        model.data.unsorted.n_1[(1:((2^nvars) * nrows)), 1:(2 * nvars)],
        mu.hat.n_1
      )
      if (nvars == 2) {
        mu.hat.n_1 <- mu.hat.n_1[order(-mu.hat.n_1$wi, -mu.hat.n_1$yj), ]
      }
      if (nvars == 3) {
        mu.hat.n_1 <- mu.hat.n_1[order(
          -mu.hat.n_1$wi, -mu.hat.n_1$yj, -mu.hat.n_1$zk,
          mu.hat.n_1$W, mu.hat.n_1$Y, mu.hat.n_1$Z
        ), ]
        if (any(pair == "WZ" | pair == "ZW" | pair == "wz" | pair == "zw")) {
          mu.hat.n_1 <- mu.hat.n_1[order(
            -mu.hat.n_1$yj, -mu.hat.n_1$wi, -mu.hat.n_1$zk,
            mu.hat.n_1$Y, mu.hat.n_1$W, mu.hat.n_1$Z
          ), ]
          model.data <- model.data[order(
            -model.data$wi, -model.data$yj, -model.data$zk,
            model.data$W, model.data$Y, model.data$Z
          ), ]
          mu.hat.star <- cbind(t(mu.hat.star), model.data[, 1:6])
          mu.hat.star <- mu.hat.star[order(
            -mu.hat.star$yj, -mu.hat.star$wi, -mu.hat.star$zk,
            mu.hat.star$Y, mu.hat.star$W, mu.hat.star$Z
          ), ]
          mu.hat.star <- t(mu.hat.star[, 1:(ncol(mu.hat.star) - 6)])
        }
        if (any(pair == "WY" | pair == "YW" | pair == "wy" | pair == "yw")) {
          mu.hat.n_1 <- mu.hat.n_1[order(
            -mu.hat.n_1$zk, -mu.hat.n_1$wi, -mu.hat.n_1$yj,
            mu.hat.n_1$Z, mu.hat.n_1$W, mu.hat.n_1$Y
          ), ]
          model.data <- model.data[order(
            -model.data$wi, -model.data$yj, -model.data$zk,
            model.data$W, model.data$Y, model.data$Z
          ), ]
          mu.hat.star <- cbind(t(mu.hat.star), model.data[, 1:6])
          mu.hat.star <- mu.hat.star[order(
            -mu.hat.star$zk, -mu.hat.star$wi, -mu.hat.star$yj,
            mu.hat.star$Z, mu.hat.star$W, mu.hat.star$Y
          ), ]
          mu.hat.star <- t(mu.hat.star[, 1:(ncol(mu.hat.star) - 6)])
        }
      }
      mu.hat.n_1 <- mu.hat.n_1[, -c(1:(2 * nvars))]
      # Calculate ORs for bootstrap resamples and jackknife samples
      if (nvars == 2) {
        or.star <- matrix(data = NA, nrow = B.use, ncol = nrows)
        or.est <- matrix(data = NA, nrow = 1, ncol = nrows)
        or.est.jack <- matrix(data = NA, nrow = n, ncol = nrows)
      }
      if (nvars == 3) {
        or.star <- matrix(data = NA, nrow = B.use, ncol = 2 * nrows)
        or.est <- matrix(data = NA, nrow = 1, ncol = 2 * nrows)
        or.est.jack <- matrix(data = NA, nrow = n, ncol = 2 * nrows)
      }
      for (i in 1:nrows) {
        or.star[, i] <- exp(log(mu.hat.star[, i]) - log(mu.hat.star[, (nrows + i)])
          - log(mu.hat.star[, (2 * nrows + i)])
          + log(mu.hat.star[, (3 * nrows + i)]))
        or.star[, i] <- or.star[do.call("order", as.data.frame(or.star[, i])), i]
        or.est.jack[, i] <- apply(
          X = mu.hat.n_1, MARGIN = 2, FUN = est.jack, i = i,
          I = I, J = J, K = K
        )
      }
      if (nvars == 3) {
        for (i in (4 * nrows + 1):(5 * nrows)) {
          or.star[, (i - 3 * nrows)] <- exp(log(mu.hat.star[, i]) - log(mu.hat.star[, (nrows + i)])
            - log(mu.hat.star[, (2 * nrows + i)])
            + log(mu.hat.star[, (3 * nrows + i)]))
          or.star[, (i - 3 * nrows)] <- or.star[do.call("order", as.data.frame(or.star[, (i - 3 * nrows)])), (i - 3 * nrows)]
          or.est.jack[, (i - 3 * nrows)] <- apply(
            X = mu.hat.n_1, MARGIN = 2,
            FUN = est.jack, i = i,
            I = I, J = J, K = K
          )
        }
      }
      # Make matrices conformable for subsequent calculations
      if (nvars == 2) {
        or.est.rep.n <- merge(t(exp(logOR.est)[, 1]), 1:n,
          by = NULL,
          all = TRUE
        )[, 1:nrows]
        or.est.rep.B <- merge(t(exp(logOR.est)[, 1]), 1:B.use,
          by = NULL,
          all = TRUE
        )[, 1:nrows]
      }
      if (nvars == 3) {
        or.est.rep.n <- merge(t(exp(logOR.est)[, 1]), 1:n,
          by = NULL,
          all = TRUE
        )[, 1:(2 * nrows)]
        or.est.rep.B <- merge(t(exp(logOR.est)[, 1]), 1:B.use,
          by = NULL,
          all = TRUE
        )[, 1:(2 * nrows)]
      }
      # Compute BCa confidence intervals
      l.jack <- (n - 1) * (or.est.rep.n - or.est.jack)
      a.jack <- 1 / 6 * colSums(l.jack^3) / colSums(l.jack^2)^(3 / 2)
      w.p <- colSums(or.star <= or.est.rep.B) / (B.use + 1)
      w <- qnorm(p = w.p)
      alpha.bca <- c(alpha / 2, 1 - alpha / 2)
      z.tilde.1 <- w + qnorm(p = alpha.bca[1])
      z.tilde.2 <- w + qnorm(p = alpha.bca[2])
      alpha.tilde.1 <- pnorm(q = w + z.tilde.1 / (1 - a.jack * z.tilde.1))
      alpha.tilde.2 <- pnorm(q = w + z.tilde.2 / (1 - a.jack * z.tilde.2))
      if (nvars == 2) {
        bca.ci.lower <- matrix(data = NA, nrow = 1, ncol = nrows)
        bca.ci.upper <- matrix(data = NA, nrow = 1, ncol = nrows)
      }
      if (nvars == 3) {
        bca.ci.lower <- matrix(data = NA, nrow = 1, ncol = 2 * nrows)
        bca.ci.upper <- matrix(data = NA, nrow = 1, ncol = 2 * nrows)
      }
      for (i in 1:nrows) {
        bca.ci.lower[1, i] <- quantile(or.star[, i], alpha.tilde.1[i])
        bca.ci.upper[1, i] <- quantile(or.star[, i], alpha.tilde.2[i])
      }
      if (nvars == 3) {
        for (i in (nrows + 1):(2 * nrows)) {
          bca.ci.lower[1, i] <- quantile(or.star[, i], alpha.tilde.1[i])
          bca.ci.upper[1, i] <- quantile(or.star[, i], alpha.tilde.2[i])
        }
        for (i in 1:(2 * nrows)) {
          if (abs(exp(logOR.est[i, 1]) - 1) < .00000001) {
            bca.ci.lower[1, i] <- 1
            bca.ci.upper[1, i] <- 1
          }
        }
      }
    }
    OR.model.BCa <- as.matrix(data.frame(
      exp(logOR.est)[, 1],
      t(bca.ci.lower), t(bca.ci.upper)
    ))
    colnames(OR.model.BCa) <- c("OR", "lower.bound", "upper.bound")
    rownames(OR.model.BCa) <- rowlabels
    boot.results <- list(
      B.use = B.use, B.discard = B.discard, OR.model.BCa =
        OR.model.BCa
    )
    output <- list(
      original.arg = original.arg, OR.obs = OR.obs, OR.model.asymp =
        OR.model.asymp, boot.results = boot.results
    )
  }
  class(output) <- "predict.genloglin"
  output
}
