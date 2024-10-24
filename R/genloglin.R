#' @name    genloglin
#' @aliases genloglin
#' @title  Model the Association Among Two or Three MRCVs
#' @description
#' The \code{genloglin} function uses a generalized loglinear modeling approach
#' to estimate the association among two or three MRCVs.  Standard errors are
#' adjusted using a second-order Rao-Scott approach.
#' 
#' @details
#' The \code{genloglin} function first converts the raw data into a form that
#' can be used for estimation.  For the two MRCV case, the reformatted data
#' frame contains 2Ix2J rows and 5 columns generically named \code{W},
#' \code{Y}, \code{wi}, \code{yj}, and \code{count}.  For the three MRCV case,
#' the reformatted data frame contains 2Ix2Jx2K rows and 7 columns generically
#' named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, \code{zk}, and
#' \code{count}.  Then, the model of interest is estimated by calling the
#' \code{\link{glm}} function where the \code{family} argument is specified as
#' \code{poisson}.  For all predictor variables, the first level is the
#' reference group (i.e., 1 is the reference for variables \code{W}, \code{Y},
#' and \code{Z}, and 0 is the reference for variables \code{wi}, \code{yj}, and
#' \code{zj}).
#'
#' The \code{boot} argument must equal \code{TRUE} in order to obtain bootstrap
#' results with the \code{genloglin} method functions.
#'
#' @param data A data frame containing the raw data where rows correspond to
#' the individual item response vectors, and columns correspond to the binary
#' items, W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
#' @param I The number of items corresponding to row variable W.
#' @param J The number of items corresponding to column variable Y.
#' @param K The number of items corresponding to strata variable Z.
#' @param model For the two MRCV case, a character string specifying one of the
#' following models: \code{"spmi"} (the complete independence model),
#' \code{"homogeneous"} (the homogeneous association model), \code{"w.main"}
#' (the w-main effects model), \code{"y.main"} (the y-main effects model),
#' \code{"wy.main"} (the w-main and y-main effects model), or
#' \code{"saturated"}.  Alternatively, a user-supplied formula can be
#' specified, where the formula is limited to the generic variables \code{W},
#' \code{Y}, \code{wi}, \code{yj}, \code{count}, \code{W1},\ldots, \code{WI},
#' and \code{Y1},\ldots, \code{YJ}.  For the three MRCV case, only
#' user-supplied formulas are accepted. In addition to the generic variables
#' defined for two MRCVs, the formula may include the generic variables
#' \code{Z}, \code{zk}, and \code{Z1},\ldots, \code{ZK}.
#' @param add.constant A positive constant to be added to all zero marginal
#' cell counts.
#' @param boot A logical value indicating whether bootstrap resamples should be
#' taken.
#' @param B The desired number of bootstrap resamples.
#' @param B.max The maximum number of bootstrap resamples.  Resamples for which
#' at least one item has all positive or negative responses are thrown out;
#' \code{genloglin} uses the first \code{B} valid resamples or all valid
#' resamples if that number is less than \code{B}.
#' @param print.status A logical value indicating whether progress updates
#' should be provided.  When \code{print.status = TRUE}, the status of the IPF
#' algorithm is printed after every 5 iterations.  Upon completion of the IPF
#' algorithm, a progress bar appears that documents progress of the bootstrap.
#' @return --- \code{genloglin} returns an object of class \code{'genloglin'}.
#' The object is a list containing at least the following objects:
#' \code{original.arg}, \code{mod.fit}, \code{sum.fit}, and \code{rs.results}.
#'
#' \code{original.arg} is a list containing the following objects: \itemize{
#' \item\code{data}: The original data frame supplied to the \code{data}
#' argument. \item\code{I}: The original value supplied to the \code{I}
#' argument. \item\code{J}: The original value supplied to the \code{J}
#' argument. \item\code{K}: The original value supplied to the \code{K}
#' argument. \item\code{nvars}: The number of MRCVs. \item\code{model}: The
#' original value supplied to the \code{model} argument.
#' \item\code{add.constant}: The original value supplied to the
#' \code{add.constant} argument. \item\code{boot}: The original value supplied
#' to the \code{boot} argument. }
#'
#' \code{mod.fit} is a list containing the same objects returned by
#' \code{\link{glm}} with a few modifications as described in
#' \code{\link{summary.genloglin}}.
#'
#' \code{sum.fit} is a list containing the same objects returned by the
#' \code{\link{summary}} method for class \code{"glm"} with a few modifications
#' as described in \code{\link{summary.genloglin}}.
#'
#' \code{rs.results} is a list containing the following objects (see Appendix A
#' of Bilder and Loughin, 2007, for more detail): \itemize{ \item\code{cov.mu}:
#' The covariance matrix for the estimated cell counts. \item\code{E}: The
#' covariance matrix for the residuals. \item\code{gamma}: Eigenvalues used in
#' computing second-order Rao-Scott adjusted statistics. }
#'
#' --- For \code{boot = TRUE}, the primary list additionally includes
#' \code{boot.results}, a list containing the following objects: \itemize{
#' \item\code{B.use}: The number of bootstrap resamples used.
#' \item\code{B.discard}: The number of bootstrap resamples discarded due to
#' having at least one item with all positive or negative responses.
#' \item\code{model.data.star}: For the two MRCV case, a numeric matrix
#' containing 2Ix2J rows and \code{B.use}+4 columns, where the first 4 columns
#' correspond to the model variables \code{W}, \code{Y}, \code{wi}, and
#' \code{yj}, and the last \code{B.use} columns correspond to the observed
#' counts for each resample.  For the three MRCV case, a numeric matrix
#' containing 2Ix2Jx2K rows and \code{B.use+6} columns, where the first 6
#' columns correspond to the model variables \code{W}, \code{Y}, \code{Z},
#' \code{wi}, \code{yj}, and \code{zk}, and the last \code{B.use} columns
#' correspond to the observed counts for each resample.
#' \item\code{mod.fit.star}: For the two MRCV case, a numeric matrix containing
#' \code{B.use} rows and 2Ix2J +1 columns, where the first 2Ix2J columns
#' correspond to the model-predicted counts for each resample, and the last
#' column corresponds to the residual deviance for each resample.  For the
#' three MRCV case, a numeric matrix containing \code{B.use} rows and
#' 2Ix2Jx2K+1 columns, where the first 2Ix2Jx2K columns correspond to the
#' model-predicted counts for each resample, and the last column corresponds to
#' the residual deviance for each resample. \item\code{chisq.star}: A numeric
#' vector of length \code{B.use} containing the Pearson statistics (comparing
#' \code{model} to the saturated model) calculated for each resample.
#' \item\code{lrt.star}: A numeric vector of length \code{B.use} containing the
#' LRT statistics calculated for each resample. \item\code{residual.star}: A
#' numeric matrix with 2Ix2J rows (or 2Ix2Jx2K rows for the three MRCV case)
#' and \code{B.use} columns containing the residuals calculated for each
#' resample. }
#' @seealso The \code{genloglin} methods \code{\link{summary.genloglin}},
#' \code{\link{residuals.genloglin}}, \code{\link{anova.genloglin}}, and
#' \code{\link{predict.genloglin}}, and the corresponding generic functions
#' \code{\link{summary}}, \code{\link{residuals}}, \code{\link{anova}}, and
#' \code{\link{predict}}.
#'
#' The \code{\link{glm}} function for fitting generalized linear models.
#'
#' The \code{\link{MI.test}} function for testing for MMI (one MRCV case) or
#' SPMI (two MRCV case).
#' @references Bilder, C. and Loughin, T. (2007) Modeling association between
#' two or more categorical variables that allow for multiple category choices.
#' \emph{Communications in Statistics--Theory and Methods}, \bold{36},
#' 433--451.
#'
#'
#' @examples
#'
#' # Estimate the y-main effects model for 2 MRCVs
#' mod.fit <- genloglin(data = farmer2, I = 3, J = 4, model = "y.main", boot = FALSE)
#' # Summarize model fit information
#' summary(mod.fit)
#' # Examine standardized Pearson residuals
#' residuals(mod.fit)
#' # Compare the y-main effects model to the saturated model
#' anova(mod.fit, model.HA = "saturated", type = "rs2")
#' # Obtain observed and model-predicted odds ratios
#' predict(mod.fit)
#'
#' # Estimate a model that is not one of the named models
#' # Note that this was the final model chosen by Bilder and Loughin (2007)
#' mod.fit.other <- genloglin(
#'   data = farmer2, I = 3, J = 4, model = count ~ -1 + W:Y +
#'     wi %in% W:Y + yj %in% W:Y + wi:yj + wi:yj %in% Y + wi:yj %in% W3:Y1, boot =
#'     FALSE
#' )
#' # Compare this model to the y-main effects model
#' anova(mod.fit, model.HA = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y + wi:yj +
#'   wi:yj %in% Y + wi:yj %in% W3:Y1, type = "rs2", gof = TRUE)
#'
#' # To obtain bootstrap results from the method functions the genloglin() boot
#' # argument must be specified as TRUE (the default)
#' # A small B is used for demonstration purposes; normally, a larger B should be used
#' \dontrun{
#' mod.fit <- genloglin(
#'   data = farmer2, I = 3, J = 4, model = "y.main", boot = TRUE,
#'   B = 99
#' )
#' residuals(mod.fit)
#' anova(mod.fit, model.HA = "saturated", type = "all")
#' predict(mod.fit)
#' }
#'
#' # Estimate a model for 3 MRCVs
#' \dontrun{
#' mod.fit.three <- genloglin(data = farmer3, I = 3, J = 4, K = 5, model = count ~
#'   -1 + W:Y:Z + wi %in% W:Y:Z + yj %in% W:Y:Z + zk %in% W:Y:Z + wi:yj +
#'   wi:yj %in% Y + wi:yj %in% W + wi:yj %in% Y:W + yj:zk + yj:zk %in% Z +
#'   yj:zk %in% Y + yj:zk %in% Z:Y, boot = TRUE, B = 99)
#' residuals(mod.fit.three)
#' anova(mod.fit.three, model.HA = "saturated", type = "all")
#' predict(mod.fit.three, pair = "WY")
#' }
#'
#' @export
genloglin <- function(data, I, J, K = NULL, model, add.constant = .5, boot = TRUE,
                      B = 1999, B.max = B, print.status = TRUE) {
  op <- options()
  on.exit(options(op))
  options(warn = 1)
  # if((class(data)!="data.frame")&(class(data)!="matrix")) {
  #   stop("The \"data\" argument requires an object of class \"data.frame\".")
  # }
  if (!(is.data.frame(data) || is.matrix(data))) {
    stop("The \"data\" argument requires an object of class \"data.frame\".")
  }

  data <- as.data.frame(data)
  if ((I == 0) | (J == 0)) {
    stop("\"I\" and \"J\" must be greater than 0.")
  }
  if (!is.numeric(I) | !is.numeric(J)) {
    stop("\"I\" and \"J\" must be numeric.")
  }
  if ((I %% 1 != 0) | (J %% 1 != 0)) {
    stop("\"I\" and \"J\" must be integers.")
  }
  I <- as.integer(I)
  J <- as.integer(J)
  if (!is.null(K)) {
    if (!is.numeric(K)) {
      stop("\"K\" must be numeric or NULL.")
    }
    if (K %% 1 != 0) {
      stop("\"K\" must be an integer.")
    }
    K <- as.integer(K)
    if (K == 0) {
      K <- NULL
    }
  }
  if (!is.null(K)) {
    # if (class(model)!="formula") {
    if (!inherits(model, "formula")) {
      stop("For the 3 MRCV case, only user-supplied formulas are accepted by the \"model\" argument.")
    }
  }
  # if ((class(model)!="formula")&(model!="spmi")&(model!="homogeneous")&(model!="w.main")&(model!="y.main")&(model!="wy.main")&(model!="saturated")) {
  if (!inherits(model, "formula") & (model != "spmi") & (model != "homogeneous") & (model != "w.main") & (model != "y.main") & (model != "wy.main") & (model != "saturated")) {
    stop("The \"model\" argument requires a formula or one of 6 recognized character strings. \n  See help(genloglin) for details.")
  }

  if (!is.numeric(B) | !is.numeric(B.max)) {
    warning("\"B\" and \"B.max\" must be numeric. \n  The input values have been changed to the default value of 1999.")
    B <- 1999
    B.max <- B
  }
  if (B.max < B) {
    warning("\"B.max\" must be greater than or equal to \"B\". \n  \"B.max\" has been set equal to \"B\".")
    B.max <- B
  }
  if ((B %% 1 != 0) | (B.max %% 1 != 0)) {
    warning("\"B\" and \"B.max\" must be integers. \n  The input values have been rounded up to the nearest integer.")
  }
  B <- as.integer(ceiling(B))
  B.max <- as.integer(ceiling(B.max))
  if (!is.numeric(add.constant)) {
    if (add.constant == FALSE) {
      add.constant <- 0
    }
    if (add.constant != FALSE) {
      warning("The \"add.constant\" argument only accepts numeric values. \n  The input value has been changed to the default value of .5.")
      add.constant <- .5
    }
  }
  if (add.constant < 0) {
    warning("The \"add.constant\" argument cannot be negative. \n  The input value has been changed to the default value of .5.")
    add.constant <- .5
  }

  # if ((class(boot)!="logical")&(boot!=1)&(boot!=0)) {
  #  warning("The \"boot\" argument requires an object of class \"logical\". \n  The input value has been changed to the default value of TRUE.")
  #  boot<-TRUE
  # }
  if (!is.logical(boot) & (boot != 1) & (boot != 0)) {
    warning("The \"boot\" argument requires an object of class \"logical\". \n  The input value has been changed to the default value of TRUE.")
    boot <- TRUE
  }

  # if ((class(print.status)!="logical")&(print.status!=1)&(print.status!=0)) {
  #  warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
  #  print.status<-TRUE
  # }
  if (!is.logical(print.status) & (print.status != 1) & (print.status != 0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status <- TRUE
  }


  boot <- as.logical(boot)
  print.status <- as.logical(print.status)
  MRCV_globals$B <- B
  MRCV_globals$print.status <- print.status
  n <- nrow(data)
  if (is.numeric(K)) {
    if (K == 0) {
      K <- NULL
    }
  }
  nvars <- 2 + is.numeric(K)
  # Reformat raw data
  model.data.unsorted <- data.format(
    data = data, I = I, J = J, K = K, nvars = nvars,
    add.constant = add.constant
  )
  # Need to sort data for later calculations
  if (nvars == 2) {
    model.data <- model.data.unsorted[order(-model.data.unsorted$wi, -model.data.unsorted$yj), ]
  }
  if (nvars == 3) {
    model.data <- model.data.unsorted[order(
      -model.data.unsorted$wi, -model.data.unsorted$yj,
      -model.data.unsorted$zk
    ), ]
  }
  # Add additional variables specified in model formula
  # if (class(model)=="formula") {
  if (inherits(model, "formula")) {
    for (i in 1:I) {
      parm <- paste("W", i, sep = "")
      if (length(agrep(parm, model, max.distance = 0)) > 0) {
        model.data <- data.frame(
          model.data,
          as.numeric((model.data[, 1] == names(data)[i]))
        )
        colnames(model.data)[ncol(model.data)] <- parm
      }
    }
    for (j in 1:J) {
      parm <- paste("Y", j, sep = "")
      if (length(agrep(parm, model, max.distance = 0)) > 0) {
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
        if (length(agrep(parm, model, max.distance = 0)) > 0) {
          model.data <- data.frame(
            model.data,
            as.numeric((model.data[, 3] == names(data)[(I + J + k)]))
          )
          colnames(model.data)[ncol(model.data)] <- parm
        }
      }
    }
  }
  # Estimate model of interest
  mod.fit <- genloglin.fit(data = model.data, model = model, nvars = nvars)
  model.data <- model.data[, 1:(2 * nvars + 1)]
  mu.hat <- mod.fit$fitted.values
  sum.fit <- summary(mod.fit)

  # RS Calculations (Appendix A)
  W.counts <- as.data.frame(table(data[, 1:I])) # Get all possible combos of W's
  cols <- c(1:I) # Need to order W's in ascending order starting with W1
  W.counts <- W.counts[do.call("order", as.data.frame(W.counts[, cols])), ]
  Y.counts <- as.data.frame(table(data[, (I + 1):(I + J)])) # All possible Y's
  cols <- c(1:J) # Need to order Y's in ascending order starting with Y1
  Y.counts <- Y.counts[do.call("order", as.data.frame(Y.counts[, cols])), ]
  if (nvars == 3) {
    Z.counts <- as.data.frame(table(data[, (I + J + 1):(I + J + K)])) # All possible Z's
    cols <- c(1:K) # Need to order Z's in ascending order starting with Z1
    Z.counts <- Z.counts[do.call("order", as.data.frame(Z.counts[, cols])), ]
  }
  n.counts <- as.data.frame(table(data)) # Get all possible combos of W's, Y's, and Z's
  cols <- c(1:(ncol(data))) # Need to order W's, Y's, Z's in ascending order W1, W2, ...
  n.counts <- n.counts[do.call("order", as.data.frame(n.counts[, cols])), ]
  G <- t(data.matrix(W.counts[, 1:I]) - 1) # rx2^r matrix
  H <- t(data.matrix(Y.counts[, 1:J]) - 1) # cx2^c matrix
  if (nvars == 3) {
    L <- t(data.matrix(Z.counts[, 1:K]) - 1)
  }
  tau <- n.counts[, ncol(n.counts)] / n # Vector of multinomial probabilities
  Jr <- matrix(data = 1, nrow = I, ncol = 2^I)
  Jc <- matrix(data = 1, nrow = J, ncol = 2^J)
  if (nvars == 3) {
    Jq <- matrix(data = 1, nrow = K, ncol = 2^K)
  }
  if (nvars == 2) {
    B.matrix <- rbind(
      kronecker(G, H), kronecker(G, (Jc - H)), kronecker((Jr - G), H),
      kronecker((Jr - G), (Jc - H))
    )
  }
  if (nvars == 3) {
    B.matrix <- rbind(
      kronecker(kronecker(G, H), L), kronecker(kronecker(G, H), (Jq - L)),
      kronecker(kronecker(G, (Jc - H)), L), kronecker(kronecker(G, (Jc - H)), (Jq - L)),
      kronecker(kronecker((Jr - G), H), L), kronecker(kronecker((Jr - G), H), (Jq - L)),
      kronecker(kronecker((Jr - G), (Jc - H)), L), kronecker(kronecker((Jr - G), (Jc - H)), (Jq - L))
    )
  }
  V <- n * B.matrix %*% tcrossprod((diag(tau) - tcrossprod(tau)), B.matrix) # Asymptotic variance for m
  X <- model.matrix(mod.fit) # Design matrix
  # Covariance matrix for B.hat
  d.mu.hat <- diag(mu.hat)
  XP.mu.X.inv <- solve(crossprod(X, d.mu.hat) %*% X)
  sigma <- tcrossprod(XP.mu.X.inv, X) %*% V %*% X %*% XP.mu.X.inv
  rs.se <- sqrt(diag(sigma)) # RS2 standard errors for B.hats
  cov.mu <- d.mu.hat %*% X %*% tcrossprod(sigma, X) %*% d.mu.hat # Cov matrix for mu.hat
  i.matrix <- diag(nrow(mod.fit$data))
  # Covariance matrix for the residuals (m - mu.hat)
  E.part <- i.matrix - d.mu.hat %*% X %*% tcrossprod(XP.mu.X.inv, X)
  E <- E.part %*% tcrossprod(V, E.part)
  gamma <- Re(eigen(diag(1 / mu.hat) %*% E)$values) # Eigenvalues (only use real part)
  # Swap default output with RS2 results
  sum.fit$coefficients[, 2] <- rs.se
  sum.fit$coefficients[, 3] <- sum.fit$coefficients[, 1] / sum.fit$coefficients[, 2]
  sum.fit$coefficients[, 4] <- 2 * (1 - pnorm(abs(sum.fit$coefficients[, 3])))
  rnames <- names(mod.fit$coefficients)
  cnames <- c("Estimate", "RS SE", "z value", "Pr(>|z|)")
  dimnames(sum.fit$coefficients) <- list(rnames, cnames)
  sum.fit$cov.unscaled <- sigma
  sum.fit$cov.scaled <- sigma

  if (boot) {
    # Use the "Gange bootstrap" algorithm
    # Get observed pairwise counts for W's
    counter <- 0
    w.m <- data.frame(matrix(data = NA, nrow = 4 * choose(I, 2), ncol = 5))
    for (i in 1:(I - 1)) {
      for (j in (i + 1):I) {
        counter <- counter + 4
        w.m[(counter - 3), ] <- c(
          names(data)[(i)], names(data)[(j)], 0, 0,
          table(data[, i], data[, j])[1, 1]
        )
        w.m[(counter - 2), ] <- c(
          names(data)[(i)], names(data)[(j)], 0, 1,
          table(data[, i], data[, j])[1, 2]
        )
        w.m[(counter - 1), ] <- c(
          names(data)[(i)], names(data)[(j)], 1, 0,
          table(data[, i], data[, j])[2, 1]
        )
        w.m[(counter), ] <- c(
          names(data)[(i)], names(data)[(j)], 1, 1,
          table(data[, i], data[, j])[2, 2]
        )
      }
    }
    # Get observed pairwise counts for Y's
    counter <- 0
    y.m <- data.frame(matrix(data = NA, nrow = 4 * choose(J, 2), ncol = 5))
    for (i in 1:(J - 1)) {
      for (j in (i + 1):J) {
        counter <- counter + 4
        y.m[(counter - 3), ] <- c(
          names(data)[(i + I)], names(data)[(j + I)], 0, 0,
          table(data[, (i + I)], data[, (j + I)])[1, 1]
        )
        y.m[(counter - 2), ] <- c(
          names(data)[(i + I)], names(data)[(j + I)], 0, 1,
          table(data[, (i + I)], data[, (j + I)])[1, 2]
        )
        y.m[(counter - 1), ] <- c(
          names(data)[(i + I)], names(data)[(j + I)], 1, 0,
          table(data[, (i + I)], data[, (j + I)])[2, 1]
        )
        y.m[(counter), ] <- c(
          names(data)[(i + I)], names(data)[(j + I)], 1, 1,
          table(data[, (i + I)], data[, (j + I)])[2, 2]
        )
      }
    }
    if (nvars == 3) {
      # Get observed pairwise counts for Z's
      counter <- 0
      z.m <- data.frame(matrix(data = NA, nrow = 4 * choose(K, 2), ncol = 5))
      for (i in 1:(K - 1)) {
        for (j in (i + 1):K) {
          counter <- counter + 4
          z.m[(counter - 3), ] <- c(
            names(data)[(i + I + J)], names(data)[(j + I + J)], 0, 0,
            table(data[, (i + I + J)], data[, (j + I + J)])[1, 1]
          )
          z.m[(counter - 2), ] <- c(
            names(data)[(i + I + J)], names(data)[(j + I + J)], 0, 1,
            table(data[, (i + I + J)], data[, (j + I + J)])[1, 2]
          )
          z.m[(counter - 1), ] <- c(
            names(data)[(i + I + J)], names(data)[(j + I + J)], 1, 0,
            table(data[, (i + I + J)], data[, (j + I + J)])[2, 1]
          )
          z.m[(counter), ] <- c(
            names(data)[(i + I + J)], names(data)[(j + I + J)], 1, 1,
            table(data[, (i + I + J)], data[, (j + I + J)])[2, 2]
          )
        }
      }
    }
    # Get model estimated counts for each pair
    nrows <- (2^nvars) * I * J * max(1, K)
    est.m <- data.frame(matrix(data = NA, nrow = nrows, ncol = (2 * nvars + 1)))
    est.m[, 1:(ncol(est.m) - 1)] <- model.data[, 1:(ncol(est.m) - 1)]
    est.m[, ncol(est.m)] <- mu.hat
    # Initialize multinomial probability matrix
    p <- n.counts
    p[, ncol(p)] <- .5
    colnames(p) <- c(names(data), "p")
    if (nvars == 2) {
      x.theta.2 <- rbind(w.m, y.m, est.m)
    }
    if (nvars == 3) {
      x.theta.2 <- rbind(w.m, y.m, z.m)
    }
    x.theta.2[, 1:2] <- lapply(x.theta.2[, 1:2], as.factor)
    x.theta.2[, 3:5] <- lapply(x.theta.2[, 3:5], as.numeric)
    x.theta.2[, 5] <- apply(
      X = as.matrix(x.theta.2[, 5]), MARGIN = 1, FUN = check.zero,
      add.constant = add.constant
    )
    x.theta.2[, 5] <- x.theta.2[, 5] / n
    p.theta.2 <- x.theta.2
    p.theta.2[, 5] <- 0
    p.theta.3 <- NULL
    x.theta.3 <- NULL
    if (nvars == 3) {
      est.m[, 7] <- apply(
        X = as.matrix(est.m[, 7]), MARGIN = 1, FUN = check.zero,
        add.constant = add.constant
      )
      est.m[, 7] <- est.m[, 7] / n
      p.theta.3 <- est.m
      p.theta.3[, 7] <- 0
      x.theta.3 <- est.m
    }
    tol <- 0.00000001
    save.p <- 1
    counter <- 1
    # Use the iterative proportional fitting algorithm
    while (max(abs(save.p - p[, ncol(p)])) > tol) {
      save.p <- p[, ncol(p)]
      p <- ipf.genloglin(
        data = data, I = I, J = J, K = K, nvars = nvars, p = p,
        p.theta.2 = p.theta.2, p.theta.3 = p.theta.3, x.theta.2 = x.theta.2,
        x.theta.3 = x.theta.3
      )
      if (print.status) {
        if (((counter %% 5) == 0) | (counter == 1)) {
          cat(counter, "iterations of the iterative proportional fitting algorithm completed", "\n")
        }
      }
      counter <- counter + 1
    }
    if (print.status) {
      cat(counter, "total iterations required", "\n")
    }
    # Generate bootstrap resamples under the null hypothesis
    boot.sample <- cbind(p[, 1:(ncol(p) - 1)], rmultinom(
      n = B.max, size = n,
      prob = p[, ncol(p)]
    ))
    model.data.star.unsorted <- data.frame(matrix(
      data = NA, nrow = nrows,
      ncol = (2 * nvars + B.max)
    ))
    model.data.star.unsorted[, 1:(2 * nvars)] <- model.data.unsorted[, 1:(2 * nvars)]
    # Create progress bar for bootstrapping
    if (print.status) {
      cat("Bootstrap Progress:", "\n")
      weight.B.max <- ((5 + 2 / 3) * B) / B.max
      MRCV_globals$weight.B.max <- weight.B.max
      MRCV_globals$pb <- txtProgressBar(min = 0, max = (B.max * weight.B.max + B), style = 3)
      MRCV_globals$pb.counter <- 0
    }
    # Reformat the bootstrap resamples into model form
    model.data.star.unsorted[, (2 * nvars + 1):(2 * nvars + B.max)] <- apply(
      X =
        as.matrix(boot.sample[, ((ncol(p)):(ncol(boot.sample)))]), MARGIN = 2,
      FUN = data.format, I = I, J = J, K = K, nvars = nvars,
      add.constant = add.constant, model.vars = boot.sample[, (1:(ncol(p) - 1))]
    )
    if (nvars == 2) {
      colnames(model.data.star.unsorted) <- c("W", "Y", "wi", "yj", rep("count", B.max))
      model.data.star <- model.data.star.unsorted[order(
        -model.data.star.unsorted$wi,
        -model.data.star.unsorted$yj
      ), ]
    }
    if (nvars == 3) {
      colnames(model.data.star.unsorted) <- c(
        "W", "Y", "Z", "wi", "yj", "zk",
        rep("count", B.max)
      )
      model.data.star <- model.data.star.unsorted[order(
        -model.data.star.unsorted$wi,
        -model.data.star.unsorted$yj,
        -model.data.star.unsorted$zk
      ), ]
    }
    # Only keep valid resamples (no items with all positive or negative counts)
    keep <- apply(
      X = as.matrix(model.data.star[, (2 * nvars + 1):(2 * nvars + B.max)]), MARGIN = 2,
      FUN = check.margins, I = I, J = J, K = K, nvars = nvars,
      model.vars = model.data.star[, 1:(2 * nvars)], item.names = names(data)
    )
    model.data.star <- model.data.star[, c(rep(TRUE, 2 * nvars), keep)]
    B.discard <- B.max - ncol(model.data.star) + 2 * nvars
    # Only keep B resamples (or all valid resamples if value < B)
    model.data.star <- model.data.star[, 1:(min((B + 2 * nvars), (ncol(model.data.star) + 2 * nvars)))]
    B.use <- ncol(model.data.star) - 2 * nvars
    MRCV_globals$B.use <- B.use
    # Add additional variables specified in model formula
    # if (class(model)=="formula") {
    if (inherits(model, "formula")) {
      for (i in 1:I) {
        parm <- paste("W", i, sep = "")
        if (length(agrep(parm, model, max.distance = 0)) > 0) {
          model.data.star <- data.frame(
            model.data.star,
            as.numeric((model.data.star[, 1] == names(data)[i]))
          )
          colnames(model.data.star)[ncol(model.data.star)] <- parm
        }
      }
      for (j in 1:J) {
        parm <- paste("Y", j, sep = "")
        if (length(agrep(parm, model, max.distance = 0)) > 0) {
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
          if (length(agrep(parm, model, max.distance = 0)) > 0) {
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
    if (ncol(model.data.star) > (B.use + 2 * nvars)) {
      var.cols <- c(1:(2 * nvars), (B.use + 2 * nvars + 1):ncol(model.data.star))
    }
    # Estimate null model for each resample and keep mu.hat and deviance output
    mod.fit.star <- t(apply(
      X = as.matrix(model.data.star[, (2 * nvars + 1):(B.use + 2 * nvars)]),
      MARGIN = 2, FUN = genloglin.fit, model = model, nvars = nvars,
      limit.output = TRUE, model.vars = model.data.star[, var.cols]
    ))
    if (print.status) {
      setTxtProgressBar(MRCV_globals$pb, (B.max * weight.B.max + B))
      close(MRCV_globals$pb)
    }
    model.data.star <- model.data.star[, 1:(B.use + 2 * nvars)]
    mu.hat.star <- mod.fit.star[, (1:nrows)]
    deviance.star <- mod.fit.star[, ncol(mod.fit.star)]
    # Assume alternative model is saturated (user can specify a different alternative
    #   model in anova.genloglin())
    mu.hat.star.HA <- t(model.data.star[, (2 * nvars + 1):(B.use + 2 * nvars)])
    # Get t*'s
    chisq.star <- rowSums(((mu.hat.star.HA - mu.hat.star)^2) / mu.hat.star)
    lrt.star <- deviance.star
    residual.star <- mu.hat.star.HA - mu.hat.star

    boot.results <- list(
      B.use = B.use, B.discard = B.discard, model.data.star =
        model.data.star, mod.fit.star = mod.fit.star, chisq.star =
        chisq.star, lrt.star = lrt.star, residual.star =
        residual.star
    )
  }
  mod.fit <- mod.fit[-11]
  sum.fit <- sum.fit[-5]
  original.arg <- list(
    data = data, I = I, J = J, K = K, nvars = nvars, model = model,
    add.constant = add.constant, boot = boot
  )
  rs.results <- list(cov.mu = cov.mu, E = E, gamma = gamma)
  if (!boot) {
    output <- list(
      original.arg = original.arg, mod.fit = mod.fit, sum.fit = sum.fit,
      rs.results = rs.results
    )
  }
  if (boot) {
    output <- list(
      original.arg = original.arg, mod.fit = mod.fit, sum.fit = sum.fit,
      rs.results = rs.results, boot.results = boot.results
    )
  }
  class(output) <- "genloglin"
  output
}
