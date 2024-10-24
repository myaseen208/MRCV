#' @name    print.predict.genloglin
#' @aliases print.predict.genloglin
#' @title Control Printed Display of MRCV predict Modeling Objects
#' @description
#' A function used to control the display of output provided by \code{predict.genloglin}.
#' @export
#' @rdname print
print.predict.genloglin <- function(x, ...) {
  op <- options()
  on.exit(options(op))
  options(scipen = 5)
  type <- names(x)
  data <- x$original.arg$data
  nvars <- x$original.arg$nvars
  I <- x$original.arg$I
  J <- x$original.arg$J
  K <- x$original.arg$K
  alpha <- x$original.arg$alpha
  coverage <- paste(((1 - alpha) * 100), "%", sep = "")
  obs <- round(x$OR.obs, 2)
  model.asymp <- round(x$OR.model.asymp, 2)
  if (nvars == 2) {
    OR.obs <- matrix(data = NA, nrow = I, ncol = J)
    OR.model.asymp <- matrix(data = NA, nrow = I, ncol = J)
    for (i in 1:I) {
      for (j in 1:J) {
        OR.obs[i, j] <- paste(obs[((i - 1) * J + j), 1], paste("(", paste(obs[((i - 1) * J + j), 2],
          obs[((i - 1) * J + j), 3],
          sep = ", "
        ), ")", sep = ""))
        OR.model.asymp[i, j] <- paste(
          model.asymp[((i - 1) * J + j), 1],
          paste("(", paste(model.asymp[((i - 1) * J + j), 2],
            model.asymp[((i - 1) * J + j), 3],
            sep = ", "
          ), ")",
          sep = ""
          )
        )
      }
    }
    rownames(OR.obs) <- names(data)[1:I]
    colnames(OR.obs) <- names(data)[(I + 1):(I + J)]
    rownames(OR.model.asymp) <- names(data)[1:I]
    colnames(OR.model.asymp) <- names(data)[(I + 1):(I + J)]
    cat("Observed odds ratios with", coverage, "asymptotic confidence intervals", "\n")
    print.default(OR.obs, quote = FALSE)
    cat("\n")
    cat("Model-predicted odds ratios with", coverage, "asymptotic confidence intervals", "\n")
    print.default(OR.model.asymp, quote = FALSE)
    cat("\n")
  }
  if (nvars == 3) {
    cat("Observed odds ratios with", coverage, "asymptotic confidence intervals", "\n")
    print.default(obs, quote = FALSE)
    cat("\n")
    cat("Model-predicted odds ratios with", coverage, "asymptotic confidence intervals", "\n")
    print.default(model.asymp, quote = FALSE)
    cat("\n")
  }

  if (any(type == "boot.results")) {
    B.use <- x$boot.results$B.use
    B.discard <- x$boot.results$B.discard
    model.BCa <- round(x$boot.results$OR.model.BCa, 2)
    if (nvars == 2) {
      OR.model.BCa <- matrix(data = NA, nrow = I, ncol = J)
      for (i in 1:I) {
        for (j in 1:J) {
          OR.model.BCa[i, j] <- paste(
            model.BCa[((i - 1) * J + j), 1],
            paste("(", paste(model.BCa[((i - 1) * J + j), 2],
              model.BCa[((i - 1) * J + j), 3],
              sep = ", "
            ), ")",
            sep = ""
            )
          )
        }
      }
      rownames(OR.model.BCa) <- names(data)[1:I]
      colnames(OR.model.BCa) <- names(data)[(I + 1):(I + J)]
      cat("Bootstrap Results:", "\n")
      if (B.discard > 0) {
        cat(B.discard, "resamples were removed from the analysis due to")
        cat(" not having all rows or columns represented in a 2x2 table", "\n")
      }
      cat("Final results based on", B.use, "resamples", "\n")
      cat("Model-predicted odds ratios with", coverage, "bootstrap BCa confidence intervals", "\n")
      print.default(OR.model.BCa, quote = FALSE)
      cat("\n")
    }
    if (nvars == 3) {
      cat("Bootstrap Results:", "\n")
      if (B.discard > 0) {
        cat(B.discard, "resamples were removed from the analysis due to")
        cat(" not having all rows or columns represented in a 2x2 table", "\n")
      }
      cat("Final results based on", B.use, "resamples", "\n")
      cat("Model-predicted odds ratios with", coverage, "bootstrap BCa confidence intervals", "\n")
      print.default(model.BCa, quote = FALSE)
      cat("\n")
    }
  }
  invisible(x)
}
