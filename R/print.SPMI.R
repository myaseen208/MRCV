#' @name    print.SPMI
#' @aliases print.SPMI
#' @title Control Printed Display of Objects of Class "SPMI"
#' @description
#' The \code{print.SPMI} method function controls the printed display of
#' objects of class \code{SPMI}.
#'
#'
#' @param x An object of class \code{"SPMI"} produced by the
#' \code{\link{MI.test}} function.
#' @param ...  Additional arguments passed to or from other methods.
#' @return A print out of selected results from \code{\link{MI.test}}.
#' @export
#' @rdname print
print.SPMI <- function(x, ...) {
  op <- options()
  on.exit(options(op))
  options(scipen = 10)
  type <- names(x)
  data <- x$general$data
  I <- x$general$I
  J <- x$general$J
  summary.data <- x$general$summary.data
  X.sq.S <- round(x$general$X.sq.S, 2)
  X.sq.S.ij <- round(x$general$X.sq.S.ij, 2)
  if (summary.data) {
    rownames(X.sq.S.ij) <- levels(data[, 1])
    colnames(X.sq.S.ij) <- levels(data[, 2])
  }
  if (!summary.data) {
    rownames(X.sq.S.ij) <- names(data)[1:I]
    colnames(X.sq.S.ij) <- names(data)[(I + 1):(I + J)]
  }
  cat("Test for Simultaneous Pairwise Marginal Independence (SPMI)", "\n", "\n")
  cat("Unadjusted Pearson Chi-Square Tests for Independence:", "\n")
  cat("X^2_S =", X.sq.S, "\n")
  cat("X^2_S.ij =", "\n")
  print.default(X.sq.S.ij)
  cat("\n")
  if (any(type == "boot")) {
    B.discard <- x$boot$B.discard
    B.use <- x$boot$B.use
    p.value.boot <- as.name(paste("=", round(x$boot$p.value.boot, 4)))
    if (x$boot$p.value.boot == 0) {
      p.value.boot <- as.name(paste("<", round(1 / B.use, 4)))
    }
    p.combo.prod <- as.name(paste("=", round(x$boot$p.combo.prod.boot, 4)))
    if (x$boot$p.combo.prod.boot == 0) {
      p.combo.prod <- as.name(paste("<", round(1 / B.use, 4)))
    }
    p.combo.min <- as.name(paste("=", round(x$boot$p.combo.min.boot, 4)))
    if (x$boot$p.combo.min.boot == 0) {
      p.combo.min <- as.name(paste("<", round(1 / B.use, 4)))
    }
    cat("Bootstrap Results:", "\n")
    if (B.discard > 0) {
      cat(B.discard, "resamples were removed from the analysis due to")
      cat(" not having all rows or columns represented in a 2x2 table", "\n")
    }
    cat("Final results based on", B.use, "resamples", "\n")
    cat("p.boot", p.value.boot, "\n")
    cat("p.combo.prod", p.combo.prod, "\n")
    cat("p.combo.min", p.combo.min, "\n", "\n")
  }
  if (any(type == "rs2")) {
    X.sq.S.rs2 <- round(x$rs2$X.sq.S.rs2, 2)
    df.rs2 <- round(x$rs2$df.rs2, 2)
    p.value.rs2 <- as.name(paste("=", round(x$rs2$p.value.rs2, 4)))
    if (round(x$rs2$p.value.rs2, 4) < .0001) {
      p.value.rs2 <- as.name(paste("<", .0001))
    }
    cat("Second-Order Rao-Scott Adjusted Results:", "\n")
    cat("X^2_S.adj =", X.sq.S.rs2, "\n")
    cat("df.adj =", df.rs2, "\n")
    cat("p.adj", p.value.rs2, "\n", "\n")
  }
  if (any(type == "bon")) {
    p.value.bon <- as.name(paste("=", round(x$bon$p.value.bon, 4)))
    if (round(x$bon$p.value.bon, 4) < .0001) {
      p.value.bon <- as.name(paste("<", .0001))
    }
    X.sq.S.ij.p.bon <- format(round(x$bon$X.sq.S.ij.p.bon, 4),
      digits = 4
    )
    rownames(X.sq.S.ij.p.bon) <- names(data)[1:I]
    colnames(X.sq.S.ij.p.bon) <- names(data)[(I + 1):(I + J)]
    if (summary.data) {
      rownames(X.sq.S.ij.p.bon) <- levels(data[, 1])
      colnames(X.sq.S.ij.p.bon) <- levels(data[, 2])
    }
    cat("Bonferroni Adjusted Results:", "\n")
    cat("p.adj", p.value.bon, "\n")
    cat("p.ij.adj =", "\n")
    print.default(X.sq.S.ij.p.bon, quote = FALSE)
    cat("\n")
  }
  invisible(x)
}
