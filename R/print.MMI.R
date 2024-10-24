#' @name    print.MMI
#' @aliases print.MMI
#' @title Control Printed Display of Objects of Class "MMI"
#' @description
#' The \code{print.MMI} method function controls the printed display of objects
#' of class \code{MMI}.
#'
#'
#' @param x An object of class \code{"MMI"} produced by the
#' \code{\link{MI.test}} function.
#' @param ...  Additional arguments passed to or from other methods.
#' @return A print out of selected results from \code{\link{MI.test}}.
#' @export
#' @rdname print
print.MMI <- function(x, ...) {
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
  srcv <- ifelse(I == 1, 1, I + 1)
  mrcv <- ifelse(I == 1, 2, 1)
  c <- ifelse(I == 1, J, I)
  if (summary.data) {
    r <- length(levels(as.factor(data[, 1])))
    rownames(X.sq.S.ij) <- ""
    colnames(X.sq.S.ij) <- levels(data[, 2])
  }
  if (!summary.data) {
    r <- length(levels(as.factor(data[, srcv])))
    rownames(X.sq.S.ij) <- ""
    colnames(X.sq.S.ij) <- names(data)[mrcv:(mrcv + c - 1)]
  }
  cat("Test for Multiple Marginal Independence (MMI)", "\n", "\n")
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
      cat(" not having all rows or columns represented in an rx2 table", "\n")
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
    if (summary.data) {
      rownames(X.sq.S.ij.p.bon) <- ""
      colnames(X.sq.S.ij.p.bon) <- levels(data[, 2])
    }
    if (!summary.data) {
      rownames(X.sq.S.ij.p.bon) <- ""
      colnames(X.sq.S.ij.p.bon) <- names(data)[mrcv:(mrcv + c - 1)]
    }
    cat("Bonferroni Adjusted Results:", "\n")
    cat("p.adj", p.value.bon, "\n")
    cat("p.ij.adj =", "\n")
    print.default(X.sq.S.ij.p.bon, quote = FALSE)
    cat("\n")
  }
  invisible(x)
}
