#' @name    print.anova.genloglin
#' @aliases print.anova.genloglin
#' @title Control Printed Display of MRCV ANOVA Modeling Objects
#' @description
#' A function used to control the display of output provided by \code{anova.genloglin}.
# 
#' @export
#' @rdname print
print.anova.genloglin <- function(x, ...) {
  op <- options()
  on.exit(options(op))
  options(scipen = 10)
  type <- names(x)
  model <- x$original.arg$model
  # if (class(model)=="formula") {
  if (inherits(model, "formula")) {
    model <- Reduce(paste, deparse(model, width.cutoff = 499L))
  }
  model.HA <- x$original.arg$model.HA
  # if (class(model.HA)=="formula") {
  if (inherits(model.HA, "formula")) {
    model.HA <- Reduce(paste, deparse(model.HA, width.cutoff = 499L))
  }
  gof <- x$original.arg$gof
  Pearson.chisq <- round(x$test.statistics$Pearson.chisq, 2)
  lrt <- round(x$test.statistics$lrt, 2)
  if (model == "saturated") {
    chisq.p.rs <- x$rs.results$Pearson.chisq.rs$p.value
    lrt.p.rs <- x$rs.results$lrt.rs$p.value
    chisq.p.boot <- x$boot.results$Pearson.chisq.boot$p.value
    lrt.p.boot <- x$boot.results$lrt.boot$p.value
    cat("\n")
    cat("Model comparison statistics for", "\n")
    cat("H0 =", model, "\n")
    cat("HA =", model.HA, "\n", "\n")
    cat("Pearson chi-square statistic =", Pearson.chisq, "\n")
    cat("LRT statistic =", lrt, "\n", "\n")
    cat("Second-Order Rao-Scott Adjusted Results:", "\n")
    cat("Pearson chi-square p-value =", chisq.p.rs, "\n")
    cat("LRT p-value =", lrt.p.rs, "\n", "\n")
    cat("Bootstrap Results:", "\n")
    cat("Pearson chi-square p-value =", chisq.p.boot, "\n")
    cat("LRT p-value =", lrt.p.boot, "\n", "\n")
  }
  if (model != "saturated") {
    cat("\n")
    cat("Model comparison statistics for", "\n")
    cat("H0 =", model, "\n")
    cat("HA =", model.HA, "\n", "\n")
    cat("Pearson chi-square statistic =", Pearson.chisq, "\n")
    cat("LRT statistic =", lrt, "\n")
    if (any(type == "rs.results")) {
      Pearson.chisq.rs <- round(x$rs.results$Pearson.chisq.rs$Pearson.chisq.rs, 2)
      df.rs <- round(x$rs.results$Pearson.chisq.rs$df, 2)
      chisq.p.rs <- as.name(paste("=", round(x$rs.results$Pearson.chisq.rs$p.value, 4)))
      if (round(x$rs.results$Pearson.chisq.rs$p.value, 4) < .0001) {
        chisq.p.rs <- as.name(paste("<", .0001))
      }
      lrt.rs <- round(x$rs.results$lrt.rs$lrt.rs, 2)
      lrt.p.rs <- as.name(paste("=", round(x$rs.results$lrt.rs$p.value, 4)))
      if (round(x$rs.results$lrt.rs$p.value, 4) < .0001) {
        lrt.p.rs <- as.name(paste("<", .0001))
      }
      cat("\n")
      cat("Second-Order Rao-Scott Adjusted Results:", "\n")
      cat("Rao-Scott Pearson chi-square statistic =", paste(Pearson.chisq.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", chisq.p.rs, "\n")
      cat("Rao-Scott LRT statistic =", paste(lrt.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", lrt.p.rs, "\n")
      if (any(type == "boot.results")) {
        B.use <- x$boot.results$B.use
        B.discard <- x$boot.results$B.discard
        p.chisq.boot <- as.name(paste("=", round(x$boot.results$p.chisq.boot, 4)))
        if (x$boot.results$p.chisq.boot == 0) {
          p.chisq.boot <- as.name(paste("<", round(1 / B.use, 4)))
        }
        p.lrt.boot <- as.name(paste("=", round(x$boot.results$p.lrt.boot, 4)))
        if (x$boot.results$p.lrt.boot == 0) {
          p.lrt.boot <- as.name(paste("<", round(1 / B.use, 4)))
        }
        cat("\n")
        cat("Bootstrap Results:", "\n")
        if (B.discard > 0) {
          cat(B.discard, "resamples were removed from the analysis due to")
          cat(" not having all rows or columns represented in a 2x2 table", "\n")
        }
        cat("Final results based on", B.use, "resamples", "\n")
        cat("Pearson chi-square p-value", p.chisq.boot, "\n")
        cat("LRT p-value", p.lrt.boot, "\n")
      }
    }
    if (gof & (model.HA != "saturated")) {
      cat("\n")
      cat("-------------------------------------------------------------------------------------")
      cat("\n", "\n")
      Pearson.chisq.gof <- round(x$test.statistics$Pearson.chisq.gof, 2)
      lrt.gof <- round(x$test.statistics$lrt.gof, 2)
      cat("Goodness of fit statistics for", "\n")
      cat("H0 =", model, "\n", "\n")
      cat("Pearson chi-square GOF statistic =", Pearson.chisq.gof, "\n")
      cat("LRT GOF statistic =", lrt.gof, "\n")
      if (any(type == "rs.results")) {
        Pearson.chisq.gof.rs <- (round(x$rs.results$Pearson.chisq.gof.rs$Pearson.chisq.gof.rs, 2))
        df.rs <- round(x$rs.results$Pearson.chisq.gof.rs$df, 2)
        chisq.gof.p.rs <- (as.name(paste("=", round(x$rs.results$Pearson.chisq.gof.rs$p.value, 4))))
        if (round(x$rs.results$Pearson.chisq.gof.rs$p.value, 4) < .0001) {
          chisq.gof.p.rs <- as.name(paste("<", .0001))
        }
        lrt.gof.rs <- round(x$rs.results$lrt.gof.rs$lrt.gof.rs, 2)
        lrt.gof.p.rs <- as.name(paste("=", round(x$rs.results$lrt.gof.rs$p.value, 4)))
        if (round(x$rs.results$lrt.gof.rs$p.value, 4) < .0001) {
          lrt.gof.p.rs <- as.name(paste("<", .0001))
        }
        cat("\n")
        cat("Second-Order Rao-Scott Adjusted Results:", "\n")
        cat("Rao-Scott Pearson chi-square GOF statistic =", paste(Pearson.chisq.gof.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", chisq.gof.p.rs, "\n")
        cat("Rao-Scott LRT GOF statistic =", paste(lrt.gof.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", lrt.gof.p.rs, "\n")
      }
      if (any(type == "boot.results")) {
        p.chisq.gof.boot <- as.name(paste("=", round(x$boot.results$p.chisq.gof.boot, 4)))
        if (x$boot.results$p.chisq.gof.boot == 0) {
          p.chisq.gof.boot <- as.name(paste("<", round(1 / B.use, 4)))
        }
        p.lrt.gof.boot <- as.name(paste("=", round(x$boot.results$p.lrt.gof.boot, 4)))
        if (x$boot.results$p.lrt.gof.boot == 0) {
          p.lrt.gof.boot <- as.name(paste("<", round(1 / B.use, 4)))
        }
        cat("\n")
        cat("Bootstrap Results:", "\n")
        cat("Pearson chi-square GOF p-value", p.chisq.gof.boot, "\n")
        cat("LRT GOF p-value", p.lrt.gof.boot, "\n")
      }
    }
    cat("\n")
  }
  invisible(x)
}
