# MMI.test()
# A function that performs bootstrap testing, the second-order Rao-Scott approach,
#   and Bonferroni adjustments to test for MMI (1 MRCV case)
# Called by the general MI.test() function

MMI.test <- function(data, I, J, type, B, B.max, summary.data, add.constant, plot.hist,
                     print.status) {
  n <- nrow(data)
  srcv <- ifelse(I == 1, 1, I + 1)
  mrcv <- ifelse(I == 1, 2, 1)
  c <- ifelse(I == 1, J, I)
  if (summary.data) {
    r <- length(levels(as.factor(data[, 1])))
  }
  if (!summary.data) {
    r <- length(levels(as.factor(data[, srcv])))
  }
  # Observed statistics
  observed <- MI.stat(
    data = data, I = I, J = J,
    summary.data = summary.data, add.constant = add.constant
  )
  p.value.ij <- 1 - pchisq(q = observed$X.sq.S.ij, df = (r - 1))
  p.value.min <- min(p.value.ij)
  p.value.prod <- prod(p.value.ij)

  # Bootstrap
  if (type == "boot" | type == "all") {
    X.sq.S.star <- numeric(length(B.max))
    p.value.b.min <- numeric(length(B.max))
    p.value.b.prod <- numeric(length(B.max))
    X.sq.S.ij.star <- matrix(data = NA, nrow = c, ncol = B.max)
    discard <- 0 # Count number of resamples with incorrect dimensions
    counter <- 0
    b <- 0
    filled <- 0
    # create progress bar
    if (print.status) {
      cat("Bootstrap Progress:", "\n")
      pb <- txtProgressBar(min = 0, max = B.max, style = 3)
    }
    while (((counter >= B) + (b >= B.max)) < 1) {
      b <- b + 1
      n_iplus <- as.matrix(table(data[, srcv]))
      W <- rep(levels(as.factor(data[, srcv])), times = n_iplus)
      Y <- sample(x = 1:n, size = n, replace = TRUE)
      data.star <- cbind(W, data[Y, mrcv:(mrcv + c - 1)])
      stat.star <- MI.stat(
        data = data.star, I = 1, J = c,
        summary.data = FALSE, add.constant = add.constant
      )
      discard <- discard + (stat.star$valid.margins < c) # Discard if any table < rx2
      drop <- 1
      if (stat.star$valid.margins == c) {
        counter <- counter + 1
        X.sq.S.star[counter] <- stat.star$X.sq.S
        X.sq.S.ij.star[, counter] <- stat.star$X.sq.S.ij
        p.value.ij <- 1 - pchisq(q = stat.star$X.sq.S.ij, df = (r - 1))
        p.value.b.min[counter] <- min(p.value.ij)
        p.value.b.prod[counter] <- prod(p.value.ij)
        drop <- 0
      }
      # update progress bar
      if (B == B.max) {
        if (print.status) {
          setTxtProgressBar(pb, b)
        }
      }
      if (B != B.max) {
        if (print.status) {
          left <- B.max - filled
          expect <- (1 - drop) * min(left, (B - counter)) + drop * max(left, (B.max - b))
          filled <- filled + left / expect
          setTxtProgressBar(pb, filled)
        }
      }
    }
    if (print.status) {
      setTxtProgressBar(pb, B.max)
      close(pb)
    }
    B.use <- min(B, (B.max - discard)) # Only use desired number of resamples
    X.sq.S.star <- X.sq.S.star[1:B.use]
    X.sq.S.ij.star <- X.sq.S.ij.star[, 1:B.use]
    p.value.b.min <- p.value.b.min[1:B.use]
    p.value.b.prod <- p.value.b.prod[1:B.use]
    p.value.boot <- mean(X.sq.S.star >= observed$X.sq.S)
    p.combo.min <- list(
      p = p.value.min, p.star = p.value.b.min,
      overall.p = mean(p.value.b.min <= p.value.min)
    )
    p.combo.prod <- list(
      p = p.value.prod, p.star = p.value.b.prod,
      overall.p = mean(p.value.b.prod <= p.value.prod)
    )
    # Histograms
    if (plot.hist) {
      layout.m <- matrix(c(1, 0, 1, 3, 2, 3, 2, 0), nrow = 2, ncol = 4)
      layout(layout.m)
      hist(
        x = X.sq.S.star, main = "Test using sum statistic", xlab =
          expression(X[S]^{
            "2*"
          }), xlim = c(
          min(X.sq.S.star, observed$X.sq.S),
          max(X.sq.S.star, observed$X.sq.S)
        )
      )
      abline(v = observed$X.sq.S, col = "darkgreen", lwd = 5)
      hist(
        x = p.value.b.prod, main = "Test using product of p-values", xlab =
          expression(tilde(p)[prod]^{
            "*"
          }), xlim = c(min(
          p.value.b.prod,
          p.value.prod
        ), max(p.value.b.prod, p.value.prod))
      )
      abline(v = p.value.prod, col = "darkgreen", lwd = 5)
      hist(
        x = p.value.b.min, main = "Test using minimum of p-values", xlab =
          expression(tilde(p)[min]^{
            "*"
          }), xlim = c(min(
          p.value.b.min,
          p.value.min
        ), max(p.value.b.min, p.value.min))
      )
      abline(v = p.value.min, col = "darkgreen", lwd = 5)
    }
    output <- list(general = list(
      data = data, I = I, J = J, summary.data =
        summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij =
        observed$X.sq.S.ij
    ), boot = list(
      B.use = B.use,
      B.discard = discard, p.value.boot = p.value.boot,
      p.combo.min.boot = p.combo.min$overall.p, p.combo.prod.boot =
        p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, X.sq.S.ij.star =
        X.sq.S.ij.star, p.combo.min.star = p.combo.min$p.star,
      p.combo.prod.star = p.combo.prod$p.star
    ))
    output.boot <- list(
      B.use = B.use, B.discard = discard,
      p.value.boot = p.value.boot, p.combo.min.boot =
        p.combo.min$overall.p, p.combo.prod.boot =
        p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star,
      X.sq.S.ij.star = X.sq.S.ij.star, p.combo.min.star =
        p.combo.min$p.star, p.combo.prod.star = p.combo.prod$p.star
    )
  }

  # Second-order Rao-Scott adjustment
  if (type == "rs2" | type == "all") {
    # Appendix A calculations
    n.counts <- as.data.frame(table(data))
    cols <- c(srcv, (mrcv:(mrcv + c - 1)))
    n.counts.ord <- n.counts[do.call("order", as.data.frame(n.counts[, cols])), ]
    Y.counts <- as.data.frame(table(data[, (mrcv:(mrcv + c - 1))]))
    cols <- c(1:c)
    Y.counts.ord <- Y.counts[do.call("order", as.data.frame(Y.counts[, cols])), ]
    n_iplus <- as.matrix(table(data[, srcv]))
    tau <- n.counts.ord[, ncol(n.counts.ord)] / rep(n_iplus, each = 2^c)
    G.tilde <- t(data.matrix(Y.counts.ord[, 1:c]) - 1)
    I.r <- diag(r)
    G <- kronecker(I.r, G.tilde)
    pi <- G %*% tau
    m <- pi * rep(n_iplus, each = c)
    a.i <- n_iplus / n
    pi.not.j <- (1 / n) * kronecker(matrix(rep(1, r), 1, r), diag(c)) %*% m
    j.r <- matrix(1, r, 1)
    pi.not <- kronecker(j.r, pi.not.j)
    I.rc <- diag(r * c)
    I.c <- diag(c)
    J.rr <- matrix(1, r, r)
    A <- diag(as.vector(a.i))
    H <- I.rc - kronecker(J.rr %*% A, I.c)
    D <- kronecker(diag(as.vector(n / n_iplus)), diag(as.vector(pi.not.j * (1 - pi.not.j))))
    V <- matrix(data = 0, nrow = r * 2^c, ncol = r * 2^c)
    for (i in 1:r) {
      V[((i - 1) * 2^c + 1):((i - 1) * 2^c + 2^c), ((i - 1) * 2^c + 1):((i - 1) * 2^c + 2^c)] <- ((1 / a.i[i])
      * (diag(as.vector(tau[((i - 1) * 2^c + 1):((i - 1) * 2^c + 2^c)]))
        - tcrossprod(tau[((i - 1) * 2^c + 1):((i - 1) * 2^c + 2^c)])))
    }
    Di.HGVGH <- diag(1 / diag(D), dim(D)) %*% H %*% G %*% tcrossprod(tcrossprod(V, G), H)
    Di.HGVGH.eigen <- Re(eigen(Di.HGVGH)$values)
    sum.Di.HGVGH.eigen.sq <- sum(Di.HGVGH.eigen^2)
    X.sq.S.rs2 <- (r - 1) * c * observed$X.sq.S / sum.Di.HGVGH.eigen.sq
    df.rs2 <- (r - 1)^2 * c^2 / sum.Di.HGVGH.eigen.sq
    X.sq.S.p.value.rs2 <- 1 - pchisq(q = X.sq.S.rs2, df = df.rs2)
    output <- list(general = list(
      data = data, I = I, J = J, summary.data =
        summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij =
        observed$X.sq.S.ij
    ), rs2 = list(
      X.sq.S.rs2 = X.sq.S.rs2,
      df.rs2 = df.rs2, p.value.rs2 = X.sq.S.p.value.rs2
    ))
    output.rs2 <- list(
      X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = df.rs2,
      p.value.rs2 = X.sq.S.p.value.rs2
    )
  }

  # Bonferroni
  if (type == "bon" | type == "all") {
    p.value.bon <- min(p.value.min * c, 1)
    X.sq.S.ij.p.bon <- apply(X = c * (1 - pchisq(
      q = observed$X.sq.S.ij,
      df = (r - 1)
    )), MARGIN = c(1, 2), FUN = check.min)
    output <- list(general = list(
      data = data, I = I, J = J, summary.data =
        summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij =
        observed$X.sq.S.ij
    ), bon = list(
      p.value.bon = p.value.bon,
      X.sq.S.ij.p.bon = X.sq.S.ij.p.bon
    ))
    output.bon <- list(p.value.bon = p.value.bon, X.sq.S.ij.p.bon = X.sq.S.ij.p.bon)
  }
  if (type == "all") {
    output <- list(
      general = list(
        data = data, I = I, J = J, summary.data =
          summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij =
          observed$X.sq.S.ij
      ), boot = output.boot, rs2 = output.rs2,
      bon = output.bon
    )
  }
  class(output) <- "MMI"
  output
}
