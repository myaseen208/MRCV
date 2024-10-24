# SPMI.test()
# A function that performs bootstrap testing, the second-order Rao-Scott approach,
#   and Bonferroni adjustments to test for SPMI (2 MRCV case)
# Called by the general MI.test() function

SPMI.test <- function(data, I, J, type, B, B.max, summary.data, add.constant, plot.hist,
                      print.status) {
  n <- nrow(data)
  # Observed statistics
  observed <- MI.stat(
    data = data, I = I, J = J,
    summary.data = summary.data, add.constant = add.constant
  )
  p.value.ij <- 1 - pchisq(q = observed$X.sq.S.ij, df = 1)
  p.value.min <- min(p.value.ij)
  p.value.prod <- prod(p.value.ij)

  # Bootstrap
  if (type == "boot" | type == "all") {
    X.sq.S.star <- numeric(length(B.max))
    p.value.b.min <- numeric(length(B.max))
    p.value.b.prod <- numeric(length(B.max))
    X.sq.S.ij.star <- matrix(data = NA, nrow = I * J, ncol = B.max)
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
      # Resample W_s and Y_s independently
      W <- sample(x = 1:n, size = n, replace = TRUE)
      Y <- sample(x = 1:n, size = n, replace = TRUE)
      data.star <- cbind(data[W, 1:I], data[Y, (I + 1):(I + J)])
      stat.star <- MI.stat(
        data = data.star, I = I, J = J,
        summary.data = FALSE, add.constant = add.constant
      )
      discard <- discard + (stat.star$valid.margins < I * J) # Discard if any table < 2x2
      drop <- 1
      if (stat.star$valid.margins == I * J) {
        counter <- counter + 1
        X.sq.S.star[counter] <- stat.star$X.sq.S
        X.sq.S.ij.star[, counter] <- stat.star$X.sq.S.ij
        p.value.ij <- 1 - pchisq(q = stat.star$X.sq.S.ij, df = 1)
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
    W.counts <- as.data.frame(table(data[, 1:I])) # Get all possible combos of W's
    cols <- c(1:I) # Need to order W's in ascending order starting with W1
    W.counts.ord <- W.counts[do.call("order", as.data.frame(W.counts[, cols])), ]
    Y.counts <- as.data.frame(table(data[, (I + 1):(I + J)])) # All possible Y's
    cols <- c(1:J) # Need to order Y's in ascending order starting with Y1
    Y.counts.ord <- Y.counts[do.call("order", as.data.frame(Y.counts[, cols])), ]
    n.counts <- as.data.frame(table(data)) # Get all possible combos of W's and Y's
    cols <- c(1:(I + J)) # Need to order W's and Y's in ascending order W1, W2, ...
    n.counts.ord <- n.counts[do.call("order", as.data.frame(n.counts[, cols])), ]
    G <- t(data.matrix(W.counts.ord[, 1:I]) - 1) # rx2^r matrix
    H <- t(data.matrix(Y.counts.ord[, 1:J]) - 1) # cx2^c matrix
    tau <- n.counts.ord[, (I + J + 1)] / n # Vector of multinomial probabilities
    m.row <- G %*% W.counts.ord[, (I + 1)] # Vector of W marginal counts
    m.col <- H %*% Y.counts.ord[, (J + 1)] # Vector of Y marginal counts
    GH <- kronecker(G, H)
    m <- GH %*% n.counts.ord[, (I + J + 1)] # Marginal table counts (row 1, row 2, ...)
    pi <- m / n # pi_ij
    pi.row <- m.row / n # pi_i.
    pi.col <- m.col / n # pi_.j
    j.2r <- matrix(data = 1, nrow = 2^I, ncol = 1) # 2^r vector of 1's
    i.2r <- diag(2^I) # (2^r)x(2^r) identity matrix
    j.2c <- matrix(data = 1, nrow = 2^J, ncol = 1) # 2^c vector of 1's
    i.2c <- diag(2^J) # (2^c)x(2^c) identity matrix
    G.ij <- G %*% kronecker(i.2r, t(j.2c))
    H.ji <- H %*% kronecker(t(j.2r), i.2c)
    F <- GH - kronecker(pi.row, H.ji) - kronecker(G.ij, pi.col)
    mult.cov <- diag(tau) - tcrossprod(tau)
    sigma <- F %*% tcrossprod(mult.cov, F)
    D <- diag(as.vector(kronecker(pi.row, pi.col) * kronecker(1 - pi.row, 1 - pi.col)))
    Di.sigma <- diag(1 / diag(D), dim(D)) %*% sigma
    Di.sigma.eigen <- Re(eigen(Di.sigma)$values) # Only use real part of eigenvalues
    sum.Di.sigma.eigen.sq <- sum(Di.sigma.eigen^2)
    X.sq.S.rs2 <- I * J * observed$X.sq.S / sum.Di.sigma.eigen.sq
    df.rs2 <- I^2 * J^2 / sum.Di.sigma.eigen.sq
    X.sq.S.p.value.rs2 <- 1 - pchisq(q = X.sq.S.rs2, df = df.rs2)
    output <- list(general = list(
      data = data, I = I, J = J, summary.data =
        summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij =
        observed$X.sq.S.ij
    ), rs2 = list(
      X.sq.S.rs2 = X.sq.S.rs2, df.rs2 =
        df.rs2, p.value.rs2 = X.sq.S.p.value.rs2
    ))
    output.rs2 <- list(
      X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = df.rs2,
      p.value.rs2 = X.sq.S.p.value.rs2
    )
  }

  # Bonferroni
  if (type == "bon" | type == "all") {
    p.value.bon <- min(p.value.min * I * J, 1) # p-value should not be greater than 1
    X.sq.S.ij.p.bon <- apply(X = I * J * (1 - pchisq(
      q = observed$X.sq.S.ij,
      df = 1
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
  class(output) <- "SPMI"
  output
}
