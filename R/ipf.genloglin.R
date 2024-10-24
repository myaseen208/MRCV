# ipf.genloglin()
# A function that uses the iterative proportional fitting algorithm to obtain
#   multinomial probabilities under the null hypothesis

ipf.genloglin <- function(data, I, J, K, nvars, p, p.theta.2, p.theta.3, x.theta.2,
                          x.theta.3) {
  p.theta.2.new <- p.theta.2
  p.theta.3.new <- p.theta.3
  if (nvars == 2) {
    for (i in 1:(I + J - 1)) {
      for (j in (i + 1):(I + J)) {
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 1)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 1)),
          ncol(p)
        ])

        p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 1)), 5])
        p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 1)), 5])

        p.theta.2 <- p.theta.2.new
      }
    }
  }
  if (nvars == 3) {
    for (i in 1:I) {
      for (j in (I + 1):(I + J)) {
        for (k in (I + J + 1):(I + J + K)) {
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 0) &
            (p.theta.3[, 6] == 0)), 7] <- sum(p[((p[, i] == 0) & (p[, j] == 0) & (p[, k] == 0)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 0) &
            (p.theta.3[, 6] == 1)), 7] <- sum(p[((p[, i] == 0) & (p[, j] == 0) & (p[, k] == 1)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 1) &
            (p.theta.3[, 6] == 0)), 7] <- sum(p[((p[, i] == 0) & (p[, j] == 1) & (p[, k] == 0)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 1) &
            (p.theta.3[, 6] == 1)), 7] <- sum(p[((p[, i] == 0) & (p[, j] == 1) & (p[, k] == 1)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 0) &
            (p.theta.3[, 6] == 0)), 7] <- sum(p[((p[, i] == 1) & (p[, j] == 0) & (p[, k] == 0)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 0) &
            (p.theta.3[, 6] == 1)), 7] <- sum(p[((p[, i] == 1) & (p[, j] == 0) & (p[, k] == 1)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 1) &
            (p.theta.3[, 6] == 0)), 7] <- sum(p[((p[, i] == 1) & (p[, j] == 1) & (p[, k] == 0)), ncol(p)])
          p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
            (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 1) &
            (p.theta.3[, 6] == 1)), 7] <- sum(p[((p[, i] == 1) & (p[, j] == 1) & (p[, k] == 1)), ncol(p)])

          p[((p[, i] == 0) & (p[, j] == 0) & (p[, k] == 0)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 0) & (p[, k] == 0)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 0) &
              (p.theta.3[, 6] == 0)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 0) & (x.theta.3[, 5] == 0) & (x.theta.3[, 6] == 0)), 7])
          p[((p[, i] == 0) & (p[, j] == 0) & (p[, k] == 1)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 0) & (p[, k] == 1)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 0) &
              (p.theta.3[, 6] == 1)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 0) & (x.theta.3[, 5] == 0) & (x.theta.3[, 6] == 1)), 7])
          p[((p[, i] == 0) & (p[, j] == 1) & (p[, k] == 0)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 1) & (p[, k] == 0)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 1) &
              (p.theta.3[, 6] == 0)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 0) & (x.theta.3[, 5] == 1) & (x.theta.3[, 6] == 0)), 7])
          p[((p[, i] == 0) & (p[, j] == 1) & (p[, k] == 1)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 1) & (p[, k] == 1)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 0) & (p.theta.3[, 5] == 1) &
              (p.theta.3[, 6] == 1)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 0) & (x.theta.3[, 5] == 1) & (x.theta.3[, 6] == 1)), 7])
          p[((p[, i] == 1) & (p[, j] == 0) & (p[, k] == 0)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 0) & (p[, k] == 0)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 0) &
              (p.theta.3[, 6] == 0)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 1) & (x.theta.3[, 5] == 0) & (x.theta.3[, 6] == 0)), 7])
          p[((p[, i] == 1) & (p[, j] == 0) & (p[, k] == 1)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 0) & (p[, k] == 1)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 0) &
              (p.theta.3[, 6] == 1)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 1) & (x.theta.3[, 5] == 0) & (x.theta.3[, 6] == 1)), 7])
          p[((p[, i] == 1) & (p[, j] == 1) & (p[, k] == 0)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 1) & (p[, k] == 0)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 1) &
              (p.theta.3[, 6] == 0)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 1) & (x.theta.3[, 5] == 1) & (x.theta.3[, 6] == 0)), 7])
          p[((p[, i] == 1) & (p[, j] == 1) & (p[, k] == 1)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 1) & (p[, k] == 1)), ncol(p)]
          / p.theta.3[((p.theta.3[, 1] == names(data)[i]) & (p.theta.3[, 2] == names(data)[j]) &
              (p.theta.3[, 3] == names(data)[k]) & (p.theta.3[, 4] == 1) & (p.theta.3[, 5] == 1) &
              (p.theta.3[, 6] == 1)), 7] * x.theta.3[((x.theta.3[, 1] == names(data)[i]) &
              (x.theta.3[, 2] == names(data)[j]) & (x.theta.3[, 3] == names(data)[k]) &
              (x.theta.3[, 4] == 1) & (x.theta.3[, 5] == 1) & (x.theta.3[, 6] == 1)), 7])

          p.theta.3 <- p.theta.3.new
        }
      }
    }
    for (i in 1:(I - 1)) {
      for (j in (i + 1):I) {
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 1)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 1)),
          ncol(p)
        ])

        p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 1)), 5])
        p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 1)), 5])

        p.theta.2 <- p.theta.2.new
      }
    }
    for (i in (I + 1):(I + J - 1)) {
      for (j in (i + 1):(I + J)) {
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 1)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 1)),
          ncol(p)
        ])

        p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 1)), 5])
        p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 1)), 5])

        p.theta.2 <- p.theta.2.new
      }
    }
    for (i in (I + J + 1):(I + J + K - 1)) {
      for (j in (i + 1):(I + J + K)) {
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 0) & (p[, j] == 1)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 0)),
          ncol(p)
        ])
        p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
          (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5] <- sum(p[
          ((p[, i] == 1) & (p[, j] == 1)),
          ncol(p)
        ])

        p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 0) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 0) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 0) & (x.theta.2[, 4] == 1)), 5])
        p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 0)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 0)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 0)), 5])
        p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)] <- (p[((p[, i] == 1) & (p[, j] == 1)), ncol(p)]
        / p.theta.2[((p.theta.2[, 1] == names(data)[i]) & (p.theta.2[, 2] == names(data)[j]) &
            (p.theta.2[, 3] == 1) & (p.theta.2[, 4] == 1)), 5]
          * x.theta.2[((x.theta.2[, 1] == names(data)[i]) & (x.theta.2[, 2] == names(data)[j]) &
            (x.theta.2[, 3] == 1) & (x.theta.2[, 4] == 1)), 5])

        p.theta.2 <- p.theta.2.new
      }
    }
  }
  p
}
