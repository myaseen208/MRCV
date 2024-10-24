# SPMI.stat()
# A function that calculates X^2_S and X^2_S.ij to test for SPMI (2 MRCV case)
# Called by the general MI.stat() function

SPMI.stat <- function(data, I, J, summary.data, add.constant) {
  op <- options()
  on.exit(options(op))
  # For inputted data set that is a summary file
  if (summary.data) {
    data[, 1:2] <- lapply(data[, 1:2], factor)
    X.sq.S.ij <- matrix(data = NA, nrow = I, ncol = J)
    counter <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        table.11 <- data[((data[, 1] == levels(data[, 1])[i]) &
          (data[, 2] == levels(data[, 2])[j]) & (data[, 3] == 0) & (data[, 4] == 0)), 5]
        table.12 <- data[((data[, 1] == levels(data[, 1])[i]) &
          (data[, 2] == levels(data[, 2])[j]) & (data[, 3] == 0) & (data[, 4] == 1)), 5]
        table.21 <- data[((data[, 1] == levels(data[, 1])[i]) &
          (data[, 2] == levels(data[, 2])[j]) & (data[, 3] == 1) & (data[, 4] == 0)), 5]
        table.22 <- data[((data[, 1] == levels(data[, 1])[i]) &
          (data[, 2] == levels(data[, 2])[j]) & (data[, 3] == 1) & (data[, 4] == 1)), 5]
        n.table <- matrix(
          data = c(table.11, table.21, table.12, table.22),
          nrow = 2, ncol = 2
        )
        # Only calculate statistics for valid 2x2 tables (no 00 in row or col)
        if (min(c(rowSums(n.table), colSums(n.table))) > 0) {
          counter <- counter + 1
          # Add .5 to 0 cell counts
          n.table <- apply(
            X = n.table, MARGIN = c(1, 2), FUN = check.zero,
            add.constant = add.constant
          )
          ####################
          # options(warn = -1)
          suppressWarnings(X.sq.S.ij[i, j] <- chisq.test(n.table, correct = F)$statistic)
          ###################
          # options(warn = 0)
        }
      }
    }
    rownames(X.sq.S.ij) <- levels(data[, 1])
    colnames(X.sq.S.ij) <- levels(data[, 2])
    output <- list(
      X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins =
        counter
    )
  }

  # For inputted data set that is a raw file
  if (!summary.data) {
    X.sq.S.ij <- matrix(data = NA, nrow = I, ncol = J)
    counter <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        # Only calculate statistics for valid tables (correct dim = 2x2)
        if (sum(dim(table(data[, i], data[, (I + j)]))) == 4) {
          counter <- counter + 1
          n.table <- table(data[, i], data[, (I + j)])
          # Add .5 to 0 cell counts
          n.table <- apply(
            X = n.table, MARGIN = c(1, 2), FUN = check.zero,
            add.constant = add.constant
          )
          ##################
          # options(warn = -1)
          suppressWarnings(X.sq.S.ij[i, j] <- chisq.test(n.table, correct = F)$statistic)
          ###################
          # options(warn = 0)
        }
      }
    }
    rownames(X.sq.S.ij) <- names(data)[1:I]
    colnames(X.sq.S.ij) <- names(data)[(I + 1):(I + J)]
    output <- list(
      X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins =
        counter
    )
  }
  output
}
