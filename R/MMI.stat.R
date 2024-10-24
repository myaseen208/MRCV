# MMI.stat()
# A function that calculates X^2_S and X^2_S.ij to test for MMI (1 MRCV case)
# Called by the general MI.stat() function

MMI.stat <- function(data, I, J, summary.data, add.constant) {
  op <- options()
  on.exit(options(op))
  srcv <- ifelse(I == 1, 1, I + 1)
  mrcv <- ifelse(I == 1, 2, 1)
  c <- ifelse(I == 1, J, I)
  # For inputted data set that is a summary file
  if (summary.data) {
    r <- length(levels(as.factor(data[, 1])))
    data[, 1:2] <- lapply(data[, 1:2], factor)
    X.sq.S.ij <- matrix(data = NA, nrow = 1, ncol = c)
    counter <- 0
    for (i in 1:c) {
      table.11 <- data[((data[, 2] == levels(data[, 2])[i]) & (data[, 3] == 0)), 4]
      table.12 <- data[((data[, 2] == levels(data[, 2])[i]) & (data[, 3] == 1)), 4]
      n.table <- matrix(data = c(table.11, table.12), nrow = r, ncol = 2)
      # Only calculate statistics for valid rx2 tables
      if (min(c(rowSums(n.table), colSums(n.table))) > 0) {
        counter <- counter + 1
        # Add .5 to 0 cell counts
        n.table <- apply(
          X = n.table, MARGIN = c(1, 2), FUN = check.zero,
          add.constant = add.constant
        )
        ####################
        # options(warn = -1)
        suppressWarnings(X.sq.S.ij[1, i] <- chisq.test(n.table, correct = F)$statistic)
        ###################
        # options(warn = 0)
      }
    }
    colnames(X.sq.S.ij) <- levels(data[, 2])
    output <- list(
      X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins =
        counter
    )
  }

  # For inputted data set that is a raw file
  if (!summary.data) {
    r <- length(levels(as.factor(data[, srcv])))
    X.sq.S.ij <- matrix(data = NA, nrow = 1, ncol = c)
    counter <- 0
    for (i in 1:c) {
      # Only calculate statistics for valid tables (correct dim = rx2)
      if (sum(dim(table(data[, srcv], data[, (mrcv + i - 1)]))) == (r + 2)) {
        counter <- counter + 1
        n.table <- table(data[, srcv], data[, (mrcv + i - 1)])
        # Add .5 to 0 cell counts
        n.table <- apply(
          X = n.table, MARGIN = c(1, 2), FUN = check.zero,
          add.constant = add.constant
        )
        ####################
        # options(warn = -1)
        suppressWarnings(X.sq.S.ij[1, i] <- chisq.test(n.table, correct = F)$statistic)
        ###################
        # options(warn = 0)
      }
    }
    colnames(X.sq.S.ij) <- names(data)[mrcv:(mrcv + c - 1)]
    output <- list(
      X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins =
        counter
    )
  }
  output
}
