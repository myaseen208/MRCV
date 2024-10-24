# check.margins()
# A function that checks whether resamples have all positive or negative responses
# for an item

check.margins <- function(data, I, J, K, nvars, model.vars, item.names) {
  nrows <- (2^nvars) * I * J * max(1, K)
  model.data <- data.frame(matrix(data = NA, nrow = nrows, ncol = (2 * nvars + 1)))
  model.data[, 1:(2 * nvars)] <- model.vars
  model.data[, (2 * nvars + 1)] <- data
  if (nvars == 2) {
    pos.count <- numeric(I + J)
    neg.count <- numeric(I + J)
    for (i in 1:I) {
      pos.count[i] <- sum(model.data[((model.data[, 1] == item.names[i]) &
        (model.data[, 3] == 1)), 5])
      neg.count[i] <- sum(model.data[((model.data[, 1] == item.names[i]) &
        (model.data[, 3] == 0)), 5])
    }
    for (j in 1:J) {
      pos.count[(I + j)] <- sum(model.data[((model.data[, 2] == item.names[(I + j)]) &
        (model.data[, 4] == 1)), 5])
      neg.count[(I + j)] <- sum(model.data[((model.data[, 2] == item.names[(I + j)]) &
        (model.data[, 4] == 0)), 5])
    }
  }
  if (nvars == 3) {
    pos.count <- numeric(I + J + K)
    neg.count <- numeric(I + J + K)
    for (i in 1:I) {
      pos.count[i] <- sum(model.data[((model.data[, 1] == item.names[i]) &
        (model.data[, 4] == 1)), 7])
      neg.count[i] <- sum(model.data[((model.data[, 1] == item.names[i]) &
        (model.data[, 4] == 0)), 7])
    }
    for (j in 1:J) {
      pos.count[(I + j)] <- sum(model.data[((model.data[, 2] == item.names[(I + j)]) &
        (model.data[, 5] == 1)), 7])
      neg.count[(I + j)] <- sum(model.data[((model.data[, 2] == item.names[(I + j)]) &
        (model.data[, 5] == 0)), 7])
    }
    for (k in 1:K) {
      pos.count[(I + J + k)] <- sum(model.data[((model.data[, 3] == item.names[(I + J + k)]) &
        (model.data[, 6] == 1)), 7])
      neg.count[(I + J + k)] <- sum(model.data[((model.data[, 3] == item.names[(I + J + k)]) &
        (model.data[, 6] == 0)), 7])
    }
  }
  ((min(pos.count) > 0) & (min(neg.count) > 0))
}
