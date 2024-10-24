# data.format()
# A function that reformats the raw data or bootstrap resample data into model form

data.format <- function(data, I, J, K, nvars, add.constant = .5, model.vars = NULL,
                        predict.func = FALSE) {
  nrows <- (2^nvars) * I * J * max(1, K)
  if (predict.func) {
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter <- MRCV_globals$pb.counter + 1.7
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
  }
  if (is.null(model.vars)) {
    model.data.unsorted <- data.frame(matrix(data = NA, nrow = nrows, ncol = (2 * nvars + 1)))
    counter <- 0
    if (nvars == 2) {
      for (i in 1:I) {
        for (j in 1:J) {
          counter <- counter + 2^nvars
          table.count <- table(data[, i], data[, (I + j)])
          model.data.unsorted[(counter - 3), ] <- c(
            names(data)[i], names(data)[(I + j)], 0, 0,
            as.numeric(table.count[1, 1])
          )
          model.data.unsorted[(counter - 2), ] <- c(
            names(data)[i], names(data)[(I + j)], 0, 1,
            as.numeric(table.count[1, 2])
          )
          model.data.unsorted[(counter - 1), ] <- c(
            names(data)[i], names(data)[(I + j)], 1, 0,
            as.numeric(table.count[2, 1])
          )
          model.data.unsorted[(counter), ] <- c(
            names(data)[i], names(data)[(I + j)], 1, 1,
            as.numeric(table.count[2, 2])
          )
        }
      }
      colnames(model.data.unsorted) <- c("W", "Y", "wi", "yj", "count")
    }
    if (nvars == 3) {
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            counter <- counter + 2^nvars
            table.count <- table(data[, i], data[, (I + j)], data[, (I + J + k)])
            model.data.unsorted[(counter - 7), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 0, 0, 0,
              as.numeric(table.count[1, 1, 1])
            )
            model.data.unsorted[(counter - 6), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 0, 0, 1,
              as.numeric(table.count[1, 1, 2])
            )
            model.data.unsorted[(counter - 5), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 0, 1, 0,
              as.numeric(table.count[1, 2, 1])
            )
            model.data.unsorted[(counter - 4), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 0, 1, 1,
              as.numeric(table.count[1, 2, 2])
            )
            model.data.unsorted[(counter - 3), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 1, 0, 0,
              as.numeric(table.count[2, 1, 1])
            )
            model.data.unsorted[(counter - 2), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 1, 0, 1,
              as.numeric(table.count[2, 1, 2])
            )
            model.data.unsorted[(counter - 1), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 1, 1, 0,
              as.numeric(table.count[2, 2, 1])
            )
            model.data.unsorted[(counter), ] <- c(
              names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)], 1, 1, 1,
              as.numeric(table.count[2, 2, 2])
            )
          }
        }
      }
      colnames(model.data.unsorted) <- c("W", "Y", "Z", "wi", "yj", "zk", "count")
    }
    model.data.unsorted[, 1:nvars] <- lapply(model.data.unsorted[, 1:nvars], as.factor)
    model.data.unsorted[, ((nvars + 1):(ncol(model.data.unsorted)))] <- (lapply(model.data.unsorted[, ((nvars + 1):(ncol(model.data.unsorted)))], as.numeric))
    model.data.unsorted[, ncol(model.data.unsorted)] <- (apply(
      X = as.matrix(model.data.unsorted[, ncol(model.data.unsorted)]), MARGIN = 1,
      FUN = check.zero, add.constant = add.constant
    ))
  }
  if (!is.null(model.vars)) {
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter <- MRCV_globals$pb.counter + MRCV_globals$weight.B.max
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
    boot.sample <- cbind(model.vars, data)
    model.data.unsorted <- data.frame(matrix(data = NA, nrow = nrows, ncol = 1))
    counter <- 0
    if (nvars == 2) {
      for (i in 1:I) {
        for (j in 1:J) {
          counter <- counter + 2^nvars
          model.data.unsorted[(counter - 3), 1] <- sum(boot.sample[((boot.sample[, i] == 0) &
            (boot.sample[, (I + j)] == 0)), (I + J + 1)])
          model.data.unsorted[(counter - 2), 1] <- sum(boot.sample[((boot.sample[, i] == 0) &
            (boot.sample[, (I + j)] == 1)), (I + J + 1)])
          model.data.unsorted[(counter - 1), 1] <- sum(boot.sample[((boot.sample[, i] == 1) &
            (boot.sample[, (I + j)] == 0)), (I + J + 1)])
          model.data.unsorted[(counter), 1] <- sum(boot.sample[((boot.sample[, i] == 1) &
            (boot.sample[, (I + j)] == 1)), (I + J + 1)])
        }
      }
    }
    if (nvars == 3) {
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            counter <- counter + 2^nvars
            model.data.unsorted[(counter - 7), 1] <- sum(boot.sample[((boot.sample[, i] == 0) &
              (boot.sample[, (I + j)] == 0) &
              (boot.sample[, (I + J + k)] == 0)), (I + J + K + 1)])
            model.data.unsorted[(counter - 6), 1] <- sum(boot.sample[((boot.sample[, i] == 0) &
              (boot.sample[, (I + j)] == 0) &
              (boot.sample[, (I + J + k)] == 1)), (I + J + K + 1)])
            model.data.unsorted[(counter - 5), 1] <- sum(boot.sample[((boot.sample[, i] == 0) &
              (boot.sample[, (I + j)] == 1) &
              (boot.sample[, (I + J + k)] == 0)), (I + J + K + 1)])
            model.data.unsorted[(counter - 4), 1] <- sum(boot.sample[((boot.sample[, i] == 0) &
              (boot.sample[, (I + j)] == 1) &
              (boot.sample[, (I + J + k)] == 1)), (I + J + K + 1)])
            model.data.unsorted[(counter - 3), 1] <- sum(boot.sample[((boot.sample[, i] == 1) &
              (boot.sample[, (I + j)] == 0) &
              (boot.sample[, (I + J + k)] == 0)), (I + J + K + 1)])
            model.data.unsorted[(counter - 2), 1] <- sum(boot.sample[((boot.sample[, i] == 1) &
              (boot.sample[, (I + j)] == 0) &
              (boot.sample[, (I + J + k)] == 1)), (I + J + K + 1)])
            model.data.unsorted[(counter - 1), 1] <- sum(boot.sample[((boot.sample[, i] == 1) &
              (boot.sample[, (I + j)] == 1) &
              (boot.sample[, (I + J + k)] == 0)), (I + J + K + 1)])
            model.data.unsorted[(counter), 1] <- sum(boot.sample[((boot.sample[, i] == 1) &
              (boot.sample[, (I + j)] == 1) &
              (boot.sample[, (I + J + k)] == 1)), (I + J + K + 1)])
          }
        }
      }
    }
    model.data.unsorted[, 1] <- apply(
      X = as.matrix(model.data.unsorted[, 1]), MARGIN = 1,
      FUN = check.zero, add.constant = add.constant
    )
    model.data.unsorted <- data.matrix(model.data.unsorted)
  }
  model.data.unsorted
}
