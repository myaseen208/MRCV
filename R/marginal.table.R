#' @name    marginal.table
#' @aliases marginal.table
#' @title Create a Marginal Table
#' @description
#' The \code{marginal.table} function is used to summarize only the positive
#' responses (joint positive responses for multiple MRCVs) for data arising
#' from one, two, or three MRCVs. This function essentially provides a subset
#' of counts from the \code{item.response.table} function.
#'
#'
#' @param data A data frame containing the raw data where rows correspond to
#' the individual item response vectors, and columns correspond to the items
#' W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
#' @param I The number of items corresponding to row variable W.  I = 1 for the
#' one MRCV case.
#' @param J The number of items corresponding to column variable Y.
#' @param K The number of items corresponding to strata variable Z.
#' @return The \code{marginal.table} function uses the
#' \code{\link[tables:tabular]{tables:tabular()}} function to produce a table
#' containing positive counts.
#' @seealso The \code{\link{item.response.table}} function for creating an item
#' response table or data frame that summarizes both the positive and negative
#' responses for each (Wi, Yj) pair (conditional on the response for each Zk in
#' the case of three MRCVs).
#' 
#' @examples
#'
#' # Create a marginal table for 1 MRCV
#' # farmer.mtable.one <- marginal.table(data = farmer1, I = 1, J = 5)
#' # farmer.mtable.one
#'
#' # Create a marginal table for 2 MRCVs
#' farmer.mtable.two <- marginal.table(data = farmer2, I = 3, J = 4)
#' farmer.mtable.two
#'
#' # Create a marginal table for 3 MRCVs
#' farmer.mtable.three <- marginal.table(data = farmer3, I = 3, J = 4, K = 5)
#' farmer.mtable.three
#'
#' @export
marginal.table <- function(data, I, J, K = NULL) {
  op <- options()
  on.exit(options(op))
  options(warn = 1)
  if (!(is.data.frame(data) || is.matrix(data))) {
    stop("The \"data\" argument requires an object of class \"data.frame\".")
  }
  data <- as.data.frame(data)
  if ((I == 0) | (J == 0)) {
    stop("\"I\" and \"J\" must be greater than 0.")
  }
  if (!is.numeric(I) | !is.numeric(J)) {
    stop("\"I\" and \"J\" must be numeric.")
  }
  if ((I %% 1 != 0) | (J %% 1 != 0)) {
    stop("\"I\" and \"J\" must be integers.")
  }
  I <- as.integer(I)
  J <- as.integer(J)
  if (!is.null(K)) {
    if (!is.numeric(K)) {
      stop("\"K\" must be numeric or NULL.")
    }
    if (K %% 1 != 0) {
      stop("\"K\" must be an integer.")
    }
    K <- as.integer(K)
    if (K == 0) {
      K <- NULL
    }
  }
  nvars <- 1 + ((I > 1) & (J > 1)) + is.numeric(K)
  if (nvars == 1) {
    srcv <- ifelse(I == 1, 1, I + 1)
    mrcv <- ifelse(I == 1, 2, 1)
    r <- length(levels(as.factor(data[, srcv])))
    c <- ifelse(I == 1, J, I)
    wy.count <- as.data.frame(matrix(data = NA, nrow = r, ncol = 4))
    merge.matrix <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 4))
    counter <- 0
    for (i in 1:c) {
      wy.count[1:r, 1] <- levels(as.factor(data[, srcv]))
      wy.count[1:r, 2] <- rep(names(data)[(mrcv + i - 1)], r)
      wy.count[1:r, 3] <- as.matrix(table(data[, srcv], data[, (mrcv + i - 1)]))[, 2]
      wy.count[1:r, 4] <- as.matrix(table(data[, srcv]))
      if (i == 1) {
        for (j in 1:r) {
          rep.matrix <- matrix(
            as.vector(rep(
              as.vector(wy.count[j, ]),
              wy.count[j, 4] - (c - 1)
            )), (wy.count[j, 4] - (c - 1)), 4,
            byrow = TRUE
          )
          merge.matrix <- rbind(merge.matrix, rep.matrix)
        }
        n.iplus <- as.matrix(table(data[, srcv]))
      }
      if (i > 1) {
        merge.matrix <- rbind(merge.matrix, wy.count)
      }
      counter <- counter + 1
    }
    merge.matrix <- merge.matrix[2:nrow(merge.matrix), ]
    merge.matrix[, 1] <- factor(unlist(merge.matrix[, 1]))
    merge.matrix[, 2] <- factor(unlist(merge.matrix[, 2]), names(data)[mrcv:(mrcv + c - 1)])
    merge.matrix[, 3] <- as.numeric(unlist(merge.matrix[, 3]))
    merge.matrix[, 4] <- as.numeric(unlist(merge.matrix[, 4]))
    colnames(merge.matrix) <- c("W", "Y", "count", "n")
    MRCV_globals$n.iplus <- n.iplus
    MRCV_globals$c <- c
    MRCV_globals$counter <- -1
    calc.percent <- function(x) {
      MRCV_globals$counter <- MRCV_globals$counter + 1
      perc <- round(100 * mean(x) / MRCV_globals$n.iplus[(trunc(MRCV_globals$counter / MRCV_globals$c) + 1)], 2)
      if (perc == 0) {
        perc <- "0.00"
      }
      perc
    }
    marginal.table <- tabular(Heading() * count * Heading() * W
    ~ Justify(l) * (n <- 1) + Heading() * Y * (Heading("count") * mean +
        Heading("  %     ") * calc.percent), data = merge.matrix)
  }
  n <- nrow(data)
  nrows <- I * J * max(1, K)
  calc.percent <- function(x) {
    perc <- round(100 * mean(x) / n, 2)
    if (perc == 0) {
      perc <- "0.00"
    }
    perc
  }
  if (nvars == 2) {
    wy.count <- as.data.frame(matrix(data = NA, nrow = nrows, ncol = (nvars + 1)))
    counter <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        counter <- counter + 1
        cell.count <- table(data[, i], data[, (I + j)])
        wy.count[counter, ] <- c(names(data)[i], names(data)[(I + j)], cell.count[2, 2])
      }
    }
    wy.count[, 1] <- factor(wy.count[, 1], names(data)[1:I])
    wy.count[, 2] <- factor(wy.count[, 2], names(data)[(I + 1):(I + J)])
    wy.count[, (nvars + 1)] <- as.numeric(wy.count[, (nvars + 1)])
    colnames(wy.count) <- c("W", "Y", "count")
    marginal.table <- tabular(
      count * W ~ Heading() * Y * (Heading("count") * mean +
        Heading("  %     ") * calc.percent),
      data = wy.count,
      suppressLabels = 2
    )
  }
  if (nvars == 3) {
    wyz.count <- as.data.frame(matrix(data = NA, nrow = nrows, ncol = (nvars + 1)))
    counter <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          counter <- counter + 1
          cell.count <- table(data[, i], data[, (I + j)], data[, (I + J + k)])
          wyz.count[counter, ] <- c(
            names(data)[i], names(data)[(I + j)], names(data)[(I + J + k)],
            cell.count[2, 2, 2]
          )
        }
      }
    }
    wyz.count[, 1] <- factor(wyz.count[, 1], names(data)[1:I])
    wyz.count[, 2] <- factor(wyz.count[, 2], names(data)[(I + 1):(I + J)])
    wyz.count[, 3] <- factor(wyz.count[, 3], names(data)[(I + J + 1):(I + J + K)])
    wyz.count[, (nvars + 1)] <- as.numeric(wyz.count[, (nvars + 1)])
    colnames(wyz.count) <- c("W", "Y", "Z", "count")
    marginal.table <- vector("list", K)
    for (k in 1:K) {
      marginal.table[k] <- list(tabular(
        count * W ~ Heading() * Y * (Heading("count") * mean +
          Heading("  %     ") * calc.percent),
        data = wyz.count[(wyz.count[, 3] == names(data)[I + J + k]), ],
        suppressLabels = 2
      ))
    }
    names(marginal.table) <- names(data[(I + J + 1):(I + J + K)])
  }
  marginal.table
}
