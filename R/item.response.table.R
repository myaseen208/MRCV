#' @name    item.response.table
#' @aliases item.response.table
#' @title Create an Item Response Table or Data Frame
#' @description
#' The \code{item.response.table} function is used to summarize data arising
#' from one, two, or three MRCVs.  For the one and two MRCV cases, a
#' cross-tabulation of the positive and negative responses for each (Wi, Yj)
#' pair is presented as a table or data frame (where Wi = W for the one MRCV
#' case).  For the three MRCV case, a cross-tabulation of the positive and
#' negative responses for each (Wi, Yj) pair is presented conditional on the
#' response for each Zk.
#'
#'
#' @param data A data frame containing the raw data where rows correspond to
#' the individual item response vectors, and columns correspond to the items
#' W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
#' @param I The number of items corresponding to row variable W.  I = 1 for the
#' one MRCV case.
#' @param J The number of items corresponding to column variable Y.
#' @param K The number of items corresponding to strata variable Z.
#' @param create.dataframe A logical value indicating whether the results
#' should be presented as a data frame instead of a table.
#' @return For \code{create.dataframe = FALSE}, \code{item.response.table} uses
#' the \code{\link[tables:tabular]{tables:tabular()}} function to produce
#' tables of marginal counts.
#'
#' For \code{create.dataframe = TRUE}, \code{item.response.table} returns the
#' same information as above but presents it as a data frame.  For the one MRCV
#' case, the data frame contains rx2J rows and 4 columns generically named
#' \code{W}, \code{Y}, \code{yj}, and \code{count}.  For the two MRCV case, the
#' data frame contains 2Ix2J rows and 5 columns named \code{W}, \code{Y},
#' \code{wi}, \code{yj}, and \code{count}.  For the three MRCV case, the data
#' frame contains 2Ix2Jx2K rows and 7 columns named \code{W}, \code{Y},
#' \code{Z}, \code{wi}, \code{yj}, \code{zk}, and \code{count}.
#' @seealso The \code{\link{marginal.table}} function for creating a marginal
#' table that summarizes only the positive responses for each pair.
#' @examples
#'
#' # Create an item response table for 1 SRCV and 1 MRCV
#' farmer.irtable.one <- item.response.table(data = farmer1, I = 1, J = 5)
#' farmer.irtable.one
#'
#' # Create an item response data frame for 1 SRCV and 1 MRCV
#' farmer.irdataframe.one <- item.response.table(
#'   data = farmer1, I = 1, J = 5,
#'   create.dataframe = TRUE
#' )
#' farmer.irdataframe.one
#'
#' # Create an item response table for 2 MRCVs
#' farmer.irtable.two <- item.response.table(data = farmer2, I = 3, J = 4)
#' farmer.irtable.two
#'
#' # Create an item response table for 3 MRCVs
#' farmer.irtable.three <- item.response.table(data = farmer3, I = 3, J = 4, K = 5)
#' farmer.irtable.three
#'
#' @export
item.response.table <- function(data, I, J, K = NULL, create.dataframe = FALSE) {
  op <- options()
  on.exit(options(op))
  options(warn = 1)
  # if((class(data)!="data.frame")&(class(data)!="matrix")) {
  #   stop("The \"data\" argument requires an object of class \"data.frame\".")
  # }
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
  # if ((class(create.dataframe)!="logical")&(create.dataframe!=1)&(create.dataframe!=0)) {
  #  warning("The \"create.dataframe\" argument requires a logical value. \n  The input value has been changed to the default value of FALSE.")
  #  create.dataframe<-FALSE
  # }
  if (!is.logical(create.dataframe) & (create.dataframe != 1) & (create.dataframe != 0)) {
    warning("The \"create.dataframe\" argument requires a logical value. \n  The input value has been changed to the default value of FALSE.")
    create.dataframe <- FALSE
  }

  create.dataframe <- as.logical(create.dataframe)
  nvars <- 1 + ((I > 1) & (J > 1)) + is.numeric(K)
  if (nvars == 1) {
    srcv <- ifelse(I == 1, 1, I + 1)
    mrcv <- ifelse(I == 1, 2, 1)
    r <- length(levels(as.factor(data[, srcv])))
    c <- ifelse(I == 1, J, I)
    wy.count <- as.data.frame(matrix(data = NA, nrow = 2 * r * c, ncol = 4))
    counter <- 0
    for (i in 1:c) {
      wy.count[((counter * r + 1):(counter * r + r)), 1] <- levels(as.factor(data[, srcv]))
      wy.count[((counter * r + 1):(counter * r + r)), 2] <- rep(names(data)[(mrcv + i - 1)], r)
      wy.count[((counter * r + 1):(counter * r + r)), 3] <- rep(0, r)
      wy.count[((counter * r + 1):(counter * r + r)), 4] <- (as.matrix(table(
        data[, srcv],
        data[, (mrcv + i - 1)]
      ))[, 1])
      wy.count[((counter * r + r + 1):(counter * r + 2 * r)), 1] <- levels(as.factor(data[, srcv]))
      wy.count[((counter * r + r + 1):(counter * r + 2 * r)), 2] <- rep(names(data)[(mrcv + i - 1)], r)
      wy.count[((counter * r + r + 1):(counter * r + 2 * r)), 3] <- rep(1, r)
      wy.count[((counter * r + r + 1):(counter * r + 2 * r)), 4] <- (as.matrix(table(
        data[, srcv],
        data[, (mrcv + i - 1)]
      ))[, 2])
      counter <- counter + 2
    }
    wy.count[, 1] <- factor(wy.count[, 1])
    wy.count[, 2] <- factor(wy.count[, 2], names(data)[mrcv:(mrcv + c - 1)])
    wy.count[, 3] <- as.numeric(wy.count[, 3])
    wy.count[, 4] <- as.numeric(wy.count[, 4])
    colnames(wy.count) <- c("W", "Y", "yj", "count")
    if (!create.dataframe) {
      output <- tabular(
        Heading() * count * Heading() * W
        ~ Heading() * Y * Heading() * Factor(yj, yj, c("0", "1   ")) * Heading() * mean,
        data = wy.count
      )
    }
    if (create.dataframe) {
      output <- wy.count
    }
  }
  nrows <- (2^nvars) * I * J * max(1, K)
  if (nvars == 2) {
    wy.count <- as.data.frame(matrix(data = NA, nrow = nrows, ncol = (2 * nvars + 1)))
    counter <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        counter <- counter + 2^nvars
        cell.count <- table(data[, i], data[, I + j])
        wy.count[(counter - 3), ] <- c(
          names(data)[i], names(data)[(I + j)], 0, 0,
          cell.count[1, 1]
        )
        wy.count[(counter - 2), ] <- c(
          names(data)[i], names(data)[(I + j)], 0, 1,
          cell.count[1, 2]
        )
        wy.count[(counter - 1), ] <- c(
          names(data)[i], names(data)[(I + j)], 1, 0,
          cell.count[2, 1]
        )
        wy.count[(counter), ] <- c(
          names(data)[i], names(data)[(I + j)], 1, 1,
          cell.count[2, 2]
        )
      }
    }
    wy.count[, 1] <- factor(wy.count[, 1], names(data)[1:I])
    wy.count[, 2] <- factor(wy.count[, 2], names(data)[(I + 1):(I + J)])
    wy.count[, (nvars + 1):(2 * nvars + 1)] <- lapply(
      wy.count[, (nvars + 1):(2 * nvars + 1)],
      as.numeric
    )
    colnames(wy.count) <- c("W", "Y", "wi", "yj", "count")
    if (!create.dataframe) {
      output <- tabular(
        count * (mean) * W * Heading() * Factor(wi, wi, c("0", "1 "))
        ~ Heading() * Y * Heading() * Factor(yj, yj, c("0", "1   ")),
        data = wy.count,
        suppressLabels = 3
      )
    }
    if (create.dataframe) {
      output <- wy.count
    }
  }
  if (nvars == 3) {
    wyz.count <- as.data.frame(matrix(data = NA, nrow = nrows, ncol = (2 * nvars + 1)))
    counter <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          counter <- counter + 2^nvars
          cell.count <- table(data[, i], data[, (I + j)], data[, (I + J + k)])
          wyz.count[(counter - 7), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 0, 0, 0, cell.count[1, 1, 1]
          )
          wyz.count[(counter - 6), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 0, 0, 1, cell.count[1, 1, 2]
          )
          wyz.count[(counter - 5), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 0, 1, 0, cell.count[1, 2, 1]
          )
          wyz.count[(counter - 4), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 0, 1, 1, cell.count[1, 2, 2]
          )
          wyz.count[(counter - 3), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 1, 0, 0, cell.count[2, 1, 1]
          )
          wyz.count[(counter - 2), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 1, 0, 1, cell.count[2, 1, 2]
          )
          wyz.count[(counter - 1), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 1, 1, 0, cell.count[2, 2, 1]
          )
          wyz.count[(counter), ] <- c(
            names(data)[i], names(data)[(I + j)],
            names(data)[(I + J + k)], 1, 1, 1, cell.count[2, 2, 2]
          )
        }
      }
    }
    wyz.count[, 1] <- factor(wyz.count[, 1], names(data)[1:I])
    wyz.count[, 2] <- factor(wyz.count[, 2], names(data)[(I + 1):(I + J)])
    wyz.count[, 3] <- factor(wyz.count[, 3], names(data)[(I + J + 1):(I + J + K)])
    wyz.count[, (nvars + 1):(2 * nvars + 1)] <- lapply(
      wyz.count[, (nvars + 1):(2 * nvars + 1)],
      as.numeric
    )
    colnames(wyz.count) <- c("W", "Y", "Z", "wi", "yj", "zk", "count")
    if (!create.dataframe) {
      output <- vector("list", K * 2)
      output.names <- matrix(data = NA, nrow = 1, ncol = K * 2)
      counter <- 0
      for (k in 1:K) {
        for (c in 0:1) {
          counter <- counter + 1
          output[counter] <- list(tabular(
            count * (mean) * W * Heading() * Factor(wi, wi, c("0", "1 "))
            ~ Heading() * Y * Heading() * Factor(yj, yj, c("0", "1   ")),
            data = wyz.count[((wyz.count[, 3] == names(data)[I + J + k]) & (wyz.count[, 6] == c)), ],
            suppressLabels = 3
          ))
          output.names[1, counter] <- paste(names(data[(I + J + k)]), "=", c)
        }
      }
      names(output) <- output.names
    }
    if (create.dataframe) {
      output <- wyz.count
    }
  }
  output
}
