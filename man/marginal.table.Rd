\name{marginal.table}
\alias{marginal.table}
\title{
Create a Marginal Table
}
\description{
The \code{marginal.table} function is used to summarize only the positive responses (joint positive responses for multiple MRCVs) for data arising from one, two, or three MRCVs. This function essentially provides a subset of counts from the \code{item.response.table} function.
}
\usage{
marginal.table(data, I, J, K = NULL)
}
\arguments{
  \item{data}{
A data frame containing the raw data where rows correspond to the individual item response vectors, and columns correspond to the items W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
}
  \item{I}{
The number of items corresponding to row variable W.  I = 1 for the one MRCV case.
}
  \item{J}{
The number of items corresponding to column variable Y.
}
  \item{K}{
The number of items corresponding to strata variable Z.
}
}
\value{
The \code{marginal.table} function uses the \code{\link[tables:tabular]{tables:tabular()}} function to produce a table containing positive counts.
}
\seealso{
The \code{\link{item.response.table}} function for creating an item response table or data frame that summarizes both the positive and negative responses for each (Wi, Yj) pair (conditional on the response for each Zk in the case of three MRCVs).
}
\examples{
# Create a marginal table for 1 MRCV
farmer.mtable.one <- marginal.table(data = farmer1, I = 1, J = 5) 
farmer.mtable.one

# Create a marginal table for 2 MRCVs
farmer.mtable.two <- marginal.table(data = farmer2, I = 3, J = 4) 
farmer.mtable.two

# Create a marginal table for 3 MRCVs
farmer.mtable.three <- marginal.table(data = farmer3, I = 3, J = 4, K = 5) 
farmer.mtable.three
}