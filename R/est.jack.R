# est.jack()
# A function that calculates model estimated odds ratios based on n-1 observations

est.jack <- function(mu.hat, i, I, J, K) {
  nrows <- I * J * max(1, K)
  output <- exp(log(mu.hat[i]) - log(mu.hat[(nrows + i)]) - log(mu.hat[(2 * nrows + i)])
    + log(mu.hat[(3 * nrows + i)]))
}
