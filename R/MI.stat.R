# MI.stat()
# A function that calculates X^2_S and X^2_S.ij

MI.stat <- function(data, I, J, summary.data = FALSE, add.constant = .5) {
  nvars <- 1 + ((I > 1) & (J > 1))
  if (nvars == 1) {
    output <- MMI.stat(
      data = data, I = I, J = J, summary.data = summary.data,
      add.constant = add.constant
    )
  }
  if (nvars == 2) {
    output <- SPMI.stat(
      data = data, I = I, J = J, summary.data = summary.data,
      add.constant = add.constant
    )
  }
  output
}
