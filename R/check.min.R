# check.min()
# Function used in Bonferroni adjustments (p-value should not be greater than 1)

check.min <- function(value) {
  min(value, 1.00000)
}
