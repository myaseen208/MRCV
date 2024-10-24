# check.zero()
# Function that adds a small constant to 0 cell counts

check.zero <- function(value, add.constant = .5) {
  if (value == 0) newvalue <- add.constant else newvalue <- value
  newvalue
}
