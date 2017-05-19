# FDRcorrectedSCCA  Copyright (C) 2017  Alexej Gossmann

simulate_MVN_data <- function(n_row, n_col_X, n_col_Y, Sigma_chol) {

  n_col <- n_col_X + n_col_Y
  if (n_col != ncol(Sigma_chol)) stop("Dimension mismatch!")

  XY <- matrix(rnorm(n_row * n_col), n_row, n_col) %*% Sigma_chol
  X <- XY[ , 1:n_col_X]
  Y <- XY[ , (n_col_X+1):n_col]

  return(list("X" = X, "Y" = Y))
}
