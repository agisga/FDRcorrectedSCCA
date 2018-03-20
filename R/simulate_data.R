# FDRcorrectedSCCA  Copyright (C) 2017  Alexej Gossmann

#' Simulate multivariate normal data
#'
#' @export
simulate_MVN_data <- function(n_row, n_col_X, n_col_Y, Sigma_chol) {

  n_col <- n_col_X + n_col_Y
  if (n_col != ncol(Sigma_chol)) stop("Dimension mismatch!")

  XY <- matrix(rnorm(n_row * n_col), n_row, n_col) %*% Sigma_chol
  X <- XY[ , 1:n_col_X]
  Y <- XY[ , (n_col_X+1):n_col]

  return(list("X" = X, "Y" = Y))
}

#' Simulate data for CCA based on a latent variable model
#'
#' @export
simulate_latent_variable_data <- function(n_row, n_col_X, n_col_Y, n_signif_X,
                                          n_signif_Y, noise_var = 1,
                                          shared_var_wthn = 0,
                                          shared_var_btwn = 1) {

  # generate matrices with independent Gaussian entries

  X <- matrix(rnorm(n_row * n_col_X, sd = sqrt(noise_var)), n_row, n_col_X)
  Y <- matrix(rnorm(n_row * n_col_Y, sd = sqrt(noise_var)), n_row, n_col_Y)

  # introduce correlations between columns within each matrix

  shared_X <- if (shared_var_wthn > 0) {
    rnorm(n_row, sd = sqrt(shared_var_wthn))
  } else {
    rep(0, n_row)
  }

  X <- X + t(rep(1, n_col_X)) %x% shared_X

  shared_Y <- if (shared_var_wthn > 0) {
    rnorm(n_row, sd = sqrt(shared_var_wthn))
  } else {
    rep(0, n_row)
  }

  Y <- Y + t(rep(1, n_col_Y)) %x% shared_Y

  # introduce cross-correlations between columns of X and Y

  if (n_signif_X > 0 & n_signif_Y > 0) {
    shared_XY <- if(shared_var_btwn > 0) {
      rnorm(n_row, sd = sqrt(shared_var_btwn))
    } else {
      rep(0, n_row)
    }
    X[ , 1:n_signif_X] <- X[ , 1:n_signif_X] + t(rep(1, n_signif_X)) %x% shared_XY
    Y[ , 1:n_signif_Y] <- Y[ , 1:n_signif_Y] + t(rep(1, n_signif_Y)) %x% shared_XY
  }

  return(list("X" = X, "Y" = Y, "noise_var" = noise_var,
              "shared_var_wthn" = shared_var_wthn,
              "shared_var_btwn" = shared_var_btwn))
}
