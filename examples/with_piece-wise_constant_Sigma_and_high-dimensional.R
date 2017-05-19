library(dplyr)
library(doParallel)
# adjust the number of cores to the particular system
cores <- 4
doParallel::registerDoParallel(cores)

# load auxilliary functions from the folder ../R/
setwd("..")
devtools::load_all()

# set some parameters
set.seed(20170318)
num_iter       <- 100
n_row_vec      <- c(100, 500, 1000)
n_signif_X_vec <- c(1, 10, 20)
n_signif_Y_vec <- c(5, 30, 60)
n_col_X        <- 1500
n_col_Y        <- 3000
n_col          <- n_col_X + n_col_Y
cor_btwn       <- 0.4 # (cross-)correlation between *significant* features of X and Y
cor_wthn       <- 0.1 # correlation between any pair of features within either X or Y

# set the target FDR level equal to the supplied command line argument
cmd_args <- commandArgs(TRUE)
fdr      <- as.numeric(cmd_args[1])

# the FDR-corrected SCCA algorithm is performed num_iter times for
# each combination of considered conditions,
# using a new randomly generated dataset each time
simulation_results_df <- NULL
for (n_row in n_row_vec) {
  for (n_signif_X in n_signif_X_vec) {
    for (n_signif_Y in n_signif_Y_vec) {

      # correlation matrix of matrix [X Y]
      Sigma <- matrix(0, n_col, n_col)

      # fill Sigma with values:
      Sigma[1:n_col_X, 1:n_col_X] <- cor_wthn
      Sigma[(n_col_X+1):n_col, (n_col_X+1):n_col] <- cor_wthn
      Sigma[1:n_signif_X, (n_col_X+1):(n_col_X+n_signif_Y)] <- cor_btwn
      Sigma[(n_col_X+1):(n_col_X+n_signif_Y), 1:n_signif_X] <- cor_btwn
      Sigma[1:n_signif_X, 1:n_signif_X] <- Sigma[1:n_signif_X, 1:n_signif_X] + cor_btwn
      Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] <- Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] + cor_btwn
      diag(Sigma) <- 1

      # Cholesky factorization
      Sigma_chol <- chol(Sigma)

      # perform FDR-corrected SCCA num_iter times
      simulation_results <- foreach (iter = 1:num_iter, .combine = "cbind") %dopar% {

        print(paste("Iteration", iter))

        XY <- simulate_MVN_data(n_row, n_col_X, n_col_Y, Sigma_chol)

        # Divide the data into two datasets
        half_size <- ceiling(n_row / 2)
        X0 <- XY$X[1:half_size, ]
        X1 <- XY$X[half_size:n_row, ]
        Y0 <- XY$Y[1:half_size, ]
        Y1 <- XY$Y[half_size:n_row, ]

        # Pre selection via CCA on X0 and Y0
        CCA0 <- L1_CCA_with_sparsity_bound(X = X0, Y = Y0,
                                           bound_X = floor(half_size / 2),
                                           bound_Y = floor(half_size / 2),
                                           tolerance = 0.01,
                                           verbose = FALSE,
                                           max_iter = 1000)
        u0 <- CCA0$best_model$u
        v0 <- CCA0$best_model$v

        # get subsets obtained by the CCA pre selection on X0 and Y0
        selected_X <- which(u0 != 0)
        selected_Y <- which(v0 != 0)

        # restrict vectors and matrices to these subsets of variables from here on
        X0 <- X0[ , selected_X]
        Y0 <- Y0[ , selected_Y]
        X1 <- X1[ , selected_X]
        Y1 <- Y1[ , selected_Y]
        u0 <- u0[selected_X]
        v0 <- v0[selected_Y]
        colnames(X0) <- selected_X
        colnames(Y0) <- selected_Y
        colnames(X1) <- selected_X
        colnames(Y1) <- selected_Y
        names(u0) <- selected_X
        names(v0) <- selected_Y

        # response variable
        Y1tX1u0 <- t(Y1) %*% X1 %*% u0

        # estimate covariance matrix of Y1tX1u0/sqrt(n1)
        # (according to the asymptotic distribution)
        S_X0 <- cov(X0)
        S_Y0 <- cov(Y0)
        S_X0Y0 <- cov(X0, Y0)
        S_Y0X0 <- t(S_X0Y0)
        Omega <- matrix(NA, ncol(Y0), ncol(Y0))
        for (i in 1:ncol(Y0)) {
          for (j in 1:ncol(Y0)) {
            Omega[i, j] <- crossprod(S_Y0X0[i, ], u0) * crossprod(S_Y0X0[j, ], u0) +
              S_Y0[i, j] * t(u0) %*% S_X0 %*% u0
          }
        }

        # Run BHq
        y <- Y1tX1u0 / sqrt(nrow(Y1))
        y <- y / sqrt(diag(Omega))
        p_values <- 2*(1 - pnorm(abs(y))) # p-values under null hypothesis of 0 mean
        p <- length(p_values)
        cutoff <- max(c(0, which(sort(p_values) <= fdr * (1:p) / p)))
        BH_selected <- which(p_values <= fdr * cutoff / p)

        # get the indices corresponding to the original matrices
        BH_selected <- as.integer(colnames(Y1)[BH_selected])

        # record the numbers of true positives and false positives
        c("TP" = length(which(BH_selected <= n_signif_Y)), "FP" = length(which(BH_selected > n_signif_Y)))
      }

      # collect all simulation results into a tidy data frame
      sim_df <- t(simulation_results) %>% tbl_df %>%
        mutate(FDP = FP / max((TP + FP), 1)) %>%
        mutate(n_row = rep(n_row, num_iter), n_signif_X = rep(n_signif_X, num_iter),
               n_signif_Y = rep(n_signif_Y, num_iter),
               target_FDR = rep(fdr, num_iter), iter = 1:num_iter)

      # concatenate the data frames obtained under different simulation settings
      if (is.null(simulation_results_df)) {
        simulation_results_df <- sim_df
      } else {
        simulation_results_df <- bind_rows(simulation_results_df, sim_df)
      }
    }
  }
}

# save the simulation results
save(list = ls(), file = paste0("./examples/RData/with_piece-wise_constant_sigma_and_high-dimensional_FDR_", fdr*100, ".RData"))
