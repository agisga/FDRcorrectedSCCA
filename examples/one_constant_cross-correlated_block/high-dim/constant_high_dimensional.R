library(dplyr)
library(readr)
library(doParallel)

out_path <- "/lustre/project/wyp/agossman/FDRcorrectedSCCA/one_constant_cross-correlated_block/high-dim/"
out_path <- paste0(out_path, "constant_high_dimensional_FDRcorrectedSCCA/")

cmd_args <- commandArgs(TRUE)

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
doParallel::registerDoParallel(cores)

setwd("../../../..")
devtools::load_all()

set.seed(20170225)

num_iter       <- 500
n_row          <- 600
n_signif_X_vec <- c(1, 20, 40, 60, 80, 100, 120)
n_signif_Y_vec <- n_signif_X_vec
n_col_X        <- 1500
n_col_Y        <- 1500
n_col          <- n_col_X + n_col_Y

cor_btwn <- 0.4 # (cross-)correlation between *significant* features of X and Y
cor_wthn <- 0.1 # correlation between any pair of features within either X or Y
fdr <- as.numeric(cmd_args[1])
subset_divisor <- as.numeric(cmd_args[2])

simulation_results_df <- NULL

for (n_signif_X in n_signif_X_vec) {
  for (n_signif_Y in n_signif_Y_vec) {

    # generate the correlation matrix of matrix [X | Y]
    Sigma <- matrix(0, n_col, n_col)
    Sigma[1:n_col_X, 1:n_col_X] <- cor_wthn
    Sigma[(n_col_X+1):n_col, (n_col_X+1):n_col] <- cor_wthn
    Sigma[1:n_signif_X, (n_col_X+1):(n_col_X+n_signif_Y)] <- cor_btwn
    Sigma[(n_col_X+1):(n_col_X+n_signif_Y), 1:n_signif_X] <- cor_btwn
    Sigma[1:n_signif_X, 1:n_signif_X] <- Sigma[1:n_signif_X, 1:n_signif_X] + cor_btwn
    Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] <- Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] + cor_btwn
    diag(Sigma) <- 1

    # Cholesky factorization
    Sigma_chol <- chol(Sigma)

    simulation_results <- foreach (iter = 1:num_iter, .combine = "cbind") %dopar% {

      print(paste("Iteration", iter))

      XY <- simulate_MVN_data(n_row, n_col_X, n_col_Y, Sigma_chol)

      #--- Divide the data into three datasets

      size_frac <- ceiling(n_row / 3)
      ind0 <- 1:size_frac
      ind1 <- (size_frac + 1):(2 * size_frac)
      ind2 <- (2 * size_frac + 1):n_row

      X0 <- XY$X[ind0, ]
      X1 <- XY$X[ind1, ]
      X2 <- XY$X[ind2, ]

      Y0 <- XY$Y[ind0, ]
      Y1 <- XY$Y[ind1, ]
      Y2 <- XY$Y[ind2, ]

      #--- Pre selection via CCA on X0 and Y0

      CCA0 <- L1_CCA_with_sparsity_bound(X = X0, Y = Y0,
                                         bound_X = floor(size_frac / subset_divisor),
                                         bound_Y = floor(size_frac / subset_divisor),
                                         tolerance = 0.01,
                                         verbose = FALSE,
                                         max_iter = 1000)

      # Try getting better coefficients by performing more CCA iterations
      CCA0_more_iter <- PMA::CCA(x = X0, z = Y0,
                                 typex = "standard", typez = "standard",
                                 penaltyx = CCA0$best_model$penaltyx,
                                 penaltyz = CCA0$best_model$penaltyz,
                                 v = CCA0$best_model$v, K = 1,
                                 trace = FALSE, niter = 5000)
      u0 <- CCA0_more_iter$u
      v0 <- CCA0_more_iter$v

      # get subsets obtained by the CCA pre selection on X0 and Y0
      selected_X <- which(u0 != 0)
      selected_Y <- which(v0 != 0)

      # restrict vectors and matrices to these subsets of variables from here on
      X0 <- X0[ , selected_X]
      Y0 <- Y0[ , selected_Y]
      X1 <- X1[ , selected_X]
      Y1 <- Y1[ , selected_Y]
      X2 <- X2[ , selected_X]
      Y2 <- Y2[ , selected_Y]

      u0 <- u0[selected_X]
      v0 <- v0[selected_Y]

      colnames(X0) <- selected_X
      colnames(Y0) <- selected_Y
      colnames(X1) <- selected_X
      colnames(Y1) <- selected_Y
      colnames(X2) <- selected_X
      colnames(Y2) <- selected_Y
      names(u0) <- selected_X
      names(v0) <- selected_Y

      #--- Use X1 and Y1 to estimate covariance matrix of Y2tX2u0/sqrt(n2) (according to the asymptotic distribution)

      S_X1 <- cov(X1)
      S_Y1 <- cov(Y1)
      S_X1Y1 <- cov(X1, Y1)
      S_Y1X1 <- t(S_X1Y1)

      Omega1 <- get_Omega(u0, S_X1, S_Y1, S_Y1X1)

      #--- Run BHq (use u0, Omega1, X2 and Y2 for the FDR-correction step)

      # response variable
      Y2tX2u0 <- c( t(Y2) %*% X2 %*% u0 )
      y <- Y2tX2u0 / sqrt(nrow(Y2))
      y <- y / sqrt(diag(Omega1))

      # BH correction
      p_values <- 2*(1 - pnorm(abs(y))) # p-values under null hypothesis of 0 mean
      Y_selected <- BH_select(p_values, fdr)

      # get the indices corresponding to the original matrices
      Y_selected <- as.integer(colnames(Y2)[Y_selected])

      c("TP" = length(which(Y_selected <= n_signif_Y)),
        "FP" = length(which(Y_selected > n_signif_Y)))
    }

    sim_df <- t(simulation_results) %>% tbl_df %>%
      mutate(FDP = FP / max((TP + FP), 1)) %>%
      mutate(n_row = rep(n_row, num_iter),
             n_col_X = rep(n_col_X, num_iter),
             n_col_Y = rep(n_col_Y, num_iter),
             n_signif_X = rep(n_signif_X, num_iter),
             n_signif_Y = rep(n_signif_Y, num_iter),
             target_FDR = rep(fdr, num_iter),
             iter = 1:num_iter)

    if (is.null(simulation_results_df)) {
      simulation_results_df <- sim_df
    } else {
      simulation_results_df <- bind_rows(simulation_results_df, sim_df)
    }

    # save just in case
    write_csv(simulation_results_df,
      path = paste0(out_path, "constant_high_dimensional_fdr_",
                    100*fdr, "_div_", subset_divisor, ".csv"))
  }
}

out_file <- paste0(out_path, "constant_high_dimensional_fdr_",
                   100*fdr, "_div_", subset_divisor, ".RData")
save(list = ls(), file = out_file)
