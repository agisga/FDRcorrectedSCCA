library(dplyr)
library(readr)
library(doParallel)

out_path <- "/lustre/project/wyp/agossman/FDRcorrectedSCCA/one_constant_cross-correlated_block/high-dim/"
out_path <- paste0(out_path, "constant_high_dimensional_ordinary_SCCA/")

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
doParallel::registerDoParallel(cores)

setwd("../../..")
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
lambda_vec <- c(0.01, 0.3) # tuning parameter for sparse CCA

simulation_results_df <- NULL

for (lambda in lambda_vec) {
  for (n_signif_X in n_signif_X_vec) {
    for (n_signif_Y in n_signif_Y_vec) {

      # correlation matrix of matrix [X | Y]
      Sigma <- matrix(0, n_col, n_col)
      # correlations within X
      Sigma[1:n_col_X, 1:n_col_X] <- cor_wthn
      # correlations within Y
      Sigma[(n_col_X+1):n_col, (n_col_X+1):n_col] <- cor_wthn
      # correlations between X and Y
      Sigma[1:n_signif_X, (n_col_X+1):(n_col_X+n_signif_Y)] <- cor_btwn
      Sigma[(n_col_X+1):(n_col_X+n_signif_Y), 1:n_signif_X] <- cor_btwn
      # adjustment for positive definiteness
      Sigma[1:n_signif_X, 1:n_signif_X] <- Sigma[1:n_signif_X, 1:n_signif_X] + cor_btwn
      Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] <- Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] + cor_btwn
      # unit variances
      diag(Sigma) <- 1

      # Cholesky factorization
      Sigma_chol <- chol(Sigma)

      simulation_results <- foreach (iter = 1:num_iter, .combine = "cbind") %dopar% {

        print(paste("Iteration", iter))

        # Generate the data
        XY <- simulate_MVN_data(n_row, n_col_X, n_col_Y, Sigma_chol)
        X <- XY$X
        Y <- XY$Y

        # Apply SCCA on X0 and Y0
        CCA0 <- PMA::CCA(x = X, z = Y,
                         typex = "standard", typez = "standard", K = 1,
                         penaltyx = lambda, penaltyz = lambda,
                         trace = FALSE, niter = 5000)
        u0 <- CCA0$u
        v0 <- CCA0$v

        # Identify true and false discoveries in v0
        Y_selected <- which(v0 != 0)

        c("TP" = length(which(Y_selected <= n_signif_Y)),
          "FP" = length(which(Y_selected > n_signif_Y)),
          "lambda_X" = CCA0$penaltyx,
          "lambda_Y" = CCA0$penaltyz)
      }

      sim_df <- t(simulation_results) %>% tbl_df %>%
        mutate(FDP = FP / max((TP + FP), 1)) %>%
        mutate(n_row = rep(n_row, num_iter),
               n_col_X = rep(n_col_X, num_iter),
               n_col_Y = rep(n_col_Y, num_iter),
               n_signif_X = rep(n_signif_X, num_iter),
               n_signif_Y = rep(n_signif_Y, num_iter),
               iter = 1:num_iter)

      if (is.null(simulation_results_df)) {
        simulation_results_df <- sim_df
      } else {
        simulation_results_df <- bind_rows(simulation_results_df, sim_df)
      }

      # save just in case
      write_csv(simulation_results_df,
        path = paste0(out_path,
          "constant_high_dimensional_ordinary_SCCA.csv"))
    }
  }
}

out_file <- paste0(out_path,
                   "constant_high_dimensional_ordinary_SCCA.RData")
save(list = ls(), file = out_file)

print("DONE!")
