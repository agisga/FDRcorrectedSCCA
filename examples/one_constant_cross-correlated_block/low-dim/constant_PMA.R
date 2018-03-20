library(PMA)
library(readr)
library(doParallel)
library(dplyr)

out_path <- "/lustre/project/wyp/agossman/FDRcorrectedSCCA/one_constant_cross-correlated_block/low-dim/"
out_path <- paste0(out_path, "constant_low_dimensional_PMA/")

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
doParallel::registerDoParallel(cores)

setwd("../../..")
devtools::load_all()

slurm_task_id  <- Sys.getenv("SLURM_ARRAY_TASK_ID") # run this as an array job 500/10=50 times
num_iter       <- 10
n_row          <- 3000
n_col_X        <- 50
n_col_Y        <- 50
n_col          <- n_col_X + n_col_Y
n_signif_X_vec <- c(1, 10, 20, 30)
n_signif_Y_vec <- c(1, 10, 20, 30)

# ensure a different random seed in each array run
rand_seed <- as.numeric(paste0(format(Sys.time(), "%M%S"), slurm_task_id))
set.seed(rand_seed)

cor_btwn   <- 0.4 # (cross-)correlation between *significant* features of X and Y
cor_wthn   <- 0.1 # correlation between any pair of features within either X or Y

simulation_results_df <- NULL

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

      set.seed(rand_seed + iter)

      print(paste("Iteration", iter))

      # Generate the data
      XY <- simulate_MVN_data(n_row, n_col_X, n_col_Y, Sigma_chol)
      X <- XY$X
      Y <- XY$Y

      # Apply SCCA on X and Y with the sparsity parameters selected by the permutation based approach of Witten et. al. (2009)
      perm_out <- CCA.permute(x = X, z = Y, typex = "standard",
                              typez = "standard", trace = FALSE)
      CCA_out  <- CCA(x = X, z = Y, K = 1,
                      typex = "standard", typez = "standard",
                      penaltyx = perm_out$bestpenaltyx,
                      penaltyz = perm_out$bestpenaltyz,
                      trace = FALSE, v = perm_out$v.init,
                      niter = 5000)

      u0 <- CCA_out$u
      v0 <- CCA_out$v

      # Identify true and false discoveries in v0
      Y_selected <- which(v0 != 0)
      X_selected <- which(u0 != 0)

      # Return the results
      c("TP_v" = length(which(Y_selected <= n_signif_Y)),
        "FP_v" = length(which(Y_selected > n_signif_Y)),
        "TP_u" = length(which(X_selected <= n_signif_X)),
        "FP_u" = length(which(X_selected > n_signif_X)),
        "lambda_X" = CCA_out$penaltyx,
        "lambda_Y" = CCA_out$penaltyz)
    }

    sim_df <- t(simulation_results) %>% tbl_df %>%
      mutate(FDP_v = FP_v / max((TP_v + FP_v), 1),
             FDP_u = FP_u / max((TP_u + FP_u), 1)) %>%
      mutate(n_row = as.integer(n_row),
             n_col_X = as.integer(n_col_X),
             n_col_Y = as.integer(n_col_Y),
             n_signif_X = as.integer(n_signif_X),
             n_signif_Y = as.integer(n_signif_Y),
             iter = 1:num_iter)

    if (is.null(simulation_results_df)) {
      simulation_results_df <- sim_df
    } else {
      simulation_results_df <- bind_rows(simulation_results_df, sim_df)
    }

    # save just in case
    write_csv(simulation_results_df,
              path = paste0(out_path,
                            "constant_low_dimensional_PMA_",
                            slurm_task_id, ".csv"))
  }
}

save(list = ls(),
     file = paste0(out_path,
                   "constant_low_dimensional_PMA_",
                   slurm_task_id,
                   ".RData"))

print("DONE!")
