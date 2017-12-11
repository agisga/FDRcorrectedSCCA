library(dplyr)
library(readr)
library(doParallel)

out_path <- "/lustre/project/wyp/agossman/FDRcorrectedSCCA/one_constant_cross-correlated_block/low-dim/"
out_path <- paste0(out_path, "constant_low_dimensional_FDRcorrectedSCCA/")

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
doParallel::registerDoParallel(cores)

setwd("../../../..")
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

cor_btwn <- 0.4 # (cross-)correlation between *significant* features of X and Y
cor_wthn <- 0.1 # correlation between any pair of features within either X or Y
lambda   <- 0.9 # tuning parameter for sparse CCA
fdr_vec  <- c(0.05, 0.1, 0.2) # target FDR

simulation_results_df <- NULL

for (fdr in fdr_vec) {
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

        #---  Pre selection via CCA on X0 and Y0

        CCA0 <- PMA::CCA(x = X0, z = Y0,
                         typex = "standard", typez = "standard",
                         K = 1, trace = FALSE, niter = 5000,
                         penaltyx = lambda, penaltyz = lambda)
        u0 <- CCA0$u
        v0 <- CCA0$v

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

        #--- Use X1 and Y1 to estimate covariance matrix of X2tY2v0/sqrt(n2) (according to the asymptotic distribution)

        Omega1 <- get_Omega(v0, S_Y1, S_X1, S_X1Y1)

        #--- Run BHq (use v0, Omega1, Y2 and X2 for the FDR-correction step)

        # response variable
        X2tY2v0 <- c( t(X2) %*% Y2 %*% v0 )
        y <- X2tY2v0 / sqrt(nrow(X2))
        y <- y / sqrt(diag(Omega1))

        # BH correction
        p_values <- 2*(1 - pnorm(abs(y))) # p-values under null hypothesis of 0 mean
        X_selected <- BH_select(p_values, fdr)

        # get the indices corresponding to the original matrices
        X_selected <- as.integer(colnames(X2)[X_selected])

        # return the results
        c("TP_v" = length(which(Y_selected <= n_signif_Y)),
          "FP_v" = length(which(Y_selected > n_signif_Y)),
          "TP_u" = length(which(X_selected <= n_signif_X)),
          "FP_u" = length(which(X_selected > n_signif_X)))
      }

      sim_df <- t(simulation_results) %>% tbl_df %>%
        mutate(FDP_v = FP_v / max((TP_v + FP_v), 1),
               FDP_u = FP_u / max((TP_u + FP_u), 1)) %>%
        mutate(n_row = as.integer(n_row),
               n_col_X = as.integer(n_col_X),
               n_col_Y = as.integer(n_col_Y),
               n_signif_X = as.integer(n_signif_X),
               n_signif_Y = as.integer(n_signif_Y),
               target_FDR = fdr,
               iter = 1:num_iter)

      if (is.null(simulation_results_df)) {
        simulation_results_df <- sim_df
      } else {
        simulation_results_df <- bind_rows(simulation_results_df, sim_df)
      }

      # save just in case
      write_csv(simulation_results_df,
                path = paste0(out_path,
                              "constant_low_dimensional_",
                              slurm_task_id, ".csv"))
    }
  }
}

save(list = ls(),
     file = paste0(out_path,
                   "constant_low_dimensional_",
                   slurm_task_id,
                   ".RData"))

print("DONE!")
