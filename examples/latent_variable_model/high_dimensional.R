library(PMA)
library(doParallel)
library(dplyr)
library(readr)

out_path <- "/lustre/project/wyp/agossman/FDRcorrectedSCCA/"
out_path <- paste0(out_path,
                   "latent_vars_blocked_design/",
                   "high_dim/")

setwd("../../..")
devtools::load_all()

cmd_args <- commandArgs(TRUE)
fdr      <- as.numeric(cmd_args[1])

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
doParallel::registerDoParallel(cores)

slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID") # run this as an array job 500/5=100 times
num_iter      <- 5
n_row         <- 900
n_signif_vars <- 60
n_col_X       <- 900
n_col_Y       <- n_col_X
#n_signif_blocks_vec <- c(1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60)
n_signif_blocks_vec <- c(1, 2, 3, 15)

set.seed(as.numeric(paste0(201710, slurm_task_id)))

# data generation from a latent variable model
# for restrictions on the input parameters see error message definitions in the beginning of function definition
simulate_blocked_data <- function(n_row, n_col_X, n_col_Y,
                                  n_signif_X, n_signif_Y,
                                  n_blocks, btwn_noise_sd,
                                  wthn_noise_sd = 1) {
  if ((n_col_X != n_col_Y) | (n_signif_X != n_signif_Y)) {
    stop("X and Y must have the same dimensions and the same number of significant features.")
  }
  n_col <- n_col_X
  n_signif <- n_signif_X

  block_size <- n_signif / n_blocks
  if (block_size %% 1.0 != 0.0) {
    stop("n_signif_X needs to be a multiple of block_size")
  }

  X_non_signif <- matrix(rnorm(n_row * (n_col - n_signif), sd = wthn_noise_sd), nrow = n_row)
  Y_non_signif <- matrix(rnorm(n_row * (n_col - n_signif), sd = wthn_noise_sd), nrow = n_row)
  X <- X_non_signif
  Y <- Y_non_signif

  for (i in 1:n_blocks) {
    z <- rnorm(n_row) # latent variable
    ones_vec <- rep(1, block_size)

    X_i <- z %x% t(ones_vec) +
      matrix(rnorm(n_row * block_size, sd = btwn_noise_sd), nrow = n_row)
    X <- cbind(X_i, X)

    Y_i <- z %x% t(ones_vec) +
      matrix(rnorm(n_row * block_size, sd = btwn_noise_sd), nrow = n_row)
    Y <- cbind(Y_i, Y)
  }

  return(list("X" = scale(X), "Y" = scale(Y)))
}

simulation_results_df <- NULL

for (n_signif_blocks in n_signif_blocks_vec) {

  simulation_results <- foreach (iter = 1:num_iter, .combine = "cbind") %dopar% {

    print(paste("Iteration", iter))

    btwn_noise_sd <- exp(-2*n_signif_blocks)
    XY <- simulate_blocked_data(n_row = n_row, n_col_X = n_col_X,
                                n_col_Y = n_col_Y,
                                n_signif_X = n_signif_vars,
                                n_signif_Y = n_signif_vars,
                                n_blocks = n_signif_blocks,
                                btwn_noise_sd = btwn_noise_sd)

    # Divide the data into two datasets

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

    # Pre selection via CCA on X0 and Y0
    CCA0 <- L1_CCA_with_sparsity_bound(X = X0, Y = Y0,
                                       bound_X = floor(size_frac / 2),
                                       bound_Y = floor(size_frac / 2),
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

    # Use X1 and Y1 to estimate covariance matrix of Y2tX2u0/sqrt(n2) (according to the asymptotic distribution)

    S_X1 <- cov(X1)
    S_Y1 <- cov(Y1)
    S_X1Y1 <- cov(X1, Y1)
    S_Y1X1 <- t(S_X1Y1)

    Omega1 <- get_Omega(u0, S_X1, S_Y1, S_Y1X1)

    # Run BHq (use u0, Omega1, X2 and Y2 for the FDR-correction step)

    # response variable
    Y2tX2u0 <- c( t(Y2) %*% X2 %*% u0 )
    y <- Y2tX2u0 / sqrt(nrow(Y2))
    y <- y / sqrt(diag(Omega1))

    # BH correction
    p_values_Y <- 2*(1 - pnorm(abs(y))) # p-values under null hypothesis of 0 mean
    Y_selected <- BH_select(p_values_Y, fdr)

    # get the indices corresponding to the original matrices
    Y_selected <- as.integer(colnames(Y2)[Y_selected])

    c("TP" = length(which(Y_selected <= n_signif_vars)),
      "FP" = length(which(Y_selected > n_signif_vars)))
  }

  sim_df <- t(simulation_results) %>% tbl_df %>%
    mutate(FDP = FP / max((TP + FP), 1)) %>%
    mutate(n_row = as.integer(n_row),
           n_signif_X = as.integer(n_signif_vars),
           n_signif_Y = as.integer(n_signif_vars),
           n_signif_blocks = as.integer(n_signif_blocks),
           target_FDR = as.integer(fdr),
           iter = 1:num_iter)

  if (is.null(simulation_results_df)) {
    simulation_results_df <- sim_df
  } else {
    simulation_results_df <- bind_rows(simulation_results_df, sim_df)
  }

  # save csv just in case
  write_csv(simulation_results_df,
            path = paste0(out_path,
                          "high_cor_fdr_", fdr * 100,
                          "_", slurm_task_id, ".csv"))
}

save(list = ls(),
     file = paste0(out_path,
                   "high_cor_fdr_", fdr * 100,
                   "_", slurm_task_id, ".RData"))
