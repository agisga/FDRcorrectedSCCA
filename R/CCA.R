#' @importFrom foreach %dopar%
NULL

#' L1-penalized SCCA with LOOCV
#'
#' L1-penalized SCCA with hyperparameter choice via leave-one-out cross-validation
#'
#' @export
L1_CCA_loocv <- function(X, Y, lambda_X, lambda_Y, num_cores) {

  doParallel::registerDoParallel(cores = num_cores)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X Y dimensions mismatch!")

  # CCA_out is a list of lists where:
  # - CCA_out[[i]] corresponds to the i'th entry of lambda_X
  # - CCA_out[[i]][[j]] corresponds to the j'th entry of lambda_Y
  # - CCA_out[[i]][[j]][[k]] contains the CCA result with k'th observation left out

  CCA_out <- vector(mode = "list", length = length(lambda_X))

  for (i in 1:length(lambda_X)) {
    CCA_out[[i]] <- vector(mode = "list", length = length(lambda_Y))
    for (j in 1:length(lambda_Y)) {
      CCA_out[[i]][[j]] <- foreach::foreach(k = 1:n) %dopar% {
        X_sub <- X[-k , ]
        Y_sub <- Y[-k , ]
        PMA::CCA(x = X_sub, z = Y_sub, typex = "standard", typez = "standard", K = 1,
                 penaltyx = lambda_X[i], penaltyz = lambda_Y[j], trace = FALSE)
      }
    }
  }

  # compute CV scores

  get_score_ijk <- function(i, j, k) {
    u_ijk <- CCA_out[[i]][[j]][[k]]$u
    v_ijk <- CCA_out[[i]][[j]][[k]]$v
    return(as.double(crossprod(u_ijk, X[k, ]) * crossprod(Y[k, ], v_ijk)))
  }

  scores <- foreach::foreach(i = 1:length(lambda_X), .combine = "rbind") %dopar% {
    sapply(1:length(lambda_Y),
           function(j) {
             sum(sapply(1:n, function(k) get_score_ijk(i, j, k) ))
           })
  }

  rownames(scores) <- paste0("lambda_X_", 1:length(lambda_X))
  colnames(scores) <- paste0("lambda_Y_", 1:length(lambda_Y))

  # identify lambdas corresponding to the largest CV score
  ind <- which(scores == max(abs(scores)), arr.ind = TRUE)
  best_lambda_X <- lambda_X[ind[1]]
  best_lambda_Y <- lambda_Y[ind[2]]

  # return the best model
  best_CCA <- PMA::CCA(x = X, z = Y, trace = FALSE, K = 1,
                       typex = "standard", typez = "standard",
                       penaltyx = best_lambda_X, penaltyz = best_lambda_Y)

  return(list("best_model" = best_CCA,
              "CV_scores"  = scores))
}

#' L1-penalized SCCA with sparsity bounds
#'
#' L1-penalized SCCA with prespecified sparsity levels for the coefficient vectors
#'
#' Sparse CCA with L1 penalty is performed repeatedly, while the tuning
#' parameters are optimized to yield a result that is close to a
#' prespecified sparsity
#'
#' @export
L1_CCA_with_sparsity_bound <- function(X, Y, bound_X, bound_Y,
                                       tolerance, verbose = FALSE,
                                       max_iter = 100) {

  n <- nrow(X)
  if (nrow(Y) != n) stop("X Y dimensions mismatch!")

  # use a bisection method to select lambdas such that the numbers of
  # non-zero entries of u and v are bounded above by `bound_X` and `bound_Y`

  lambda_X_upper <- 1
  lambda_X_lower <- 0
  lambda_Y_upper <- 1
  lambda_Y_lower <- 0

  lambda_X <- (lambda_X_upper + lambda_X_lower) / 2
  lambda_Y <- (lambda_Y_upper + lambda_Y_lower) / 2

  lambda_X_old <- 0
  lambda_Y_old <- 0

  lambda_X_not_found = TRUE
  lambda_Y_not_found = TRUE

  iter <- 0

  while((lambda_X_not_found | lambda_Y_not_found) & iter < max_iter) {
    iter <- iter + 1

    lambda_X <- (lambda_X_upper + lambda_X_lower) / 2
    lambda_Y <- (lambda_Y_upper + lambda_Y_lower) / 2

    model <- PMA::CCA(x = X, z = Y, typex = "standard", typez = "standard", K = 1,
                      penaltyx = lambda_X, penaltyz = lambda_Y, trace = FALSE)

    if (verbose) {
      print(paste("iteration", iter, "lambda_X =", lambda_X,
                  "|", "lambda_Y =", lambda_Y, "|",
                  "cor(Xu, Yv) =", cor(X %*% model$u, Y %*% model$v)))
    }

    # stop when the number of non-zeros is equal to `bound`,
    # otherwise adjust lower and upper bounds on lambdas

    if (sum(model$u != 0) < bound_X) {
      lambda_X_lower <- lambda_X
    } else if (sum(model$u != 0) > bound_X) {
      lambda_X_upper <- lambda_X
    } else {
      lambda_X_not_found <- FALSE
    }

    if (sum(model$v != 0) < bound_Y) {
      lambda_Y_lower <- lambda_Y
    } else if (sum(model$v != 0) > bound_Y) {
      lambda_Y_upper <- lambda_Y
    } else {
      lambda_Y_not_found <- FALSE
    }

    # alternatively stop when upper and lower bounds get too close
    if ((lambda_X_upper - lambda_X_lower) < tolerance &
        (lambda_Y_upper - lambda_Y_lower) < tolerance) {
      lambda_X_not_found <- FALSE
      lambda_Y_not_found <- FALSE
    }

    lambda_X_old <- lambda_X
    lambda_Y_old <- lambda_Y
  }

  return(list("best_model" = model,
              "lambda_X"   = lambda_X,
              "lambda_Y"   = lambda_Y))
}

#' L1-penalized SCCA with CV
#'
#' L1-penalized SCCA with hyperparameter choice via k-fold cross-validation
#'
#' The considered search space is effectively one-dimensional by considering
#' penalty pairs (lambda_X[1], lambda_Y[1]), (lambda_X[2], lambda_Y[2]), etc.
#' The number of cross-validation folds can be specified.
#'
#' @export
L1_CCA_CV <- function(X, Y, lambda_X, lambda_Y, n_folds = 5) {

  n <- nrow(X)
  if (nrow(Y) != n) stop("X Y dimensions mismatch!")

  # considered are penalty pairs (lambda_X[1], lambda_Y[1]), (lambda_X[2], lambda_Y[2]), etc. (one-dimensional search)
  if (length(lambda_X) != length(lambda_Y)) stop("lambda_X and lambda_Y need to be of the same length!")

  # CCA_out is a list of lists where:
  # - CCA_out[[i]] corresponds to lambda_X[i], lambda_Y[i]
  # - CCA_out[[i]][[j]] contains the SCCA run on the j'th CV fold combination

  fold_ind <- caret::createFolds(rnorm(n), k = n_folds, list = TRUE,
                                 returnTrain = FALSE)

  CCA_out <- vector(mode = "list", length = length(lambda_X))

  for (i in 1:length(lambda_X)) {
    for (j in 1:length(fold_ind)) {
      X_sub <- X[-fold_ind[[j]], ]
      Y_sub <- Y[-fold_ind[[j]], ]
      CCA_out[[i]][[j]] <- PMA::CCA(x = X_sub, z = Y_sub,
                                    typex = "standard",
                                    typez = "standard",
                                    K = 1, trace = FALSE,
                                    penaltyx = lambda_X[i],
                                    penaltyz = lambda_Y[i])
    }
  }

  # compute CV scores on test folds

  get_score_ij <- function(i, j) {
    u_ij <- CCA_out[[i]][[j]]$u
    v_ij <- CCA_out[[i]][[j]]$v
    return(abs(as.double(t(X[fold_ind[[j]], ] %*% u_ij) %*% (Y[fold_ind[[j]], ] %*% v_ij))))
  }

  scores <- sapply(1:length(lambda_X),
    function(i) {
      mean(sapply(1:n_folds, function(j) get_score_ij(i, j) ))
    })

  # identify lambdas corresponding to the largest CV score
  ind <- min(which(scores == max(abs(scores))))
  best_lambda_X <- lambda_X[ind]
  best_lambda_Y <- lambda_Y[ind]

  # return the best model
  best_CCA <- PMA::CCA(x = X, z = Y, trace = FALSE, K = 1, niter = 1000,
                       typex = "standard", typez = "standard",
                       penaltyx = best_lambda_X, penaltyz = best_lambda_Y)

  return(best_CCA)
}
