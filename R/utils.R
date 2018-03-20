#' Benjamini-Hochberg procedure
#'
#' @export
BH_select <- function(p_values, fdr) {
  p <- length(p_values)
  below_threshold <- sort(p_values) <= fdr * (1:p) / p
  all_below_threshold <- sapply(1:length(below_threshold),
                                function(i) { all(below_threshold[1:i]) } )
  cutoff <- max(c(0, which(all_below_threshold)))
  if (cutoff == 0) {
    selected_ind <- c()
  } else {
    selected_ind <- sort.int(p_values, index.return = TRUE)$ix[1:cutoff]
  }

  return(selected_ind)
}

#' Covariance matrix estimation
#'
#' Covariance matrix from the asymptotic distribution derived in Gossmann et. al. (2018)
#'
#' Get the covariance matrix defined in Theorem 1 of Gossmann et. al. (2018)
#' "FDR-Corrected Sparse Canonical Correlation Analysis with Applications to Imaging Genomics"
#' (https://arxiv.org/pdf/1705.04312).
#'
#' @export
get_Omega <- function(u0, S_X1, S_Y1, S_Y1X1) {
  Omega1 <- matrix(NA, ncol(S_Y1), ncol(S_Y1))

  for (i in 1:ncol(S_Y1)) {
    for (j in 1:ncol(S_Y1)) {
      c1 <- crossprod(S_Y1X1[i, ], u0)
      c2 <- crossprod(S_Y1X1[j, ], u0)
      c3 <- S_Y1[i, j] * t(u0) %*% S_X1 %*% u0
      Omega1[i, j] <- c1 * c2 + c3
    }
  }

  return(Omega1)
}
