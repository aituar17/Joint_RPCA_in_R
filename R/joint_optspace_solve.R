#' Joint OptSpace Optimization Across Multiple Train/Test Splits
#'
#' Internal function that performs joint matrix factorization using OptSpace across a set of paired train/test compositional tables.
#' Stacks training matrices horizontally, applies low-rank optimization, splits feature loadings per table,
#' and evaluates projection error on test data via Frobenius norm.
#'
#' @param train_test_pairs A list of paired matrices where each element is a two-item list: \code{[[test, train]]}.
#' @param n_components Integer specifying number of components to retain in the OptSpace model.
#' @param max_iter Maximum number of optimization iterations. Default is 50.
#' @param verbose Logical; whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{U}{Shared sample embedding matrix across all input tables.}
#'   \item{S}{Singular values matrix from OptSpace decomposition.}
#'   \item{V_list}{List of per-table feature loading matrices.}
#'   \item{dists}{Matrix of reconstruction errors (rows: error type, columns: tables).}
#' }
#'
#' @keywords internal

.joint_optspace_solve <- function(train_test_pairs, n_components,
                                  max_iter = 50, verbose = TRUE) {
  #prepare lists to hold training matrices and dimensions
  train_matrices <- list()
  test_matrices <- list()
  dims <- list()
  
  for (pair in train_test_pairs) {
    test_mat <- pair[[1]]
    train_mat <- pair[[2]]
    train_matrices <- append(train_matrices, list(train_mat))
    test_matrices <- append(test_matrices, list(test_mat))
    dims <- append(dims, list(dim(train_mat)))
  }
  
  #stack training matrices horizontally 
  train_stacked <- do.call(cbind, train_matrices)
  
  #apply OptSpace to stacked matrix
  if (verbose) message("Running optspace() on stacked training data...")
  fit <- .optspace(train_stacked, ropt = n_components, niter = max_iter, tol = 1e-5, verbose = FALSE)
  
  #extract sample loadings
  U_shared <- fit$X
  S_shared <- fit$S
  
  #split V back into per-table pieces
  feat_indices <- cumsum(sapply(dims, function(d) d[2]))
  feat_starts <- c(1, head(feat_indices, -1) + 1)
  V_list <- Map(function(start, end) fit$Y[start:end, , drop = FALSE],
                feat_starts, feat_indices)
  
  dists <- matrix(0, nrow = 2, ncol = max_iter)
  for (i in seq_along(test_matrices)) {
    V_k <- V_list[[i]]
    test_mat <- test_matrices[[i]]
    
    #project test samples: U_test = test × V
    U_test <- as.matrix(test_mat) %*% V_k
    U_test <- sweep(U_test, 2, diag(S_shared), "/")
    recon_test <- U_test %*% S_shared %*% t(V_k)
    
    #center for consistency
    recon_test <- scale(recon_test, center = TRUE, scale = FALSE)
    recon_test <- t(scale(t(recon_test), center = TRUE, scale = FALSE))
    
    error <- test_mat - recon_test
    error[is.na(error)] <- 0
    error_val <- norm(error, "F") / sqrt(sum(!is.na(test_mat)))
    
    dists[1, i] <- error_val
    dists[2, i] <- 0  
  }
  
  list(U = U_shared, S = S_shared, V_list = V_list, dists = dists)
}