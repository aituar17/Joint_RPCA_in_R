#' Joint RPCA Ordination Across Multiple Compositional Tables
#'
#' Internal function that performs Robust PCA via joint OptSpace decomposition across multiple compositional tables.
#' It splits each table into train/test sets, applies joint factorization, reconstructs sample and feature embeddings,
#' optionally projects test samples, computes a sample distance matrix, and returns cross-validation error statistics.
#'
#' @param tables A list of compositional matrices or data frames with features as rows and samples as columns.
#' @param n_components Number of principal components to compute.
#' @param max_iterations Maximum number of optimization iterations for OptSpace.
#' @param test_samples Character vector of sample IDs to be projected into the ordination space.
#' @param train_samples Character vector of sample IDs used to fit the ordination.
#'
#' @return A list with:
#' \describe{
#'   \item{ord_res}{An \code{OrdinationResults} object containing embeddings, loadings, and variance explained.}
#'   \item{dist}{A \code{DistanceMatrix} object for sample embeddings.}
#'   \item{cv_stats}{A data frame summarizing reconstruction error across iterations and tables.}
#' }
#'
#' @keywords internal

.joint_optspace_helper <- function(tables,
                                   n_components,
                                   max_iterations,
                                   test_samples,
                                   train_samples) {
  #split and transpose training/test data per table
  tables_split <- lapply(tables, function(tbl) {
    list(t(tbl[, test_samples, drop = FALSE]),
         t(tbl[, train_samples, drop = FALSE]))
  })
  
  #format input for solver
  tables_for_solver <- lapply(tables_split, function(pair) {
    lapply(pair, as.matrix)
  })
  
  #run joint OptSpace solver
  opt_result <- .joint_optspace_solve(tables_for_solver,
                                      n_components = n_components,
                                      max_iter = max_iterations)
  
  U <- opt_result$U
  S <- opt_result$S
  V_list <- opt_result$V_list
  dists <- opt_result$dists
  
  #assign row/column names to loadings
  pc_names <- paste0("PC", seq_len(n_components))
  
  #combine feature loadings with table-derived row names
  vjoint <- do.call(rbind, Map(function(tbl, V) {
    rownames(V) <- rownames(tbl)
    colnames(V) <- pc_names
    V
  }, tables, V_list))
  
  U <- U[seq_along(train_samples), , drop = FALSE]
  rownames(U) <- train_samples
  colnames(U) <- pc_names
  
  #recenter & re-factor via SVD
  X <- U %*% S %*% t(vjoint)
  X <- sweep(X, 2, colMeans(X))
  X <- sweep(X, 1, rowMeans(X))
  svd_res <- svd(X)
  u <- svd_res$u[, seq_len(n_components), drop = FALSE]
  v <- svd_res$v[, seq_len(n_components), drop = FALSE]
  s_eig <- svd_res$d[seq_len(n_components)]
  
  rownames(u) <- train_samples
  rownames(v) <- rownames(vjoint)
  colnames(u) <- colnames(v) <- pc_names
  
  #create ordination object
  prop_exp <- s_eig^2 / sum(s_eig^2)
  ord_res <- .OrdinationResults(
    method = "rpca",
    eigvals = setNames(s_eig, pc_names),
    samples = u,
    features = v,
    proportion_explained = setNames(prop_exp, pc_names)
  )
  
  #project test samples
  if (length(test_samples) > 0) {
    test_matrices <- lapply(tables, function(tbl) {
      tbl[, test_samples, drop = FALSE]
    })
    ord_res <- .transform(ord_res, test_matrices, apply_rclr = FALSE)
  }
  
  #compute distance matrix and CV error summary
  dist_mat <- as.matrix(dist(ord_res$samples))
  dist_res <- .DistanceMatrix(dist_mat, ids = rownames(ord_res$samples))
  
  cv_dist <- data.frame(t(dists))
  colnames(cv_dist) <- c("mean_CV", "std_CV")
  cv_dist$run <- sprintf("tables_%d.n_components_%d.max_iterations_%d.n_test_%d",
                         length(tables), n_components, max_iterations, length(test_samples))
  cv_dist$iteration <- seq_len(nrow(cv_dist))
  rownames(cv_dist) <- seq_len(nrow(cv_dist))
  
  list(ord_res = ord_res, dist = dist_res, cv_stats = cv_dist)
}