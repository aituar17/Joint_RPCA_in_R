#' OptSpace-Based Dimensionality Reduction and RPCA Biplot Generation
#'
#' Internal function that fits an OptSpace model to a rCLR-transformed compositional table,
#' reconstructs the low-rank matrix, applies PCA, and constructs an ordination result
#' capturing sample embeddings, feature loadings, and explained variance. A distance matrix
#' is also generated using Aitchison geometry.
#'
#' @param rclr_table A numeric matrix representing rCLR-transformed compositional data.
#' @param feature_ids Character vector of feature names (used for row labeling of loadings).
#' @param subject_ids Character vector of sample names (used for row labeling of embeddings).
#' @param n_components Integer specifying number of principal components to retain. Default is 3.
#' @param max_iterations Maximum number of iterations to run OptSpace optimization. Default is 5.
#'
#' @return A list with:
#' \describe{
#'   \item{ord_res}{An \code{OrdinationResults} object containing PCA scores, loadings, and metadata.}
#'   \item{dist}{A sample-by-sample \code{DistanceMatrix} object using Aitchison geometry.}
#'   \item{opt_fit}{The raw OptSpace fit result containing matrices \code{X}, \code{Y}, and \code{S}.}
#' }
#'
#' @keywords internal

.optspace_helper <- function(rclr_table,
                             feature_ids,
                             subject_ids,
                             n_components = 3,
                             max_iterations = 5) {
  
  #fit OptSpace
  opt_result <- .optspace(rclr_table, ropt = n_components, niter = max_iterations, tol = 1e-5, verbose = FALSE)
  
  #update n_components
  n_components <- ncol(opt_result$S)
  
  #reconstruct and re-center matrix
  X_hat <- opt_result$X %*% opt_result$S %*% t(opt_result$Y)
  X_hat <- scale(X_hat, center = TRUE, scale = FALSE)
  X_hat <- t(scale(t(X_hat), center = TRUE, scale = FALSE))
  
  #PCA
  svd_out <- svd(X_hat)
  u <- svd_out$u[, 1:n_components, drop = FALSE]
  s <- svd_out$d[1:n_components]
  v <- svd_out$v[, 1:n_components, drop = FALSE]
  
  #label loadings
  rename_cols <- paste0("PC", seq_len(n_components))
  sample_scores <- data.frame(u, row.names = subject_ids)
  feature_scores <- data.frame(v, row.names = feature_ids)
  colnames(sample_scores) <- rename_cols
  colnames(feature_scores) <- rename_cols
  
  #proportion explained
  prop_var <- s^2 / sum(svd_out$d^2)
  names(prop_var) <- rename_cols
  names(s) <- rename_cols
  
  #add PC3 for 2D case
  if (n_components == 2) {
    sample_scores$PC3 <- 0
    feature_scores$PC3 <- 0
    s <- c(s, PC3 = 0)
    prop_var <- c(prop_var, PC3 = 0)
    rename_cols <- c(rename_cols, "PC3")
  }
  
  #compute distance
  dist_matrix_raw <- as.matrix(dist(u))
  rownames(dist_matrix_raw) <- subject_ids
  colnames(dist_matrix_raw) <- subject_ids
  
  #wrap with DistanceMatrix
  dist_res <- .DistanceMatrix(dist_matrix_raw, ids = subject_ids, method = "aitchison")
  
  #build OrdinationResults object
  ord_res <- .OrdinationResults(
    method = "rpca_biplot",
    eigvals = s,
    samples = sample_scores,
    features = feature_scores,
    proportion_explained = prop_var,
    dist = dist_matrix_raw,  
    metadata = list(
      long_method_name = "(Robust Aitchison) RPCA Biplot",
      run_id = sprintf("optspace_helper_n_components_%d.max_iterations_%d", 
                       n_components, max_iterations)
    )
  )
  
  return(list(
    ord_res = ord_res,
    dist = dist_res,  
    opt_fit = opt_result
  ))
}