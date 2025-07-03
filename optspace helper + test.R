#install all necessary packages
install.packages("Matrix")
install.packages("testthat")
install.packages("vegan")

#load all necessary packages
library(Matrix)
library(testthat)
library(vegan)

#ORDINATION_RESULTS FUNCTION

#constructor function
OrdinationResults <- function(method, eigvals, samples, features,
                              proportion_explained, dist = NULL, metadata = list()) {
  structure(list(
    method = method,
    eigvals = eigvals,
    samples = samples,
    features = features,
    proportion_explained = proportion_explained,
    dist = dist,
    metadata = metadata
  ), class = "OrdinationResults")
}

#print method
print.OrdinationResults <- function(x, ...) {
  cat("OrdinationResults (method:", x$method, ")\n")
  cat("Number of components:", length(x$eigvals), "\n")
  cat("Variance explained:\n")
  print(round(x$proportion_explained, 3))
  invisible(x)
}

#summary method
summary.OrdinationResults <- function(object, ...) {
  print(object)
  cat("\nSample scores (first few rows):\n")
  print(head(object$samples))
  cat("\nFeature loadings (first few rows):\n")
  print(head(object$features))
  invisible(object)
}

#plot method
plot.OrdinationResults <- function(x, comps = c(1, 2), ...) {
  if (length(comps) != 2) stop("Please select two components to plot.")
  plot(x$samples[, comps], col = "blue", pch = 19,
       xlab = paste0("PC", comps[1]),
       ylab = paste0("PC", comps[2]),
       main = paste("Ordination (", x$method, ")", sep = ""))
  points(x$features[, comps], col = "red", pch = 4)
  legend("topright", legend = c("Samples", "Features"),
         col = c("blue", "red"), pch = c(19, 4))
}

#RCLR TRANSFORM FUNCTION

rclr_transform <- function(df) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for rclr_transform.")
  }
  vegan::decostand(df, method = "rclr")
}

#DISTANCE MATRIX FUNCTION

DistanceMatrix <- function(matrix, ids = NULL, method = "euclidean") {
  if (!is.matrix(matrix)) stop("Input must be a matrix.")
  if (!isSymmetric(matrix)) stop("Distance matrix must be symmetric.")
  if (!is.null(ids)) {
    if (length(ids) != nrow(matrix)) stop("Length of 'ids' must match matrix dimensions.")
    rownames(matrix) <- ids
    colnames(matrix) <- ids
  }
  structure(list(
    data = matrix,
    ids = rownames(matrix),
    method = method
  ), class = "DistanceMatrix")
}

print.DistanceMatrix <- function(x, ...) {
  cat("DistanceMatrix (", x$method, ")\n", sep = "")
  cat("Number of objects:", length(x$ids), "\n")
  print(head(x$data, 6))  # Show only top part
  invisible(x)
}

summary.DistanceMatrix <- function(object, ...) {
  cat("Summary of DistanceMatrix\n")
  cat("Method:", object$method, "\n")
  cat("Size:", nrow(object$data), "x", ncol(object$data), "\n")
  cat("IDs:\n")
  print(head(object$ids, 6))
  cat("\nDistance Summary Stats:\n")
  print(summary(as.vector(object$data[upper.tri(object$data)])))
  invisible(object)
}

#FUNCTION FOR OPTSPACE HELPER

optspace_helper <- function(rclr_table,
                            feature_ids,
                            subject_ids,
                            n_components = 3,
                            max_iterations = 5) {
  
  #fit OptSpace
  opt_result <- optspace(rclr_table, ropt = n_components, niter = max_iterations, tol = 1e-5, verbose = FALSE)
  
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
  dist_res <- DistanceMatrix(dist_matrix_raw, ids = subject_ids, method = "aitchison")
  
  #build OrdinationResults object
  ord_res <- OrdinationResults(
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

#TESTING THE OPTSPACE HELPER FUNCTION

#set seed
set.seed(123)

#generate synthetic data
n_features <- 20
n_samples <- 10
table <- matrix(rlnorm(n_features * n_samples, meanlog = 1), nrow = n_features)
rownames(table) <- paste0("Feature", 1:n_features)
colnames(table) <- paste0("Sample", 1:n_samples)

#introduce some missing values
table[sample(length(table), 10)] <- NA

#apply rclr transform
table_rclr <- rclr_transform(table)

#define IDs
feature_ids <- rownames(table)
subject_ids <- colnames(table)

#run optspace_helper
result <- optspace_helper(
  rclr_table = t(table_rclr), 
  feature_ids = feature_ids,
  subject_ids = subject_ids,
  n_components = 3,
  max_iterations = 50
)

#view ordination summary
print(result$ord_res)
plot(result$ord_res)
