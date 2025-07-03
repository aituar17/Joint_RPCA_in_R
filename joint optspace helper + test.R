#install all necessary packages
install.packages("Matrix")
install.packages("testthat")
install.packages("vegan")
install.packages("ggplot2")

#load all necessary packages
library(Matrix)
library(testthat)
library(vegan)
library(ggplot2)

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

#TRANSFORM FUNCTION

#rclr_transform function
rclr_transform <- function(df) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for rclr_transform.")
  }
  vegan::decostand(df, method = "rclr")
}

#helper to normalize and project new data into ordination space
transform_helper <- function(Udf, Vdf, s_eig, table_rclr_project) {
  #align features
  table_rclr_project <- table_rclr_project[rownames(Vdf), , drop = FALSE]
  
  #transpose to get samples as rows
  M_project <- t(as.matrix(table_rclr_project))
  
  #deduplicate sample names by stripping suffixes
  sample_names_clean <- sub("_\\d+$", "", rownames(M_project))
  M_project <- cbind(Sample = sample_names_clean, M_project)
  M_project <- aggregate(. ~ Sample, data = as.data.frame(M_project), FUN = function(x) {
    if (all(is.na(x))) NA else mean(as.numeric(x), na.rm = TRUE)
  })
  rownames(M_project) <- M_project$Sample
  M_project$Sample <- NULL
  M_project <- as.matrix(M_project)
  
  #double-centering
  M_project <- sweep(M_project, 1, rowMeans(M_project, na.rm = TRUE), "-")
  M_project <- sweep(M_project, 2, colMeans(M_project, na.rm = TRUE), "-")
  
  #projection
  U_projected <- M_project %*% as.matrix(Vdf)
  U_projected <- U_projected / sqrt(sum(s_eig^2))
  colnames(U_projected) <- colnames(Udf)
  
  #merge with training samples
  U_combined <- rbind(Udf[setdiff(rownames(Udf), rownames(U_projected)), , drop = FALSE],
                      U_projected)
  return(U_combined)
}

#main transform function
transform <- function(ordination, tables,
                      subset_tables = TRUE,
                      apply_rclr = TRUE) {
  
  Udf <- ordination$samples
  project_ids <- unique(unlist(lapply(tables, colnames)))
  Udf <- Udf[!rownames(Udf) %in% project_ids, , drop = FALSE]
  Vdf <- ordination$features
  s_eig <- ordination$eigvals
  all_features <- rownames(Vdf)
  
  #apply rclr transformation if needed
  rclr_tables <- lapply(tables, function(tab) {
    mat <- as.matrix(tab)
    if (apply_rclr) {
      rclr_transform(mat)
    } else {
      mat
    }
  })
  
  #pad each table to match training feature set
  rclr_tables <- lapply(rclr_tables, function(mat) {
    missing_feats <- setdiff(all_features, rownames(mat))
    if (length(missing_feats) > 0) {
      pad <- matrix(0, nrow = length(missing_feats), ncol = ncol(mat),
                    dimnames = list(missing_feats, colnames(mat)))
      mat <- rbind(mat, pad)
    }
    mat[all_features, , drop = FALSE]  
  })
  
  #combine all tables column-wise (samples)
  rclr_combined <- do.call(cbind, rclr_tables)
  colnames(rclr_combined) <- make.unique(colnames(rclr_combined), sep = "_")
  
  if (subset_tables) {
    unshared_N <- nrow(rclr_combined) - length(all_features)
    if (unshared_N > 0) {
      warning(sprintf("Removing %d feature(s) in table(s) but not the ordination.", unshared_N))
    }
    rclr_combined <- rclr_combined[all_features, , drop = FALSE]
  } else {
    stop("Features in the input tables do not match the ordination. Set subset_tables = TRUE to proceed.")
  }
  
  ordination$samples <- transform_helper(Udf, Vdf, s_eig, rclr_combined)
  return(ordination)
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

#JOINT SOLVE FUNCTION

joint_optspace_solve <- function(train_test_pairs, n_components,
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
  fit <- optspace(train_stacked, ropt = n_components, niter = max_iter, tol = 1e-5, verbose = FALSE)
  
  #extract sample loadings 
  U_shared <- fit$X
  S_shared <- fit$S
  
  #split V (feature loadings) back into per-table pieces
  feat_indices <- cumsum(sapply(dims, function(d) d[2]))
  feat_starts <- c(1, head(feat_indices, -1) + 1)
  V_list <- Map(function(start, end) fit$Y[start:end, , drop = FALSE],
                feat_starts, feat_indices)
  
  #compute test reconstruction error
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

#FUNCTION FOR JOINT OPTSPACE HELPER

joint_optspace_helper <- function(tables,
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
  opt_result <- joint_optspace_solve(tables_for_solver,
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
  ord_res <- OrdinationResults(
    method = "rpca",
    eigvals = setNames(s_eig, pc_names),
    samples = u,
    features = v,
    proportion_explained = setNames(prop_exp, pc_names)
  )
  
  #project test samples if provided
  if (length(test_samples) > 0) {
    test_matrices <- lapply(tables, function(tbl) {
      tbl[, test_samples, drop = FALSE]
    })
    ord_res <- transform(ord_res, test_matrices, apply_rclr = FALSE)
  }
  
  #compute distance matrix and CV error summary
  dist_mat <- as.matrix(dist(ord_res$samples))
  dist_res <- DistanceMatrix(dist_mat, ids = rownames(ord_res$samples))
  
  cv_dist <- data.frame(t(dists))
  colnames(cv_dist) <- c("mean_CV", "std_CV")
  cv_dist$run <- sprintf("tables_%d.n_components_%d.max_iterations_%d.n_test_%d",
                         length(tables), n_components, max_iterations, length(test_samples))
  cv_dist$iteration <- seq_len(nrow(cv_dist))
  rownames(cv_dist) <- seq_len(nrow(cv_dist))
  
  list(ord_res = ord_res, dist = dist_res, cv_stats = cv_dist)
}

#TESTING THE JOINT OPTSPACE HELPER FUNCTION

#generate two synthetic tables with overlapping samples
set.seed(77)
samples <- paste0("Sample", 1:6)
features1 <- paste0("F1_", 1:4)
features2 <- paste0("F2_", 1:5)

table1 <- matrix(rnorm(6 * 4), nrow = 6, dimnames = list(samples, features1))
table2 <- matrix(rnorm(6 * 5), nrow = 6, dimnames = list(samples, features2))

#introduce light missingness
table1[sample(length(table1), 2)] <- NA
table2[sample(length(table2), 3)] <- NA

tables <- list(t(table1), t(table2))

#split into training and test sets
train_samples <- samples[1:4]
test_samples  <- samples[5:6]

#run joint optspace helper
result <- joint_optspace_helper(
  tables = tables,
  n_components = 2,
  max_iterations = 10,
  test_samples = test_samples,
  train_samples = train_samples
)

#inspect results
summary(result$ord_res)
summary(result$dist)
print(result$cv_stats)

#extract ordination results
ord <- result$ord_res

#prepare sample data for plotting
sample_df <- as.data.frame(ord$samples)
sample_df$Label <- rownames(sample_df)
sample_df$Set <- ifelse(rownames(sample_df) %in% c(test_samples), "Test", "Train")

#prepare feature vectors
feature_df <- as.data.frame(ord$features)
feature_df$Label <- rownames(feature_df)

#build ggplot
ggplot() +
  geom_point(data = sample_df, aes(x = PC1, y = PC2, color = Set), size = 3) +
  geom_text(data = sample_df, aes(x = PC1, y = PC2, label = Label), vjust = -1.2, size = 3) +
  geom_segment(data = feature_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "gray50") +
  geom_text(data = feature_df, aes(x = PC1, y = PC2, label = Label), 
            color = "gray30", vjust = 1.2, size = 3) +
  theme_minimal() +
  labs(title = "Joint RPCA Ordination", x = "PC1", y = "PC2") +
  coord_equal()
