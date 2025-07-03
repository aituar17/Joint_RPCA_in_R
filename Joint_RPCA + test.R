#load all necessary packages
library(Matrix)
library(irlba)
library(RSpectra)
library(pracma)
library(optimx)
library(glmnet)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(factoextra)
library(vegan)
library(randomForest)
library(ROCR)
library(vegan)
library(ggforce)
library(concaveman)

#FUNCTION FOR RPCA TABLE PROCESSING

rpca_table_processing <- function(table,
                                  min_sample_count = 0,
                                  min_feature_count = 0,
                                  min_feature_frequency = 0) {
  #ensure the input is a matrix
  if (is.data.frame(table)) {
    table <- as.matrix(table)
  }
  
  n_features <- nrow(table)
  n_samples  <- ncol(table)
  
  #filter features by total count
  if (!is.null(min_feature_count)) {
    feature_totals <- rowSums(table, na.rm = TRUE)
    keep_features <- feature_totals > min_feature_count
    table <- table[keep_features, , drop = FALSE]
  }
  
  #filter features by frequency across samples
  if (!is.null(min_feature_frequency)) {
    freq_threshold <- min_feature_frequency / 100
    feature_freq <- rowMeans(table > 0, na.rm = TRUE)
    keep_features <- feature_freq > freq_threshold
    table <- table[keep_features, , drop = FALSE]
  }
  
  #filter samples by total count
  if (!is.null(min_sample_count)) {
    sample_totals <- colSums(table, na.rm = TRUE)
    keep_samples <- sample_totals > min_sample_count
    table <- table[, keep_samples, drop = FALSE]
  }
  
  #check for duplicate IDs
  if (any(duplicated(colnames(table)))) {
    stop("Data table contains duplicate sample (column) IDs.")
  }
  if (any(duplicated(rownames(table)))) {
    stop("Data table contains duplicate feature (row) IDs.")
  }
  
  #remove empty rows and columns if sample filtering applied
  if (!is.null(min_sample_count)) {
    nonzero_features <- rowSums(table, na.rm = TRUE) > 0
    nonzero_samples  <- colSums(table, na.rm = TRUE) > 0
    table <- table[nonzero_features, nonzero_samples, drop = FALSE]
  }
  
  return(table)
}

#FUNCTION FOR MASK VALUE ONLY

mask_value_only <- function(mat) {
  #ensure matrix is at least 2D
  if (is.vector(mat)) {
    mat <- matrix(mat, nrow = 1)
  }
  
  #ensure matrix is not more than 2D
  if (length(dim(mat)) > 2) {
    stop("Input matrix can only have two dimensions or less")
  }
  
  #generate logical mask: TRUE where values are missing
  mask <- !is.finite(mat)  
  
  #create masked matrix
  masked_mat <- mat
  masked_mat[!is.finite(mat)] <- NA
  
  #return as a masked matrix
  return(structure(list(
    data = masked_mat,
    mask = mask
  ), class = "MaskedMatrix"))
}

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
  
  #project test samples
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

#FUNCTION FOR JOINT RPCA

joint_rpca <- function(tables,
                       n_test_samples = 10,
                       sample_metadata = NULL,
                       train_test_column = NULL,
                       n_components = 3,
                       rclr_transform_tables = TRUE,
                       min_sample_count = 0,
                       min_feature_count = 0,
                       min_feature_frequency = 0,
                       max_iterations = 5) {
  
  if (n_components < 2) stop("n_components must be at least 2.")
  if (max_iterations < 1) stop("max_iterations must be at least 1.")
  
  #filter each table
  if (rclr_transform_tables) {
    tables <- lapply(tables, function(tbl) {
      rpca_table_processing(tbl,
                            min_sample_count = min_sample_count,
                            min_feature_count = min_feature_count,
                            min_feature_frequency = min_feature_frequency)
    })
  }
  
  #find shared samples
  sample_sets <- lapply(tables, colnames)
  shared_all_samples <- Reduce(intersect, sample_sets)
  if (length(shared_all_samples) == 0) {
    stop("No samples overlap between all tables. If using pre-transformed tables, set rclr_transform_tables = FALSE.")
  }
  unshared_samples <- setdiff(unique(unlist(sample_sets)), shared_all_samples)
  if (length(unshared_samples) > 0) {
    warning(sprintf("Removing %d sample(s) that do not overlap in tables.", length(unshared_samples)))
  }
  
  #re-filter tables to shared samples
  tables <- lapply(tables, function(tbl) {
    tbl <- tbl[, shared_all_samples, drop = FALSE]
    if (rclr_transform_tables) {
      rpca_table_processing(tbl,
                            min_sample_count = min_sample_count,
                            min_feature_count = min_feature_count,
                            min_feature_frequency = min_feature_frequency)
    } else {
      tbl
    }
  })
  shared_all_samples <- Reduce(intersect, lapply(tables, colnames))
  
  #transform tables
  rclr_tables <- lapply(tables, function(tbl) {
    mat <- as.matrix(tbl)
    if (rclr_transform_tables) {
      rclr_transform(mat)
    } else {
      mask_value_only(mat)$data
    }
  })
  
  #determine train/test split
  if (!is.null(sample_metadata) && !is.null(train_test_column)) {
    md <- as.data.frame(sample_metadata)
    md <- md[shared_all_samples, , drop = FALSE]
    train_samples <- rownames(md)[md[[train_test_column]] == "train"]
    test_samples  <- rownames(md)[md[[train_test_column]] == "test"]
  } else {
    ord_tmp <- optspace_helper(
      rclr_table = t(rclr_tables[[1]]),
      feature_ids = rownames(rclr_tables[[1]]),
      subject_ids = colnames(rclr_tables[[1]]),
      n_components = n_components,
      max_iterations = max_iterations
    )$ord_res
    sorted_ids <- rownames(ord_tmp$samples[order(ord_tmp$samples[, 1]), ])
    idx <- round(seq(1, length(sorted_ids), length.out = n_test_samples))
    test_samples <- sorted_ids[idx]
    train_samples <- setdiff(shared_all_samples, test_samples)
  }
  
  #run joint RPCA via helper
  result <- joint_optspace_helper(
    tables = rclr_tables,
    n_components = n_components,
    max_iterations = max_iterations,
    test_samples = test_samples,
    train_samples = train_samples
  )
  
  return(result)
}


#TESTING THE JOINT RPCA FUNCTION

#TEST 1

#set up synthetic data
set.seed(123)

#create shared samples
samples <- paste0("Sample", 1:10)
features1 <- paste0("F1_", 1:8)
features2 <- paste0("F2_", 1:10)

#generate two synthetic tables
table1 <- matrix(rpois(80, lambda = 25), nrow = length(features1),
                 dimnames = list(features1, samples))
table2 <- matrix(rpois(100, lambda = 30), nrow = length(features2),
                 dimnames = list(features2, samples))

#introduce a few NAs
table1[sample(length(table1), 5)] <- NA
table2[sample(length(table2), 6)] <- NA

#sample metadata with train/test split
metadata <- data.frame(Set = c(rep("train", 7), rep("test", 3)),
                       row.names = samples)

#bundle the tables
tables <- list(table1, table2)

#run joint_rpca
result <- joint_rpca(
  tables = tables,
  n_components = 2,
  n_test_samples = 3,
  sample_metadata = metadata,
  train_test_column = "Set",
  rclr_transform_tables = TRUE,
  min_sample_count = 0,
  min_feature_count = 0,
  min_feature_frequency = 0,
  max_iterations = 5
)

#view results
print(result$ord_res$samples)
print(result$ord_res$features)
print(result$cv_stats)

#extract from ordination result
samples_df <- as.data.frame(result$ord_res$samples)
samples_df$Label <- rownames(samples_df)
test_samples <- c("Sample8", "Sample9", "Sample10") 
samples_df$Set <- ifelse(samples_df$Label %in% test_samples, "Test", "Train")

features_df <- as.data.frame(result$ord_res$features)
features_df$Label <- rownames(features_df)

#plot
ggplot(samples_df, aes(x = PC1, y = PC2, color = Set)) +
  geom_point(size = 3) +
  geom_text(aes(label = Label), vjust = -1.2) +
  theme_minimal() +
  scale_color_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(title = "Joint RPCA Sample Ordination", x = "PC1", y = "PC2")

ggplot(samples_df, aes(x = PC1, y = PC2, color = Set)) +
  geom_point(size = 3) +
  geom_text(aes(label = Label), vjust = -1.2) +
  geom_segment(data = features_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray50") +
  geom_text(data = features_df, aes(x = PC1, y = PC2, label = Label),
            color = "gray30", vjust = 1.2, size = 3) +
  theme_minimal() +
  coord_equal() +
  scale_color_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(title = "Joint RPCA Biplot", x = "PC1", y = "PC2")

#extract loadings matrix
loadings <- result$ord_res$features

#for each PC, rank features by absolute contribution
ranked_features <- lapply(colnames(loadings), function(pc) {
  df <- data.frame(Feature = rownames(loadings),
                   Loading = loadings[, pc],
                   AbsLoading = abs(loadings[, pc]))
  df <- df[order(-df$AbsLoading), ]
  rownames(df) <- NULL
  df
})
names(ranked_features) <- colnames(loadings)

head(ranked_features$PC1, 5)  #top 5 for PC1
head(ranked_features$PC2, 5)  #top 5 for PC2

top_PC1 <- head(ranked_features$PC1, 10)

ggplot(top_PC1, aes(x = reorder(Feature, AbsLoading), y = AbsLoading)) +
  geom_col(fill = "darkslateblue") +
  coord_flip() +
  labs(title = "Top Features Driving PC1",
       x = NULL, y = "Absolute Loading")

#TEST 2 (3 TABLES)

#set up synthetic data
set.seed(123)

#create shared samples
samples <- paste0("Sample", 1:10)
features1 <- paste0("F1_", 1:8)
features2 <- paste0("F2_", 1:10)
features3 <- paste0("F3_", 1:6)

#generate three synthetic tables
table1 <- matrix(rpois(80, lambda = 25), nrow = length(features1),
                 dimnames = list(features1, samples))
table2 <- matrix(rpois(100, lambda = 30), nrow = length(features2),
                 dimnames = list(features2, samples))
table3 <- matrix(rpois(60, lambda = 35), nrow = length(features3),
                 dimnames = list(features3, samples))

#introduce a few NAs
table1[sample(length(table1), 5)] <- NA
table2[sample(length(table2), 6)] <- NA
table3[sample(length(table3), 4)] <- NA

#sample metadata with train/test split
metadata <- data.frame(Set = c(rep("train", 7), rep("test", 3)),
                       row.names = samples)

#bundle the tables
tables <- list(table1, table2, table3)

#run joint_rpca
result <- joint_rpca(
  tables = tables,
  n_components = 2,
  n_test_samples = 3,
  sample_metadata = metadata,
  train_test_column = "Set",
  rclr_transform_tables = TRUE,
  min_sample_count = 0,
  min_feature_count = 0,
  min_feature_frequency = 0,
  max_iterations = 5
)

#view results
print(result$ord_res$samples)
print(result$ord_res$features)
print(result$cv_stats)

#extract from ordination result
samples_df <- as.data.frame(result$ord_res$samples)
samples_df$Label <- rownames(samples_df)
test_samples <- c("Sample8", "Sample9", "Sample10") 
samples_df$Set <- ifelse(samples_df$Label %in% test_samples, "Test", "Train")

features_df <- as.data.frame(result$ord_res$features)
features_df$Label <- rownames(features_df)

#plot
ggplot(samples_df, aes(x = PC1, y = PC2, color = Set)) +
  geom_point(size = 3) +
  geom_text(aes(label = Label), vjust = -1.2) +
  theme_minimal() +
  scale_color_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(title = "Joint RPCA Sample Ordination", x = "PC1", y = "PC2")

ggplot(samples_df, aes(x = PC1, y = PC2, color = Set)) +
  geom_point(size = 3) +
  geom_text(aes(label = Label), vjust = -1.2) +
  geom_segment(data = features_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray50") +
  geom_text(data = features_df, aes(x = PC1, y = PC2, label = Label),
            color = "gray30", vjust = 1.2, size = 3) +
  theme_minimal() +
  coord_equal() +
  scale_color_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(title = "Joint RPCA Biplot", x = "PC1", y = "PC2")

#extract loadings matrix
loadings <- result$ord_res$features

#for each PC, rank features by absolute contribution
ranked_features <- lapply(colnames(loadings), function(pc) {
  df <- data.frame(Feature = rownames(loadings),
                   Loading = loadings[, pc],
                   AbsLoading = abs(loadings[, pc]))
  df <- df[order(-df$AbsLoading), ]
  rownames(df) <- NULL
  df
})
names(ranked_features) <- colnames(loadings)

head(ranked_features$PC1, 5)  #top 5 for PC1
head(ranked_features$PC2, 5)  #top 5 for PC2

top_PC1 <- head(ranked_features$PC1, 10)

ggplot(top_PC1, aes(x = reorder(Feature, AbsLoading), y = AbsLoading)) +
  geom_col(fill = "darkslateblue") +
  coord_flip() +
  labs(title = "Top Features Driving PC1",
       x = NULL, y = "Absolute Loading")

#COVARIANCE MATRIX

numeric_loadings <- as.matrix(features_df[, c("PC1", "PC2")])

#cross-product of feature loadings
cov_matrix <- tcrossprod(numeric_loadings)

print(cov_matrix)