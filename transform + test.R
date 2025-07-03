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
  table_rclr_project <- table_rclr_project[rownames(Vdf), ]
  M_project <- t(as.matrix(table_rclr_project))
  M_project <- sweep(M_project, 1, rowMeans(M_project, na.rm = TRUE))
  M_project <- sweep(M_project, 2, colMeans(M_project, na.rm = TRUE))
  U_projected <- t(M_project %*% as.matrix(Vdf))
  U_projected <- t(U_projected) / sqrt(sum(s_eig^2))
  rownames(U_projected) <- colnames(table_rclr_project)
  colnames(U_projected) <- colnames(Udf)
  rbind(Udf, U_projected)
}

#main transform function
transform <- function(ordination, tables,
                          subset_tables = TRUE,
                          apply_rclr = TRUE) {
  
  Udf <- ordination$samples
  Vdf <- ordination$features
  s_eig <- ordination$eigvals
  
  #apply rclr transformation if needed
  rclr_tables <- lapply(tables, function(tab) {
    mat <- as.matrix(tab)
    if (apply_rclr) {
      rclr_transform(mat)
    } else {
      mat
    }
  })
  
  #combine all tables row-wise
  rclr_combined <- do.call(rbind, rclr_tables)
  
  shared_features <- intersect(rownames(rclr_combined), rownames(Vdf))
  
  if (length(shared_features) < nrow(Vdf)) {
    stop("The input tables do not contain all the features in the ordination.")
  } else if (subset_tables) {
    unshared_N <- nrow(rclr_combined) - length(shared_features)
    if (unshared_N > 0) {
      warning(sprintf("Removing %d feature(s) in table(s) but not the ordination.", unshared_N))
    }
    rclr_combined <- rclr_combined[shared_features, ]
  } else {
    stop("Features in the input tables do not match the ordination. Set subset_tables = TRUE to proceed.")
  }
  
  ordination$samples <- transform_helper(Udf, Vdf, s_eig, rclr_combined)
  return(ordination)
}

#TESTING THE TRANSFORM FUNCTION

#define toy ordination object (mock eigenvectors, scores, loadings)
eigvals <- c(2.5, 1.5)
samples <- data.frame(PC1 = c(0.1, 0.2), PC2 = c(-0.1, 0.3),
                      row.names = c("s1", "s2"))
features <- data.frame(PC1 = c(0.4, -0.3, 0.5), PC2 = c(0.2, 0.6, -0.1),
                       row.names = c("f1", "f2", "f3"))
prop_exp <- eigvals / sum(eigvals)

ordination <- OrdinationResults(method = "TestPCA",
                                eigvals = eigvals,
                                samples = samples,
                                features = features,
                                proportion_explained = prop_exp)

#create synthetic test table (as list of matrices)
table1 <- matrix(runif(3 * 3, 1, 10), nrow = 3,
                 dimnames = list(c("f1", "f2", "f3"),
                                 c("t1", "t2", "t3")))
table2 <- matrix(runif(3 * 3, 1, 10), nrow = 3,
                 dimnames = list(c("f1", "f2", "f3"),
                                 c("t4", "t5", "t6")))

#wrap in list
tables <- list(table1, table2)

#apply transformation
updated_ord <- transform(ordination, tables)

#inspect result
print(updated_ord)

summary(updated_ord)
