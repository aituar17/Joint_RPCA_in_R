#install all necessary packages
install.packages("Matrix")
install.packages("testthat")
install.packages("vegan")

#load all necessary packages
library(Matrix)
library(testthat)
library(vegan)

#DEFINING THE FUNCTION

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

#TESTING THE FUNCTION

set.seed(42)

#generate a synthetic table: 10 features × 8 samples
features <- paste0("F", 1:10)
samples  <- paste0("Sample", 1:8)
mat <- matrix(rpois(80, lambda = 20), nrow = 10,
              dimnames = list(features, samples))

#introduce some low-abundance features and samples
mat[1:2, ] <- 0                     #feature F1, F2: zero counts
mat[3, 1:6] <- 0                    #feature F3: low frequency
mat[4, ] <- c(rep(0, 7), 200)       #feature F4: very sparse
mat[, 8] <- 0                       #sample 8: empty

#add a duplicate sample name to trigger an error
colnames(mat)[7] <- "Sample6"       #now Sample6 appears twice

#inspect raw matrix
print(mat)

#try rpca_table_processing (expecting an error due to duplicate IDs)
processed <- tryCatch({
  rpca_table_processing(mat,
                        min_sample_count = 30,
                        min_feature_count = 50,
                        min_feature_frequency = 20)
}, error = function(e) e)

#inspect the result
print(processed)

#fix duplicate sample name
colnames(mat) <- paste0("Sample", 1:8)

#optionally inspect to confirm
print(colnames(mat))         #should all be unique
any(duplicated(colnames(mat)))  #should be FALSE

processed <- rpca_table_processing(
  mat,
  min_sample_count = 30,
  min_feature_count = 50,
  min_feature_frequency = 20
)

print(processed)
