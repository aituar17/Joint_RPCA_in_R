#' RPCA Table Filtering and Preprocessing
#'
#' Internal function that performs filtering and cleanup on a compositional data table
#' prior to Robust PCA analysis. Removes low-count features/samples, enforces non-zero frequency thresholds,
#' checks for ID duplication, and returns a matrix suitable for transformation and ordination.
#'
#' @param table A matrix or data frame with features as rows and samples as columns.
#' @param min_sample_count Minimum total count required for a sample to be retained. Default is 0.
#' @param min_feature_count Minimum total count required for a feature to be retained. Default is 0.
#' @param min_feature_frequency Minimum percentage (0–100) of samples in which a feature must be non-zero. Default is 0.
#'
#' @return A filtered numeric matrix containing non-empty features and samples.
#' @keywords internal

.rpca_table_processing <- function(table,
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