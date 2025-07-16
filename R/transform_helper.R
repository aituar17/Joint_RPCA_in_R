#' Project New Data into Existing Ordination Space
#'
#' Internal function to align rCLR-transformed samples to an existing RPCA ordination space.
#' Handles feature alignment, deduplication of sample names, double-centering normalization,
#' and projection into low-rank space using previously learned components.
#'
#' @param Udf Matrix of training sample embeddings (samples × components).
#' @param Vdf Matrix of feature loadings (features × components).
#' @param s_eig Singular values from the RPCA decomposition.
#' @param table_rclr_project New rCLR-transformed table for projection (features × samples).
#'
#' @return A combined matrix of training and projected samples (samples × components).
#' @keywords internal

.transform_helper <- function(Udf, Vdf, s_eig, table_rclr_project) {
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