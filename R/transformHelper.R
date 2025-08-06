#' Project New Data into Existing Ordination Space
#'
#' Internal function to align rCLR-transformed samples to an existing RPCA ordination space.
#' Handles feature alignment, deduplication of sample names, double-centering normalization,
#' and projection into low-rank space using previously learned components.
#'
#' @param Udf Matrix of training sample embeddings (samples × components).
#' @param Vdf Matrix of feature loadings (features × components).
#' @param s.eig Singular values from the RPCA decomposition.
#' @param table.rclr.project New rCLR-transformed table for projection (features × samples).
#'
#' @return A combined matrix of training and projected samples (samples × components).
#' @keywords internal

.transform_helper <- function(Udf, Vdf, s.eig, table.rclr.project) {
    #align features
    table.rclr.project <- table.rclr.project[rownames(Vdf), , drop = FALSE]
  
    #transpose to get samples as rows
    M.project <- t(as.matrix(table.rclr.project))
  
    #deduplicate sample names by stripping suffixes
    sample.names.clean <- sub("_\\d+$", "", rownames(M.project))
    M.project <- cbind(Sample = sample.names.clean, M.project)
    M.project <- aggregate(. ~ Sample, data = as.data.frame(M.project), FUN = function(x) {
        if (all(is.na(x))) NA else mean(as.numeric(x), na.rm = TRUE)
    })
    rownames(M.project) <- M.project$Sample
    M.project$Sample <- NULL
    M.project <- as.matrix(M.project)
  
    #double-centering
    M.project <- sweep(M.project, 1, rowMeans(M.project, na.rm = TRUE), "-")
    M.project <- sweep(M.project, 2, colMeans(M.project, na.rm = TRUE), "-")
  
    #projection
    U.projected <- M.project %*% as.matrix(Vdf)
    U.projected <- U.projected / sqrt(sum(s.eig^2))
    colnames(U.projected) <- colnames(Udf)
  
    #merge with training samples
    U.combined <- rbind(Udf[setdiff(rownames(Udf), rownames(U.projected)), , drop = FALSE],
                        U.projected)
    return(U.combined)
}