#' Apply Projection of New Compositional Tables to Existing Ordination
#'
#' Internal function that transforms and projects new sample tables into an existing Joint RPCA ordination space.
#' It handles feature alignment, optional rCLR preprocessing, padding of missing features, and merges projected samples with existing ones.
#'
#' @param ordination A list containing previous ordination results: `samples`, `features`, and `eigvals`.
#' @param tables A list of new compositional tables (matrices or data frames) with features as rows and samples as columns.
#' @param subset_tables Logical; if `TRUE`, removes unshared features in new tables not present in training. If `FALSE`, stops with an error if features mismatch.
#' @param apply_rclr Logical; whether to apply rCLR transformation to new input tables before projection. Default is `TRUE`.
#'
#' @return An updated ordination list with the `samples` matrix extended to include new projected samples.
#' @keywords internal

.transform <- function(ordination, tables,
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
      .rclr_transform(mat)
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
  
  ordination$samples <- .transform_helper(Udf, Vdf, s_eig, rclr_combined)
  return(ordination)
}