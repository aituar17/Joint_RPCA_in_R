#' Apply Projection of New Compositional Tables to Existing Ordination
#'
#' Internal function that transforms and projects new sample tables into an existing Joint RPCA ordination space.
#' It handles feature alignment, optional rCLR preprocessing, padding of missing features, and merges projected samples with existing ones.
#'
#' @param ordination A list containing previous ordination results: `samples`, `features`, and `eigvals`.
#' @param tables A list of new compositional tables (matrices or data frames) with features as rows and samples as columns.
#' @param subset.tables Logical; if `TRUE`, removes unshared features in new tables not present in training. If `FALSE`, stops with an error if features mismatch.
#' @param apply.rclr Logical; whether to apply rCLR transformation to new input tables before projection. Default is `TRUE`.
#'
#' @return An updated ordination list with the `samples` matrix extended to include new projected samples.
#' @keywords internal

.transform <- function(ordination, tables,
                       subset.tables = TRUE,
                       apply.rclr = TRUE) {
  
    Udf <- ordination$samples
    project.ids <- unique(unlist(lapply(tables, colnames)))
    Udf <- Udf[!rownames(Udf) %in% project.ids, , drop = FALSE]
    Vdf <- ordination$features
    s.eig <- ordination$eigvals
    all.features <- rownames(Vdf)
  
    #apply rclr transformation if needed
    rclr.tables <- lapply(tables, function(tab) {
        mat <- as.matrix(tab)
        if (apply.rclr) {
            vegan::decostand(mat, method = "rclr")
        } else {
            mat
        }
    })
  
    #pad each table to match training feature set
    rclr.tables <- lapply(rclr.tables, function(mat) {
        missing.feats <- setdiff(all.features, rownames(mat))
        if (length(missing.feats) > 0) {
            pad <- matrix(0, nrow = length(missing.feats), ncol = ncol(mat),
                        dimnames = list(missing.feats, colnames(mat)))
            mat <- rbind(mat, pad)
        }
        mat[all.features, , drop = FALSE]  
    })
  
    #combine all tables column-wise (samples)
    rclr.combined <- do.call(cbind, rclr.tables)
    colnames(rclr.combined) <- make.unique(colnames(rclr.combined), sep = "_")
  
    if (subset.tables) {
        unshared.N <- nrow(rclr.combined) - length(all.features)
        if (unshared.N > 0) {
            warning(sprintf("Removing %d feature(s) in table(s) but not the ordination.", unshared.N))
        }
        rclr.combined <- rclr.combined[all.features, , drop = FALSE]
    } else {
        stop("Features in the input tables do not match the ordination. Set subset.tables = TRUE to proceed.")
    }
  
    ordination$samples <- .transform_helper(Udf, Vdf, s.eig, rclr.combined)
    return(ordination)
}