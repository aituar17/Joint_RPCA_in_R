#' Joint RPCA wrapper for MultiAssayExperiment
#'
#' @param mae A MultiAssayExperiment object containing compositional tables
#' @param assays Character vector of assay names to extract. If NULL, use all assays.
#' @param ... Additional arguments passed to jointRPCA()
#'
#' @return Output from jointRPCA()
#' @export
jointRPCAmae <- function(mae, assays = NULL, ...) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
        stop("Package 'MultiAssayExperiment' required for MAE input.")
    }
  
    if (is.null(assays)) assays <- names(experiments(mae))
  
    tables <- lapply(assays, function(a) {
        assay(mae, a)
    })
  
    jointRPCA(tables = tables, ...)
}
