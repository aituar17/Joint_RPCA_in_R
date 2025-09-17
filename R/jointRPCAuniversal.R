#' Universal Joint RPCA Wrapper
#'
#' @param data Input object: MultiAssayExperiment, TreeSummarizedExperiment, SummarizedExperiment, list of matrices, or single matrix.
#' @param assays Character vector of assay names to extract (only used for MultiAssayExperiment). If NULL, use all assays.
#' @param ... Additional arguments passed to jointRPCA()
#'
#' @return Output from jointRPCA()
#' @export
jointRPCAuniversal <- function(data, assays = NULL, ...) {
    #MultiAssayExperiment
    if (inherits(data, "MultiAssayExperiment")) {
        if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
            stop("Package 'MultiAssayExperiment' required for MAE input.")
        }
        if (is.null(assays)) assays <- names(experiments(data))
        tables <- lapply(assays, function(a) assay(data, a))
    
    #TreeSummarizedExperiment or SummarizedExperiment
    } else if (inherits(data, "TreeSummarizedExperiment") || inherits(data, "SummarizedExperiment")) {
        tables <- list(assay(data))
    
    #list of matrices
    } else if (is.list(data) && all(sapply(data, is.matrix))) {
        tables <- data
    
    #single matrix
    } else if (is.matrix(data)) {
        tables <- list(data)
    
    #unsupported format
    } else {
        stop("Unsupported input type. Provide a MultiAssayExperiment, TreeSummarizedExperiment, SummarizedExperiment, list of matrices, or a single matrix.")
    }
  
    jointRPCA(tables = tables, ...)
}