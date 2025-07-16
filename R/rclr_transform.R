#' Robust Centered Log Ratio Transformation
#'
#' Internal wrapper that applies robust centered log ratio (rclr) transformation
#' to compositional data using the \code{vegan::decostand()} function.
#'
#' @param df A numeric matrix or data frame with features as rows and samples as columns.
#'
#' @return Transformed matrix with rclr-applied values.
#'
#' @importFrom vegan decostand
#' @keywords internal

.rclr_transform <- function(df) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for rclr_transform.")
  }
  vegan::decostand(df, method = "rclr")
}