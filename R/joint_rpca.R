#load all necessary packages
suppressPackageStartupMessages({
library(Matrix)
library(irlba)
library(RSpectra)
library(pracma)
library(optimx)
library(glmnet)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(factoextra)
library(vegan)
library(randomForest)
library(ROCR)
library(vegan)
library(ggforce)
library(concaveman)
})

#FUNCTION FOR JOINT RPCA

#' Joint Robust PCA on Multiple Compositional Tables
#'
#' Performs Joint Robust Principal Component Analysis (RPCA) using OptSpace on multiple compositional tables.
#' Automatically handles shared sample alignment, optional train/test splitting, and compositional preprocessing.
#'
#' @param tables A list of compositional data tables (matrices or data frames) with features as rows and samples as columns.
#' @param n_test_samples Integer specifying the number of samples to hold out for testing (only used if `sample_metadata` is NULL). Default is 10.
#' @param sample_metadata Optional data frame containing sample-level metadata. Must include a column indicating train/test labels.
#' @param train_test_column The name of the column in `sample_metadata` that defines training vs test samples.
#' @param n_components Integer specifying the number of principal components to compute. Must be at least 2.
#' @param rclr_transform_tables Logical; whether to apply rclr transformation to each input table before ordination. Default is TRUE.
#' @param min_sample_count Minimum total count required for a sample to be retained during filtering. Default is 0 (no filtering).
#' @param min_feature_count Minimum total count required for a feature to be retained. Default is 0 (no filtering).
#' @param min_feature_frequency Minimum percentage (0–100) of samples in which a feature must be non-zero to be retained. Default is 0.
#' @param max_iterations Maximum number of optimization iterations. Must be at least 1.
#'
#' @return A list with:
#' \describe{
#'   \item{ord_res}{An \code{OrdinationResults} object containing sample scores, feature loadings, and variance explained.}
#'   \item{U_dist_res}{A distance matrix (samples × samples) based on the learned sample embeddings.}
#'   \item{cv_dist}{A data frame of cross-validation error statistics across iterations.}
#' }
#'
#' @examples
#' \dontrun{
#' # See examples/joint_rpca_example.qmd for a full reproducible test case.
#' }
#'
#' @export

joint_rpca <- function(tables,
                       n_test_samples = 10,
                       sample_metadata = NULL,
                       train_test_column = NULL,
                       n_components = 3,
                       rclr_transform_tables = TRUE,
                       min_sample_count = 0,
                       min_feature_count = 0,
                       min_feature_frequency = 0,
                       max_iterations = 5) {
  
  if (n_components < 2) stop("n_components must be at least 2.")
  if (max_iterations < 1) stop("max_iterations must be at least 1.")
  
  #filter each table
  if (rclr_transform_tables) {
    tables <- lapply(tables, function(tbl) {
      .rpca_table_processing(tbl,
                            min_sample_count = min_sample_count,
                            min_feature_count = min_feature_count,
                            min_feature_frequency = min_feature_frequency)
    })
  }
  
  #find shared samples
  sample_sets <- lapply(tables, colnames)
  shared_all_samples <- Reduce(intersect, sample_sets)
  if (length(shared_all_samples) == 0) {
    stop("No samples overlap between all tables. If using pre-transformed tables, set rclr_transform_tables = FALSE.")
  }
  unshared_samples <- setdiff(unique(unlist(sample_sets)), shared_all_samples)
  if (length(unshared_samples) > 0) {
    warning(sprintf("Removing %d sample(s) that do not overlap in tables.", length(unshared_samples)))
  }
  
  #re-filter tables to shared samples
  tables <- lapply(tables, function(tbl) {
    tbl <- tbl[, shared_all_samples, drop = FALSE]
    if (rclr_transform_tables) {
      .rpca_table_processing(tbl,
                            min_sample_count = min_sample_count,
                            min_feature_count = min_feature_count,
                            min_feature_frequency = min_feature_frequency)
    } else {
      tbl
    }
  })
  shared_all_samples <- Reduce(intersect, lapply(tables, colnames))
  
  #transform tables
  rclr_tables <- lapply(tables, function(tbl) {
    mat <- as.matrix(tbl)
    if (rclr_transform_tables) {
      .rclr_transform(mat)
    } else {
      .mask_value_only(mat)$data
    }
  })
  
  #determine train/test split
  if (!is.null(sample_metadata) && !is.null(train_test_column)) {
    md <- as.data.frame(sample_metadata)
    md <- md[shared_all_samples, , drop = FALSE]
    train_samples <- rownames(md)[md[[train_test_column]] == "train"]
    test_samples  <- rownames(md)[md[[train_test_column]] == "test"]
  } else {
    ord_tmp <- .optspace_helper(
      rclr_table = t(rclr_tables[[1]]),
      feature_ids = rownames(rclr_tables[[1]]),
      subject_ids = colnames(rclr_tables[[1]]),
      n_components = n_components,
      max_iterations = max_iterations
    )$ord_res
    sorted_ids <- rownames(ord_tmp$samples[order(ord_tmp$samples[, 1]), ])
    idx <- round(seq(1, length(sorted_ids), length.out = n_test_samples))
    test_samples <- sorted_ids[idx]
    train_samples <- setdiff(shared_all_samples, test_samples)
  }
  
  #run joint RPCA via helper
  result <- .joint_optspace_helper(
    tables = rclr_tables,
    n_components = n_components,
    max_iterations = max_iterations,
    test_samples = test_samples,
    train_samples = train_samples
  )
  
  return(result)
}