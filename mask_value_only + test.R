#install all necessary packages
install.packages("Matrix")
install.packages("testthat")
install.packages("vegan")

#load all necessary packages
library(Matrix)
library(testthat)
library(vegan)

#DEFINING THE FUNCTION

mask_value_only <- function(mat) {
  #ensure matrix is at least 2D
  if (is.vector(mat)) {
    mat <- matrix(mat, nrow = 1)
  }
  
  #ensure matrix is not more than 2D
  if (length(dim(mat)) > 2) {
    stop("Input matrix can only have two dimensions or less")
  }
  
  #generate logical mask: TRUE where values are missing
  mask <- !is.finite(mat)  
  
  #create masked matrix (using NA for invalid entries)
  masked_mat <- mat
  masked_mat[!is.finite(mat)] <- NA
  
  #return as a masked matrix
  return(structure(list(
    data = masked_mat,
    mask = mask
  ), class = "MaskedMatrix"))
}

#TESTING THE FUNCTION

test_matrix <- matrix(c(
  1, 2, NA,
  4, Inf, 6,
  7, 8, NaN
), nrow = 3, byrow = TRUE)

masked_result <- mask_value_only(test_matrix)

print(masked_result$data)  # Should replace non-finite values with NA
print(masked_result$mask)  # Should be TRUE where values are non-finite
