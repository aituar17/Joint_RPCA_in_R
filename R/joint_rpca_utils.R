#OPTSPACE FUNCTION AND ITS AUXILIARY FUNCTIONS

.optspace  <- function(x, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE){
    
    ## Preprocessing : x : partially revealed matrix
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    
    idxna <- is.na(x)
    
    if (!is.matrix(x)) {
      stop("* optspace : input 'x' should be a matrix")
    }
    if (!any(idxna)) {
      x # no NAs to be filled in and can be returned immediately
    }  
    if (any(is.infinite(x))) {
      stop("* optspace : infinite values are not allowed in 'x'")
    }
    
    m_e <- array(0, c(nrow(x), ncol(x)))
    m_e[!idxna] <- x[!idxna]
    
    ## Preprocessing : size information
    n <- nrow(x)
    m <- ncol(x)
    
    ## Preprocessing : other sparse-related concepts
    nnz_e <- sum(!idxna)
    E <- array(0, c(nrow(x), ncol(x)))
    E[!idxna] <- 1
    eps <- nnz_e / sqrt(m * n)
    
    ## Preprocessing : ropt  : implied rank
    if (ropt) {
      r <- round(ropt)
      if (!is.numeric(ropt) || (!is.numeric(r)) || (r < 1) || (r > m) || (r > n)) {
        stop("* optspace: value of argument 'ropt' should be integer
            in [1, min(nrow(x), ncol(x))]")
      }
    } else {
      r <- min(max(round(.guess_rank(m_e, nnz_e)), 2), m - 1)
      if (verbose) {
        message(paste0("* optspace: Guessing an implicit rank: Estimated rank 'ropt': ", r))
      }
    }
    
    ## Preprocessing : niter : maximum number of iterations
    if ((is.infinite(niter)) || (niter <= 1) || (!is.numeric(niter))) {
      stop("* optspace: invalid number provided for argument 'niter'")
    }
    niter <- round(niter)
    rho <-  eps * n
    
    ## Main Computation
    rescal_param <- sqrt(nnz_e * r / (norm(m_e, 'f')^2))
    m_e <- m_e * rescal_param
    
    # 1. SVD
    if (verbose) {
      message("* optspace: Step 2: SVD ...")
    }
    
    svdEt <- svd(m_e)
    X0 <- matrix(svdEt$u[, seq_len(r)], ncol=r)
    X0 <- matrix(X0[, rev(seq_len(ncol(X0)))], ncol=r)
    S0 <- diag(rev(svdEt$d[seq_len(r)]))
    Y0 <- matrix(svdEt$v[, seq_len(r)], ncol=r)
    Y0 <- matrix(Y0[, rev(seq_len(ncol(Y0)))], ncol=r)
    
    # 3. Initial Guess
    if (verbose) {
      message("* optspace: Step 3: Initial Guess ...")
    }
    X0 <- X0 * sqrt(n)
    Y0 <- Y0 * sqrt(m)
    S0 <- S0 / eps
    
    # 4. Gradient Descent
    if (verbose) {
      message("* optspace: Step 4: Gradient Descent ...")
    }
    X <- X0
    Y <- Y0
    S <- .aux_getoptS(X, Y, m_e, E)
    # initialize
    dist <- array(0, c(1, (niter + 1)))
    dist[1] <- norm((m_e - (X %*% S %*% t(Y))) * E, 'f') / sqrt(nnz_e)
    # Resolution/regularization for gradient calculation
    m0 <- 10000
    
    for (i in seq_len(niter)) {
      # compute the gradient
      tmpgrad <- .aux_gradF_t(X, Y, S, m_e, E, m0, rho)
      W <- tmpgrad$W
      Z <- tmpgrad$Z
      # line search for the optimum jump length
      t <- .aux_getoptT(X, W, Y, Z, S, m_e, E, m0, rho)
      X <- X + t * W
      Y <- Y + t * Z
      S <- .aux_getoptS(X, Y, m_e, E)
      # compute the distortion
      dist[i + 1] <- norm(((m_e - X %*% S %*% t(Y)) * E), 'f') / sqrt(nnz_e)    
      if (dist[i + 1] < tol) {
        dist <- dist[seq_len(i + 1)]
        break
      }
    }
    S <- S / rescal_param  
    # Return Results
    out <- list()
    
    # re-order Optspace may change order during iters
    index_order <- order(diag(S), decreasing = TRUE)
    X <- matrix(X[, index_order], ncol=length(index_order))
    Y <- matrix(Y[, index_order], ncol=length(index_order))
    S <- matrix(S[index_order, index_order], ncol=length(index_order))
    out$X <- X
    out$S <- S
    out$Y <- Y
    out$dist <- dist
    if (verbose) {
      message('* optspace: estimation finished.')
    }
    
    # -------------------------------------------
    
    # This part is not in the Python / Gemelli implementation
    # but has been added in R to provide more direct access
    # to the imputed matrix.
    
    # Reconstruct the matrix
    M <- X %*% S %*% t(Y)
    
    # Centering is common operation supporting output visualization
    # Center cols to 0
    M <- as.matrix(scale(M, center = TRUE, scale = FALSE))
    # Center rows to 0
    M <- as.matrix(t(scale(t(M), center = TRUE, scale = FALSE)))
    
    # Add imputed matrix to the output
    out$M <- M
    
    # -------------------------------------------
    
    out
  }


# Estimate Matrix Rank 
#
# x: numeric matrix for which the rank is to be estimated
#
# nnz: estimated number of non-zero entries to derive noise threshold
#
# keywords internal
.guess_rank <- function(x, nnz)
{
  maxiter <- 10000
  n <- nrow(x)
  m <- ncol(x)
  epsilon <- nnz / sqrt(m * n)
  svdX <- svd(x)
  S0 <- svdX$d
  
  nsval0 <- length(S0)
  S1 <- S0[seq_len(nsval0 - 1)] - S0[seq(2, nsval0)]  
  nsval1 <- length(S1)
  if (nsval1 > 10) {
    S1_ <- S1 / mean(S1[seq((nsval1 - 10), nsval1)])
  } else {
    S1_ <- S1 / mean(S1[seq_len(nsval1)])
  }
  r1 <- 0
  lam <- 0.05
  
  itcounter <- 0
  while (r1 <= 0) {
    itcounter <- itcounter + 1
    cost <- array(0, c(1, length(S1_)))
    for (idx in seq_len(length(S1_))) {
      cost[idx] <- lam * max(S1_[seq(idx, length(S1_))]) + idx
    }
    v2 <- min(cost)
    i2 <- which(cost == v2)
    if (length(i2) == 1) {
      r1 <- i2 - 1
    } else {
      r1 <- max(i2) - 1
    }
    lam <- lam + 0.05
    if (itcounter > maxiter) {
      break
    }
  }
  
  if (itcounter <= maxiter) {
    cost2 <- array(0, c(1, (length(S0) - 1)))
    for (idx in seq_len(length(S0) - 1)) {
      cost2[idx] <- (S0[idx + 1] + sqrt(idx * epsilon) * S0[1] / epsilon) / S0[idx]
    }
    v2 <- min(cost2)
    i2 <- which(cost2 == v2)
    if (length(i2) == 1) {
      r2 <- i2
    } else {
      r2 <- max(i2)
    }
    
    if (r1 > r2) {
      r <- r1
    } else {
      r <- r2
    }
  } else {
    r <- min(nrow(x), ncol(x))
  }
  r
}


# Auxiliary Gradient Contribution Function
#
# Compute distortion. Computes nonlinear transformation of the squared
# row norms of a matrix, with thresholding based on scaled values.
# Typically used as part of a gradient computation or optimization routine,
# where it selectively activates rows based on their scaled magnitude.
#
# x: numeric matrix; the function computes the squared norm of each row.
#
# m0 positive scalar that acts as a scaling parameter in the nonlinear regularization term
#
# r: numeric scalar; another scaling factor applied alongside m0 in the normalization of row norms.
#' Its effect is similar to m0, the two are often coupled in optimization problems.
#
# keywords internal
.aux_G <- function(x, m0, r)
{
  z <- rowSums(x^2) / (2 * m0 * r)
  y <- exp((z - 1)^2) - 1
  idxfind <- (z < 1)
  y[idxfind] <- 0
  out <- sum(y)
  out
}

# Total Loss Function with Regularization
#
# Composite loss function combining a weighted squared error term and 
# nonlinear regularization penalties based on the row norms of the input matrices. 
# This function is commonly used in matrix factorization or low-rank modeling 
# frameworks where latent matrices are learned under structural constraints.
#
# x numeric matrix, typically representing latent factors or components.
#
# y numeric matrix, typically representing latent factors or components.
#
# s numeric matrix used for linear transformation between \code{x} and \code{y}.
#
# m_e numeric matrix of size representing observed or target values to approximate.
#
# e numeric matrix of the same dimension as \code{m_e}, used as a weight or mask matrix. Values typically range from 0 to 1.
#
# m0 positive scalar that acts as a scaling parameter in the nonlinear regularization term
#
# rho positive scalar controlling the strength of the regularization terms
#
# keywords internal
.aux_F_t <- function(x, y, s, m_e, e, m0, rho)
{
  n <- nrow(x)
  r <- ncol(x)
  out1 <- sum((((x %*% s %*% t(y)) - m_e) * e)^2) / 2
  out2 <- rho * .aux_G(y, m0, r)
  out3 <- rho * .aux_G(x, m0, r)
  out  <- out1 + out2 + out3
  out
}



# Gradient of Auxiliary Regularization Term
#
# Computes the gradient of the nonlinear regularization function \code{.aux_G} 
# with respect to its matrix input \code{x}. 
#
# x numeric matrix of size \code{n x r}. Each row is treated as a vector
# whose squared norm determines the activation of the gradient.
#
# m0 positive scalar that acts as a scaling parameter in the nonlinear regularization term
#
# r: numeric scalar; another scaling factor applied alongside m0 in the normalization of row norms.
# Its effect is similar to m0, the two are often coupled in optimization problems.
#
# keywords internal
.aux_Gp <- function(x, m0, r)
{
  z <- rowSums(x^2) / (2 * m0 * r)
  z <- 2 * exp((z - 1)^2) / (z - 1)
  idxfind <- (z < 0)
  z[idxfind] <- 0  
  out <- x * matrix(z, nrow = nrow(x), ncol = ncol(x), byrow = FALSE) / (m0 * r)
}


# Gradient of Composite Loss Function
#
# Computes the gradient of the total loss. The result includes both the data 
# reconstruction gradient and the nonlinear regularization gradient.
#
# x numeric matrix of size \code{n x r}, representing one set of latent 
# variables or factor matrix.
#
# y numeric matrix of size \code{m x r}, representing the counterpart latent 
# variable matrix.
#
# s numeric matrix used in the bilinear transformation between \code{x} and \code{y}.
#
# m_e A numeric matrix representing the observed data or target matrix 
#
# e numeric matrix used as a weighting or masking matrix.
#
# m0 positive scalar for gradient regularization 
#
# rho weight applied to the regularization gradient terms.
#
# @keywords internal
.aux_gradF_t <- function(x, y, s, m_e, e, m0, rho)
{
  n <- nrow(x)
  r <- ncol(x)
  m <- nrow(y)
  if (ncol(y) != r) {
    stop("dimension error from the internal function .aux_gradF_t")
  }
  
  XS  <- x %*% s
  YS  <- y %*% t(s)
  XSY <- XS %*% t(y)
  
  Qx <- t(x) %*% ((m_e - XSY) * e) %*% YS / n
  Qy <- t(y) %*% t((m_e - XSY) * e) %*% XS / m
  
  W <- ((XSY - m_e) * e) %*% YS  + (x %*% Qx) + rho * .aux_Gp(x, m0, r)
  Z <- t((XSY - m_e) * e) %*% XS + (y %*% Qy) + rho * .aux_Gp(y, m0, r)
  
  resgrad <- list()
  resgrad$W <- W
  resgrad$Z <- Z
  resgrad
  
}


# Solve for Optimal Transformation Matrix S
#
# Computes the optimal transformation matrix that minimizes the squared 
# error, optionally weighted by importance matrix \code{e}.
# This function forms and solves a linear least-squares problem.
#
# x A numeric matrix representing one side of the bilinear transformation.
#
# y A numeric matrix representing the other side of the transformation.
#
# m_e A numeric matrix representing the observed data or target matrix to approximate.
#
# e numeric matrix used as a weighting or masking matrix.
#
# # @keywords internal
.aux_getoptS <- function(x, y, m_e, e)
{
  n <- nrow(x)
  r <- ncol(x)  
  C <- t(x) %*% (m_e) %*% y
  C <- matrix(as.vector(C))  
  nnrow <- ncol(x) * ncol(y)
  A <- matrix(NA, nrow = nnrow, ncol = (r^2))
  
  for (i in seq_len(r)) {
    for (j in seq_len(r)) {
      ind <- (j - 1) * r + i
      tmp <- t(x) %*% (outer(x[, i], y[, j]) * e) %*% y      
      A[, ind] <- as.vector(tmp)
    }
  }
  
  S <- solve(A, C)
  out <- matrix(S, nrow = r)
  out
}



# Backtracking Line Search for Step Size Optimization
#
# Optimal line search.
# Determines an optimal step size \code{t} for descending along a direction 
# defined by perturbations \code{w} and \code{z} for the matrices \code{x} and \code{y}. 
#
# x numeric matrix representing current values of one latent factor.
#
# w numeric matrix representing the search direction (gradient or descent vector) for \code{x}.
#
# y numeric matrix representing current values of the second latent factor.
#
# z numeric matrix representing the search direction for \code{y}.
#
# s numeric matrix used in the bilinear product in \code{.aux_F_t}.
#
# m_e numeric the observed or target matrix.
#
# e numeric mask or weighting matrix.
#
# m0 positive scalar regularization parameter.
#
# rho positive scalar weighting the regularization term in the total loss function.
#
# keywords internal
.aux_getoptT <- function(x, w, y, z, s, m_e, e, m0, rho)
{
  norm2WZ <- norm(w, 'f')^2 + norm(z, 'f')^2
  f <- array(0, c(1, 21))
  f[1] <- .aux_F_t(x, y, s, m_e, e, m0, rho)
  t <- -1e-1
  for (i in seq_len(21)) {
    f[i + 1] <- .aux_F_t(x + t * w, y + t * z, s, m_e, e, m0, rho)
    if ((f[i + 1] - f[1]) <= 0.5 * t * norm2WZ) {
      out <- t
      break
    }
    t <- t / 2
  }
  out <- t
  t
}

#FUNCTION FOR RPCA TABLE PROCESSING

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

#FUNCTION FOR MASK VALUE ONLY

.mask_value_only <- function(mat) {
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
  
  #create masked matrix
  masked_mat <- mat
  masked_mat[!is.finite(mat)] <- NA
  
  #return as a masked matrix
  return(structure(list(
    data = masked_mat,
    mask = mask
  ), class = "MaskedMatrix"))
}

#ORDINATION_RESULTS FUNCTION

#constructor function
.OrdinationResults <- function(method, eigvals, samples, features,
                              proportion_explained, dist = NULL, metadata = list()) {
  structure(list(
    method = method,
    eigvals = eigvals,
    samples = samples,
    features = features,
    proportion_explained = proportion_explained,
    dist = dist,
    metadata = metadata
  ), class = "OrdinationResults")
}

#print method
.print.OrdinationResults <- function(x, ...) {
  cat("OrdinationResults (method:", x$method, ")\n")
  cat("Number of components:", length(x$eigvals), "\n")
  cat("Variance explained:\n")
  print(round(x$proportion_explained, 3))
  invisible(x)
}

#summary method
.summary.OrdinationResults <- function(object, ...) {
  print(object)
  cat("\nSample scores (first few rows):\n")
  print(head(object$samples))
  cat("\nFeature loadings (first few rows):\n")
  print(head(object$features))
  invisible(object)
}

#plot method
.plot.OrdinationResults <- function(x, comps = c(1, 2), ...) {
  if (length(comps) != 2) stop("Please select two components to plot.")
  plot(x$samples[, comps], col = "blue", pch = 19,
       xlab = paste0("PC", comps[1]),
       ylab = paste0("PC", comps[2]),
       main = paste("Ordination (", x$method, ")", sep = ""))
  points(x$features[, comps], col = "red", pch = 4)
  legend("topright", legend = c("Samples", "Features"),
         col = c("blue", "red"), pch = c(19, 4))
}

#TRANSFORM FUNCTION

#rclr_transform function
.rclr_transform <- function(df) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for rclr_transform.")
  }
  vegan::decostand(df, method = "rclr")
}

#helper to normalize and project new data into ordination space
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

#main transform function
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

#DISTANCE MATRIX FUNCTION

.DistanceMatrix <- function(matrix, ids = NULL, method = "euclidean") {
  if (!is.matrix(matrix)) stop("Input must be a matrix.")
  if (!isSymmetric(matrix)) stop("Distance matrix must be symmetric.")
  if (!is.null(ids)) {
    if (length(ids) != nrow(matrix)) stop("Length of 'ids' must match matrix dimensions.")
    rownames(matrix) <- ids
    colnames(matrix) <- ids
  }
  structure(list(
    data = matrix,
    ids = rownames(matrix),
    method = method
  ), class = "DistanceMatrix")
}

.print.DistanceMatrix <- function(x, ...) {
  cat("DistanceMatrix (", x$method, ")\n", sep = "")
  cat("Number of objects:", length(x$ids), "\n")
  print(head(x$data, 6))  # Show only top part
  invisible(x)
}

.summary.DistanceMatrix <- function(object, ...) {
  cat("Summary of DistanceMatrix\n")
  cat("Method:", object$method, "\n")
  cat("Size:", nrow(object$data), "x", ncol(object$data), "\n")
  cat("IDs:\n")
  print(head(object$ids, 6))
  cat("\nDistance Summary Stats:\n")
  print(summary(as.vector(object$data[upper.tri(object$data)])))
  invisible(object)
}

#FUNCTION FOR OPTSPACE HELPER

.optspace_helper <- function(rclr_table,
                            feature_ids,
                            subject_ids,
                            n_components = 3,
                            max_iterations = 5) {
  
  #fit OptSpace
  opt_result <- .optspace(rclr_table, ropt = n_components, niter = max_iterations, tol = 1e-5, verbose = FALSE)
  
  #update n_components
  n_components <- ncol(opt_result$S)
  
  #reconstruct and re-center matrix
  X_hat <- opt_result$X %*% opt_result$S %*% t(opt_result$Y)
  X_hat <- scale(X_hat, center = TRUE, scale = FALSE)
  X_hat <- t(scale(t(X_hat), center = TRUE, scale = FALSE))
  
  #PCA
  svd_out <- svd(X_hat)
  u <- svd_out$u[, 1:n_components, drop = FALSE]
  s <- svd_out$d[1:n_components]
  v <- svd_out$v[, 1:n_components, drop = FALSE]
  
  #label loadings
  rename_cols <- paste0("PC", seq_len(n_components))
  sample_scores <- data.frame(u, row.names = subject_ids)
  feature_scores <- data.frame(v, row.names = feature_ids)
  colnames(sample_scores) <- rename_cols
  colnames(feature_scores) <- rename_cols
  
  #proportion explained
  prop_var <- s^2 / sum(svd_out$d^2)
  names(prop_var) <- rename_cols
  names(s) <- rename_cols
  
  #add PC3 for 2D case
  if (n_components == 2) {
    sample_scores$PC3 <- 0
    feature_scores$PC3 <- 0
    s <- c(s, PC3 = 0)
    prop_var <- c(prop_var, PC3 = 0)
    rename_cols <- c(rename_cols, "PC3")
  }
  
  #compute distance
  dist_matrix_raw <- as.matrix(dist(u))
  rownames(dist_matrix_raw) <- subject_ids
  colnames(dist_matrix_raw) <- subject_ids
  
  #wrap with DistanceMatrix
  dist_res <- .DistanceMatrix(dist_matrix_raw, ids = subject_ids, method = "aitchison")
  
  #build OrdinationResults object
  ord_res <- .OrdinationResults(
    method = "rpca_biplot",
    eigvals = s,
    samples = sample_scores,
    features = feature_scores,
    proportion_explained = prop_var,
    dist = dist_matrix_raw,  
    metadata = list(
      long_method_name = "(Robust Aitchison) RPCA Biplot",
      run_id = sprintf("optspace_helper_n_components_%d.max_iterations_%d", 
                       n_components, max_iterations)
    )
  )
  
  return(list(
    ord_res = ord_res,
    dist = dist_res,  
    opt_fit = opt_result
  ))
}

#JOINT SOLVE FUNCTION

.joint_optspace_solve <- function(train_test_pairs, n_components,
                                 max_iter = 50, verbose = TRUE) {
  #prepare lists to hold training matrices and dimensions
  train_matrices <- list()
  test_matrices <- list()
  dims <- list()
  
  for (pair in train_test_pairs) {
    test_mat <- pair[[1]]
    train_mat <- pair[[2]]
    train_matrices <- append(train_matrices, list(train_mat))
    test_matrices <- append(test_matrices, list(test_mat))
    dims <- append(dims, list(dim(train_mat)))
  }
  
  #stack training matrices horizontally 
  train_stacked <- do.call(cbind, train_matrices)
  
  #apply OptSpace to stacked matrix
  if (verbose) message("Running optspace() on stacked training data...")
  fit <- .optspace(train_stacked, ropt = n_components, niter = max_iter, tol = 1e-5, verbose = FALSE)
  
  #extract sample loadings
  U_shared <- fit$X
  S_shared <- fit$S
  
  #split V back into per-table pieces
  feat_indices <- cumsum(sapply(dims, function(d) d[2]))
  feat_starts <- c(1, head(feat_indices, -1) + 1)
  V_list <- Map(function(start, end) fit$Y[start:end, , drop = FALSE],
                feat_starts, feat_indices)
  
  dists <- matrix(0, nrow = 2, ncol = max_iter)
  for (i in seq_along(test_matrices)) {
    V_k <- V_list[[i]]
    test_mat <- test_matrices[[i]]
    
    #project test samples: U_test = test × V
    U_test <- as.matrix(test_mat) %*% V_k
    U_test <- sweep(U_test, 2, diag(S_shared), "/")
    recon_test <- U_test %*% S_shared %*% t(V_k)
    
    #center for consistency
    recon_test <- scale(recon_test, center = TRUE, scale = FALSE)
    recon_test <- t(scale(t(recon_test), center = TRUE, scale = FALSE))
    
    error <- test_mat - recon_test
    error[is.na(error)] <- 0
    error_val <- norm(error, "F") / sqrt(sum(!is.na(test_mat)))
    
    dists[1, i] <- error_val
    dists[2, i] <- 0  
  }
  
  list(U = U_shared, S = S_shared, V_list = V_list, dists = dists)
}

#FUNCTION FOR JOINT OPTSPACE HELPER

.joint_optspace_helper <- function(tables,
                                  n_components,
                                  max_iterations,
                                  test_samples,
                                  train_samples) {
  #split and transpose training/test data per table
  tables_split <- lapply(tables, function(tbl) {
    list(t(tbl[, test_samples, drop = FALSE]),
         t(tbl[, train_samples, drop = FALSE]))
  })
  
  #format input for solver
  tables_for_solver <- lapply(tables_split, function(pair) {
    lapply(pair, as.matrix)
  })
  
  #run joint OptSpace solver
  opt_result <- .joint_optspace_solve(tables_for_solver,
                                     n_components = n_components,
                                     max_iter = max_iterations)
  
  U <- opt_result$U
  S <- opt_result$S
  V_list <- opt_result$V_list
  dists <- opt_result$dists
  
  #assign row/column names to loadings
  pc_names <- paste0("PC", seq_len(n_components))
  
  #combine feature loadings with table-derived row names
  vjoint <- do.call(rbind, Map(function(tbl, V) {
    rownames(V) <- rownames(tbl)
    colnames(V) <- pc_names
    V
  }, tables, V_list))
  
  U <- U[seq_along(train_samples), , drop = FALSE]
  rownames(U) <- train_samples
  colnames(U) <- pc_names
  
  #recenter & re-factor via SVD
  X <- U %*% S %*% t(vjoint)
  X <- sweep(X, 2, colMeans(X))
  X <- sweep(X, 1, rowMeans(X))
  svd_res <- svd(X)
  u <- svd_res$u[, seq_len(n_components), drop = FALSE]
  v <- svd_res$v[, seq_len(n_components), drop = FALSE]
  s_eig <- svd_res$d[seq_len(n_components)]
  
  rownames(u) <- train_samples
  rownames(v) <- rownames(vjoint)
  colnames(u) <- colnames(v) <- pc_names
  
  #create ordination object
  prop_exp <- s_eig^2 / sum(s_eig^2)
  ord_res <- .OrdinationResults(
    method = "rpca",
    eigvals = setNames(s_eig, pc_names),
    samples = u,
    features = v,
    proportion_explained = setNames(prop_exp, pc_names)
  )
  
  #project test samples
  if (length(test_samples) > 0) {
    test_matrices <- lapply(tables, function(tbl) {
      tbl[, test_samples, drop = FALSE]
    })
    ord_res <- .transform(ord_res, test_matrices, apply_rclr = FALSE)
  }
  
  #compute distance matrix and CV error summary
  dist_mat <- as.matrix(dist(ord_res$samples))
  dist_res <- .DistanceMatrix(dist_mat, ids = rownames(ord_res$samples))
  
  cv_dist <- data.frame(t(dists))
  colnames(cv_dist) <- c("mean_CV", "std_CV")
  cv_dist$run <- sprintf("tables_%d.n_components_%d.max_iterations_%d.n_test_%d",
                         length(tables), n_components, max_iterations, length(test_samples))
  cv_dist$iteration <- seq_len(nrow(cv_dist))
  rownames(cv_dist) <- seq_len(nrow(cv_dist))
  
  list(ord_res = ord_res, dist = dist_res, cv_stats = cv_dist)
}
