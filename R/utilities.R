require(copula)

make_closest_pos_semidef_corr_matrix <- function(M) {
  # Compute eigenvalue decomposition
  eigen_decomp <- eigen(M)
  
  # Get the eigenvalues and eigenvectors
  eig_values <- eigen_decomp$values
  eig_vectors <- eigen_decomp$vectors
  
  # Set negative eigenvalues to zero
  eig_values[eig_values < 0] <- 0
  
  # Reconstruct the positive semidefinite matrix
  M_pos_semidef <- eig_vectors %*% diag(sqrt(eig_values)) %*% t(eig_vectors)
  
  # Normalize to ensure it represents a correlation matrix
  M_corr <- M_pos_semidef / sqrt(diag(M_pos_semidef) %*% t(diag(M_pos_semidef)))
  
  return(M_corr)
}

make_corr_mtx <- function(ncovar, unobs.cor) {
  corr.mtx <- diag(ncovar + 1)
  corr.mtx[1,ncovar + 1] = corr.mtx[ncovar + 1,1] = unobs.cor
  # multiple covariates are correlated with the unobserved
  if (ncovar > 1) {
    for (i in 2:ncovar) {
      corr.mtx[i, ncovar + 1] = corr.mtx[ncovar + 1, i] = 0.8
      corr.mtx[1, i] = corr.mtx[i, 1] = unobs.cor
    }
    # force it to be a positive semidefinite correlation matrix
    # by ignoring evecs w negative evals and normalizing
    corr.mtx <- make_closest_pos_semidef_corr_matrix(corr.mtx)
  }
  return(corr.mtx)
}

generate_correlated_betas <- function(n, alphas, betas, corr.mtx=NULL) {
  ndim <- ifelse(is.null(dim(corr.mtx)), 1, nrow(corr.mtx))
  if (ndim > 1) {
    # Create a Gaussian copula with the specified correlation
    copula_obj <- normalCopula(param=P2p(corr.mtx), dim = ndim, dispstr="un")
    
    # Generate uniform samples from the copula
    u <- rCopula(n, copula_obj) 
  } else {
      u <- matrix(runif(n), ncol=1)
  }
  
  # Transform uniform samples to (2*beta-1)-distributed variables
  Xs <- do.call(cbind, lapply(1:ndim, function(k) {
    2*sapply(1:n, function(i) {
      qbeta(u[i,k], alphas[k], betas[k])
    }) - 1
  }))
  return(matrix(Xs, ncol=ndim))
}

generate_correlated_betas_pergroup <- function(batches, alphas1, betas1, alphas2, betas2, corr.mtx=NULL) {
  n = length(batches)
  ndim <- ifelse(is.null(dim(corr.mtx)), 1, nrow(corr.mtx))
  if (ndim > 1) {
    # Create a Gaussian copula with the specified correlation
    copula_obj <- normalCopula(param=P2p(corr.mtx), dim = ndim, dispstr="un")
    
    # Generate uniform samples from the copula
    u <- rCopula(n, copula_obj)  } else {
      u <- matrix(runif(n), ncol=1)
    }
  
  # Transform uniform samples to (2*beta-1)-distributed variables
  Xs <- do.call(cbind, lapply(1:ndim, function(k) {
    2*sapply(1:length(batches), function(i) {
      ifelse(batches[i] == 0, qbeta(u[i,k], alphas1[k], betas1[k]), 
             qbeta(u[i,k], alphas2[k], betas2[k]))
    }) - 1
  }))
  return(matrix(Xs, ncol=ndim))
}