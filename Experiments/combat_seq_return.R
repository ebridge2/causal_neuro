ComBat_seq_return_model <- function (counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE, 
                                     shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) 
{
  batch <- as.factor(batch)
  counts <- t(counts)
  if (any(table(batch) <= 1)) {
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  keep_lst <- lapply(levels(batch), function(b) {
    which(apply(counts[, batch == b], 1, function(x) {
      !all(x == 0)
    }))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  dge_obj <- edgeR:::DGEList(counts = counts)
  n_batch <- nlevels(batch)
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  n_batches <- sapply(batches_ind, length)
  n_sample <- sum(n_batches)
  cat("Found", n_batch, "batches\n")
  batchmod <- model.matrix(~-1 + batch)
  group <- as.factor(group)
  if (full_mod & nlevels(group) > 1) {
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }
  else {
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data = as.data.frame(t(counts)))
  }
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), 
                                         function(i) {
                                           model.matrix(~covar_mod[, i])
                                         }))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x) {
      all(x == 1)
    })]
  }
  mod <- cbind(mod, covar_mod)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")
    }
    if (ncol(design) > (n_batch + 1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, 
                                                          -c(1:n_batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")
      }
    }
  }
  NAs = any(is.na(counts))
  if (NAs) {
    cat(c("Found", sum(is.na(counts)), "Missing Data Values\n"), 
        sep = " ")
  }
  cat("Estimating dispersions\n")
  disp_common <- sapply(1:n_batch, function(i) {
    if ((n_batches[i] <= ncol(design) - ncol(batchmod) + 
         1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)) {
      return(edgeR:::estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                                   design = NULL, subset = nrow(counts)))
    }
    else {
      return(edgeR:::estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                                   design = mod[batches_ind[[i]], ], subset = nrow(counts)))
    }
  })
  genewise_disp_lst <- lapply(1:n_batch, function(j) {
    if ((n_batches[j] <= ncol(design) - ncol(batchmod) + 
         1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)) {
      return(rep(disp_common[j], nrow(counts)))
    }
    else {
      return(edgeR:::estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], 
                                    design = mod[batches_ind[[j]], ], dispersion = disp_common[j], 
                                    prior.df = 0))
    }
  })
  names(genewise_disp_lst) <- paste0("batch", levels(batch))
  phi_matrix <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (k in 1:n_batch) {
    phi_matrix[, batches_ind[[k]]] <- sva:::vec2mat(genewise_disp_lst[[k]], 
                                              n_batches[k])
  }
  cat("Fitting the GLM model\n")
  glm_f <- edgeR:::glmFit(dge_obj, design = design, dispersion = phi_matrix, 
                  prior.count = 1e-04)
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample)
  new_offset <- t(sva:::vec2mat(edgeR:::getOffset(dge_obj), nrow(counts))) + 
    sva:::vec2mat(alpha_g, ncol(counts))
  glm_f2 <- edgeR:::glmFit.default(dge_obj$counts, design = design, 
                           dispersion = phi_matrix, offset = new_offset, prior.count = 1e-04)
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  if (shrink) {
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- sva:::monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(dat = counts[, batches_ind[[ii]]], 
                           mu = mu_hat[, batches_ind[[ii]]], gamma = gamma_hat[, 
                                                                               ii], phi = phi_hat[, ii], gene.subset.n = gene.subset.n)
      }
      else {
        invisible(capture.output(mcres <- mcint_fun(dat = counts[, 
                                                                 batches_ind[[ii]]], mu = mu_hat[, batches_ind[[ii]]], 
                                                    gamma = gamma_hat[, ii], phi = phi_hat[, ii], 
                                                    gene.subset.n = gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0("batch", levels(batch))
    gamma_star_mat <- lapply(monte_carlo_res, function(res) {
      res$gamma_star
    })
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res) {
      res$phi_star
    })
    phi_star_mat <- do.call(cbind, phi_star_mat)
    if (!shrink.disp) {
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }
  else {
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  mu_star <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (jj in 1:n_batch) {
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]]) - 
                                          sva:::vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  
  return(list(glm=glm_f2, phi_hat=phi_hat, mu_star=mu_star, phi_star=phi_star, keep=keep, rm=rm))
}

ComBat_seq_apply_model <- function(counts, batch, model) {
  counts <- t(counts)
  countsOri <- counts
  batch <- as.factor(batch)
  n_batch <- nlevels(batch)
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  dge_obj <- edgeR:::DGEList(counts = counts)
  
  
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (kk in 1:n_batch) {
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- model$mu_hat[, batches_ind[[kk]]]
    old_phi <- model$phi_hat[, kk]
    new_mu <- model$mu_star[, batches_ind[[kk]]]
    new_phi <- model$phi_star
    adjust_counts[, batches_ind[[kk]]] <- sva:::match_quantiles(counts_sub = counts_sub, 
                                                          old_mu = old_mu, old_phi = old_phi, new_mu = new_mu, 
                                                          new_phi = new_phi)
  }
  adjust_counts_whole <- matrix(NA, nrow = nrow(countsOri), 
                                ncol = ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[model$keep, ] <- adjust_counts
  adjust_counts_whole[model$rm, ] <- countsOri[model$rm, ]
  return(t(adjust_counts_whole))
}

combat_seq_mod_fit <- function(Ys, Ts, Xs) {
  
  Ns.samp <- apply(Ys, 1, sum)
  mods <- list()
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys[,k], Batch=Ts, Covar=Xs, N=log(Ns.samp))
    mods[[k]] <- glm.nb(Y ~ factor(Batch) + Covar + offset(N), data=dat)
  }
  
  return(mod)
}

combat_seq_mod_apply <- function(Ys, Ts, Xs, mods) {
  Ns.samp <- apply(Ys, 1, sum)
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys[,k], Batch=Ts, Covar=Xs, N=log(Ns.samp))
    muhats <- predict(mods[[k]], dat, type="response")
    mustar <- muhats/exp(mods[[k]]$coefficients["factor(Batch)1"])
  }
  
}