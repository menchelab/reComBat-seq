# Based on ComBat_seq from the sva package (Zhang et al. 2020).
# Extended by Zhasmina Stoyanova, 2026
# Changes: replaced standard NB GLM with elastic net regularized NB GLM;
#          added lambda.reg, alpha.reg, num.threads parameters throughout.
#' Adjust for batch effects using an empirical Bayes framework in RNA-seq raw counts
#'
#' reComBat_seq is an extension to the ComBat_seq method using regularized Negative Binomial model.
#'
#' @param counts Raw count matrix from genomic studies (dimensions gene x sample)
#' @param batch vector containing batch assignment of samples
#' @param wanted.variation a data.frame containing the covariates you want to preserve in the data
#' @param shrink Boolean, whether to apply empirical Bayes estimation on parameters
#' @param shrink.disp Boolean, whether to apply empirical Bayes estimation on dispersion
#' @param gene.subset.n Number of genes to use in empirical Bayes estimation, only useful when shrink = TRUE
#' @param lambda.reg Regularization strength
#' @param alpha.reg Elastic Net Tuner, 1 for pure LASSO, 0 for pure ridge
#' @param num.threads Threads for Paralellisation
#' @return data A probe x sample count matrix, adjusted for batch effects.
#'
#' @examples
#'
#'
#'
#' @export
#'

reComBat.seq <- function(
  counts, 
  batch, 
  wanted.variation=NULL,
  shrink=FALSE, 
  shrink.disp=FALSE, 
  gene.subset.n=NULL,
  lambda.reg=0.8, 
  alpha.reg=0.3, 
  num.threads=1
){
  ########  Preparation  ########
  counts <- as.matrix(counts)

  ## Does not support 1 sample per batch yet
  batch <- as.factor(batch)
  if(any(table(batch)<=1)){
    stop("reComBat-seq doesn't support 1 sample per batch yet")
  }

  message("Using THREADS: " , num.threads)

  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=counts)

  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch
  n_batches <- sapply(batches_ind, length)
  n_sample <- sum(n_batches)
  message("Found ",n_batch,' batches')

  ## Make design matrix
  design <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariates to preserve
  if(!is.null(wanted.variation)) {
    n_covariates <- ncol(wanted.variation)
    wanted.variation[] <- lapply(wanted.variation, as.factor)
    wanted.variation <- do.call(
      cbind, 
      lapply(
        1:ncol(wanted.variation), 
        function(i){model.matrix(~wanted.variation[,i])}
      )
    )
    # this is to emulate the mod variable which is used further down
    intercept <- data.frame(intercept = rep(1, dim(wanted.variation)[1]))
    colnames(intercept) <- c('(Intercept)')
    mod <- cbind(intercept, wanted.variation[, !apply(wanted.variation, 2, function(x){all(x==1)})])
    message("Adjusting for ", n_covariates, ' covariates with a total of ', ncol(wanted.variation), ' covariate level(s)')
    design <- cbind(design, wanted.variation)
  } else {
    mod <- model.matrix(~1, data = as.data.frame(t(counts)))
    message('No covariates to adjust for')
  }

  ## Check for intercept in covariates, and drop if present
  is.intercept <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!is.intercept])

  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)) {
    if(ncol(design)==(n_batch+1)) {stop("The covariate is confounded with batch!")}
    if(ncol(design)>(n_batch+1)) {
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){
        message('The covariates are confounded!\n')
      } else {
        message("At least one covariate is confounded with batch!")
      }
    }
  }

  ## Check for missing values in count matrix
  NAs = any(is.na(counts))
  if(NAs){message(c('Found',sum(is.na(counts)),'Missing Data Values'),sep=' ')}


  ########  Estimate gene-wise dispersions within each batch  ########
  message("Estimating common dispersions per batch")
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- simplify2array(
    mclapply(
      1:n_batch, 
      function(i) {
        if((n_batches[i] <= ncol(design)-n_batch+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){
          # not enough residual degree of freedom
          return(
            estimateGLMCommonDisp(
              counts[, batches_ind[[i]]], 
              design=NULL, 
              subset=nrow(counts)
            )
          )
        } else {
          return(
            estimateGLMCommonDisp(
              counts[, batches_ind[[i]]], 
              design=mod[batches_ind[[i]], ], 
              subset=nrow(counts), 
              lambda_reg=lambda.reg, 
              alpha_reg=alpha.reg
            )
          )
        }
      },
      mc.cores = num.threads
    )
  )

  ## Estimate gene-wise dispersion within each batch
  message("Estimating gene-wise dispersions per batch")
  genewise_disp_lst <- mclapply(
    1:n_batch, 
    function(j) {
      if((n_batches[j] <= ncol(design)-n_batch+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
        # not enough residual degrees of freedom - use the common dispersion
        return(rep(disp_common[j], nrow(counts)))
      } else {
        return(
          estimateGLMTagwiseDisp(
            counts[, batches_ind[[j]]], 
            design=mod[batches_ind[[j]], ],
            dispersion=disp_common[j], 
            prior.df=0, 
            lambda_reg=lambda.reg, 
            alpha_reg=alpha.reg
          )
        )
      }
    },
    mc.cores = num.threads
  )
  names(genewise_disp_lst) <- paste0('batch', levels(batch))

  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch) {
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k])
  }

  ########  Estimate parameters from NB GLM  ########
  message("Fitting the GLM model")
  # conform with multithreading logic
  if(num.threads == 1) {
    use_threads <- 0
  } else {
    use_threads <- num.threads
  }
  # no intercept - nonEstimable; compute offset (library sizes) within function
  glm_f <- glmFit(
    dge_obj, 
    design=design, 
    dispersion=phi_matrix, 
    prior.count=1e-4, 
    lambda_reg=lambda.reg, 
    alpha_reg=alpha.reg, 
    num_threads=use_threads
  )
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- (
    t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size
  )
  glm_f2 <- glmFit.default(
    dge_obj$counts, 
    design=design, 
    dispersion=phi_matrix, 
    offset=new_offset, 
    prior.count=1e-4,
    maxit=51, 
    lambda_reg=lambda.reg, 
    alpha_reg=alpha.reg,  
    num_threads=use_threads
  )
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)

  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########
  if(shrink){
    message("Apply shrinkage - computing posterior estimates for parameters")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- mclapply(
      1:n_batch, 
      function(ii) {
        if(ii==1) {
          mcres <- mcint_fun(
            dat=counts[, batches_ind[[ii]]], 
            mu=mu_hat[, batches_ind[[ii]]],
            gamma=gamma_hat[, ii], 
            phi=phi_hat[, ii], 
            gene.subset.n=gene.subset.n)
        } else {
          invisible(
            capture.output(
              mcres <- mcint_fun(
                dat=counts[, batches_ind[[ii]]], 
                mu=mu_hat[, batches_ind[[ii]]],
                gamma=gamma_hat[, ii], 
                phi=phi_hat[, ii], 
                gene.subset.n=gene.subset.n
              )
            )
          )
        }
        return(mcres)
      },
      mc.cores = num.threads
    )
    names(monte_carlo_res) <- paste0('batch', levels(batch))

    gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
    phi_star_mat <- do.call(cbind, phi_star_mat)

    if(!shrink.disp) {
      message("Apply shrinkage to mean only")
      phi_star_mat <- phi_hat
    }
  } else {
    message("Shrinkage off - using GLM estimates for parameters")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }

  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(jj in 1:n_batch) {
    mu_star[, batches_ind[[jj]]] <- exp(
      log(mu_hat[, batches_ind[[jj]]]) - 
      vec2mat(gamma_star_mat[, jj], n_batches[jj])
    )
  }
  phi_star <- rowMeans(phi_star_mat)


  ########  Adjust the data  ########
  message("Adjusting the data")
  adjusted_counts <- do.call(
    cbind,
    mclapply(
      1:n_batch,
      function(kk) {
        return(
          match_quantiles(
            counts_sub=counts[, batches_ind[[kk]]],
            old_mu=mu_hat[, batches_ind[[kk]]], 
            old_phi=phi_hat[, kk],
            new_mu=mu_star[, batches_ind[[kk]]], 
            new_phi=phi_star
          ) 
        )
      },
      mc.cores = num.threads
    )
  )

  dimnames(adjusted_counts) <- dimnames(counts[,do.call(c, batches_ind)])
  return(adjusted_counts[,colnames(counts)])
}

