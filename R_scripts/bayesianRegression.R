# TODO: Add comment
# 
# Author: Christoph
# Modified to predict degradation rates by Kostya
###############################################################################

BBSR <- function(X, Y, clr.mat, nS, no.pr.val, weights.mat, cores, X.dr, ts.conds, use.deg.rates, dynamics.prior.mat, dynamics.prior.weights, scale.desres) {
  
  # Scale and permute design and response matrix
  # t(scale(t(X))) makes every row of X a unit vector (in Euclidean norm) with mean 0.
  if(scale.desres) {
    X <- t(scale(t(X))) 
    Y <- t(scale(t(Y)))
  }
    
  #Added by Kostya: First we scale X.dr, then we zero the steady-state conditions:
  if(use.deg.rates) {
    if(scale.desres) {
      X.dr <- t(scale(t(X.dr)))
    }
    X.dr[,setdiff(colnames(X.dr), ts.conds)]=0
  }
  
  G <- nrow(Y)  # number of genes
  K <- nrow(X)  # max number of possible predictors (number of TFs)
  
  pp <- matrix(FALSE, G, K, dimnames=list(rownames(Y), rownames(X)))  # predictors that will be used in the regression
  
  # keep all predictors that we have priors for
  pp[(weights.mat != no.pr.val) & !is.na(clr.mat)] <- TRUE
  
  # for each gene, add the top nS predictors of the list to possible predictors
  clr.mat[clr.mat == 0] <- NA
  for (ind in 1:G) {
    clr.order <- order(clr.mat[ind, ], decreasing=TRUE, na.last=NA)
    pp[ind, clr.order[1:min(K, nS, length(clr.order))]] <- TRUE
  }
  diag(pp) <- FALSE
  
  # NEW: remove predictors that are constant or essentially duplicates of others
  const <- which(apply(is.na(X), 1, sum) > 0)
  cat('The following predictors are constant and will not be used:\n')
  print(const)
  pp[, const] <- FALSE
  exp.mat <- X
  done <- FALSE
  while (!done) {
    cor.mat <- cor(t(exp.mat))
    diag(cor.mat) <- 0
    cor.mat[is.na(cor.mat)] <- 0
    todo <- which(apply(cor.mat > 0.99, 1, sum) > 0)
    if (length(todo) > 0) {
      remove <- which(cor.mat[todo[1], ] > 0.99)[1]
      cat('Removing predictor', rownames(exp.mat)[remove], 'because it has correlation > 0.99 with', rownames(exp.mat)[todo[1]], '\n')
      pp[, rownames(exp.mat)[todo[1]]] <- pp[, rownames(exp.mat)[todo[1]]] | pp[, rownames(exp.mat)[remove]]
      pp[, rownames(exp.mat)[remove]] <- FALSE
      exp.mat <- exp.mat[-remove, ]
    } else {
      done <- TRUE
    }
  }
  
  out.list <- mclapply(1:G, BBSRforOneGene, X, Y, pp, weights.mat, nS, X.dr, use.deg.rates, dynamics.prior.mat, dynamics.prior.weights, mc.cores=cores)
  return(out.list)
}

BBSRforOneGene <- function(ind, X, Y, pp, weights.mat, nS, X.dr, use.deg.rates, dynamics.prior.mat, dynamics.prior.weights) {
  if (ind %% 100 == 0) {
    cat('Progress: BBSR for gene', ind, '\n')
  }
  
  pp.i <- pp[ind, ]
  
  if (sum(pp.i) == 0) {
    return(list(ind=ind, pp=rep(TRUE, length(pp.i)), betas=0, betas.resc=0))
  }

  # create BestSubsetRegression input  
  y <- as.vector(Y[ind, ], mode="numeric")
  x <- t(matrix(X[pp.i, ], ncol=ncol(X)))
  g <- matrix(weights.mat[ind, pp.i], ncol=sum(pp.i))
  if(use.deg.rates) {
    beta.priors <- dynamics.prior.mat[ind, c(TRUE, pp.i)] #added by Kostya
    x.ind <- X.dr[ind, ]
    g <- cbind(dynamics.prior.weights[ind], g)
  } else { beta.priors <- dynamics.prior.mat[ind, pp.i] }
  
  # experimental stuff
  spp <- ReduceNumberOfPredictors(y, x, g, nS, beta.priors, use.deg.rates, x.ind)
  pp.i[pp.i == TRUE] <- spp
  x <- t(matrix(X[pp.i, ], ncol=ncol(X)))
  g <- matrix(weights.mat[ind, pp.i], ncol=sum(pp.i))
  if(use.deg.rates) {
    beta.priors <- beta.priors[c(TRUE, spp)]
    g <- cbind(dynamics.prior.weights[ind], g)
  } else { beta.priors <- beta.priors[spp] }
  
  #Kostya: here I will modify BestSubsetRegression to take x.ind = X.dr[ind,] as an argument
  betas <- BestSubsetRegression(y, x, g, x.ind, use.deg.rates, beta.priors)
  if(use.deg.rates) { # added by Kostya because now betas has length nS+1
    betas.resc <- PredErrRed(y, cbind(x.ind, x), betas)
  } else {
    betas.resc <- PredErrRed(y, x, betas)
  }
  
  return(list(ind=ind, pp=pp.i, betas=betas, betas.resc=betas.resc))
}


ReduceNumberOfPredictors <- function(y, x, g, n, beta.priors, use.deg.rates, x.ind) {
  K <- ncol(x) + use.deg.rates
  spp <- 
  if (K <= n) {
    return(rep(TRUE, K - use.deg.rates))
  }

  combos <- cbind(diag(K - use.deg.rates) == 1, CombCols(diag(K - use.deg.rates)))
  if(use.deg.rates) {
    combos <- rbind(TRUE, combos)
    x <- cbind(x.ind, x)
  }
  bics <- ExpBICforAllCombos(y, x, g, combos, beta.priors)
  bics.avg <- apply(t(t(combos) * bics), 1, sum)
  if(use.deg.rates) { bics.avg = bics.avg[-1] }
  ret <- rep(FALSE, K - use.deg.rates)
  ret[order(bics.avg)[1:n]] <- TRUE
  
  return(ret)
}

BayesianModelAveraging <- function(y, x, g) {
  
  mprior.size <- g
  bms.out <- bms(cbind(y,x), nmodel=0, mprior='pip', 
                 mprior.size=mprior.size, user.int=F)
  #bms.out <- bms(cbind(y,x), nmodel=0, user.int=F)
  
  tmp <- coef(bms.out)
  tmp <- tmp[order(tmp[, 'Idx']), ]
  
  ret <- tmp[, 'Post Mean']
  ret[tmp[, 'PIP'] < 0.9] <- 0
  #print(sum(tmp[, 'PIP'] > 0.5))
  #return(as.numeric(tmp[, 'PIP'] > 0.9)[order(tmp[, 'Idx'])])
  return(ret)
}

BestSubsetRegression <- function(y, x, g, x.ind, use.deg.rates, beta.priors) {
  # Do best subset regression by using all possible combinations of columns of
  # x as predictors of y. Model selection criterion is BIC using results of
  # Bayesian regression with Zellner's g-prior.
  #
  # Args:
  #   y: dependent variable
  #   x: independent variable
  #   g: value for Zellner's g-prior; can be single value or vector
  #
  # Returns:
  #   Beta vector of best model

  K <- ncol(x) + use.deg.rates
  N <- nrow(x)
  ret <- c()

  combos <- AllCombinations(K - use.deg.rates)
  if(use.deg.rates) {
    combos <- rbind(TRUE, combos)
    x <- cbind(x.ind, x)
  }
  bics <- ExpBICforAllCombos(y, x, g, combos, beta.priors)
  #UNDER CONSTRUCTION
  #if(use.deg.rates) { #added by Kostya, in case we want to include deg.rates in model selection step
  #  bics <- ExpBICforAllCombos(y, x, g, combos)
    
  
  not.done <- TRUE
  while (not.done) {
    best <- which.min(bics)
    #cat('best is', best, '\n')
  
    # For the return value, re-compute beta ignoring g-prior.
    betas <- rep(0, K)
    if (best > 1 | use.deg.rates) { #Kostya: because combos[,1] is the null model with all F's.
      x.tmp <- matrix(x[,combos[, best]], N)
      
      tryCatch({
        bhat <- solve(crossprod(x.tmp), crossprod(x.tmp, y))
        #if(use.deg.rates) {
        #  if(bhat[1] > 0) { bhat[1] = 0 } #because deg rate can't be negative
        #  betas[c(T, combos[, best])] <- bhat
        #} else {
        betas[combos[, best]] <- bhat
        #}
        not.done <- FALSE
      }, error = function(e) {
        if (any(grepl('solve.default', e$call)) & grepl('singular', e$message)) {
          # error in solve - system is computationally singular
          cat(bics[best], 'at', best, 'replaced\n')
          bics[best] <<- Inf
        } else {
          stop(e)
        }
      })
    }
    else {
      not.done <- FALSE
    }
  }
  return(betas)
}


AllCombinations <- function(k) {
  # Create a boolean matrix with all possible combinations of 1:k.
  # Output has k rows and 2^k columns where each column is one combination.
  # Note that the first column is all FALSE and corresponds to the null model.
  if (k < 1) {
    stop("No combinations for k < 1")
  }
  
  N <- 2^k
  out <- matrix(FALSE, k, N)
  out[1, 2] <- TRUE
  
  row <- 2
  col <- 3
  while (col < N) {
    out[row, col] <- TRUE
    
    for (i in 1:(col-2)) {
      out[, col + i] <- out[, col] | out[, i + 1]
    }
    
    row <- row + 1
    col <- col * 2 - 1
  }
  
  return(out)
}


CombCols <- function(m) {
  K <- ncol(m)
  ret <- matrix(TRUE, nrow(m), K * (K - 1) / 2)
  ret.col <- 1
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      ret[, ret.col] <- m[, i] | m[, j]
      ret.col <- ret.col + 1
    }
  }
  return(ret)
}


ExpBICforAllCombos <- function(y, x, g, combos, beta.priors) {
  # For a list of combinations of predictors do Bayesian linear regression,
  # more specifically calculate the parametrization of the inverse gamma 
  # distribution that underlies sigma squared using Zellner's g-prior method.
  # Parameter g can be a vector. The expected value of the log of sigma squared
  # is used to compute expected values of BIC.
  # Returns list of expected BIC values, one for each model.
  K <- ncol(x)
  N <- nrow(x)
  
  C <- ncol(combos)
  bics <- rep(0, C)
  
  # is the first combination the null model?
  first.combo <- 1
  if (sum(combos[, 1]) == 0) {  
    bics[1] <- N * log(var(y))
    first.combo <- 2
  }
  
  # shape parameter for the inverse gamma sigma squared would be drawn from  
  shape <- N / 2
  
  # compute digamma of shape here, so we can re-use it later
  dig.shape <- digamma(shape)

  # pre-compute the crossproducts that we will need to solve for beta  
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
  
  # In Zellner's formulation there is a factor in the calculation of the rate 
  # parameter: 1 / (g + 1)
  # Here we replace the factor with the approriate matrix since g is a vector
  # now.
  var.mult <- matrix(sqrt(1 / (g + 1)), K, K)
  var.mult <- var.mult * t(var.mult)
      
  for (i in first.combo:C){
    comb <- combos[, i]
    x.tmp <- matrix(x[, comb], N)
    k <- sum(comb)
    
    tryCatch({
      # this is faster than calling lm
      bhat <- solve(xtx[comb, comb], xty[comb])
      bnot = beta.priors[comb]
      
      ssr <- sum((y - x.tmp %*% bhat)^2)  # sum of squares of residuals
    
        # rate parameter for the inverse gamma sigma squared would be drawn from
        # our guess on the regression vector beta is all 0 for sparse models
        rate <- (ssr + 
        #        (0 - t(bhat)) %*%
                (bnot - t(bhat)) %*% #bnot business added by Kostya
                (xtx[comb, comb] * var.mult[comb, comb]) %*% 
                t(bnot - t(bhat))) / 2
        #        t(0 - t(bhat))) / 2
        
        # the expected value of the log of sigma squared based on the 
        # parametrization of the inverse gamma by rate and shape
        exp.log.sigma2 <- log(rate) - dig.shape
        
        # expected value of BIC
        bics[i] <- N * exp.log.sigma2 + k * log(N)
      
    }, error = function(e) {
      if (any(grepl('solve.default', e$call)) & grepl('singular', e$message)) {
        # error in solve - system is computationally singular
        bics[i] <<- Inf
      } else {
        stop(e)
      }
    })
    
  }
  
  return(bics)
}



PredErrRed <- function(y, x, beta) {
  # Calculates the error reduction (measured by variance of residuals) of each
  # predictor - compare full model to model without that predictor
  N <- nrow(x)
  K <- ncol(x)
  pred <- beta != 0
  P <- sum(pred)

  # compute sigma^2 for full model
  residuals <- y - x %*% beta
  sigma.sq.full <- var(residuals)

  # this will be the output
  err.red <- rep(0, K)

  # special case if there is only one predictor
  if (P == 1) {
    err.red[pred] <- 1 - sigma.sq.full / var(y)
    return(err.red)
  }

  # one by one leave out each predictor and re-compute the model with the
  # remaining ones
  for (i in (1:K)[pred]) {
    pred.tmp <- pred
    pred.tmp[i] <- FALSE
    x.tmp <- matrix(x[,pred.tmp], N, P-1)
    #bhat <- solve(t(x.tmp) %*% x.tmp) %*% t(x.tmp) %*% y
    tryCatch({
      bhat <- solve(crossprod(x.tmp), crossprod(x.tmp, y))
    }, error = function(e) {
      if (any(grepl('solve.default', e$call)) & grepl('singular', e$message)) {
        # error in solve - system is computationally singular
        cat('PredErrRed: error in solve - system is computationally singular\n')
        browser()
        stop(e)
      } else {
        stop(e)
      }
    })
    residuals <- y - x.tmp %*% bhat
    sigma.sq <- var(residuals)

    err.red[i] <- 1 - sigma.sq.full / sigma.sq
  }
  return(err.red)
}


