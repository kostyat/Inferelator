## Bonneau lab 
## NYU - Center for Genomics and Systems Biology

# Call this script with a job config file as arguments
# Example call: Rscript inferelator.R jobs/dream4_cfg.R


require(Matrix)

rm(list=ls())
gc()

source('R_scripts/utils.R')
source('R_scripts/design_and_response.R')
source('R_scripts/priors.R')
source('R_scripts/mi_and_clr.R')
source('R_scripts/bayesianRegression.R')
source('R_scripts/men.R')
source('R_scripts/evaluate.R')
source('R_scripts/tfa.R')


date.time.str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
print(date.time.str)
start.proc.time <- proc.time()

# default job parameters
PARS <- list()

PARS$input.dir <- 'input/yeast/gpl90'

PARS$exp.mat.file <- 'expression_svd.tsv'
PARS$tf.names.file <- 'tf_names.tsv'
PARS$meta.data.file <- 'metadata_gpl90s.tsv'
PARS$priors.file <- 'gold_standard_signed.tsv'
PARS$gold.standard.file <- 'gold_standard.tsv'
PARS$deg.rates.file <- 'sun_dr_norm_KOcnames.tsv' # can also be a file of Taus.
PARS$leave.out.file <- NULL

PARS$job.seed <- 42  # set to NULL if a random seed should be used
PARS$job.seed.prior <- 42
PARS$save.to.dir <- NULL
PARS$num.boots <- 2
PARS$max.preds <- 10 #This is the maximum number of predictors for each gene.
PARS$mi.bins <- 10
PARS$cores <- 8

PARS$delT.max <- 110
PARS$delT.min <- 0
PARS$tau <- 45

PARS$perc.tp <- 0
PARS$perm.tp <- 1
PARS$perc.fp <- 0
PARS$perm.fp <- 1
PARS$pr.sel.mode <- 'random'  # prior selection mode: 'random' or 'tf'

PARS$eval.on.subset <- FALSE

PARS$method <- 'BBSR'  # 'BBSR' or 'MEN'
PARS$prior.weight <- 3

PARS$use.tfa <- TRUE

PARS$deg.rates <- TRUE #added by Kostya; whether or not to use degradation rates
PARS$scale.desres <- TRUE #added by Kostya
PARS$set.deg.rates <- "none" #added by Kostya, 3 options: "taus", "decay", and "none".
#In case of "taus"/"decay" it is assumed that taus/deg.rates are given by PARS$deg.rates.file.

# some of the elastic net parameters that are essentially constants;
# only override in config script if you know what you are doing
PARS$enet.sparseModels <- TRUE    # sparser models
PARS$enet.nCv <- 10               # number of cross-validations
PARS$enet.lambda <- c(0, 1, 100)  # l2 weights
PARS$enet.verbose <- FALSE        # print progress to screen
PARS$enet.plot.it <- FALSE        # generate cross-validation plots
PARS$enet.plot.file.name <- NULL  # file name for plots


# input argument is the job config script which overrides the default parameters
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0) {
  job.cfg <- args[1]
} else {
  #job.cfg <- 'jobs/bsubtilis_eu_201310_stfa_bbsr_22.R'
  #job.cfg <- 'jobs/bsubtilis_eu_201310_stfa_bbsr_1_tp0_fp0.R'
  #job.cfg <- '/home/ch1421/Projects/Rice/inferelator_jobs/ALL_htseq_intersection-strict_noprior.R'
  #job.cfg <- 'jobs/bsubtilis_eu_201310_stfa_bbsr_tf_11_tp50_fp0.R'
  #job.cfg <- 'jobs/bsubtilis_us_201310_stfa_bbsr_11_tp100_fp0.R'
  #job.cfg <- 'jobs/bsubtilis_201408_us_noKO_weight_110_tp_100_fp_000.R'
  #job.cfg <- 'jobs/yeast_kemmer.R'
  #job.cfg <- 'jobs/yeast_affy.R'#_sgd.R'#_timeseries_noDR.R'
  #job.cfg <- 'jobs/yeast_affy_hclust_k5_1.R'
  #job.cfg <- 'jobs/bsubtilis_kostya.R'
  job.cfg <- 'jobs/yeast_affy_PS_tp50_deg_rates.R'
  #job.cfg <- 'jobs/yeast_affy_timeseries_dw1.R'
  #job.cfg <- 'jobs/yeast_affy_timeseries_dw1.4.R'
  #job.cfg <- 'jobs/yeast_affy_timeseries_dw0.5.R'
  #job.cfg <- 'jobs/bsubtilis_us_test.R'
}

# load job specific parameters from input config file
if (!is.null(job.cfg)) {
  source(job.cfg)
}

if (length(args) > 1) {
  #PARS$priors.file <- args[2]
  #PARS$gold.standard.file <- args[3]
  #gs.name <- args[4]
  PARS$exp.mat.file <- args[2]
  PARS$meta.data.file <- args[3]
  #expr.name <- args[4]
  #nclusts <- as.numeric(args[4])
  #PARS$tf.names.file <- args[2]
  PARS$deg.rates.file <- args[4]
  clustn <- as.numeric(args[5])
  #PARS$tau <- as.numeric(args[5])
  #PARS$job.seed <- as.numeric(args[6])
  #PARS$job.seed.prior <- as.numeric(args[7])
  save.to.dir = strsplit(PARS$save.to.dir, '/')[[1]]
  std.last = strsplit(save.to.dir[length(save.to.dir)], '_')[[1]]
  #std.last.new=paste(c(std.last[1:4], sprintf("tau%03.0f",PARS$tau), std.last[5:6]), collapse='_')
  std.last.new=paste(c(std.last[1:4], paste("clust",clustn,sep=''), std.last[5:6]), collapse='_')
  PARS$save.to.dir <- paste(paste(save.to.dir[-length(save.to.dir)], collapse='/'), std.last.new, sep='/')
}

# read input data
IN <- read.input(PARS$input.dir, PARS$exp.mat.file, PARS$tf.names.file, 
                 PARS$meta.data.file, PARS$priors.file, PARS$gold.standard.file,
                 PARS$leave.out.file, PARS$deg.rates.file)

# keep only TFs that are part of the expression data
IN$tf.names <- IN$tf.names[IN$tf.names %in% rownames(IN$exp.mat)]

# order genes so that TFs come before the other genes
gene.order <- rownames(IN$exp.mat)
gene.order <- c(gene.order[match(IN$tf.names, gene.order)],
                gene.order[which(!(gene.order %in% IN$tf.names))])

IN$exp.mat <- IN$exp.mat[gene.order, ]
if (!is.null(IN$priors.mat)) {
  IN$priors.mat <- IN$priors.mat[gene.order, IN$tf.names]
}
if (!is.null(IN$gs.mat)) {
  IN$gs.mat <- IN$gs.mat[gene.order, IN$tf.names]
}

if(PARS$reduce.expr == TRUE) {
  is.gs.mat.row.0 <- rowSums(IN$gs.mat!=0)>0
}


# no meta data given - assume all steady state measurements
if (is.null(IN$meta.data)) {
  IN$meta.data <- trivial.meta.data(colnames(IN$exp.mat))
}

#keep deg. rates for genes only in the expression matrix:
IN$deg.mat <- IN$deg.mat[rownames(IN$deg.mat) %in% rownames(IN$exp.mat),]

# create the dynamics priors matrix
dynamics.prior.mat = matrix(0, nrow = nrow(IN$exp.mat), ncol = length(IN$tf.names))
rownames(dynamics.prior.mat) = rownames(IN$exp.mat)
if (PARS$deg.rates) {
  deg.rate.mat = IN$deg.mat[which(rowSums(!is.na(IN$deg.mat)) > 2),]
  deg.rates = apply(deg.rate.mat,1,median, na.rm=TRUE)
  dynamics.prior.mat = cbind(NA, dynamics.prior.mat)
  dynamics.prior.mat[names(deg.rates),1] = -deg.rates #known deg rates
  dynamics.prior.mat[which(is.na(dynamics.prior.mat[,1])),1] = -median(deg.rates, na.rm=TRUE) #unknown deg rates are set to the mean known deg rate
  #setting the weights of these priors
  deg.rate.sd = apply(deg.rate.mat, 1, sd, na.rm=TRUE)
  dynamics.prior.weights = rep(PARS$dr.weight, nrow(IN$exp.mat))
  names(dynamics.prior.weights) = rownames(IN$exp.mat)
  dynamics.prior.weights[names(deg.rate.sd)] = deg.rate.sd * (PARS$dr.weight/median(deg.rate.sd, na.rm=TRUE))
}

if (PARS$set.deg.rates == "taus" | PARS$set.deg.rates == "decay") {
  if(!is.null(dim(IN$deg.mat))) {
    deg.rate.mat = IN$deg.mat[which(rowSums(!is.na(IN$deg.mat)) > 0),]
    deg.rates = apply(deg.rate.mat,1,median, na.rm=TRUE)
  } else {
    deg.rates = IN$deg.mat
  }
  deg.rates.all = rep(NA, nrow(IN$exp.mat))
  names(deg.rates.all) = rownames(IN$exp.mat)
  if(PARS$set.deg.rates == "taus") {
    const = 1
  } else if(PARS$set.deg.rates == "decay") {
    const = -1
  }
  deg.rates.all[names(deg.rates)] = const * deg.rates
  #deg.rates.all[which(is.na(deg.rates.all))] = const * mean(deg.rates, na.rm=TRUE) #unknown deg rates are set to the mean known deg rate
  #deg.rates.all[which(is.na(deg.rates.all))] = const * exp(mean(log(deg.rates), na.rm=TRUE))
  deg.rates.all[which(is.na(deg.rates.all))] = const * median(deg.rates, na.rm=TRUE)
  IN$decay.or.taus = deg.rates.all
}

# create dummy clusterStack - a real clusterStack is only needed when inferring 
# on bi-clusters
clusterStack <- trivial.cluster.stack(IN$exp.mat)


# set the random seed
if(!is.null(PARS$job.seed)) {
  set.seed(PARS$job.seed, "Mersenne-Twister", "Inversion")
  cat("RNG seed has been set to ", PARS$job.seed, "\n")
} else {
  ignore <- runif(1)
}
SEED <- .Random.seed

if(is.null(PARS$save.to.dir)) {
  PARS$save.to.dir <- file.path(PARS$input.dir, date.time.str)
}
cat("Output dir:", PARS$save.to.dir, "\n")
if (!file.exists(PARS$save.to.dir)){
  dir.create(PARS$save.to.dir, recursive=TRUE)
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# create design and response matrix
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

#des.res <- design.and.response(IN$meta.data, IN$exp.mat, PARS$delT.min, 
#                               PARS$delT.max, PARS$tau, clusterStack, 
#                               IN$tf.names, time_delayed=T, all_intervals=F, 
#                               use_t0_as_steady_state=F, 
#                               use_delt_bigger_than_cutoff_as_steady_state=T)
                               
des.res <- design.and.response(IN$meta.data, IN$exp.mat, PARS$delT.min, 
                              PARS$delT.max, PARS$tau, PARS$deg.rates,
                              PARS$set.deg.rates, IN$decay.or.taus,
                              tfa.corr = 1, PARS$avg.diffs, PARS$cores)

IN$final_response_matrix <- des.res$final_response_matrix
IN$final_design_matrix <- des.res$final_design_matrix
resp.idx <- des.res$resp.idx
                                
if (!all(apply(resp.idx, 1, identical, resp.idx[1,]))) {
    stop('This version of the Inferelator does not support biclusters. Sorry.')
}
    

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# set up the bootstrap permutations
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

IN$bs.pi <- matrix(0, nrow=PARS$num.boots, ncol=ncol(resp.idx))
if (PARS$num.boots == 1) {
  IN$bs.pi[1, ] <- resp.idx[1, ]
} else {
  for (bootstrap in 1:PARS$num.boots) {
    IN$bs.pi[bootstrap, ] <- resp.idx[1, sample(ncol(resp.idx), replace=TRUE)]
  }
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# parse priors parameters and set up priors list
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

IN$priors <- getPriors(IN$exp.mat, IN$tf.names, IN$priors.mat, IN$gs.mat, 
                       PARS$eval.priors.on.subset, PARS$job.seed.prior, PARS$perc.tp,
                       PARS$perm.tp, PARS$perc.fp, PARS$perm.fp, 
                       PARS$pr.sel.mode)


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# TFA specific initialization
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

if(PARS$use.tfa) {
  
  IN$tf.activities <- list()
  
  des.res <- design.and.response(IN$meta.data, IN$exp.mat, PARS$delT.min, 
                                 PARS$delT.max, PARS$tau, PARS$deg.rates,
                                 PARS$set.deg.rates, IN$decay.or.taus,
                                 tfa.corr=0.5, PARS$avg.diffs, PARS$cores)

  IN$final_response_matrix_halftau <- des.res$final_response_matrix
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# main loop
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

for (prior.name in names(IN$priors)) {
  cat('Method:', PARS$method, '\nWeight:', PARS$prior.weight, '\nPriors:', 
      prior.name, '\n')
  prior <- as.matrix(IN$priors[[prior.name]])
  
  # estimate transcription factor activities
  if(PARS$use.tfa) {
    IN$tf.activities[[prior.name]] <- tfa(prior, IN$final_design_matrix, IN$final_response_matrix_halftau)
  }
    
  # set the prior weights matrix
  no.pr.weight <- 1
  if (sum(prior != 0) > 0) {
    if (PARS$prior.weight == no.pr.weight) {
      warning(paste('Priors present, but they will not be used, because \
                    PARS$prior.weight is set to ', no.pr.weight, '.', sep=''), 
              immediate. = TRUE)
    }
    if (PARS$method == 'BBSR') {
      no.pr.weight <- 1 / PARS$prior.weight
    }
  }
  weights.mat <- matrix(no.pr.weight, nrow(IN$exp.mat), length(IN$tf.names))
  weights.mat[prior != 0] <- PARS$prior.weight
  
  if(PARS$reduce.expr == TRUE) {
    weights.mat <- weights.mat[is.gs.mat.row.0,]
  }
  
  betas <- list()
  betas.resc <- list()
  for (bootstrap in 1:PARS$num.boots) {
    cat("Bootstrap", bootstrap, "of", PARS$num.boots, "\n")
    
    # set up bootstrap specific design and response
    X <- IN$final_design_matrix[, IN$bs.pi[bootstrap, ]]
    Y <- IN$final_response_matrix[, IN$bs.pi[bootstrap, ]]
    
    if(PARS$reduce.expr == TRUE) {
      Y <- Y[is.gs.mat.row.0,]
    }
      

    if (nrow(X) > 6000) {
      X <- X[IN$tf.names, ]  # speeds up MI calculation for large datasets
    }
    
    if(PARS$use.tfa) {
      X <- IN$tf.activities[[prior.name]][, IN$bs.pi[bootstrap, ]]
    }

    # fill mutual information matrices
    cat("Calculating MI\n")	
    Ms <- mi(t(Y), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
    diag(Ms) <- 0
    cat("Calculating Background MI\n")
    Ms_bg <- mi(t(X), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
    diag(Ms_bg) <- 0
    
    # get CLR matrix
    cat("Calculating CLR Matrix\n")
    clr.mat = mixedCLR(Ms_bg,Ms)
    dimnames(clr.mat) <- list(rownames(Y), rownames(X))
    clr.mat <- clr.mat[, IN$tf.names]
    
    #the following is added by Kostya:
    #Became mostly unnecessary when I changed one line in design.and.response function
    if(PARS$deg.rates) {
      X.dr <- IN$final_design_matrix[, IN$bs.pi[bootstrap, ]]
      ts.conds <- as.vector(IN$meta.data[which((IN$meta.data$isTs == T) 
        & (((IN$meta.data$del.t >= PARS$delT.min) & (IN$meta.data$del.t <= PARS$delT.max)) 
          | (IN$meta.data$is1stLast == "f"))),"condName"])
    #  X.to.subtract.from.Y <- X.dr
    #  X.to.subtract.from.Y[,setdiff(colnames(X.to.subtract.from.Y), ts.conds)] = 0
    #  Y <- Y - X.to.subtract.from.Y # now for steady-state conditions, Y is still the same (and the same as X and IN$exp.mat), but for time-series conditions, it is tau/del.t*(x_j-x_i)
    }
    
    # DREAM8 induced change:
    #for (tf1 in IN$tf.names) {
    #  for (tf2 in IN$tf.names) {
    #    if (tf1 != tf2) {
    #      #if (clr.mat[tf1, tf2] > clr.mat[tf2, tf1]) {
    #      if (Ms[tf1, tf2] > Ms[tf2, tf1]) {
    #        clr.mat[tf2, tf1] <- min(clr.mat)
    #      } else if (Ms[tf1, tf2] < Ms[tf2, tf1]) {
    #        clr.mat[tf1, tf2] <- min(clr.mat)
    #      }
    #    }
    #  }
    #}
    
    # get the sparse ODE models
    X <- X[IN$tf.names, ]
    cat('Calculating sparse ODE models\n')
    if (PARS$method == 'BBSR') {
      x <- BBSR(X, Y, clr.mat, PARS$max.preds, no.pr.weight, weights.mat, PARS$cores,
    X.dr, ts.conds, use.deg.rates = PARS$deg.rates, dynamics.prior.mat, dynamics.prior.weights, PARS$scale.desres)
    }
    if (PARS$method == 'MEN' ) {
      x <- mclapply(1:nrow(Y), callMEN, Xs=X, Y=Y, 
                    clr.mat=clr.mat, nS=PARS$max.preds, nCv=PARS$enet.nCv,
                    lambda=PARS$enet.lambda, verbose=PARS$enet.verbose, 
                    plot.it=PARS$enet.plot.it, 
                    plot.file.name=PARS$enet.plot.file.name, 
                    weights.mat=weights.mat, no.pr.val=no.pr.weight, 
                    mc.cores=PARS$cores)
    }
    cat('\n')
    
    # our output will be a list holding two matrices: betas and betas.resc
    if(PARS$deg.rates) { # added by Kostya
      bs.betas <- Matrix(0, nrow(Y), nrow(X)+1, 
                       dimnames=list(rownames(Y), c("deg.rate",rownames(X))))
      bs.betas.resc <- Matrix(0, nrow(Y), nrow(X)+1, 
                            dimnames=list(rownames(Y), c("deg.rate",rownames(X))))
      for (res in x) {
        bs.betas[res$ind, c(T,res$pp)] <- res$betas
        bs.betas.resc[res$ind, c(T,res$pp)] <- res$betas.resc
      }
    }
    else {
      bs.betas <- Matrix(0, nrow(Y), nrow(X), 
                       dimnames=list(rownames(Y), rownames(X)))
      bs.betas.resc <- Matrix(0, nrow(Y), nrow(X), 
                            dimnames=list(rownames(Y), rownames(X)))
      for (res in x) {
        bs.betas[res$ind, res$pp] <- res$betas
        bs.betas.resc[res$ind, res$pp] <- res$betas.resc
      }
    }
    betas[[bootstrap]] <- bs.betas
    betas.resc[[bootstrap]] <- bs.betas.resc
    
  }  # end bootstrap for loop
  
  res.file <- paste(PARS$save.to.dir, "/betas_", prior.name, "_", PARS$prior.weight, ".RData", sep="")
  save(betas, betas.resc, file = res.file)
  
  # rank-combine the rescaled betas (confidence scores) of the bootstraps
  confs.file <- sub('/betas_', '/combinedconf_', res.file)
  comb.confs <- Matrix(0, nrow(betas.resc[[1]]), ncol(betas.resc[[1]]), 
                       dimnames=dimnames(betas.resc[[1]]))
  for (beta.resc in betas.resc) {
    comb.confs <- comb.confs + rank(as.matrix(beta.resc), ties.method='average')
  }
  save(comb.confs, file=confs.file)
  
}  # end prior.name loop

if(PARS$reduce.expr == TRUE) {
  IN$priors[[prior.name]] <- IN$priors[[prior.name]][is.gs.mat.row.0,]
  IN$gs.mat <- IN$gs.mat[is.gs.mat.row.0,]
}

PROCTIME <- proc.time() - start.proc.time
save(PARS, IN, SEED, PROCTIME, file = paste(PARS$save.to.dir, "/params_and_input.RData", sep=""))

if (!is.null(IN$gs.mat)) {
  cat('Using gold standard to evaluate results. Evaluate on subset is set to', PARS$eval.on.subset, '. \n')
  summarizeResults(PARS$save.to.dir, PARS$eval.on.subset)
}

