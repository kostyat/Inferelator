##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input: cM - colMap
#			r - expression matrix (assumed already normalized)
#			param - a vector of 2 params:
#						1-  c1 - cutoff1: represent the maximal time interval allowed for time series (ts) data
#						2- 'trivial' or 'time_delayed' (param to choose the type of design matrix)
#						3- 'consecutive' or 'all_intervals' (determine if consecutive time measurements [no longer than c1], 
#							 or all permutations of time measurements [up to c1] respectively)
#						4- 'TRUE' or 'FALSE' use_t0_as_steady_state
#						5- 'TRUE' or 'FALSE' use_delt_bigger_than_cutoff_as_steady_state
#
# output: 
#			 steadyStateDesignMat, timeSeriesDesignMat

library(parallel)

get_usr_chosen_design <- function(cM, r, delT_min, delT_max, time_delayed, 
                                  all_intervals, use_t0_as_steady_state, 
                                  use_delt_bigger_than_cutoff_as_steady_state) {
  
  delT_min_vec = rep(delT_min, nrow(cM))
  delT_max_vec = rep(delT_max, nrow(cM))
  
  delT_vec = cM$del.t
  isTs_vec = cM$isTs
  eq_idx   = which(!isTs_vec)
  ts_idx   = which(isTs_vec)
  
  delT_vec 		 = delT_vec		 [ts_idx]
  delT_min_vec = delT_min_vec[ts_idx]
  delT_max_vec = delT_max_vec[ts_idx]
  # set delT_vec: 
  #		0 - last time measurement in ts
  #		>0 - first and middle time measurements in ts
  # following line make 0s indicate first time measurement in ts
  delT_vec[which(is.na(delT_vec))] = 0

  delT_vec_trivial = delT_vec
  # following 2 lines make 0s indicate last time measurement in ts
  delT_vec[-length(delT_vec)] = delT_vec[-1]
  delT_vec[length(delT_vec)] = 0
  
  # data for steady state
  rSS = r[,eq_idx]
  # if we have time series experiments
  if(length(ts_idx)>0) { 
    # data for time series
    if( length(eq_idx) > 0 ){ #simple check to see if there are any equlibrium conditions 
      rTS = r[,-eq_idx]
    }else{
      rTS = r
    }
    
    eq_idx_pseudo = numeric()
    
    # get ts starting conditions, we treat these as equilibrium
    if (use_t0_as_steady_state)
      eq_idx_pseudo = which(delT_vec_trivial == 0)
    
    # get ts conditions with larger than c1 delt, we treat these as equilibrium
    if (use_delt_bigger_than_cutoff_as_steady_state)
      eq_idx_pseudo = c(eq_idx_pseudo, which(delT_vec_trivial > max(delT_max_vec)))
    
    # create design matrix for steady state
    if( length(eq_idx) > 0 ){
      DesignMatSS = cbind(rSS, rTS[,eq_idx_pseudo])
    }else{
      DesignMatSS = rTS[,eq_idx_pseudo,drop=F]
    }
    
    if (all_intervals) { # all permutations time series
      #x 					 = get_all_perms(delT_vec, delT_min_vec[1], delT_max_vec[1]) #assumes delT_min & max are scalars
      x 					 = get_all_perms_vec(delT_vec, delT_min_vec, delT_max_vec) #assumes delT_min & max are vectors
      init_ind 		 = x[[1]]
      boundary_ind = x[[2]]
      DesignMatTS  = rTS[,init_ind]
    } else { # consecutive time series 
      if(time_delayed) {
        DesignMatTS = rTS[,which(delT_vec != 0 & delT_vec <= delT_max_vec)]
      } else {
        DesignMatTS = rTS[,which(delT_vec_trivial != 0 & delT_vec_trivial <= delT_max_vec)]
      }
    }
  } else {
    DesignMatSS = rSS
    DesignMatTS = matrix(0, nrow(DesignMatSS), 0)
  }
  return(list(DesignMatSS,as.matrix(DesignMatTS)))
}


##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input: cM - colMap
#			r - ratios matrix (assumed already normalized)
#			param - a vector of 5 params:
#						1-  c1 - cutoff1: represent the maximal time interval allowed for time series (ts) data
#						2- 'trivial','time_difference','rate','inf_1', or 'inf_1_all_intervals'
#						3- tau
#						4- TRUE or FALSE use_t0_as_steady_state
#						5- TRUE or FALSE use_delt_bigger_than_cutoff_as_steady_state
#
# output:
#			 steadyStateResponseMat timeSeriesResponseMat

get_usr_chosen_response <- function(cM, r, delT_min, delT_max, method, tau, 
                                    use_t0_as_steady_state, 
                                    use_delt_bigger_than_cutoff_as_steady_state) {

  delT_min_vec = rep(delT_min, nrow(cM))
  delT_max_vec = rep(delT_max, nrow(cM))
  
  delT_vec = cM$del.t
  isTs_vec = cM$isTs
  eq_idx   = which(!isTs_vec)
  ts_idx   = which(isTs_vec)
  
  delT_vec     = delT_vec[ts_idx]
  delT_min_vec = delT_min_vec[ts_idx]
  delT_max_vec = delT_max_vec[ts_idx]
  
# set delT_vec: 
#		0 - last time measurement in ts
#		>0 - first and middle time measurements in ts
# following line make 0's indicate first time measurement in ts
  delT_vec[which(is.na(delT_vec))] = 0
  delT_vec_trivial = delT_vec
# following 2 lines make 0's indicate last time measurement in ts
  delT_vec[-length(delT_vec)] = delT_vec[-1]
  delT_vec[length(delT_vec)] = 0
  
  rSS = r[,eq_idx]
  # if we have time series experiments
  if(length(ts_idx)>0) { 
    if(any(eq_idx)){
      rTS = r[,-eq_idx]
    }else{
      rTS=r
    }
    
    eq_idx_pseudo = numeric()
    # get ts starting conditions, we treat these as equilibrium
    if (use_t0_as_steady_state)
      eq_idx_pseudo = which(delT_vec_trivial == 0)
    
    # get ts conditions with larger than c1 delt, we treat these as equilibrium
    if (use_delt_bigger_than_cutoff_as_steady_state)
      eq_idx_pseudo = c(eq_idx_pseudo, which(delT_vec_trivial > max(delT_max_vec)))
    
    if(any(eq_idx)){
      response_matrixSS = cbind(rSS, rTS[,eq_idx_pseudo])
    }else{
      if(any(eq_idx_pseudo)){
        response_matrixSS = rTS[,eq_idx_pseudo,drop=F]
      }else{
        response_matrixSS = rSS
      }
    }
    
    init_ind = which(delT_vec != 0 & delT_vec <= delT_max_vec)
    boundary_ind = init_ind+1
    
    # finished response matrices now go on to response
    if (method == 'trivial') {
      response_matrixTS = rTS[,boundary_ind]
    } 	else if (method == 'time_difference') {
      response_matrixTS = (rTS[,boundary_ind] - rTS[,init_ind])
    }	else if (method == 'rate') {
      response_matrixTS = t(1/delT_vec[init_ind] * t(rTS[,boundary_ind] - rTS[,init_ind]))
    }	else if (method == 'inf_1') {
      response_matrixTS = t(tau/delT_vec[init_ind] * t(rTS[,boundary_ind] - rTS[,init_ind])) + (rTS[,init_ind])
    } else if (method == 'inf_1_all_intervals') {
      #x 					 		  = get_all_perms    (delT_vec, delT_min_vec[1], delT_max_vec[1]) #assumes delT_min and delT_max are scalars
      x 					 		  = get_all_perms_vec(delT_vec, delT_min_vec, delT_max_vec) #assumes delT_min and delT_max are vectors
      init_ind 		      = x[[1]]
      boundary_ind      = x[[2]]
      response_matrixTS = t(tau/delT_vec[init_ind] * t(rTS[,boundary_ind] - rTS[,init_ind])) + (rTS[,init_ind])
    } else {
      stop("unknown response read")
    }
  } else {
    response_matrixSS = rSS
    response_matrixTS = matrix(0, nrow(response_matrixSS), 0)		
  }
  return(list(response_matrixSS ,response_matrixTS))
}



##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input:
#	1- dMSS: design matrix steady state
#	2- dMTS: design matrix time series
#	3- rMSS: response matrix steady state
#	4- dMTS: response matrix time series
#	5- param: what final design matrix? choose from all, ts, or ss

# output:
#	resopnse and corresponding design matrices
make_final_design_and_response_matrix <- function(dMSS, dMTS, rMSS, rMTS, param, cS, r, tf.names, make.des.red.exp = F) {
  if (param == 'all') {
    final_response_matrix = cbind(rMSS, rMTS) 
    final_design_matrix = cbind(dMSS, dMTS)
  } else if (param == 'ts') {
    final_response_matrix = rMTS 
    final_design_matrix = dMTS
  } else if (param == 'ss') {
    final_response_matrix = rMSS		
    final_design_matrix = dMSS
  } else {
    stop("unknown final design or response matrices read")
  }
  
  #check if we have biclusters
  #if we do, make final_response_martix to be the 
  #median of the genes in the bicluster, and add TFs to the end of that
  
  cs.num.rows <- as.vector(unlist(lapply(cS, function(i) i$nrows)))
  if(any(cs.num.rows > 1) & (nrow(r) != length(cS))){
    cat("We have biclusters--generating redExp for response and design--\n")
    #finding if we have any single genesin the clusterStack
    #also adding any genes that are not in clusterStack as singletons
    singletons        <- unlist(lapply( cS[ which(cs.num.rows == 1) ], function(i) i$rows))		
    ag                <- unique(c(unlist(lapply(cS, function(i) i$rows)), tf.names))
    
    singletons.to.add <- unique(c(singletons, rownames(r)[which(!rownames(r) %in% ag)]))
    #tf.names #hack to keep full list of tfs, used in making the design matrix
    if(!is.null(singletons.to.add) & any(singletons.to.add %in% tf.names)) 
      singletons.to.add <- singletons.to.add[-which(singletons.to.add %in% tf.names)]

    #intialize red.exp, which contains the reduced expression (median expression)
    #for reach gene in the cluster
    red.exp <- matrix( NA, length(tf.names) + length(cs.num.rows) + length(singletons.to.add), ncol(final_response_matrix))
    colnames(red.exp) <- colnames(final_response_matrix)
  
    des.red.exp <- matrix( NA, length(tf.names) + length(cs.num.rows) + length(singletons.to.add), ncol(final_response_matrix))
    #initialize resp_idx, 
    #index of i,j contains a 0 if condition i is not in bicluster
    #or the index of j in colnames(final_response_matrix) if condition j is in bicluster i
    resp_idx <- matrix( 0, length(tf.names) + length(cs.num.rows) + length(singletons.to.add), ncol(final_response_matrix))
    
    cs.start.idx <- length(tf.names) + 1
    cs.end.idx   <- cs.start.idx + length(cS) - 1
    mat.idx      <- cs.start.idx:cs.end.idx
    for(i in 1:length(cS)){ 
      #get the rows for the bicluster
      cs.rows <- cS[[i]]$rows
      #get the cols for the bicluster...this if statement is because some biclusters
      #have $cols and some have $conds...WTF
      if( !is.null(cS[[i]]$cols) ){
        cs.cols <- cS[[i]]$cols
      }else{
        cs.cols <- cS[[i]]$conds
      }
      
      #calculate the reduced expression for this bicluster, and put it in red.exp
      #NOTE: some biclusters contain the first time-series...this can't be in the 
      #response matrix, because we use a time lag, so the first time-series is removed
      #THUS, we only take the columns which are in the response matrix
      cs.cols     <- colnames(final_response_matrix)[which(colnames(final_response_matrix) %in% cs.cols)]
      cs.cols.idx <- which(colnames(final_response_matrix) %in% cs.cols)
      
      cs.red.exp <- apply(final_response_matrix[cs.rows, cs.cols, drop=F], 2, mean)
      red.exp[mat.idx[i], cs.cols.idx] <- cs.red.exp
      
      #create the indexing vector for this bicluster and store it in resp_idx
      all.idx <- rep(0, ncol(final_response_matrix))
      all.idx[ which(colnames(final_response_matrix) %in% cs.cols) ] <- which(colnames(final_response_matrix) %in% cs.cols) 
      resp_idx[mat.idx[i], ] <- all.idx
    
      #getting the red_exp for the design matrix
      if(make.des.red.exp)
         des.red.exp[mat.idx[i],] <- apply( final_design_matrix[cs.rows, ,drop=F],2,mean)
    }
    #add the genes which are not in the any clusters to the clusterStack
    sing.start.idx <- length(tf.names) + 1 + length(cS)
    sing.stop.idx  <- sing.start.idx + length(singletons.to.add) - 1
    red.exp[sing.start.idx:sing.stop.idx, ] <- final_response_matrix[singletons.to.add, ]
    
    #add the TFs to red.exp	
    tfs.start.idx = 1
    tfs.stop.idx  = tfs.start.idx + length(tf.names) - 1
    red.exp[tfs.start.idx:tfs.stop.idx, ] <- final_response_matrix[tf.names,]
    #create rownames: its either the bicluster, or the name of the TF, if the bicluster is a singleton then we use the name of the gene
    if (is.null(names(cS))) {
      r.names <- c(paste(sep="","BC_",seq(1:length(cS))))
      r.names <- c(tf.names, r.names, singletons.to.add)#r.names[ which(cs.num.rows == 1) ] <- singletons
    } else {
      r.names <- c(tf.names,names(cS), singletons.to.add)
    }
    rownames(red.exp) <- r.names
    
    #add the indexing vectors for the TFs and singletons, 
    #this is trivial as for the TFs they have expression in every condition
    resp_idx[c(tfs.start.idx:tfs.stop.idx, sing.start.idx:sing.stop.idx) , ] <- t(matrix(rep(1:ncol(final_response_matrix),length(tf.names) + length(singletons.to.add)), 
                                                                                    nrow = ncol(final_response_matrix), ncol= length(tf.names) + length(singletons.to.add)))
    final_response_matrix <- red.exp
    
    #finishing the final_design_matrix
    if(make.des.red.exp){
      #des.red.exp[(length(cS)+1):nrow(des.red.exp), ] <- final_design_matrix[tf.names,]
      des.red.exp[sing.start.idx:sing.stop.idx, ] <- final_design_matrix[singletons.to.add, ]
      #add the TFs to the design red.exp
      des.red.exp[tfs.start.idx:tfs.stop.idx, ] <- final_design_matrix[tf.names, ]
      
      if (is.null(names(cS))) {
        r.names <- c(paste(sep="","BC_",seq(1:length(cS))))
        r.names <- c(tf.names, r.names, singletons.to.add)
        #r.names[ which(cs.num.rows == 1) ] <- singletons
      } else {
        r.names <- c(tf.names, names(cS), singletons.to.add)
      }
      rownames(des.red.exp) <- r.names
      final_design_matrix   <- des.red.exp
    } else {
      final_design_matrix <- final_design_matrix[tf.names, ]
    }
  }else{
    resp_idx <- t(matrix(rep(c(1:ncol(final_response_matrix)), nrow(final_response_matrix)), ncol = nrow(final_response_matrix), nrow=ncol(final_response_matrix)))
  }
    
  return (list(final_response_matrix, final_design_matrix, resp_idx))
}



################################################################################
# CH Jun 3, 2014
# All the functions in this file should be replaced with something sane, but
# I don't have time for this right now.

design.and.response.old <- function(meta.data, exp.mat, delT.min, delT.max, tau,
                                clusterStack, tf.names, time_delayed=T, 
                                all_intervals=F, use_t0_as_steady_state=F, 
                                use_delt_bigger_than_cutoff_as_steady_state=T) {
  # get design matrix
  x <- get_usr_chosen_design(meta.data, exp.mat, delT.min, 
                             delT.max, time_delayed=T, all_intervals=F, 
                             use_t0_as_steady_state=F, 
                             use_delt_bigger_than_cutoff_as_steady_state=T)
  design_matrix_steady_state <- x[[1]]
  design_matrix_time_series <- x[[2]]

  # get response matrix
  x <- get_usr_chosen_response(meta.data, exp.mat, delT.min, 
                               delT.max, 'inf_1', tau, 
                               use_t0_as_steady_state=F, 
                               use_delt_bigger_than_cutoff_as_steady_state=T)
  response_matrix_steady_state <- x[[1]]
  response_matrix_time_series <- x[[2]]

  # make final design/response matrices
  x <- make_final_design_and_response_matrix(design_matrix_steady_state,
                                             design_matrix_time_series,
                                             response_matrix_steady_state,
                                             response_matrix_time_series,
                                             'all', clusterStack, exp.mat, 
                                             tf.names, make.des.red.exp=T)
  
  return(list(final_design_matrix=x[[2]], final_response_matrix=x[[1]], resp.idx=x[[3]]))
}


# CH Aug 15, 2014
# This replaces all of the above code and should be more robust.

# notes
# design matrix is same as exp.mat leaving out last time points
# response matrix is same as design for steady state; linear interpolation else
design.and.response <- function(meta.data, exp.mat, delT.min, delT.max, tau, use.deg.rates, set.deg.rates, deg.rates, tfa.corr, avg.diffs, cores) {
  
  cond <- as.character(meta.data$condName)
  prev <- as.character(meta.data$prevCol)
  delt <- meta.data$del.t
  
  # break time series if del.t is larger than delT.max
  prev[delt > delT.max] <- NA
  delt[delt > delT.max] <- NA
  
  # fix condition names if needed (R is picky about row and column names and 
  # might have used different condition names in expression matrix)
  not.in.mat <- setdiff(cond, colnames(exp.mat))
  cond.dup <- duplicated(cond)
  if (length(not.in.mat) > 0) {
    cond <- gsub('[/+-]', '.', cond)
    prev <- gsub('[/+-]', '.', prev)
    if (!all(cond.dup == duplicated(cond))) {
      stop('Tried to fix condition names in meta data so that they would match column names in expression matrix, but failed')
    }
  }
  # check if there are condition names missing in expression matrix
  not.in.mat <- setdiff(cond, colnames(exp.mat))
  if (length(not.in.mat) > 0) {
    print(not.in.mat)
    stop('Error when creating design and response. The conditions printed above are in the meta data, but not in the expression matrix')
  }
  
  des.mat <- matrix(0, nrow(exp.mat), 0, dimnames=list(rownames(exp.mat), NULL))
  res.mat <- matrix(0, nrow(exp.mat), 0, dimnames=list(rownames(exp.mat), NULL))
  
  # handle all the steady state conditions first
  #Replaced old code with a faster way. This is 400x faster and saves me 3.5 minutes for yeast.
  steady <- is.na(prev) & !(cond %in% prev)
  des.mat <- cbind(des.mat, exp.mat[, cond[steady]])
  colnames(des.mat) <- cond[steady]
  if(set.deg.rates == 'decay') {
    res.mat <- cbind(res.mat, exp.mat[, cond[steady]] * -deg.rates)
  } else {
    res.mat <- cbind(res.mat, exp.mat[, cond[steady]])
  }
  colnames(res.mat) <- cond[steady]
  
  # handle time series
  des.res.ts.list <- mclapply(which(!steady), function(i) {
  #for (i in which(!steady)) {
    # find the conditions that follow this one and are at least delT.min far apart
    following <- which(prev == cond[i])
    following.delt <- delt[following]
    off <- which(following.delt < delT.min)
    while (length(off) > 0) {
      off.fol <- which(prev == cond[following[off[1]]])
      off.fol.delt <- delt[off.fol]
      following <- c(following[-off[1]], off.fol)
      following.delt <- c(following.delt[-off[1]], off.fol.delt + following.delt[off[1]])
      off <- which(following.delt < delT.min)
    }
    
    # for each following condition (time point)
    n <- length(following)
    cntr <- 1
    des.mat.temp <- matrix(0, nrow(exp.mat), 0, dimnames=list(rownames(exp.mat), NULL))
    res.mat.temp <- matrix(0, nrow(exp.mat), 0, dimnames=list(rownames(exp.mat), NULL))
    cnames <- c()
    for (j in following) {
      if (n > 1) {
        this.cond <- sprintf('%s_dupl%02d', cond[i], cntr)
      } else {
        this.cond <- cond[i]
      }
      #des.mat <- cbind(des.mat, exp.mat[, cond[i]])
      #colnames(des.mat)[ncol(des.mat)] <- this.cond
      des.mat.temp <-  cbind(des.mat.temp, exp.mat[, cond[i]])
      cnames <- cbind(cnames, this.cond)
      
      dX.dt <- (exp.mat[, cond[j]] - exp.mat[, cond[i]]) / following.delt[cntr]
      
      if(avg.diffs) {
        if(!is.na(prev[i])) {
          h = which(cond == prev[i]) #this is the previous condition
          if(length(h) > 1) { #this case is most certainly a replicate or a bug
            exp.mat.h <- rowSums(exp.mat[, cond[h]])/length(h)
          } else { exp.mat.h <- exp.mat[, cond[h]] }
          dX.dt <- (dX.dt + (exp.mat[, cond[i]] - exp.mat.h) / delt[i]) / 2
        }
      }
      
      if(!use.deg.rates) { #added by Kostya
        if(set.deg.rates == "decay") {
          interp.res <- tfa.corr * dX.dt - exp.mat[, cond[i]] * deg.rates
      } else if(set.deg.rates == "taus") {
          interp.res <- tfa.corr * deg.rates * dX.dt + exp.mat[, cond[i]]
      } else {
        interp.res <- tfa.corr * tau * dX.dt + exp.mat[, cond[i]]
      } }
      else {
        interp.res <- tfa.corr * dX.dt
      }
      #res.mat <- cbind(res.mat, interp.res)
      #colnames(res.mat)[ncol(res.mat)] <- this.cond
      res.mat.temp <- cbind(res.mat.temp, interp.res)
      
      cntr <- cntr + 1
    }
    
    # special case: nothing is following this condition within delT.min
    # and it is the first of a time series --- treat as steady state
    if (n == 0 & !(cond[i] %in% prev)) {
      #des.mat <- cbind(des.mat, exp.mat[, cond[i]])
      #colnames(des.mat)[ncol(des.mat)] <- cond[i]
      this.cond <- cond[i]
      cnames <- cbind(cnames, this.cond)
      des.mat.temp <- cbind(des.mat.temp, exp.mat[, cond[i]])
      if(set.deg.rates == "decay") {
        #res.mat <- cbind(res.mat, exp.mat[, cond[i]] * -deg.rates)
        res.mat.temp <- cbind(res.mat.temp, exp.mat[, cond[i]] * -deg.rates)
      } else {
        #res.mat <- cbind(res.mat, exp.mat[, cond[i]])
        res.mat.temp <- cbind(res.mat.temp, exp.mat[, cond[i]])
      }
      #colnames(res.mat)[ncol(res.mat)] <- cond[i]
    }
    return(list(des.mat = des.mat.temp, res.mat = res.mat.temp, cnames = cnames))
  }, mc.cores = cores)
  
  #"Unpacking" the time-series list of design and response matrices
  des.ts <- lapply(des.res.ts.list, function(x) { return(x$des.mat) } )
  res.ts <- lapply(des.res.ts.list, function(x) { return(x$res.mat) } )
  des.mat.ts = Reduce(cbind, des.ts)
  res.mat.ts = Reduce(cbind, res.ts)
  cnames.ts <- unlist(lapply(des.res.ts.list, function(x) { return(x$cnames) } ))
  colnames(des.mat.ts) = cnames.ts
  colnames(res.mat.ts) = cnames.ts
  des.mat <- cbind(des.mat, des.mat.ts)
  res.mat <- cbind(res.mat, res.mat.ts)
  resp.idx <- t(matrix(1:ncol(res.mat), ncol(res.mat), nrow(exp.mat)))
  return(list(final_design_matrix=des.mat, final_response_matrix=res.mat, resp.idx=resp.idx))
}


########################################
#KMT April 9 2015
#This function removes extra rows in the meta data file
#Which is useful when deleting extra columns in the expression file:

meta.reduce <- function(meta, expr.cnames) {
  rm.cols = setdiff(meta$condName, expr.cnames)#conditions to remove
  which.rm = which(meta$condName %in% rm.cols)
  #Before we remove these, we need to look at the ones that are time series and change the table accordingly:
  which.ts = which.rm[which(meta$isTs[which.rm])]
  for(rm in which.ts) {
    if(meta$is1stLast[rm] == "f") {
      next.tps = which(meta$prevCol == as.vector(meta$condName[rm]))
      next.last = next.tps[which(meta$is1stLast[next.tps] == "l")]
      next.mid = next.tps[which(meta$is1stLast[next.tps] != "l")]
      #first we deal with the ones that end at the next time point
      #so they need to be converted to a steady state entry:
      meta[next.last, "isTs"] = FALSE
      meta[next.last, "is1stLast"] = "e"
      meta[next.last, "prevCol"] = NA
      meta[next.last,"del.t"] = 0
      #next we approach the ones that stay as time series
      meta[next.mid,"is1stLast"] = "f"
      meta[next.mid,"prevCol"] = NA
      meta[next.mid,"del.t"] = 0
    } else if(meta$is1stLast[rm] == "m") {
      next.tps = which(meta$prevCol == as.vector(meta$condName[rm]))
      meta[next.tps,"prevCol"] = meta[rm,"prevCol"]
      meta[next.tps,"del.t"] = meta[rm,"del.t"] + meta[next.tps,"del.t"]
    } else {
      prev = which(meta$condName == as.vector(meta$prevCol[rm]))
      prev.first = prev[which(meta$is1stLast[prev] == "f")]
      prev.mid = prev[which(meta$is1stLast[prev] != "f")]
      #first we need to deal with the ones that start at the previous pt
      #so they need to be converted to a steady state entry:
      meta[prev.first, "isTs"] = FALSE
      meta[prev.first, "is1stLast"] = "e"
      meta[prev.first, "prevCol"] = NA
      meta[prev.first,"del.t"] = 0
      #then we can approach the ones that stay as time series:
      meta[prev.mid, "is1stLast"] = "l"
    }
  }
  #Now we are ready to remove the extra rows:
  meta.new = meta[-which.rm,]
  return(meta.new)
}
  
