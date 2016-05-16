#This file will contain functions that deal with sorting through the time series, plotting them and their predictions based on the Inferelator-given TF network.

library(scales)
library(parallel)

#This function takes in a meta data table like IN$meta.data and breaks it up into separate time series, outputting a list of arrays with sample names in temporal order.
TS2List <- function(meta.data) {
  meta.data = as.matrix(meta.data)
  out = lapply(which(meta.data[,"is1stLast"] == "l"), function(i.last) {
    i = i.last
    cond.array = meta.data[i, "condName"]
    time.diffs = meta.data[i, "del.t"]
    while(meta.data[i, "is1stLast"] != "f") {
      i = which(meta.data[, "condName"] == meta.data[i, "prevCol"])[1]
      if(length(i) != 1) {
        stop(cat("gene", meta.data[i, "prevCol"], "is not unique in the meta data matrix"))
      }
      cond.array = c(meta.data[i, "condName"], cond.array)
      time.diffs = c(meta.data[i, "del.t"], time.diffs)
      if(meta.data[i, "is1stLast"] == "f") {break}
    }
    time = cumsum(time.diffs)
    return(list(conds = cond.array, time = time))
  } )
  names(out) = meta.data[which(meta.data[,"is1stLast"] == "l"), "condName"]
  return(out)
}

#USEFUL FOR DEBUGGING (also "GSM356527"):
#ts.list = TS2List(IN$meta.data)
#conds = ts.list[["GSM237066"]]$conds
#time = ts.list[["GSM237066"]]$time
#gene = "YAL001C"
#betas = bs.betas[gene,which(bs.betas[gene,] != 0)]
#names(betas)[1] = gene
#expr.mat = IN$final_design_matrix
#profile = expr.mat[gene, conds]
#pred.mat = expr.mat[names(betas),conds[1:(length(conds)-1)]]
#tau = PARS$tau


#this function predicts and potentially graphs all genes for all time series. act.mat can be a matrix of TF activities. bs.betas is a named matrix of interactions between TFs and targets.
predict.all.genes <- function(expr.mat, act.mat, meta.data, tau, predict.dr, network, min.cond = 3, CORES = 1) {
  ts.list = TS2List(meta.data)
  predictions = matrix(NA, nrow = nrow(expr.mat), ncol = ncol(expr.mat))
  rownames(predictions) = rownames(expr.mat)
  colnames(predictions) = colnames(expr.mat)
  if(predict.dr) { tau = 1 }
  for(ts.name in names(ts.list)) {
    conds = intersect(ts.list[[ts.name]]$conds, intersect(colnames(expr.mat), colnames(act.mat)))
    time = ts.list[[ts.name]]$time[which(ts.list[[ts.name]]$conds %in% conds)]
    #cat(ts.name)
    cat(".")
    if(length(conds) > min.cond) {
      all.preds = mclapply(rownames(expr.mat), function(gene) {
        profile = expr.mat[gene, conds]
        if(sum(diff(profile) != 0) != 0) { #making sure the target gene is not constant:
          if(sum(network[gene,]) > 0) { #we don't need to make predictions based only on decay.
            pred.names = which(network[gene,]==1)
            pred.mat = act.mat[pred.names, conds[1:(length(conds)-1)]]
            if(is.null(nrow(pred.mat))) { #this is necessary in case we are left with one row
              pred.mat = t(as.matrix(pred.mat))
              rownames(pred.mat) = names(which(network[gene,]==1))
            }
            if(predict.dr) {
              pred.mat = rbind("decay term" = expr.mat[gene, conds[1:(length(conds)-1)]], pred.mat)
            }
            if(nrow(pred.mat) >= ncol(pred.mat)) { #a bad hack to avoid comp. singular systems.
              #tryCatch({
              cors = apply(pred.mat, 1, cor, diff(profile))
              cors[which(is.na(cors))] = 0
              cors.sorted = sort(abs(cors), decreasing=T)
              pred.mat = pred.mat[names(cors.sorted)[1:(ncol(pred.mat)-1)], ]
              #}, error = function(err) {
              #  print(gene)
              #  print(ts.name)
              #})
            }
            if(det(crossprod(t(pred.mat))) > 1e-7) { #this matrix must be invertible
            #tryCatch({
              pred = predict.gene(profile, pred.mat, betas = NULL, time, tau)
            #}, error = function(err) {
            #  print(gene)
            #  print(ts.name)
              result = pred$pred.profile
            #})
            } else { result = NA }
            return(result)
          } else { result = NA }
        } else { result = NA }
      }, mc.cores = CORES)
      for(ngene in 1:length(all.preds)) {
        predictions[ngene,conds] = all.preds[[ngene]]
      }
    }
  }
  return(predictions)
}
  

graph.all.genes <- function(folder.name, predictions, expr.mat, act.mat, meta.data, network, min.cond = 3, predict.dr = TRUE) {
  ts.list = TS2List(meta.data)
  dir.create(folder.name)
  for(ts.name in names(ts.list)) {
    conds = intersect(ts.list[[ts.name]]$conds, intersect(colnames(expr.mat), colnames(act.mat)))
    time = ts.list[[ts.name]]$time[which(ts.list[[ts.name]]$conds %in% conds)]
    if(length(conds) > min.cond) {
      pdf(file = paste(paste(folder.name, ts.name, sep="/"), "pdf", sep="."), height=15, width=15, paper="a4r")
      par(mfrow = c(2,3))
      for(gene in rownames(expr.mat)) {
        profile = expr.mat[gene, conds]
        pred.profile = predictions[gene, conds]
        pred.mat = act.mat[which(network[gene,] == 1), conds[1:(length(conds)-1)]]
        if(is.null(nrow(pred.mat))) { #this is necessary in case we are left with one row
          pred.mat = t(as.matrix(pred.mat))
          rownames(pred.mat) = names(which(network[gene,]==1))
        }
        if(predict.dr) {
          pred.mat = rbind("decay term" = expr.mat[gene, conds[1:(length(conds)-1)]], pred.mat)
        }
        if(sum(!is.na(pred.profile)) > 0) {
          graph.gene.ts(gene, time, pred.profile, pred.mat, profile)
        }
      }
      dev.off()
    }
  }
}


#this function plots the predicted and observed RNA level for gene i, as well as its TFs,
#for a given gene and a given time-series.
#TODO: add an option for the degradation term to be calculated from the predicted profile, rather than from the observed profile.
predict.gene <- function(profile, pred.mat, betas, time, tau) {
  time.diffs = diff(time)
  if(is.null(betas)) {
    profile.diffs = diff(profile) * tau / time.diffs
    betas = t(solve(crossprod(t(pred.mat)), pred.mat %*% profile.diffs))
  }
  diffs = c(0, betas %*% pred.mat * time.diffs) / tau
  pred.profile = profile[1] + cumsum(diffs)
  return(list(time = time, pred.profile = pred.profile))
}  

#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jet.colors <- colorRampPalette(c("#00007F", "blue", "#FF7F00", "red"))
  
graph.gene.ts <- function(gene, time, pred.profile, pred.mat, profile, legend.bg="#ffffff88") {
  par(mar=c(5,4,4,4)+0.1) #leaves space on the right for axis 4
  plot.new()
  r = ylim=range(c(profile, pred.profile))
  plot.window(xlim = range(time), ylim = r)
  title(main=gene, xlab="Time (min)", ylab="RNA concentration")
  axis(side=1, at=time, labels=time)
  #plot the predicted gene profile:
  axis(side=2, at=pretty(r))
  lines(time, pred.profile, type="o", col="black", pch="x", lwd=3)
  legend.col <- "black"
  legend.text <- "predicted"
  legend.pch <- "x"
  legend.lwd <- 3
  #plot the observed gene profile:
  par(new=T)
  lines(time, profile, type="o", col="brown", pch="|", lwd=3)
  legend.col <- c(legend.col, "brown")
  legend.text <- c(legend.text, "observed")
  legend.pch <- c(legend.pch, "|")
  legend.lwd <- c(legend.lwd, 3)
  #plot the TF abundances on the right axis:
  par(new=T)
  lr <- c(min(pred.mat), min(pred.mat) + 2*(max(pred.mat) - min(pred.mat)))
  cols = alpha(jet.colors(nrow(pred.mat)), 0.4)
  for(i in 1:nrow(pred.mat)) {
    plot(time[1:(length(time)-1)], pred.mat[i,], axes=F, bty="n", xlab="", ylab="", type="o", pch=".", col=cols[i], ylim=lr, xlim=range(time), lwd=2)
    legend.col <- c(legend.col, cols[i])
    legend.text <- c(legend.text, rownames(pred.mat)[i])
    legend.pch <- c(legend.pch, ".")
    legend.lwd <- c(legend.lwd, 2)
    par(new = T)
  }
  axis(4, at=pretty(lr), col=cols[1], col.axis=cols[1], col.lab=cols[1])
  mtext("TF Activity", 4, 3, col=cols[1], cex=par("cex"))
  legend("topleft", legend.text, col=legend.col, lty=1, pch=legend.pch, bg=legend.bg, lwd=legend.lwd)
}


#This function creates a genes-by-conditions (grouped by time series) table of correlations between predicted and observed profiles.
corr.table <- function(predictions, expr.mat, meta.data, cor.method = "spearman", CORES = 1) {
  ts.list = TS2List(meta.data)
  cor.table = matrix(NA, nrow = nrow(expr.mat), ncol = length(ts.list), dimnames = list(rownames(expr.mat), names(ts.list)))
  for(ts.name in names(ts.list)) {
    conds = intersect(ts.list[[ts.name]]$conds, colnames(expr.mat))
    cor.table[,ts.name] = unlist(mclapply(rownames(expr.mat), function(gene) {
      profile = expr.mat[gene, conds]
      pred.profile = predictions[gene, conds]
      if(sum(!is.na(pred.profile)) > 0) {
        if(cor.method != "OLS.norm") {
          result = cor(profile, pred.profile, method = cor.method)
        } else {
          result = sum((profile - pred.profile)^2)/length(profile)
        }
      } else {
        result = NA
      }
      return(result)
    }, mc.cores = CORES))
    cat("*")
  }
  return(cor.table)
}
        
#####
#####
# Here is a different approach to refitting betas and predicting new data, and probably a faster one (May 5, 2016):

#Create the response function: rns is an array of gene names, expr.mat is the expression matrix (genes-by-conds), time.diffs are the differences of times between consecutive timepoints, and taus is an genes-by-1 matrix of taus for each gene
response.calc <- function(rns, expr.mat, time.diffs, taus) {
  expr.mat = as.data.frame(expr.mat)
  Reduce(rbind, lapply(rns, function(rn) { 
	diff(t(expr.mat[rn,]))/time.diffs*taus[rn,1] + expr.mat[rn,-dim(expr.mat)[2]] } ) )
}

#Recalculating betas, given:
#tfa.mat is the TFA matrix (TF-by-conditions), rns is an array of genes for which to calculate betas, network.mat (genes-by-TFs) is the connectivity matrix denoting which betas to calculate, resp.mat (genes-by-conditions) is the response matrix.
betas.recalc <- function(tfa.mat, rns, network.mat, resp.mat) {
  betas.mat <- matrix(0, nrow=length(rns), ncol=nrow(tfa.mat), dim=list(rns, rownames(tfa.mat)))
  for(rn in rns) {
    pp.i = which(network.mat[rn,]!=0)
    if(length(pp.i) > 0) {
      x <- t(matrix(tfa.mat[pp.i,], ncol=ncol(tfa.mat)))
      y <- as.vector(resp.mat[rn,], mode='numeric')
      xtx <- crossprod(x)
      xty <- crossprod(x, y)
      bhat <- solve(xtx, xty)
      betas.mat[rn, pp.i] = bhat
  } }
  return(betas.mat)
}

#Predict data from the time=0 expression (expr.t0, array of genes) given betas.mat (genes-by-TFs), rns a list of target genes to make predictions on, tfa.mat (TFs-by-conds) calculated from expr.mat or its response function (if available), time.diffs - a list of time differences between time points, and taus - a list of tau values for each gene.
predict.data <- function(betas.mat, rns, tfa.mat, expr.t0, time.diffs, taus, units='min') {
  add.mat = ((betas.mat[rns,] %*% tfa.mat) * time.diffs) / taus[rns,]
  pred.mat = matrix(as.vector(expr.t0), length(rns), 1, dimnames=list(rns,c()))
  for(i in 1:ncol(tfa.mat)) {
    pred.mat = cbind(pred.mat, add.mat[,i] + pred.mat[,i]*(1-time.diffs[i]/taus[rns,]))
  }
  colnames(pred.mat) = paste(c(0,cumsum(time.diffs)), units, sep='')
  #diffs.mat = cbind(0, ((betas.mat[rns,] %*% tfa.mat) - expr.mat[rns, -ncol(expr.mat)]) * time.diffs) / taus[rns,]
  #pred.mat = t(apply(diffs.mat, 1, cumsum)) + expr.mat[rns,1]
  return(pred.mat)
}


train.and.predict <- function(rns, expr.train, expr.valid, time.diffs.train, time.diffs.valid, taus, network, tfa.train.use.response = TRUE) {
  #Calculate training response matrices:
  resp = response.calc(rns, expr.train, time.diffs.train, taus)
  resp.halftau = response.calc(rns, expr.train, time.diffs.train, taus/2)
  #Calculate the TFA for the training data:
  source("R_scripts/tfa.R")
  if(tfa.train.use.response == TRUE) {
    tfa.train = tfa(as.matrix(network)[rns,], as.matrix(expr.train)[rns,-dim(expr.train)[2]], as.matrix(resp.halftau)[rns,])
  } else {
    tfa.train = tfa(as.matrix(network)[rns,], as.matrix(expr.train)[rns,-dim(expr.train)[2]], as.matrix(expr.train)[rns,-dim(expr.train)[2]])
  }
  #Train betas:
  betas.mat = betas.recalc(tfa.train, rns, network, resp)
  #Recalculate predictions:
  pred.train = predict.data(betas.mat,rns,tfa.train,expr.train[rns,1],time.diffs.train,taus)
  #Calculate TFA on unseen/validation data:
  tfa.valid = tfa(as.matrix(network)[rns,], as.matrix(expr.valid)[rns,-dim(expr.valid)[2]], as.matrix(expr.valid)[rns,-dim(expr.valid)[2]])
  #Make new data prediction:
  pred.valid = predict.data(betas.mat, rns, tfa.valid, expr.valid[rns,1], time.diffs.valid, taus)
  return(list(betas = betas.mat, pred.train = pred.train, pred.valid = pred.valid))
}


  

  
  
  





