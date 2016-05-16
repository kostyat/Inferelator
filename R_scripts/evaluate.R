require(Matrix)
library("reshape2")

ChristophsPR <- function(ord.idx, gs) {
  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx)))
  rec <- cumsum(gs[ord.idx]) / sum(gs)

  prec <- c(prec[1], prec)
  rec <- c(0, rec)

  auc <- ChristophsAUC(rec, prec)
  return(list(prec=prec, rec=rec, auc=auc))
}


TP.FP <- function(ord.idx, gs) {
  tp <- cumsum(gs[ord.idx])
  fp <- cumsum(rep(1, length(ord.idx))) - tp
  fn <- sum(gs) - tp
  
  auc <- ChristophsAUC(fp/max(fp), tp/max(tp))
  return(list(tp=tp, fp=fp, fn=fn, auc=auc))
}


ChristophsAUC <- function(x, y) {
  dx <- diff(x)
  my <- y[1:(length(y) - 1)] + diff(y) / 2
  return(sum(dx * my))
}


aupr <- function(mat, gs, eval.on.subset=FALSE) {
  rows <- rep(TRUE, nrow(gs))
  cols <- rep(TRUE, ncol(gs))
  if (eval.on.subset) {
    rows <- apply(gs, 1, sum) > 0
    cols <- apply(gs, 2, sum) > 0
  }
  return(ChristophsPR(order(mat[rows, cols], decreasing=TRUE), gs[rows, cols])$auc)
}

auroc <- function(mat, gs, eval.on.subset=FALSE) {
  rows <- rep(TRUE, nrow(gs))
  cols <- rep(TRUE, ncol(gs))
  if (eval.on.subset) {
    rows <- apply(gs, 1, sum) > 0
    cols <- apply(gs, 2, sum) > 0
  }
  return(TP.FP(order(mat[rows, cols], decreasing=TRUE), gs[rows, cols])$auc)
}


summarizeResults <- function(full.dir, eval.on.subset) {
  params.and.input <- paste(full.dir, 'params_and_input.RData', sep='/')
  if (!file.exists(params.and.input)) {
    cat('No params_and_input.RData - skipping', full.dir, '\n')
    return()
  }
  load(params.and.input)
  files <- list.files(full.dir, "combinedconf_.+\\.RData$")
  print(file.path(full.dir, files))
  
  gs <- IN$gs.mat
  
  out <- matrix('', 0, 2)
  out.roc <- matrix('', 0, 2)
  for (res.file in files) {
    load(file.path(full.dir, res.file))
    if(PARS$deg.rates) {
      comb.confs = comb.confs[,-1]
    }
    aupr.tot <- aupr(comb.confs, gs, eval.on.subset)
    out <- rbind(out, c(res.file, aupr.tot))
    auroc.tot <- auroc(comb.confs, gs, eval.on.subset)
    out.roc <- rbind(out.roc, c(res.file, auroc.tot))
    aupr.rand <- sapply(1:100, function(x) {
        return(aupr(comb.confs[,sample(ncol(comb.confs))], gs, eval.on.subset))
    } )
    aupr.rand.95c <- quantile(aupr.rand, c(0.05, 0.95))
    out <- rbind(out, c(aupr.rand.95c[1], aupr.rand.95c[2]))
    auroc.rand <- sapply(1:100, function(x) {
        return(auroc(comb.confs[,sample(ncol(comb.confs))], gs, eval.on.subset))
    } )
    auroc.rand.95c <- quantile(auroc.rand, c(0.05, 0.95))
    out.roc <- rbind(out.roc, c(auroc.rand.95c[1], auroc.rand.95c[2]))
    gs.li = IN$priors[[which(res.file %in% files)]] != 0
    #mode(gs.li) = "numeric"
    #if(PARS$perc.tp < 100) {
    if(sum(gs != gs.li*gs) > 0) {
      gs.lo = gs - gs.li*gs
      #Now we need to exclude all predictions that are in the prior:
      gs.li.inv = abs(1 - gs.li)
      comb.confs = comb.confs * gs.li.inv
      aupr.lo <- aupr(comb.confs, gs.lo, eval.on.subset)
      out <- rbind(out, c("Leave-Out AUPR", aupr.lo))
      auroc.lo <- auroc(comb.confs, gs.lo, eval.on.subset)
      out.roc <- rbind(out.roc, c("Leave-Out AUROC", auroc.lo))
      aupr.lo.rand <- sapply(1:100, function(x) {
          return(aupr(comb.confs[,sample(ncol(comb.confs))], gs.lo, eval.on.subset))
      } )
      aupr.lo.rand.95c <- quantile(aupr.lo.rand, c(0.05, 0.95))
      out <- rbind(out, c(aupr.lo.rand.95c[1], aupr.lo.rand.95c[2]))
      auroc.lo.rand <- sapply(1:100, function(x) {
          return(auroc(comb.confs[,sample(ncol(comb.confs))], gs.lo, eval.on.subset))
      } )
      auroc.lo.rand.95c <- quantile(auroc.lo.rand, c(0.05, 0.95))
      out.roc <- rbind(out.roc, c(auroc.lo.rand.95c[1], auroc.lo.rand.95c[2]))
      }
  }

  # write results for this directory
  write.table(out, file=file.path(full.dir, 'auprs.tsv'), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  write.table(out.roc, file=file.path(full.dir, 'aurocs.tsv'), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
}  

#### The following functions were added by Kostya to make some things easier:

getNW <- function(gs, priors, comb.confs, cutoff) {
  gs.li = abs(priors)
  if(sum(gs != gs.li*gs) > 0) {
    gs.lo = gs - gs.li*gs
    gs.li.inv = abs(1 - gs.li)
    comb.confs.mat = as.matrix(comb.confs * gs.li.inv)
  } else {
    comb.confs.mat = as.matrix(comb.confs)
  }
  comb.confs.mat.1 = comb.confs.mat
  comb.confs.mat.1[order(comb.confs.mat, decreasing = TRUE)[1:cutoff]] = 1
  comb.confs.mat.1[order(comb.confs.mat, decreasing = TRUE)[(cutoff+1):length(comb.confs.mat)]] = 0
  comb.confs.long = melt(comb.confs.mat.1)
  comb.confs.long = comb.confs.long[which(comb.confs.long[,3] == 1),]
  return(list(wide = comb.confs.mat.1, long = comb.confs.long))
}

make.ChristophsPR <- function(gs, priors, comb.confs, eval.on.subset = TRUE) {
  gs.li = abs(priors)
  if(sum(gs != gs.li*gs) > 0) {
    gs.lo = gs - gs.li*gs
    gs.li.inv = abs(1 - gs.li)
    comb.confs = comb.confs * gs.li.inv
    if(eval.on.subset) {
      rows <- apply(gs.lo, 1, sum) > 0
      cols <- apply(gs.lo, 2, sum) > 0
    } else {
      rows <- rep(TRUE, nrow(comb.confs))
      cols <- rep(TRUE, ncol(comb.confs))
    }
    pr = ChristophsPR(order(comb.confs[rows, cols], decreasing=TRUE), gs.lo[rows, cols])
  } else {
    if(eval.on.subset) {
      rows <- apply(gs, 1, sum) > 0
      cols <- apply(gs, 2, sum) > 0
    } else {
      rows <- rep(TRUE, nrow(comb.confs))
      cols <- rep(TRUE, ncol(comb.confs))
    }
    pr = ChristophsPR(order(comb.confs[rows, cols], decreasing=TRUE), gs[rows, cols])
  }
  return(pr)
}

make.TP.FP <- function(gs, priors, comb.confs, eval.on.subset = TRUE) {
  rows <- rep(TRUE, nrow(comb.confs))
  cols <- rep(TRUE, ncol(comb.confs))
  gs.li = abs(priors)
  if(sum(gs != gs.li*gs) > 0) {
    gs.lo = gs - gs.li*gs
    gs.li.inv = abs(1 - gs.li)
    comb.confs = comb.confs * gs.li.inv
    if(eval.on.subset) {
      rows <- apply(gs.lo, 1, sum) > 0
      cols <- apply(gs.lo, 2, sum) > 0
    }
    pr = TP.FP(order(comb.confs[rows, cols], decreasing=TRUE), gs.lo[rows, cols])
  } else {
    if(eval.on.subset) {
      rows <- apply(gs, 1, sum) > 0
      cols <- apply(gs, 2, sum) > 0
    }
    pr = TP.FP(order(comb.confs[rows, cols], decreasing=TRUE), gs[rows, cols])
  }
  return(pr)
}


#this function determines the index of the last element of an array of precisions which is bigger than a certain precision value.
prec.border <- function(prec.array, prec.cutoff) {
  more.than.cutoff = cumsum(prec.array[length(prec.array):1] > prec.cutoff)
  return(length(prec.array) - min(which(more.than.cutoff > 0)) + 1)
}

#This function determines the index of the first element of an array of recalls that's bigger than a certain recall cutoff value.
rec.border <- function(rec.array, rec.cutoff) {
  return(min(which(rec.array > rec.cutoff)))
}

#This function yields how many connections are in the intersection of a list of networks generated with an array of cutoff sizes.
conn.shared <- function(sizes, comb.confs.list, gs, priors) {
  return(sapply(sizes, function(n) {
    net.list = lapply(comb.confs.list, function(comb.confs) {
      return(getNW(gs, priors, comb.confs, n)$wide)
    })
    return(sum(Reduce("*", net.list)))
  }))
}


