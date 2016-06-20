#This is to be used on original Christoph's design_and_response.R file from August 2015, not the one I have modified

source('R_scripts/design_and_response.R')
library('Matrix')

delT.max <- 4
delT.min <- 2
tau <- 2

#"isTs","is1stLast","prevCol","del.t","condName")
# TRUE	 "f"	 NA	 NA	 "ts1"
# TRUE	 "m"	 "ts1"	 1	 "ts2"
# TRUE	 "m"	 "ts2"	 2	 "ts3"
# TRUE	 "l"	 "ts3"	 5	 "ts4"
# FALSE	 "e"	 NA	 NA	 "ss"

meta = data.frame(isTs=c(TRUE,TRUE,TRUE,TRUE,FALSE),
  is1stLast = c('f','m','m','l','e'),
  prevCol = c(NA,'ts1','ts2','ts3',NA),
  del.t = c(NA, 1, 2, 5, NA),
  condName = c('ts1','ts2','ts3','ts4','ss'))

exp.mat = matrix(c(1:10),nrow=2, byrow=T, dimnames=list(paste("gene",1:2,sep=''), c(paste("ts",1:4,sep=''),'ss')))

des.res <- design.and.response(meta, exp.mat, delT.min, delT.max, tau)

#Outputs: 

#des.res$final_response_matrix
#      ts4 ss      ts1 ts2 ts3
#gene1   7  9 3.666667   5   5
#gene2   8 10 4.666667   6   6

#des.res$final_design_matrix
#      ts4 ss ts1 ts2 ts3
#gene1   7  9   1   3   5
#gene2   8 10   2   4   6

#des.res$resp.idx
#     [,1] [,2] [,3] [,4] [,5]
#[1,]    1    2    3    4    5
#[2,]    1    2    3    4    5


### Weird behavior:
# When delT.max is 5 or higher, ss column gets renamed to V1
# When delT.min is greater than 2 (and everything else is the same), ts2 disappears.
