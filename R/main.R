#' Estimate the parameters via ERCC spike-in for TASC model.
#'
#' The parameters Alpha and Beta for each cell is estimated by linear regression.
#' @param ercc A vector indicating molecule numbers of ERCC spike-in's in each cell.
#' @param y A matrix or data frame, with each row indicating the read counts of ERCC spike-in's in each cell.
#' @return A two-column data frame including the estimated Alpha and Beta for each cell.
#' @author Rong Ma, University of Pennsylvania <rongm@mail.med.upenn.edu>
#' @export
ab.est = function(ercc, y){
  y = as.matrix(y)
  est.ab = matrix(nrow=dim(y)[1],ncol=2)
  for(i in 1:dim(y)[1]){
    est.ab[i,] = lm(log(y[i,])~log(ercc))$coefficients
  }
  est.ab = as.data.frame(est.ab)
  names(est.ab) = c("Alpha","Beta")
  return(est.ab)
}


#' Recover the gene expression level using ERCC spike-in.
#'
#' Using the TASC model introduced by Jia and Li (2016), with parameters Alpha and Beta for each cell estimated from ERCC spike-in, we can recover the non-zero counts while ignoring the zero counts. The preprocessed data is readily be used for downstream cell type clustering.
#' @param X The UMI based gene expression data set (can be data.frame or matrix) with columns indexing the cell id and the rows indexing the genes.
#' @param est.ab A two-column data frame including the estimated Alpha and Beta for each cell.
#' @return A list including the recovered gene expression data matrix and its logarithmic transformed version, the reduced gene expression data (with 0's substituted by NA's) and its logarithmic transformed version, with the same structure as the input data set X.
#' @author Rong Ma, University of Pennsylvania <rongm@mail.med.upenn.edu>
#' @export
recover <- function(X, est.ab){
  if(length(est.ab$Alpha)!= dim(X)[2]) warning("The number of Alpha should match the number of cells in the data!")
  if(length(est.ab$Beta)!= dim(X)[2]) warning("The number of Beta should match the number of cells in the data!")

  recover.data = matrix(nrow = dim(X)[1], ncol = dim(X)[2])
  for(genenumber in 1:dim(X)[1]){
    for(cellnumber in 1:dim(X)[2]){
  if(X[genenumber, cellnumber] != 0) {
    recover.data[genenumber, cellnumber] = exp((log(X[genenumber, cellnumber])-est.ab$Alpha[cellnumber])/est.ab$Beta[cellnumber])
  } else{
    recover.data[genenumber, cellnumber] = 0}
  }
  if((genenumber/dim(X)[1] <0.25) &((genenumber+1)/dim(X)[1] >= 0.25) ) message("20%")
  if((genenumber/dim(X)[1] <0.5) &((genenumber+1)/dim(X)[1] >= 0.5) ) message("40%")
  if((genenumber/dim(X)[1] <0.75) &((genenumber+1)/dim(X)[1] >= 0.75) ) message("60%")
  if((genenumber/dim(X)[1] <0.98) &((genenumber+1)/dim(X)[1] >= 0.98) ) message("80%")
  }

  recover.data.log = log(recover.data)
  recover.data.log[which(recover.data.log == -Inf)] = unique(sort(recover.data.log))[2]

  reduced.data = recover.data
  for(i in 1:dim(recover.data)[1]){
    for(j in 1:dim(recover.data)[2]){
      if(recover.data[i,j] == 0) reduced.data[i,j] = NA
    }
  }
  reduced.data.log = log(reduced.data)

  message("100%")
  return.data = list(recover.data, recover.data.log, reduced.data, reduced.data.log)
  return(return.data)
}

