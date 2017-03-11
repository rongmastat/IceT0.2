#' Estimate the parameters via ERCC spike-in for TASC model.
#'
#' The parameters Alpha and Beta for each cell is estimated by linear regression.
#' @param x A vector indicating molecule numbers of ERCC spike-in's in each cell.
#' @param y A matrix or data frame, with each row indicating the read counts of ERCC spike-in's in each cell.
#' @return A two-column data frame including the estimated Alpha and Beta for each cell.
#' @author Rong Ma, University of Pennsylvania <rongm@mail.med.upenn.edu>
#' @export
ab.est = function(x, y){
  if(is.list(y)) warning("y should be a matrix or data frame, but shouldn't be a list object!")
  y = as.matrix(y)
  est.ab = matrix(nrow=dim(y)[2],ncol=2)
  for(i in 1:dim(y)[2]){
    est.ab[i,] = lm(log(y[,i]+0.5)~log(x))$coefficients
  }
  est.ab = as.data.frame(est.ab)
  names(est.ab) = c("Alpha","Beta")
  return(est.ab)
}


#' Recover the gene expression level using ERCC spike-in.
#'
#' Using the TASC model introduced by Jia and Li (2016), with parameters Alpha and Beta for each cell estimated from ERCC spike-in, we can recover the non-zero counts while ignoring the zero counts. The preprocessed data is readily be used for downstream cell type clustering.
#' @param X The gene expression data set (can be data.frame or matrix) with columns indexing the cell id and the rows indexing the genes.
#' @param est.ab A two-column data frame including the estimated Alpha and Beta for each cell.
#' @return A list including the recovered gene expression data matrix and its logarithmic transformed version, the reduced gene expression data (with 0's substituted by NA's) and its logarithmic transformed version, with the same structure as the input data set X.
#' @author Rong Ma, University of Pennsylvania <rongm@mail.med.upenn.edu>
#' @export
recover <- function(X, est.ab){
  if(length(est.ab$Alpha)!= dim(X)[2]) {warning("The number of Alpha should match the number of cells in the data!")
    break}
  if(length(est.ab$Beta)!= dim(X)[2]) {warning("The number of Beta should match the number of cells in the data!")
    break}

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
  return.data = list(recover=recover.data, recover.log=recover.data.log, reduced=reduced.data, reduced.log=reduced.data.log)
  return(return.data)
}

#' Identification of Cell Types using nonparametric density estimation
#'
#' Nonparametric density estimation is used for cell type clustering. The intermediate parameters such as bandwith and signal-to-noise ratio are chosen by build-in algorithms but can also be user-determined.
#' @param recover.data The recovered or imputated gene expression data set (can be data.frame or matrix) with columns indexing the cell id and the rows indexing the genes.
#' @param ratio Numeric; the signal-to-noise ratio thresholding the noise level.
#' @param threshold Numeric between 0 to 1; the threshold that eliminate all the genes below that quantile in terms of mean expression levels (default: 0.5)
#' @param fdr Numeric; the false discovery rate used in bandwith selection, the smaller fdr is the less genes would be selected in the end (default: 0.002)
#' @param imputation A logical indicating whether the imputated data is used (default: FALSE). If `FALSE`, the list object generated from the `recover` function in IceT package should be used.
#' @param pca logical; Whether an initial PCA step should be performed (default: TRUE)
#' @param max_iter  integer; Number of iterations in Rtsne function (default: 1000)
#' @param numeric; Perplexity parameter in Rtsne function (default: 50)
#' @param eps Reachability distance in dbscan function (default: 2)
#' @param MinPts Reachability minimum no. of points in dbscan function (default: 5)
#' @return points A matrix indicating the coordinates of the visualization of the cells on a 2 dimensional plane.
#' @return cluster A vector indicating the group index for each cell.
#' @return selected_genes A vector of the index of the informative genes selected by IceT.
#' @return num_of_clusters A numeric value indicating the number of clusters generated by IceT (including singletons as an extra cluster).
#' @author Rong Ma, University of Pennsylvania <rongm@mail.med.upenn.edu>
#' @export
icet <- function(recover.data, ratio=20, threshold = 0.5, fdr = 0.002, imputation = FALSE, pca = TRUE, max_iter = 1000, perplexity=50, eps = 2, MinPts = 5){
  #preparations
  localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
      y <- y[-1]
    }
    y
  }
  require(dplyr)
  require(Rtsne)
  require(fpc)

  if(imputation == TRUE) {

  recover.data[recover.data<0.6] = 0
  recover.data.log = log(recover.data)
  recover.data.log[which(recover.data.log == -Inf)] = unique(sort(recover.data.log))[2]

  #get reduced data
  reduced.data = recover.data
  for(i in 1:dim(recover.data)[1]){
    for(j in 1:dim(recover.data)[2]){
      if(recover.data[i,j] == 0) reduced.data[i,j] = NA
    }
  }
  reduced.data.log = log(reduced.data)

  return.data = list(recover=recover.data, recover.log=recover.data.log, reduced=reduced.data, reduced.log=reduced.data.log)
  }


  if(imputation == TRUE){
    #initial elimination
    cv = c()
    colmean = recover.data %>% apply(1,mean)
    colsd = recover.data %>% apply(1,sd)
    for(i in 1:dim(recover.data)[1]){
      cv[i] = colsd[i]/max(colmean[i],0.5)
    }
    target.index = which(colmean>=quantile(colmean,threshold))
    null.index = target.index[which(cv[target.index]<quantile(cv[target.index],0.5))]
    target.index = target.index[which(cv[target.index]>quantile(cv[target.index],0.5))]

    message("Initial gene selection succeed!")

    #banwidth selection
    bwd = seq(1,3,0.1)
    criterion = c()
    for(p in 1:length(bwd)){
      modal = c()
      for(i in 1:length(null.index)){
        if(is.na(sd(reduced.data.log[null.index[i],],na.rm=T))) {
          message("Error: Larger threshold should be used!")
          break}
        den = density(reduced.data.log[null.index[i],],
                      if(sd(reduced.data.log[null.index[i],],na.rm=T)==0) bw = 1
                      else
                        bw = sd(reduced.data.log[null.index[i],],na.rm=T)/bwd[p],na.rm=T)
        test = c()
        poiw = sort(den$y[localMaxima(den$y)],decreasing = T)
        for(j in 1:length(localMaxima(den$y))){
          test[j] = (poiw[1]/poiw[j] < ratio)
        }
        modal[i] = sum(test)
      }
      criterion[p] = sum(modal>1,na.rm=T)
    }
    bwdindex = max(which(criterion<fdr*dim(recover.data)[1]))

    message("Bandwidth selection succeed!")

    modal=c()
    for(i in 1:length(target.index)){
      den = density(recover.data.log[target.index[i],],
                    if(sd(recover.data.log[target.index[i],],na.rm=T)==0) bw = 1
                    else
                      bw = sd(recover.data.log[target.index[i],],na.rm=T)/bwd[bwdindex],na.rm=T)
      test = c()
      poiw = sort(den$y[localMaxima(den$y)],decreasing = T)
      for(j in 1:length(localMaxima(den$y))){
        test[j] = (poiw[1]/poiw[j] < ratio)
      }
      modal[i] = sum(test)
    }
    total.detection.best.log = target.index[which(modal>1)]

    clustset = recover.data.log[total.detection.best.log,]
    rtsne = Rtsne(t(clustset),  check_duplicates = F, pca = pca, max_iter = max_iter, perplexity=perplexity)

    dbscan.detect = dbscan(rtsne$Y, eps = eps, MinPts = MinPts)
    cluster = (dbscan.detect$cluster)
    num_of_clusters = length(levels(as.factor(cluster)))

    result.icet = list(data = return.data, points = rtsne$Y, cluster = cluster, selected_genes = total.detection.best.log, num_of_clusters = num_of_clusters)
  } else{
    #initial elimination
    cv = c()
    colmean = recover.data$recover %>% apply(1,mean)
    colsd = recover.data$recover %>% apply(1,sd)
    for(i in 1:dim(recover.data$recover)[1]){
      cv[i] = colsd[i]/max(colmean[i],0.5)
    }
    target.index = which(colmean>=quantile(colmean,threshold))
    null.index = target.index[which(cv[target.index]<quantile(cv[target.index],0.5))]
    target.index = target.index[which(cv[target.index]>quantile(cv[target.index],0.5))]

    message("Initial gene selection succeeded!")
    #banwidth selection
    bwd = seq(1,3,0.1)
    criterion = c()
    for(p in 1:length(bwd)){
      modal = c()
      for(i in 1:length(null.index)){
        den = density(recover.data$reduced.log[null.index[i],],
                      if(sd(recover.data$reduced.log[null.index[i],],na.rm=T)==0) bw = 1
                      else
                        bw = sd(recover.data$reduced.log[null.index[i],],na.rm=T)/bwd[p],na.rm=T)
        test = c()
        poiw = sort(den$y[localMaxima(den$y)],decreasing = T)
        for(j in 1:length(localMaxima(den$y))){
          test[j] = (poiw[1]/poiw[j] < ratio)
        }
        modal[i] = sum(test)
      }
      criterion[p] = sum(modal>1,na.rm=T)
    }
    bwdindex = max(which(criterion<fdr*dim(recover.data$recover)[1]))

    message("Bindwidth selection succeeded!")
    message(bwd[bwdindex])

    modal=c()
    for(i in 1:length(target.index)){
      den = density(recover.data$recover.log[target.index[i],],
                    if(sd(recover.data$recover.log[target.index[i],],na.rm=T)==0) bw = 1
                    else
                      bw = sd(recover.data$recover.log[target.index[i],],na.rm=T)/bwd[bwdindex],na.rm=T)
      test = c()
      poiw = sort(den$y[localMaxima(den$y)],decreasing = T)
      for(j in 1:length(localMaxima(den$y))){
        test[j] = (poiw[1]/poiw[j] < ratio)
      }
      modal[i] = sum(test)
    }
    total.detection.best.log = target.index[which(modal>1)]

    clustset = recover.data$recover.log[total.detection.best.log,]
    rtsne = Rtsne(t(clustset), check_duplicates = F, pca = pca, max_iter = max_iter, perplexity=perplexity)

    dbscan.detect = dbscan(rtsne$Y, eps = eps, MinPts = MinPts)
    cluster = (dbscan.detect$cluster)
    num_of_clusters = length(levels(as.factor(cluster)))

    result.icet = list(points = rtsne$Y, cluster = cluster, selected_genes = total.detection.best.log, num_of_clusters = num_of_clusters)
    }

  return(result.icet)
  }


#' Cell Type Visualization via IceT
#'
#' Make plots based on results of IceT.
#' @param icet The list object returned from the icet function.
#' @param pch A graphical parameter, specifying the characteristics of plotting points.
#' @param reset A logical (default: FALSE), if TRUE, the user can set other parameters provided the recover.data is included.
#' @param recover.data The recovered data generated by recover function; or imputated data. Used only when reset = TRUE.
#' @param imputation A logical (default: FALSE) indicating whether the imputated data is used. Used only when reset = TRUE.
#' @param max_iter integer; Number of iterations in Rtsne function (default: 1000). Can only be changed when reset = TRUE.
#' @param numeric; Perplexity parameter in Rtsne function (default: 50). Can only be changed when reset = TRUE.
#' @param eps Reachability distance in dbscan function (default: 2). Can only be changed when reset = TRUE.
#' @param MinPts Reachability minimum no. of points in dbscan function (default: 5). Can only be changed when reset = TRUE.
#' @param set.color Prescribe the colors for the data representations. Defaul as the colors by IceT.
#' @author Rong Ma, University of Pennsylvania <rongm@mail.med.upenn.edu>
#' @export
icet.plot <- function(icet, pch = 20, reset = FALSE, recover.data = NA, imputation = FALSE, max_iter = 1000, perplexity = 50,
                      eps = 2, MinPts = 5, set.color = NA){
if(reset == TRUE){
    clustset = recover.data$recover.log[icet$selected_genes,]
    rtsne = Rtsne(t(clustset), check_duplicates = F, max_iter = max_iter, perplexity=perplexity)
    dbscan.detect = dbscan(rtsne$Y, eps = eps, MinPts = MinPts)
    if(is.na(set.color)) {
      plot(rtsne$Y, col = rainbow(length(levels(as.factor(dbscan.detect$cluster))))[as.factor(dbscan.detect$cluster)], pch = pch)
      result.plot = list(points = rtsne$Y, cluster = dbscan.detect$cluster)
      return(result.plot)
      } else{
      plot(rtsne$Y, col = set.color, pch = pch)}
      result.plot = list(points = rtsne$Y, cluster = dbscan.detect$cluster)
      return(result.plot)
} else{

    if(is.na(set.color)) {
      plot(icet$points[,1], icet$points[,2], col = rainbow(icet$num_of_clusters)[as.factor(icet$cluster)], pch = pch)
    } else{
      plot(icet$points[,1], icet$points[,2], col = set.color, pch = pch)
    }
}

}
