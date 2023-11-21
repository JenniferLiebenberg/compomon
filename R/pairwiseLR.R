
#' Title
#'
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
plr = function(dat){
  n = nrow(dat)
  p = ncol(dat)
  m = ((1/2)*p*(p-1))
  L = log(dat)


  pair_lr = matrix(data = NA,nrow=n,ncol=m)
  col_names = rep(0,m)


  k=1
  for(i in (1:(p-1))){
    for(j in ((i+1):p)){
      pair_lr[,k] = L[,i]-L[,j]
      col_names[k] = paste0(colnames(dat)[i],"/",colnames(dat)[j])
      k=k+1
    }
  }
  colnames(pair_lr) = col_names



  Y = pair_lr
  for(i in 1:ncol(pair_lr)){
    Y[,i] = pair_lr[,i]-colMeans(pair_lr)[i]
  }
  Y = apply(pair_lr,2, function(y) y - mean(y))

  Y_svd = svd(Y)


  return(Y)
}
