#' Title
#'
#' @param dat
#' @param pcs
#'
#' @return
#' @export
#'
#' @examples
dii2 = function(dat,pcs=c(1,2)){

  dat_trans = dat
  n=nrow(dat_trans)
  p=ncol(dat_trans)

  d = matrix(data=0, nrow=n,ncol=n)

  for(i in 1:n){
    for(j in 1:n){
      d[i,j] = sqrt(sum((dat_trans[i,]-dat_trans[j,])^2))
      }
  }

  row_points = princomp(dat_trans,cor=FALSE,scores=TRUE)$scores[,pcs]
  d_est = matrix(data=0, nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      d_est[i,j] = eucl_dist(row_points[i,],row_points[j,])
    }
  }


  out = matrix(data=0,nrow=n,ncol=n)
  out[upper.tri(out)] = d[upper.tri(d)]
  out[lower.tri(out)] = d_est[lower.tri(d_est)]

  return(out)
}




