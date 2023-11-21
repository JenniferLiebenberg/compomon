
#' Axis Predictivity
#'
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
axis_pred = function(dat){
  dat_svd = svd(dat)
  V = (dat_svd$v)
  U = (dat_svd$u)
  D = (dat_svd$d)
  p = ncol(dat)

  tss = diag(V%*%diag(D)%*%t(V))
  ax_pred = NULL

  for(i in 1:(p-1)){
    for(j in (i+1):p){
      ss =  diag(V[,c(i,j)]%*%diag(D[c(i,j)])%*%t(V[,c(i,j)]))
      AP_row = c(i,j,ss/tss)
      ax_pred = rbind(ax_pred,AP_row)
    }}
  ax_pred = as.matrix(ax_pred)
  rownames(ax_pred)=NULL
  colnames(ax_pred) = c("PC1","PC2",colnames(dat))

  return(ax_pred)
}
