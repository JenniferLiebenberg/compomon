#' Title
#'
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
clr_double_centered = function(dat){
  n = nrow(dat)
  p = ncol(dat)

  L = log(dat)
  L_idot = colMeans(L)
  L_dotj = rowMeans(L)
  L_dotdot = mean(L_dotj)

  Z = L
  for(i in 1:n){
    for(j in 1:p){
      Z[i,j] = L[i,j]-L_idot[j]-L_dotj[i]+L_dotdot
      }
  }

  return(Z)
  }
