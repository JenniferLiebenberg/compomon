#' Title
#'
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
clr = function(dat){
  n = nrow(dat)
  p = ncol(dat)
  gx = rep(0,n)
  for(i in 1:n){
    gx[i] = (prod(dat[i,]))^(1/p)
  }
  rat = dat/gx
  log_rat = log(dat/gx)

  return(log_rat)
}
