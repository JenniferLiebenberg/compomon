#' Title
#'
#' @param dat
#' @param dat_untransformed
#'
#' @return
#' @export
#'
#' @examples
prop_relat = function(dat,dat_untransformed){
   dat_trans = plr(dat_untransformed)
   n=nrow(dat_trans)
   p=ncol(dat_trans)

   ratio_mat = matrix(NA,nrow=p,ncol=p)
   colnames(ratio_mat)=colnames(dat_trans)
   rownames(ratio_mat)=colnames(dat_trans)

   lengths = diff_ax(dat,dat_untransformed)$link_lengths

   for(i in 1:p){
     for(j in 1:p){

        ratio_mat[i,j] = lengths[i]/lengths[j]
     }
   }

   return(ratio_mat)

   }
