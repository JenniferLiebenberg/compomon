#' Title
#'
#' @param dat
#' @param dat_untransformed
#' @param pcs
#'
#' @return
#' @export
#'
#' @examples
link_cos = function(dat,dat_untransformed,pcs=c(1,2)){


  bplt = diff_ax(dat,dat_untransformed,pcs)

  links = bplt$links

  m = length(links)

  cos_angle = function (Line1, Line2,option=3)
  {
    m1 = coefficients(Line1)[2]
    m2 = coefficients(Line2)[2]

      vector1 = c(1, m1)
      vector2 = c(1, m2)
      num = (vector1[1] * vector2[1] + vector1[2] * vector2[2])
      den = sqrt(vector1[1]^2 + vector1[2]^2) * sqrt(vector2[1]^2 + vector2[2]^2)

      if(option==1){
        cos_angle = num/den
        angle=cos_angle
        angle = (360 * angle)/(2 * pi)
        return(angle)

      }

      if(option==2){
        tan_angle = (m1-m2)/(1+(m1*m2))
        angle = atan(tan_angle)
        return(angle)

      }

      if(option==3){
        inc_angle1 = atan(m1)
        inc_angle2 = atan(m2)
        angle = inc_angle1-inc_angle2
        return(angle)
      }


  }

  angle_mat = matrix(NA,nrow=m,ncol=m)
  colnames(angle_mat) = names(links)
  rownames(angle_mat) = names(links)
  angle_mat_op1 = matrix(NA,nrow=m,ncol=m)
  angle_mat_op2 = matrix(NA,nrow=m,ncol=m)

  for(i in 1:m){
    for(j in 1:m){
      angle_mat[i,j] = cos_angle(links[[i]],links[[j]],option=3)
      angle_mat_op1[i,j] = cos_angle(links[[j]],links[[i]],option=1)
      angle_mat_op2[i,j] = cos_angle(links[[j]],links[[i]],option=2)

    }
  }
  angle_mat = cos(angle_mat)
  rad_angles = acos(angle_mat)
  deg_angles = (360 * rad_angles)/(2 * pi)
  angle_mat = abs(angle_mat)

  pair_lr = plr(dat_untransformed)
  correlations = cor(pair_lr)
  corr_angle_mat = angle_mat
  corr_angle_mat[upper.tri(corr_angle_mat)] = correlations[upper.tri(correlations)]

  abs(angle_mat)-abs(correlations)

  return(list(cosine_angles = angle_mat, angle_degrees = deg_angles, corr_angle_mat=corr_angle_mat,corr=correlations))

}
