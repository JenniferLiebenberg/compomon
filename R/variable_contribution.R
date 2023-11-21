
#' Title
#'
#' @param dat
#' @param dat_new
#' @param npcs
#' @param alpha
#' @param type
#'
#' @return
#' @export
#'
#' @examples
var_contrib = function(dat,dat_new,npcs=c(1,2),alpha = 0.01,type="T2"){
  dat = as.matrix(dat)
  S = (1/(n-1)) * t(dat) %*% (dat)


  n = nrow(dat)
  p = ncol(dat)
  r = length(npcs)

  S_svd = svd(S)

  D_svd = diag((S_svd$d)[npcs])
  U_svd = (S_svd$u)[,npcs]
  V_svd = (S_svd$v)[,npcs]
  V_svd_rem = (S_svd$v)[,-npcs]
  Lam_svd = (D_svd%*%t(D_svd))
  D_svd_rem = diag((S_svd$d)[-npcs])
  Lam_svd_rem = D_svd_rem%*%t(D_svd_rem)


  S_eig = eigen(S)
  V_eig = (S_eig$vectors)[,npcs]
  Lam_eig = diag((S_eig$values)[npcs])
  V_eig_rem = (S_eig$vectors)[,-npcs]
  Lam_eig_rem = diag((S_eig$values)[-npcs])

  V_r = V_eig
  Lam = Lam_eig
  V_rem = V_eig_rem
  Lam_rem = Lam_eig_rem


    X = as.matrix(dat_new)

  T_latent = X %*% V_r

  T_latent_rem = X %*% V_rem


  T2_crit = qchisq(p=1-alpha, df=r)
  theta1 = sum(diag(Lam_rem))
  theta2 = sum(diag(Lam_rem)^2)
  k = theta2/theta1
  df_spe = (theta1^2)/(theta2)
  SPE_crit = k*qchisq(p = 1-alpha,df = df_spe)
  phi = ((V_rem%*%t(V_rem))/(SPE_crit))+((V_r%*%solve(Lam)%*%t(V_r))/(T2_crit))

  D = V_r %*% solve(Lam) %*% t(V_r)

  xtilde = (V_rem) %*% t(V_rem) %*% t(X)

  C = V_r %*% t(V_r)
  C_tilde = (V_rem) %*% t(V_rem)

  n_fault = nrow(X)
  p_fault = ncol(X)
  Xi = diag(rep(1,p_fault))

  if(type=="T2"){
    M = D
  }
  if(type == "SPE"){
    M = C_tilde
  }
  if(type == "comb"){
    M = phi
    }

  CDC = matrix(0,n_fault,p_fault)
  for(i in 1:n_fault){
    x = cbind(X[i,])
    for(j in 1:p_fault){
      CDC[i,j] = t(x) %*% expm::sqrtm(M) %*% Xi[,j] %*% t(Xi[,j]) %*% expm::sqrtm(M) %*% x
      }

  }
  CDC = Re(CDC)

  PDC = matrix(0,n_fault,p_fault)
  for(i in 1:n_fault){
    x = cbind(X[i,])
    for(j in 1:p_fault){
      PDC[i,j] = t(x) %*% M %*% (Xi[,j]) %*% t(Xi[,j] ) %*% x
    }
  }


  RBC = matrix(0,n_fault,p_fault)
  for(i in 1:n_fault){
    x = cbind(X[i,])
    for(j in 1:p_fault){
      RBC[i,j] = ((t(Xi[,j])%*%M%*%x)^(2)) / (t(Xi[,j])%*%M%*%Xi[,j])
    }
  }
  return(list(CDC=CDC,PDC=PDC,RBC=RBC))
}
