
#' Title
#'
#' @param dat
#' @param dat_new
#' @param npcs
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
monitor = function(dat,dat_new=NULL,npcs=c(1,2),alpha = 0.01){
dat = as.matrix(dat)
S = (1/(n-1)) * t(dat) %*% (dat)

n = nrow(dat)
p = ncol(dat)
r = length(npcs)

S_svd = svd(S)

D_svd = diag((S_svd$d)[1:r])
U_svd = (S_svd$u)[,1:r]
V_svd = (S_svd$v)[,1:r]
V_svd_rem = (S_svd$v)[,(r+1):p]
Lam_svd = (D_svd%*%t(D_svd))
D_svd_rem = diag((S_svd$d)[(r+1):p])
Lam_svd_rem = D_svd_rem%*%t(D_svd_rem)


S_eig = eigen(S)
V_eig = (S_eig$vectors)[,1:r]
Lam_eig = diag((S_eig$values)[1:r])
V_eig_rem = (S_eig$vectors)[,(r+1):p]
Lam_eig_rem = diag((S_eig$values)[(r+1):p])

V_r = V_eig
Lam = Lam_eig
V_rem = V_eig_rem
Lam_rem = Lam_eig_rem

  X = as.matrix(dat_new)

T_latent = X %*% V_r

T_latent_rem = X %*% V_rem

n_fault = nrow(X)
p_fault = ncol(X)
D = (V_r) %*% solve(Lam) %*% t(V_r)
T2 = rep(0,n_fault)
T2_2 = rep(0,n_fault)
for(i in 1:n_fault){
  x = cbind(X[i,])

  T2[i] = t(T_latent[i,]) %*% solve(Lam) %*% (T_latent[i,])
  T2_2[i] = t(x) %*% D %*% x
}


T2_crit = qchisq(p=1-alpha, df=r)

T2_fault = which(T2>T2_crit)



SPE = rep(0,n_fault)
for(i in 1:n_fault){
  SPE[i] = t(T_latent_rem[i,]) %*% (T_latent_rem[i,])
}


theta1 = sum(diag(Lam_rem))
theta2 = sum(diag(Lam_rem)^2)
k = theta2/theta1
df_spe = (theta1^2)/(theta2)

SPE_crit = k*qchisq(p = 1-alpha,df = df_spe)

SPE_fault = which(SPE>SPE_crit)

phi = ((V_rem%*%t(V_rem))/(SPE_crit))+((V_r%*%solve(Lam)%*%t(V_r))/(T2_crit))
comb = rep(0,n_fault)
for(i in 1:n_fault){
  comb[i] = t(X[i,]) %*% phi %*% (X[i,])
}

g = ((r/T2_crit^2)+(theta2/SPE_crit^2))/((r/T2_crit)+(theta1/SPE_crit))
h = (((r/T2_crit)+(theta1/SPE_crit))^2)/((r/T2_crit^2)+(theta2/SPE_crit^2))


comb_crit = g*qchisq(p = 1-alpha,df = h)

comb_fault = which(comb>comb_crit)

intersect(intersect(T2_fault,SPE_fault),comb_fault)


out = list(T2=T2,T2_crit=T2_crit,T2_fault=T2_fault,SPE=SPE,SPE_crit=SPE_crit,SPE_fault=SPE_fault,combined=comb,combined_crit=comb_crit,combined_fault=comb_fault)

return(out)

}
