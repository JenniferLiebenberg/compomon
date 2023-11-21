#' Title
#'
#' @param dat_transformed
#' @param dat_untransformed
#' @param pcs
#'
#' @return
#' @export
#'
#' @examples
diff_ax = function(dat_transformed,dat_untransformed,pcs=c(1,2)){

  dat_trans = dat_transformed


  n = nrow(dat_trans)
  p = ncol(dat_trans)
  r = length(pcs)


  PCA = princomp(dat_trans,cor=FALSE,scores=TRUE)

  PCA$loadings = -(PCA$loadings)/sqrt(n-1)

  factoextra::fviz_pca_biplot(PCA)

  coords = PCA$loadings[,pcs]

  svdZ = svd(dat_trans)

  Fform = matrix(data=NA,nrow=n,ncol=r)
  Gform = t(svdZ$v[,pcs])
  for(i in 1:r){
    Fform[,i] = cbind(svdZ$d[i]*svdZ$u[,i])
    }

  Matform = Fform%*%Gform
  Matform
  t(Gform)

  Fcov = (svdZ$u[,pcs])
  Gcov = matrix(data=NA,nrow=r,ncol=p)
  for(i in 1:r){
    Gcov[i,] = t(cbind(svdZ$d[pcs[i]]*svdZ$v[,pcs[i]]))
  }
  Matcov = Fcov%*%Gcov
  Matcov
  t(Gcov)

  coordscov = t(Gcov)/sqrt(n-1)

  coordscov[,1] = -coordscov[,1]


  points = Fcov
  points[1,] = -points[1,]

  eucl_dist = function(a,b){

    return(sqrt(sum((a-b)^2)))

  }

  intersect = function (Line1, Line2)
  {
    m1 = coefficients(Line1)[2]
    m2 = coefficients(Line2)[2]
    n1 = coefficients(Line1)[1]
    n2 = coefficients(Line2)[1]
    if (m1 == m2 & m1 != "Inf" & m2 != "Inf") {
      int = "The lines do not intersect, they are parallel"
    }
    if (m1 == "Inf" | m2 != "Inf") {
      int = c(as.numeric(n1), m2 * as.numeric(n1) + n2)
    }
    if (m2 == "Inf" | m1 != "Inf") {
      int = c(as.numeric(n2), m1 * as.numeric(n2) + n1)
    }
    if (m1 == "Inf" & m2 == "Inf") {
      if (as.numeric(n1) == as.numeric(n2)) {
        int = "The lines are coincident"
      }
      else {
        int = "The lines do not intersect, they are parallel"
      }
    }
    if (m1 != "Inf" & m2 != "Inf") {
      x = (n2 - n1)/(m1 - m2)
      y = m1 * x + n1
      int = c(x, y)
      names(int) = c("X", "Y")
    }
    return(int)
  }

  project = function (P, Line)
  {
    m = coefficients(Line)[2]
    n = coefficients(Line)[1]
    v = c(-m, 1)
    P1 = P + v
    x1 = c(P[1],P1[1])
    y1 = c(P[2],P1[2])

    Line1 = lm(y1~x1)
    proj = intersect(Line, Line1)
    dist = eucl_dist(P,proj)
    out = c(proj,dist)
    names(out) = c("X", "Y","Length")

    return(out)
  }

  link_len = (p*p)/2 - p/2
  link_dist = rep(0,link_len)

  k=1
  for(i in (1:(p-1))){
    for(j in ((i+1):p)){
       link_dist[k]=eucl_dist(coordscov[i,],coordscov[j,])
       names(link_dist)[k]=paste0(rownames(coords)[i],"/",rownames(coords)[j])
       k=k+1
    }
  }

  link_line = list()
  k=1
  for(i in (1:(p-1))){
    for(j in ((i+1):p)){
      x = c(coordscov[i,1],coordscov[j,1])
      y = c(coordscov[i,2],coordscov[j,2])
      link_line[[k]]=lm(y~x)
      names(link_line)[[k]]=paste0(rownames(coords)[i],"/",rownames(coords)[j])
      k=k+1
    }
  }

  points_proj = list()
  k=1
  for(i in 1:nrow(points)){
    for(j in 1:length(link_line)){
      points_proj[[k]] = project(points[i,],link_line[[j]])
      names(points_proj)[[k]]=paste0("point",i,"/",names(link_line)[j])
      k=k+1
      }
  }

  origin_proj = list()

  k=1
  for(j in 1:length(link_line)){
    origin_proj[[k]] = project(c(0,0),link_line[[j]])
    names(origin_proj)[[k]]=paste0("(0,0) on ",names(link_line)[j])
    k=k+1
  }

  pair_lr = plr(dat_untransformed)

  standard_deviations = apply(pair_lr,2,sd)

  est_sd = link_dist
  standard_deviations[(1:5)]
  est_sd[(1:5)]
  sd_sdEst = matrix(data=0,nrow=p,ncol=p)
  colnames(sd_sdEst)=names(dat_untransformed)
  rownames(sd_sdEst)=names(dat_untransformed)
  k=1
  for(i in (1:(p-1))){
    for(j in ((i+1):p)){
      sd_sdEst[i,j]=standard_deviations[k]
      k=k+1
    }
  }
  sd_sdEst[lower.tri(sd_sdEst)] = est_sd



  return(list(sd_sdEst=sd_sdEst,link_lengths=link_dist,projections=points_proj,origin_projection=origin_proj,links=link_line,rowpoints=points,colpoints=coordscov, estimated_sd=est_sd, actual_sd=standard_deviations))
}
















