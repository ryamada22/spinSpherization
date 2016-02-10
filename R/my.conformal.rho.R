#' Transform shape with curvature information
#'
#' transfromation with curvature
#' @param vertices is a quaternion vector, coordinates of vertices
#' @param faces.v is a 3 x No.faces integer matrix
#' @param rho.v is a scaler vector of length = no. vertices, or no. faces
#' @param face is a logical indicater rho.v is vertices vs. faces
#' @keywords transformation curvature
#' @export

my.conformal.rho <- function(vertices,faces.v,rho.v,face=FALSE){
  xyz.ori <- t(as.matrix(vertices)[2:4,])
  n.v <- length(vertices) # 頂点数
  n.f <- length(faces.v[1,]) # 三角形数
  if(!face){
    rho.f <- rho.fromVtoTri(rho.v,faces.v)
  }else{
    rho.f <- rho.v
  }
  
  
  # 三角形の面積
  edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
	edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	tmp <- edge1 * edge2
  A <- Mod(Im(tmp))/2
  
  # rho の面積重みつき総和は0でないと「閉じ」ない
  s <- sum(A*rho.f)
  rho.f <- rho.f -s*A/sum(A)
  
  
  E <- my.make.E.v(vertices,faces.v,rho.f)
  E.re <- my.qMtorM(E)
  
  lambda.v <- my.inv.pow.2(E.re)[[1]]
  
  L <- my.make.L(vertices,faces.v)
  L.q <- my.make.quatList(L)
  L.re <- my.qMtorM(L.q)
  
  omega <- my.make.omega(vertices,faces.v,lambda.v)
  
  new.vertices <- as.quaternion(matrix(solve(L.re,omega),nrow=4))
  xyz.new <- t(as.matrix(new.vertices)[2:4,])
  mean.new <- apply(xyz.new,2,mean)
  xyz.new. <- t(t(xyz.new)-mean.new)
  max.new. <- max(abs(xyz.new.))
  xyz.new.st <- xyz.new./max.new.
  
  new.q <- as.quaternion(t(cbind(rep(0,n.v),xyz.new.st)))
  #ret <- xyz.new.st[,1]*Hi + xyz.new.st[,2]*Hj + xyz.new.st[,3]*Hk
  
  ret <- list(xyz.new=xyz.new.st,xyz.ori=sp.mesh$xyz,xyz.new.q=new.q,xyz.ori.q=vertices,faces.v=faces.v,E=E,L=L,lambda.v=lambda.v,omega=omega,n.psi=n.psi,rho.fx=rho.fx,rho.v=rho.v,rho.f=rho.f,sp.mesh=sp.mesh)
  ret
}