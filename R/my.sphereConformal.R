#' Conformal change from sphere
#'
#' conformal change from sphere
#' @param n.psi is fineness index of triangle mesh of sphere
#' @param rho.fx is a function giving rho value for coordinates on sphere
#' @keywords transformation curvature
#' @export
#' @examples
#' n.psi <- 30 # 球面メッシュの緯度刻み数
#' rho.fx <- function(X){
#'   ret <- sin(X[,3]*pi*2 )*3
#'   return(ret)
#' }
#' out1 <- my.sphereConformal(n.psi,rho.fx)
#' plot.sp.conformal(out1$xyz.ori,out1$faces.v,out1$sp.mesh$n.triangles,out1$rho.f)


my.sphereConformal <- function(n.psi,rho.fx){
  sp.mesh <- my.sphere.tri.mesh(n.psi)
  vertices.mat <- t(sp.mesh$xyz) # ３行行列化
  vertices <- vertices.mat[1,]*Hi + vertices.mat[2,]*Hj + vertices.mat[3,]*Hk
  edges <- sp.mesh$edge
  faces.v <- t(sp.mesh$triangles)
  rho.v <- rho.fx(sp.mesh$xyz) # 頂点における
  out <- my.conformal.rho(vertices,faces.v,rho.v)
  ret <- list(xyz.new=out$xyz.new,xyz.ori=sp.mesh$xyz,xyz.new.q=out$xyz.new.q,xyz.ori.q=vertices,faces.v=faces.v,E=out$E,L=out$L,lambda.v=out$lambda.v,omega=out$omega,n.psi=n.psi,rho.fx=rho.fx,rho.v=rho.v,rho.f=out$rho.f,sp.mesh=sp.mesh)
  ret
}