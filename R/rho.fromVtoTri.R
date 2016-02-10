#' Calculate rho of triangle from rho values of veritices
#'
#' Rho of triangle
#' @param rho.v is a scaler vector of rho of vertices
#' @param faces.v is a 3 x No.faces integer matrix
#' @keywords rho
#' @export

rho.fromVtoTri <- function(rho.v,faces.v){
  tmp.rho <- matrix(rho.v[faces.v],nrow=3)
  apply(tmp.rho,2,mean)
}