#' Deformation functions
#'
#' @export


my.deform.k.serial <- function(vertices,faces.v,k,n.iter=10){
  v.list <- rho.cot.list <- list()
  v.list[[1]] <- vertices
  rho.cot.list[[1]] <- my.curvature.cot(v.list[[1]],faces.v)
  for(i in 1:n.iter){
    
    tmp.rho.f <- rho.cot.list[[i]][[3]] * Mod(Im(rho.cot.list[[i]][[2]]))
    tmp.out <- my.conformal.rho(v.list[[i]],faces.v,k*tmp.rho.f,face=TRUE)
    v.list[[i+1]] <- tmp.out$xyz.new.q
    rho.cot.list[[i+1]] <- my.curvature.cot(v.list[[i+1]],faces.v)
  }
  return(list(v.list=v.list,rho.cot.list=rho.cot.list,k=k))
}

#' @export
my.deform.k.serial.m <- function(vertices,faces.v,k,m=2,n.iter=10){
  v.list <- rho.cot.list <- list()
  v.list[[1]] <- vertices
  rho.cot.list[[1]] <- my.curvature.cot(v.list[[1]],faces.v)
  for(i in 1:n.iter){
    
    #tmp.rho.f <- rho.cot.list[[i]][[3]] * Mod(Im(rho.cot.list[[i]][[2]]))
    tmp.rho.f <- rho.cot.list[[i]][[3]] * Mod(Im(rho.cot.list[[i]][[2]]))^m
    #sign.tmp.rho.f <- sign(tmp.rho.f)
    #tmp.rho.f.m <- sign.tmp.rho.f * abs(tmp.rho.f)^2
    #tmp.out <- my.conformal.rho(v.list[[i]],faces.v,k*tmp.rho.f,face=TRUE)
    tmp.out <- my.conformal.rho(v.list[[i]],faces.v,k*tmp.rho.f,face=TRUE)
    v.list[[i+1]] <- tmp.out$xyz.new.q
    rho.cot.list[[i+1]] <- my.curvature.cot(v.list[[i+1]],faces.v)
  }
  return(list(v.list=v.list,rho.cot.list=rho.cot.list,k=k))
}
