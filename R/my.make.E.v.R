#' Make E-matrix
#' 
#' Four Real-valued matrices representing four components of quatenion
#' Returned value is a list of 4 sparse vectors
#' @export

my.make.E.v <- function(vertices,faces.v,rho){
  # 三角形の面積
	edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
	edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	tmp <- edge1 * edge2
  A <- Mod(Im(tmp))/2
	#A <- abs(i(tmp)+j(tmp)+k(tmp))/2
	# 三角形ごとの計算用係数
	coef.a <- -1/(4*A)
	coef.b <- rho/6
	coef.c <- A*rho^2/9
	
	# Rでは四元数を要素とする行列がないので、re,i,j,kごとに正方行列を作ることにする
	E.re <- E.i <- E.j <- E.k <- sparseVector(c(0),i=c(1),length=length(vertices)^2)
  
	e.q <- list()
	e.q[[1]] <- vertices[faces.v[2,]]-vertices[faces.v[3,]]
	e.q[[2]] <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	e.q[[3]] <- vertices[faces.v[1,]]-vertices[faces.v[2,]]
  # すべての頂点ペアについて、すべての三角形ごとに処理をする
	for(i in 1:3){
		for(j in 1:3){
			tmp <- coef.a * e.q[[i]] * e.q[[j]]+ coef.b * (e.q[[j]] -e.q[[i]] ) + coef.c
			addr <- faces.v[i,] + (length(vertices)*(faces.v[j,]-1))
			
			tmp.v <- t(as.matrix(tmp))
			tmp.out <- my.vector.access(tmp.v,addr)
			E.re <- E.re + sparseVector(tmp.out[[2]][,1],tmp.out[[1]],length(vertices)^2)
			E.i <- E.i + sparseVector(tmp.out[[2]][,2],tmp.out[[1]],length(vertices)^2)
			E.j <- E.j + sparseVector(tmp.out[[2]][,3],tmp.out[[1]],length(vertices)^2)
			E.k <- E.k + sparseVector(tmp.out[[2]][,4],tmp.out[[1]],length(vertices)^2)
		}
	}

	return(list(E.re=E.re,E.i=E.i,E.j=E.j,E.k=E.k))
}
