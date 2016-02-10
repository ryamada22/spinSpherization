#' Utility function to apply functions to elements in sparse matrix
#'
#' utility function
#' @param v is a vector or matrix/array that can be vectorized 
#' @param a is a value vector corresponding to v
#' @param func funtion to be applied to the selected elelments by a
#' @param an initial value of calculation
#' @keywords sparse matrix 
#' @export
#' @examples
#' m <- matrix(c(1,2,3,4,5,6),2,3)
#' a <- c(3,4)
#' my.vector.access(m,a)
#' sum(m[a])

my.vector.access <- function(v,a,func=sum,zero=0){
  if(is.vector(v)){
		v <- matrix(v,ncol=1)
	}
	ord <- order(a)
	rle.out <- rle(a[ord])
	num.row <- length(rle.out[[1]])
	num.col <- max(rle.out[[1]])
	tmp1 <- rep(1:num.row,rle.out[[1]])
	tmp2 <- c()
	for(i in 1:num.row){
		tmp2 <- c(tmp2,1:rle.out[[1]][i])
	}
	addr <- tmp1 + num.row*(tmp2-1)
	ret.v <- matrix(0,num.row,ncol(v))
	for(i in 1:ncol(v)){
		if(zero==0){
			tmp.V <- sparseVector(v[ord,i],i=addr,length=num.row*num.col)
			M <- Matrix(tmp.V,num.row,num.col)
		}else{
			M <- matrix(zero,num.row,num.col)
			M[addr] <- v[ord,i]

		}
		ret.v[,i] <- apply(M,1,func)

	}
	return(list(A = rle.out[[2]],V = ret.v))
}