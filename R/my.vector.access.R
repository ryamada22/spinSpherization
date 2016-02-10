#' Utility function to apply functions to elements in sparse matrix
#'
#' utility function
#' @export


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