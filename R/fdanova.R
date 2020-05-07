#' Create a Fourier basis of K functions in [0, 1]
#'
#' @param K A positive integer specifying the number of eigenfunctions to generate.
#' @param pts A vector specifying the time points to evaluate the basis functions.
#' @return A K by len(pts) matrix, each column containing a basis function.
#'
#' @examples
#' basis <- create.basis(3)
#' head(basis)
#' @keywords internal
fourier.basis <- function(K,pts = seq(0, 1, length.out = 50)) 
{
    nGrid <- length(pts)
    
    stopifnot(is.numeric(K) && length(K) == 1 && K > 0)

        res <- sapply(seq_len(K), function(k)
            if (k == 1) {
                rep(1, nGrid)
            } else if (k %% 2 == 0) {
                sqrt(2) * sin(k * pi * pts)
            } else {
                sqrt(2) * cos((k - 1) * pi * pts)
            })
        
    res <- matrix(res, ncol = K) # prevent single length pts
    res
}


#' Hypothesis test for high-dimensional data
#' @description Test the mean or differences of means of high-dimensional vectors are zero or not
#' @param X a matrix (one sample) or a list of matrices (two samples or more), with observations contained in rows; equal spacing is assumed when \code{transform} is TRUE
#' @param alpha significance level
#' @param tau the decay parameter, automatically selected if set to \code{NULL}
#' @param B the number of bootstrap replicates
#' @param pairs a matrix with two columns, used when there are more than two populations, each row specifying a pair of populations whose means are compared, and if set to \code{NULL}, comparisons for all pairs are performed
#' @param transform TRUE/FALSE, whether to transform the data into frequency domain via Fourier basis
#' @param K a positive integer, the number of bais functions for transforming the data when \code{transform} is TRUE
#' @return a list of objects: \code{reject} is a boolean variable indicating whether the null is rejected or not. If the null if reject, \code{rej.pairs} optionally gives the pairs of samples that lead to rejection. \code{pvalue} is also returned.
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lopes2019+}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) synfd::dense.fd(mu=0.05*g, X=synfd::gaussian.process(), n=30, m=50))
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' res <- fdtest(X,alpha=0.05,pairs=matrix(1:4,2,2),tau=c(0.4,0.5,0.6))
#' 
#' # get p-value
#' res$pvalue 
#' 
#' # get selected tau
#' res$tau
#' @export
#' 
fdtest <- function(X,alpha=0.05,tau=NULL,B=1000,pairs=NULL,transform=T,K=50,verbose=F)
{
    if(transform)
    {
        
        if(is.matrix(X)) X <- X %*% fourier.basis(K,ncol(X)) / ncol(X)
        else
        {
            M <- ncol(X[[1]])
            h <- 1/(2*M)
            pts <- seq(h,1-h,length.out=M)
            basis <- fourier.basis(K,pts)
            X <- lapply(X,function(z){
                z %*% basis / ncol(z)
            })
        }
    }
    
    return(hdtest(X,alpha=alpha,tau=tau,B=B,pairs=pairs,verbose=verbose))
}
