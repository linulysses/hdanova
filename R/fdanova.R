#' Hypothesis Test for Densely and Regularly Observed Functional Data
#' @description Test the mean or differences in means of functional data are zero or not.
#' @param X a matrix (one-sample) or a list of matrices (multiple-samples), with each row representing observations from a function.
#' @param alpha significance level; default value: 0.05.
#' @param tau real number(s) in the interval \code{[0,1)} that specifies the decay parameter and is automatically selected if it is set to \code{NULL} or multiple values are provided; default value: \code{NULL}, which is equivalent to \code{tau=1/(1+exp(-0.8*seq(-6,5,by=1))).}
#' @param B the number of bootstrap replicates; default value: \code{ceiling(50/alpha)}.
#' @param pairs a matrix with two columns, only used when there are more than two populations, where each row specifies a pair of populations for which the SCI is constructed; default value: \code{NULL}, so that SCIs for all pairs are constructed.
#' @param transform TRUE/FALSE, whether to transform the data into frequency domain via a basis; default: TRUE.
#' @param basis basis for transformation, for which possible options are \code{'fourier'} and \code{'eigen'}; default value: \code{'eigen'}.
#' @param K a positive integer specifying the number of basis functions for transforming the data when \code{transform} is TRUE.
#' @param tau.method the method to select tau; possible values are 'MGB' (default), 'MGBA', 'WB' and 'WBA' (see \code{\link{hdsci}}).
#' @return a list that includes all objects returned by \code{\link{hdsci}} and the following additional objects:
#'      \describe{
#'          \item{\code{reject}}{a T/F value indicating whether the hypothesis is rejected.}
#'          \item{\code{accept}}{a T/F value indicating whether the hypothesis is rejected.}
#'          \item{\code{rej.paris}}{optionally gives the pairs of samples that lead to rejection.}
#'          \item{\code{pvalue}}{the p-value of the test.}
#'          }
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lopes2020}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) synfd::reg.fd(mu=0.05*g, X=synfd::gaussian.process(), n=30, m=50)$y)
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' res <- fdtest(X,alpha=0.05,pairs=matrix(1:4,2,2),tau=c(0.4,0.5,0.6))
#' 
#' # get p-value
#' res$pvalue 
#' 
#' # get selected tau
#' res$selected.tau
#' @export
#' 
fdtest <- function(X,alpha=0.05,tau=NULL,B=ceiling(50/alpha),pairs=NULL,transform=T,K=50,verbose=F,tau.method='MGB',basis='fourier')
{
    if(transform)
    {
        
        if(is.matrix(X)) X <- X %*% get.basis(X,K,ncol(X),basis) / ncol(X)
        else
        {
            M <- ncol(X[[1]])
            h <- 1/(2*M)
            pts <- seq(h,1-h,length.out=M)
            
            phi <- get.basis(X,K,M,basis)
            
            X <- lapply(X,function(z){
                z %*% phi / ncol(z)
            })
        }
    }
    
    return(hdtest(X,alpha=alpha,tau=tau,B=B,pairs=pairs,verbose=verbose))
}

get.basis <- function(X,K,M,basis)
{
    if(tolower(basis)=='fourier') fourier.basis(K,seq(0, 1, length.out = M))
    else stop(paste0(basis,': not supported'))
}


#' Create a Fourier basis of K functions in [0, 1]
#'
#' @param K A positive integer specifying the number of eigenfunctions to generate.
#' @param pts A vector specifying the time points to evaluate the basis functions.
#' @return A K by len(pts) matrix, each column containing a basis function.
#'
#' @examples
#' basis <- fourier.basis(3)
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

