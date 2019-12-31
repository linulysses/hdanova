#' Hypothesis test for high-dimensional data
#' @description Test the mean or differences of means of high-dimensional vectors are zero or not
#' @param X a matrix (one sample) or a list of matrices (two samples or more), with observations contained in rows
#' @param alpha significance level
#' @param tau the decay parameter
#' @param B the number of bootstrap replicates
#' @param pairs a vector with two columns. When number of samples is larger than two, each row specifies the pairs of samples for which the SCI are constructed
#' @param Sig a matrix (one sample) or a list of matrices, each of them is the covariance matrix of a sample
#' @return a list of objects: \code{reject} is a boolean variable indicating whether the null is rejected or not. If the null if reject, \code{rej.pairs} optionally gives the pairs of samples that lead to rejection
#' @export
hdtest <- function(X,alpha=0.05,tau=NULL,B=1000,pairs=NULL,Sigma=NULL,verbose=F)
{
    if(is.list(X)) G <- length(X)
    else if(is.matrix(X)) G <- 1
    else stop('X must be a matrix or a list of matrices')
    
    sci <- hdsci(X,alpha,'both',tau,B,pairs,Sigma,verbose)
    
    lo <- sci$sci.lower
    up <- sci$sci.upper
    if(is.list(lo))
    {
        lo <- unlist(lo)
        up <- unlist(up)
    }
    if(any(lo > 0) || any(up < 0)) reject <- TRUE
    else reject <- FALSE
    
    res <- list(reject=reject,accept=!reject)
    
    if(reject && G > 2)
    {
        rej.idx <- sapply(sci$sci.lower,function(z) any(z>0)) |
            sapply(sci$sci.upper, function(z) any(z<0))
        rej.pairs <- sci$pairs[which(rej.idx==1),]
        res$rej.pairs <- rej.pairs
    }
    
    return(res)
}