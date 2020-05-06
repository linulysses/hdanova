#' Hypothesis test for high-dimensional data
#' @description Test the mean or differences of means of high-dimensional vectors are zero or not
#' @param X a matrix (one sample) or a list of matrices (two samples or more), with observations contained in rows
#' @param alpha significance level
#' @param tau the decay parameter, automatically selected if set to \code{NULL}
#' @param B the number of bootstrap replicates
#' @param pairs a matrix with two columns, used when there are more than two populations, each row specifying a pair of populations whose means are compared, and if set to \code{NULL}, comparisons for all pairs are performed
#' @param Sig a matrix (one sample) or a list of matrices, each of them is the covariance matrix of a sample and automatically estimated if \code{NULL}
#' @return a list of objects: \code{reject} is a boolean variable indicating whether the null is rejected or not. If the null if reject, \code{rej.pairs} optionally gives the pairs of samples that lead to rejection. \code{pvalue} is also returned.
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lopes2019+}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0,10),diag((1:10)^(-0.5*g))))
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' hdtest(X,alpha=0.05,pairs=matrix(1:4,2,2))$reject
#' @export
hdtest <- function(X,alpha=0.05,tau=NULL,B=1000,pairs=NULL,Sig=NULL,verbose=F)
{
    if(is.list(X)) G <- length(X)
    else if(is.matrix(X)) G <- 1
    else stop('X must be a matrix or a list of matrices')
    
    sciobj <- hdsci(X,alpha,'both',tau,B,pairs,Sig,verbose)
    
    sci <- sciobj$sci
    
    lo <- sci$sci.lower
    up <- sci$sci.upper
    if(is.list(lo))
    {
        lo <- unlist(lo)
        up <- unlist(up)
    }
    if(any(lo > 0) || any(up < 0)) reject <- TRUE
    else reject <- FALSE
    
    res <- list(reject=reject,accept=!reject,tau=sci$tau,sci=sci,sciobj=sciobj)
    
    if(reject && G > 2)
    {
        rej.idx <- sapply(sci$sci.lower,function(z) any(z>0)) |
            sapply(sci$sci.upper, function(z) any(z<0))
        rej.pairs <- sci$pairs[which(rej.idx==1),]
        res$rej.pairs <- rej.pairs
    }
    
    # compute p-value
    
    res$pvalue <- pvalue(X,sci$pairs,sci$sigma2,sci$tau,sci$Mn.sorted,sci$Ln.sorted,B=1000)
    
    return(res)
}