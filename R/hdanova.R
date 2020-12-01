#' Hypothesis Test for High-dimensional Data
#' @description Test the mean or differences of means of high-dimensional vectors are zero or not.
#' @param X a matrix (one-sample) or a list of matrices (multiple-samples), with each row representing an observation.
#' @param alpha significance level; default value: 0.05.
#' @param side either of \code{'<=','>='} or \code{'=='}; default value: \code{'=='}.
#' @param tau real number(s) in the interval \code{[0,1)} that specifies the decay parameter and is automatically selected if it is set to \code{NULL} or multiple values are provided; default value: \code{NULL}, which is equivalent to \code{tau=1/(1+exp(-0.8*seq(-6,5,by=1))).}
#' @param B the number of bootstrap replicates; default value: \code{ceiling(50/alpha)}.
#' @param pairs a matrix with two columns, only used when there are more than two populations, where each row specifies a pair of populations for which the SCI is constructed; default value: \code{NULL}, so that SCIs for all pairs are constructed.
#' @param Sig a matrix (one-sample) or a list of matrices (multiple-samples), each of which is the covariance matrix of a sample; default value: \code{NULL}, so that it is automatically estimated from data.
#' @param verbose TRUE/FALSE, indicator of whether to output diagnostic information or report progress; default value: FALSE.
#' @param tau.method the method to select tau; possible values are 'MGB' (default), 'MGBA', 'RMGB','RMGBA', 'WB' and 'WBA' (see \code{\link{hdsci}}).
#' @param ncore the number of CPU cores to be used; default value: 1.
#' @param cuda T/F to indicate whether to use CUDA GPU implementation when the package \code{hdanova.cuda} is installed. This option takes effect only when \code{ncore=1}.
#' @param R the number of Monte Carlo replicates for estimating the empirical size; default: \code{ceiling(25/alpha)}
#' @param nblock the number of block in CUDA computation
#' @param tpb number of threads per block; the maximum number of total number of parallel GPU threads is then \code{nblock*tpb}
#' @param seed the seed for random number generator
#' @param sci T/F, indicating whether to construct SCIs or not; default: FALSE.
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
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0.3*g,10),diag((1:10)^(-0.5*g))))
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' hdtest(X,alpha=0.05,pairs=matrix(1:4,2,2),tau=c(0.4,0.5,0.6))$reject
#' @export
hdtest <- function(X,alpha=0.05,side='==',tau=1/(1+exp(-0.8*seq(-6,5,by=1))),
                   B=ceiling(50/alpha),pairs=NULL,Sig=NULL,verbose=F,
                   tau.method='MGB',R=10*ceiling(1/alpha),ncore=1,cuda=T,
                   nblock=32,tpb=64,seed=sample.int(2^30,1),sci=F)
{
    if(is.list(X)) G <- length(X)
    else if(is.matrix(X)) G <- 1
    else stop('X must be a matrix or a list of matrices')
    
    sci.side <- switch (side,
            '>=' = 'upper',
            '<=' = 'lower',
            '==' = 'both',
            '='  = 'both',
            'both' = 'both'
    )
    
    if(ncore<=1 && cuda)
    {
        if('hdanova.cuda' %in% installed.packages()[,"Package"])
            return( eval(parse(text = 'hdanova.cuda::hdtest(X,alpha,side,tau,B,pairs,Sig,verbose,tau.method,R,nblock,tpb,seed,sci)')) )
        else
            message('Package hdanova.cuda is not detected. Automatically switch to the non-CUDA version.')
    }
    
    res <- hdsci(X,alpha,sci.side,tau,B,pairs,Sig,verbose,
                 tau.method,R,ncore,cuda,nblock,tpb,seed)
    
    res$pvalue <- pvalue(X,sci.side,res$sci$pairs,res$sci$sigma2,
                         res$sci$tau,res$sci$Mn,res$sci$Ln,B=B)
    
    res$reject <- (res$pvalue < alpha)
    res$accept <- !res$reject
    
    if(G > 2)
    {
        rej.idx <- sapply(res$sci$sci.lower,function(z) any(z>0)) |
            sapply(res$sci$sci.upper, function(z) any(z<0))
        if(res$reject)
        {
            rej.pairs <- res$sci$pairs[rej.idx,]
            res$rej.pairs <- rej.pairs
        }
        
        tmp <- cbind(res$sci$pairs,rej.idx)
        colnames(tmp) <- c('g1','g2','reject')
        res$pairs <- tmp
    }
    
    
    if(sci==FALSE){
        res <- within(res,rm(sci,sci.tau))
    } 
    
    return(res)
}
