#' Hypothesis Test for High-dimensional Data
#' @description Test the mean or differences of means of high-dimensional vectors are zero or not.
#' @param X a matrix (one-sample) or a list of matrices (multiple-samples), with each row representing an observation.
#' @param alpha significance level; default value: 0.05.
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
hdtest <- function(X,alpha=0.05,tau=1/(1+exp(-0.8*seq(-6,5,by=1))),
                   B=ceiling(50/alpha),pairs=NULL,Sig=NULL,verbose=F,
                   tau.method='MGB',R=10*ceiling(1/alpha),ncore=1,cuda=T,
                   nblock=32,tpb=64,seed=sample.int(2^30,1),sci=F)
{
    if(is.list(X)) G <- length(X)
    else if(is.matrix(X)) G <- 1
    else stop('X must be a matrix or a list of matrices')
    
    if(ncore<=1 && cuda && 'hdanova.cuda' %in% installed.packages()[,"Package"])
    {
        return( hdanova.cuda::hdtest(X,alpha,tau,B,pairs,Sig,verbose,tau.method,R,nblock,tpb,seed,sci) )
    }
    
    res <- hdsci(X,alpha,'both',tau,B,pairs,Sig,verbose,
                 tau.method,R,ncore,cuda,nblock,tpb,seed)
    
    sci <- res$sci
    
    lo <- sci$sci.lower
    up <- sci$sci.upper
    if(is.list(lo))
    {
        lo <- unlist(lo)
        up <- unlist(up)
    }
    if(any(lo > 0) || any(up < 0)) reject <- TRUE
    else reject <- FALSE
    
    res$reject <- reject
    res$accept <- !reject
    
    if(G > 2)
    {
        rej.idx <- sapply(sci$sci.lower,function(z) any(z>0)) |
            sapply(sci$sci.upper, function(z) any(z<0))
        if(reject)
        {
            rej.pairs <- sci$pairs[rej.idx,]
            res$rej.pairs <- rej.pairs
        }
        
        tmp <- cbind(sci$pairs,rej.idx)
        colnames(tmp) <- c('g1','g2','reject')
        res$pairs <- tmp
    }
    
    # compute p-value
    
    res$pvalue <- pvalue(X,sci$pairs,sci$sigma2,sci$tau,sci$Mn,sci$Ln,B=B)
    
    return(res)
}

#' Empirical Size Associated with Decay Parameter
#' @description Find the empirical size associated with the decay parameter conditional on a dataset
#' @param X a matrix (one-sample) or a list of matrices (multiple-samples), with each row representing an observation.
#' @param alpha significance level; default value: 0.05.
#' @param tau real number(s) in the interval \code{[0,1)} for which the empirical size will be evaluated.
#' @param B the number of bootstrap replicates; default value: \code{ceiling(50/alpha)}.
#' @param pairs a matrix with two columns, only used when there are more than two populations, where each row specifies a pair of populations for which the SCI is constructed; default value: \code{NULL}, so that SCIs for all pairs are constructed.
#' @param verbose TRUE/FALSE, indicator of whether to output diagnostic information or report progress; default value: FALSE.
#' @param R the number of iterations; default value: \code{ceiling(25/alpha)}.
#' @param method the evaluation method tau; possible values are 'MGB' (default), 'MGBA', 'RMGB', 'RMGBA', 'WB' and 'WBA' (see \code{\link{hdsci}} for details).
#' @param ncore the number of CPU cores to be used; default value: 1.
#' @param cuda T/F to indicate whether to use CUDA GPU implementation when the package \code{hdanova.cuda} is installed. This option takes effect only when \code{ncore=1}.
#' @param nblock the number of block in CUDA computation
#' @param tpb number of threads per block; the maximum number of total number of parallel GPU threads is then \code{nblock*tpb}
#' @param seed the seed for random number generator
#' @return a vector of empirical size corresponding to \code{tau}.
#' @importFrom Rdpack reprompt
#' @importFrom utils combn installed.packages
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterExport
#' @import foreach
#' @references 
#' \insertRef{Lopes2020}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0.3*g,10),diag((1:10)^(-0.5*g))))
#' 
#' size.tau(X,tau=seq(0,1,by=0.1),alpha=0.05,pairs=matrix(1:4,2,2),R=100)
#' @export
size.tau <- function(X,tau,alpha=0.05,B=ceiling(50/alpha),pairs=NULL,verbose=F,
                     method='MGB',R=ceiling(10/alpha),ncore=1,cuda=T,
                     nblock=32,tpb=64,seed=sample.int(2^30,1))
{
    if(ncore<=1 && cuda && 'hdanova.cuda' %in% installed.packages()[,"Package"])
    {
        return( hdanova.cuda::size.tau(X,tau,alpha,B,pairs,verbose,method,R,nblock,tpb,seed) )
    }
    
    if(is.matrix(X)) X <- scale(X,scale=F)
    else
    {
        X <- lapply(X,function(x) scale(x,scale=F))
        if(is.null(pairs)) pairs <- t(combn(1:length(X),2))
    }
    
    if(method %in% c('MGB','MGBA'))
    {
        bs <- bootstrap(X,B,pairs,tau,Sig=NULL)
        sigma2 <- bs$sigma2
        Mn.sorted <- bs$Mn.sorted
        Ln.sorted <- bs$Ln.sorted

        
        if(is.matrix(X)) sigma2 <- apply(X,2,var)
        else sigma2 <- lapply(X, function(x) apply(x,2,var))
    }
    else
    {
        sigma2 <- NULL
        Mn.sorted <- NULL
        Ln.sorted <- NULL
    }
    
    test <- function(X,alpha,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method)
    {
        if(method %in% c('RMGB','RMGBA'))
        {
            if(is.matrix(X)) Y <- mgauss(X[sample.int(nrow(X),replace=T),],nrow(X),NULL)
            else Y <- lapply(X,function(x) mgauss(x[sample.int(nrow(x),replace=T),],nrow(x),NULL))
        }
        else if(method %in% c('WB','WBA'))
        {
            if(is.matrix(X)) Y <- X[sample.int(nrow(X),replace=T),]
            else Y <- lapply(X,function(x) x[sample.int(nrow(x),replace=T),])
        }
        else # MGB and MGBA
        {
            if(is.matrix(X)) Y <- mgauss(X,nrow(X),NULL)
            else Y <- lapply(X,function(x) mgauss(x,nrow(x),NULL))
        }
        
        
        pv <- pvalue(Y,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B=B)
        pv < alpha
    }
    
    if(ncore==1)
        test.result <- sapply(1:R,function(j){
            test(X,alpha,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method)
        })
    else
    {
     
        ncore <- min(parallel::detectCores(),ncore)
        
        
        cl <- parallel::makeCluster(ncore, type="SOCK")  
        parallel::clusterExport(cl, c('mgauss','pvalue'))
        
        doParallel::registerDoParallel(cl)  
        
        result <- foreach(i=1:R,
                          .packages=c('R.utils')#, .options.snow=opts
                          ) %dopar% {
                              test(X,alpha,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method)
                              
                          }
        parallel::stopCluster(cl)
        test.result <- do.call(cbind,result)
    }
    apply(test.result,1,mean)
}
    