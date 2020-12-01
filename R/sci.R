#' Construct Simultaneous Confidence Interval
#' @description Construct (1-\code{alpha}) simultaneous confidence interval (SCI)  for the mean or difference of means of high-dimensional vectors.
#' @param X a matrix (one-sample) or a list of matrices (multiple-samples), with each row representing an observation.
#' @param alpha significance level; default value: 0.05.
#' @param side either of \code{'lower','upper'} or \code{'both'}; default value: \code{'both'}.
#' @param tau real number(s) in the interval \code{[0,1)} that specifies the decay parameter and is automatically selected if it is set to \code{NULL} or multiple values are provided; default value: \code{NULL}, which is equivalent to \code{tau=1/(1+exp(-0.8*seq(-6,5,by=1))).}
#' @param B the number of bootstrap replicates; default value: \code{ceiling(50/alpha)}.
#' @param pairs a matrix with two columns, only used when there are more than two populations, where each row specifies a pair of populations for which the SCI is constructed; default value: \code{NULL}, so that SCIs for all pairs are constructed.
#' @param Sig a matrix (one-sample) or a list of matrices (multiple-samples), each of which is the covariance matrix of a sample; default value: \code{NULL}, so that it is automatically estimated from data.
#' @param verbose TRUE/FALSE, indicator of whether to output diagnostic information or report progress; default value: FALSE.
#' @param tau.method the method to select tau; possible values are 'MGB' (default), 'MGBA', 'RMGB', 'RMGBA', 'WB' and 'WBA' (see details).
#' @param ncore the number of CPU cores to be used; default value: 1.
#' @param cuda T/F to indicate whether to use CUDA GPU implementation when the package \code{hdanova.cuda} is installed. This option takes effect only when \code{ncore=1}.
#' @param R the number of Monte Carlo replicates for estimating the empirical size; default: \code{ceiling(25/alpha)}
#' @param nblock the number of block in CUDA computation
#' @param tpb number of threads per block; the maximum number of total number of parallel GPU threads is then \code{nblock*tpb}
#' @param seed the seed for random number generator
#' @return a list of the following objects: 
#'      \describe{
#'          \item{\code{sci}}{the constructed SCI, which is a list of the following objects:
#'              \describe{
#'                  \item{\code{sci.lower}}{a vector (when <= two samples) or a list of vectors (when >= 3 samples) specifying the lower bound of the SCI for the mean (one-sample) or the difference of means of each pair of samples.}
#'                  \item{\code{sci.upper}}{a vector (when <= two samples) or a list of vectors (when >= 3 samples) specifying the upper bound of the SCI.}
#'                  \item{\code{pairs}}{a matrix of two columns, each row containing the a pair of indices of samples of which the SCI of the difference in mean is constructed.}
#'                  \item{\code{tau}}{the decay parameter that is used to construct the SCI.}
#'                  \item{\code{Mn}}{the sorted (in increasing order) bootstrapped max statistic.}
#'                  \item{\code{Ln}}{the sorted (in increasing order) bootstrapped min statistic.}
#'                  \item{\code{side}}{the input \code{side}.}
#'                  \item{\code{alpha}}{the input \code{alpha}.}
#'              }
#'          }
#'          \item{\code{tau}}{a vector of candidate values of the decay parameter.}
#'          \item{\code{sci.tau}}{a list of \code{sci} objects corresponding to the candidate values in \code{tau}.}
#'          \item{\code{selected.tau}}{the selected value of the decay parameter from \code{tau}.}
#'          \item{\code{side}}{the input \code{side}.}
#'          \item{\code{alpha}}{the input \code{alpha}.}
#'          \item{\code{pairs}}{a matrix of two columns, each row containing the a pair of indices of samples of which the SCI of the difference in mean is constructed.}
#'          \item{\code{sigma2}}{a vector (for one sample) or a list (for multiple samples) of vectors containing variance for each coordinate.}
#'          }
#' @details Four methods to select the decay parameter \code{tau} are provided. Using the fact that a SCI is equivalent to a hypothesis test problem, all of them first identify a set of good candidates which give rise to test that respects the specified level \code{alpha}, and then select a candidate that minimizes the p-value. These methods differ in how to identify the good candidates.
#'     \describe{
#'         \item{\code{MGB}}{for this method, conditional on the data \code{X}, \code{R=10*ceiling(1/alpha)} i.i.d. zero-mean multivariate Gaussian samples (called MGB samples here) are drawn, where the covariance of each sample is equal to the sample covariance matrix \code{Sig} of the data \code{X}. For each candidate value in \code{tau}, 1) the empirical distribution of the corresponding max/min statistic is obtained by reusing the same bootstrapped sample, 2) the corresponding p-value is obtained, and 3) the size is estimated by applying the test to all MGB samples. The candidate values with the empirical size closest to \code{alpha} are considered as good candidates.}
#'         \item{\code{MGBA}}{an slightly more aggressive version of \code{MGB}, where the candidate values with the estimated empirical size no larger than \code{alpha} are considered good candidates.}    
#'         \item{\code{RMGB}}{this method is similar to \code{MGB}, except that for each MGB sample, the covariance matrix is the sample covariance matrix of a resampled (with replacement) data \code{X}.}
#'         \item{\code{RMGBA}}{an slightly more aggressive version of \code{RMGB}, where the candidate values with the estimated empirical size no larger than \code{alpha} are considered good candidates.}
#'         \item{\code{WB}}{for this method, conditional on \code{X}, \code{R=10*ceiling(1/alpha)} i.i.d. samples (called WB samples here) are drawn by resampling \code{X} with replacement. For each candidate value in \code{tau}, 1) the corresponding p-value is obtained, and 2) the size is estimated by applying the test to all WB samples without reusing the bootstrapped sample. The candidate values with the empirical size closest to \code{alpha} are considered as good candidates.}
#'         \item{\code{WBA}}{an slightly more aggressive version of \code{WB}, where the candidate values with the estimated empirical size no larger than \code{alpha} are considered good candidates.}
#'     }
#'     Among these methods, MGB and MGBA are recommended, since they are computationally more efficiently and often yield good performance. The MGBA might have slightly larger empirical size. The WB and WBA methods may be subject to outliers, in which case they become more conservative. The RMGB is computationally slightly slower than WB, but is less subject to outliers.
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lopes2020}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples  
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0,10),diag((1:10)^(-0.5*g))))
#' 
#' # construct SCIs for the mean vectors with pairs={(1,3),(2,4)}
#' hdsci(X,alpha=0.05,pairs=matrix(1:4,2,2))$sci
#' @export
hdsci <- function(X,alpha=0.05,side='both',tau=1/(1+exp(-0.8*seq(-6,5,by=1))),
                  B=ceiling(50/alpha),pairs=NULL,Sig=NULL,
                  verbose=F,tau.method='MGB',R=10*ceiling(1/alpha),ncore=1,cuda=T,
                  nblock=32,tpb=64,seed=sample.int(2^30,1))
{
    
    if(ncore<=1 && cuda)
    {
        if('hdanova.cuda' %in% installed.packages()[,"Package"])
            return( eval(parse(text = 'hdanova.cuda::hdsci(X,alpha,side,tau,B,pairs,Sig,verbose,tau.method,nblock,tpb,seed,R)')) )
        else message('Package hdanova.cuda is not detected. Automatically switch to the non-CUDA version.')
    }
    
    
    res <- hdanova(X,alpha,side,tau,B,pairs,Sig,verbose,ncore=1)
    
    if(length(tau) > 1){
            D <- hdsci.tau(X,alpha,side,res$pairs,res$sigma2,tau,res$Mn,res$Ln,B,tau.method,verbose,R,ncore)
            v <- which(tau==D$selected.tau)
            res$sci <- res$sci.tau[[v]]
            res$selected.tau <- D$selected.tau
    } 
    else
    {
        res$sci <- res$sci.tau[[1]]
        res$selected.tau <- res$tau
    }
    
    res <- within(res,rm(Mn,Ln))
 
    return(res)
    
}

hdanova <- function(X,alpha,side,tau,B,pairs,Sig,verbose,ncore=1)
{
    if(is.matrix(X)) # one-sample
    {
        pairs <- NULL
        K <- 1
        ns <- nrow(X)
        p <- ncol(X)
        n <- ns
    }
    else if(is.list(X)) # K-sample with K>1
    {
        K <- length(X)
        if(is.null(pairs)) pairs <- t(combn(1:K,2))
        p <- ncol(X[[1]])
        ns <- sapply(X,function(x){nrow(x)}) # size of each sample
        n <- sum(ns)
    }
    else stop('X must be matrix or a list')
    
    # bootstrapping
    bres <- bootstrap.mc(X,B,pairs,tau,Sig,ncore)
    Mn.sorted <- bres$Mn.sorted
    Ln.sorted <- bres$Ln.sorted
    sigma2 <- bres$sigma2
    
    if(K==1)
    {
        rtn <- sqrt(n)
        X.bar <- apply(X,2,mean)
        
        sci.tau <- lapply(1:length(tau),function(v){
            
            # construct SCI
            side <- tolower(side)
            if(side == 'both')
            {
                a1 <- alpha/2
                a2 <- 1 - alpha/2
            }
            else
            {
                a1 <- alpha
                a2 <- alpha
            }
            
            b1 <- max(1,round(a1*B))
            b2 <- round(a2*B)
            
            sigma <- sqrt(sigma2)^tau[v]
            
            idx <- sigma == 0
            
            sci.lower <- X.bar - Mn.sorted[[v]][b2] * sigma / rtn
            sci.lower[idx] <- 0
            sci.upper <- X.bar - Ln.sorted[[v]][b1] * sigma / rtn
            sci.upper[idx] <- 0
            
            if(side == 'upper')  sci.lower[] <- -Inf
            if(side == 'lower')  sci.upper[] <- Inf
            
            list(sci.lower=sci.lower,
                 sci.upper=sci.upper,
                 tau=tau[v],
                 side=side,
                 alpha=alpha,
                 Mn=Mn.sorted[[v]],
                 Ln=Ln.sorted[[v]])
        })
        
    }
    else
    {
        sci.tau <- lapply(1:length(tau),function(v){
            # construct SCI
            side <- tolower(side)
            if(side == 'both')
            {
                a1 <- alpha/2
                a2 <- 1 - alpha/2
            }
            else
            {
                a1 <- alpha
                a2 <- alpha
            }
            
            
            b1 <- max(1,round(a1*B))
            b2 <- round(a2*B)
            
            sci.lower <- list()
            sci.upper <- list()
            
            for(q in 1:nrow(pairs))
            {
                j <- pairs[q,1]
                k <- pairs[q,2]
                
                lamj <- sqrt(ns[k]/(ns[j]+ns[k]))
                lamk <- sqrt(ns[j]/(ns[j]+ns[k]))
                
                sig2j <- sigma2[[j]]
                sig2k <- sigma2[[k]]
                
                sigjk <- sqrt(lamj^2 * sig2j + lamk^2 * sig2k)^tau[v] 
                
                X.bar <- apply(X[[j]],2,mean)
                Y.bar <- apply(X[[k]],2,mean)
                sqrt.harm.n <- sqrt(ns[j]*ns[k]/(ns[j]+ns[k]))
                
                idx <- (sigjk==0)
                
                if(side %in% c('both','lower'))
                {
                    tmp <- (X.bar-Y.bar) - Mn.sorted[[v]][b2] * sigjk / sqrt.harm.n
                    tmp[idx] <- 0
                    sci.lower[[q]] <- tmp
                }
                if(side %in% c('both','upper'))
                {
                    tmp <- (X.bar-Y.bar) - Ln.sorted[[v]][b1] * sigjk / sqrt.harm.n
                    tmp[idx] <- 0
                    sci.upper[[q]] <- tmp
                }
                
                if(side == 'upper')  sci.lower[[q]] <- rep(-Inf,p)
                if(side == 'lower')  sci.upper[[q]] <- rep(Inf,p)
            }
            
            if(K <= 2) 
            {
                sci.lower <- sci.lower[[1]]
                sci.upper <- sci.upper[[1]]
            }
            
            list(sci.lower=sci.lower,
                 sci.upper=sci.upper,
                 tau=tau[v],
                 side=side,
                 alpha=alpha,
                 pairs=pairs,
                 Mn=Mn.sorted[[v]],
                 Ln=Ln.sorted[[v]])
        })
    }
    
    # output
    res <- list(tau=tau,
                side=side,
                alpha=alpha,
                pairs=pairs,
                sigma2=sigma2,
                sci.tau=sci.tau,
                Mn=Mn.sorted,
                Ln=Ln.sorted)
    
    return(res)
}

#' @import MASS
bs.mc.helper <- function(X,G,ns,p,Sig,sigma2,B,tau,pairs)
{
    # no sorting here
    Mn.tau <- matrix(0,B,length(tau))
    Ln.tau <- matrix(0,B,length(tau))
    
    if(is.matrix(X))
    {
        n <- ns
        rtn <- sqrt(n)
        
        # sample Gaussian vectors
        
        if(p <= 100)
        {
            W <- MASS::mvrnorm(B,rep(0,p),Sig)
        }
        else
        {
            xi <- matrix(rnorm(n*B),B,n)
            W <- (xi %*% X - apply(xi,1,sum) %*% t(apply(X,2,mean))) / sqrt(n)
        }
        
        
        for(v in 1:length(tau))
        {
            sigma <- sqrt(sigma2)^tau[v]
            idx <- (sigma <= 0)
            S <- W / matrix(sigma,B,p,byrow=T)
            S[,idx] <- 0
            Mn <- apply(S,1,max)
            Ln <- apply(S,1,min)
            
            Mn.tau[,v] <- Mn
            Ln.tau[,v] <- Ln
        }
    }
    else
    {

        n <- sum(ns) # total sample size
        
        # bootstrap max statistic
        
        if(p <= 100)
        {
            # sample a Gaussian vector for each sample
            W <- list()
            for(g in 1:G)
            {
                W[[g]] <- MASS::mvrnorm(B,rep(0,p),Sig[[g]])
            }
        }
        else
        {

            xi <- lapply(ns,function(m) matrix(rnorm(m*B),B,m))
            
            W <- lapply(1:G,function(k){
                (xi[[k]] %*% X[[k]] - apply(xi[[k]],1,sum) %*% t(apply(X[[k]],2,mean))) / sqrt(ns[k])
            })
        }
        
        
        for(v in 1:length(tau))
        {
            Mn <- matrix(0,B,nrow(pairs))
            Ln <- matrix(0,B,nrow(pairs))
            
            for(q in 1:nrow(pairs))
            {
                j <- pairs[q,1]
                k <- pairs[q,2]
                
                lamj <- sqrt(ns[k]/(ns[j]+ns[k]))
                lamk <- sqrt(ns[j]/(ns[j]+ns[k]))
                
                sig2j <- sigma2[[j]] 
                sig2k <- sigma2[[k]]
                
                sigjk <- sqrt(lamj^2 * sig2j + lamk^2 * sig2k)^tau[v]
                S <- (lamj * W[[j]] - lamk * W[[k]]) / matrix(sigjk,B,p,byrow=T)
                idx <- (sigjk==0)
                S[,idx] <- 0
                Mn[,q] <- apply(S,1,max)
                Ln[,q] <- apply(S,1,min)
            }
            
            Mn.tau[,v] <- apply(Mn,1,max)
            Ln.tau[,v] <- apply(Ln,1,min)
        }
    }
    
    list(Mn=Mn.tau,Ln=Ln.tau)
}

# bootstrap
#' @importFrom stats var
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterExport
#' @import foreach
bootstrap.mc <- function(X,B,pairs,tau,Sig,ncore)
{
    #if(B <= 1000) return(bootstrap(X,B,pairs,tau,Sig))
    
    ncore_ <- ncore
    ncore <- min(parallel::detectCores(),ncore_)
    
    if(ncore <= 1)
    {
        if(ncore_ > 1) warning('only one core is detected.')
        return(bootstrap(X,B,pairs,tau,Sig))
    } 
    
    
    
    BB <- rep(floor(B/ncore),ncore)
    BB[ncore] <- BB[ncore] + B %% ncore
    
    if(is.matrix(X)) # one-sample
    {
        ns <- nrow(X)
        p <- ncol(X)
        
        rtn <- sqrt(ns)
        if(is.null(Sig) && p <= 100){
            Sig <- var(scale(X,scale=F))
            sigma2 <- diag(Sig)
        } 
        else sigma2 <- sigma2 <- apply(X,2,var)
        
    }
    else # K-sample
    {
        # size of each sample
        ns <- sapply(X,function(x){nrow(x)})
        
        # dimension
        p <- ncol(X[[1]])
        
        G <- length(X) # number of samples
        
        # estimated covariance matrices of each sample if not provided
        if(is.null(Sig) && p <= 100)
        {
            Sig <- lapply(X,function(x){var(scale(x,center=TRUE,scale=FALSE))})
            sigma2 <- lapply(Sig,diag)
        }
        else sigma2 <- lapply(X,function(Z) apply(Z,2,var))
    }
    

    cl <- parallel::makeCluster(ncore, type="SOCK")  
    parallel::clusterExport(cl, c('bs.mc.helper'),envir = environment())
    
    doParallel::registerDoParallel(cl)  
    
    i <- NULL # a dirty trick to avoid a note: no visible binding for global variable 'i' Undefined global functions or variables: i
    result <- foreach(i=1:ncore) %dopar% {
        
        bs.mc.helper(X,G,ns,p,Sig,sigma2,BB[i],tau,pairs)
    }
    
    parallel::stopCluster(cl)
    
    # merge
    Mn.sorted <- list()
    Ln.sorted <- list()
    
    Mn <- do.call(rbind,lapply(result,'[[','Mn'))
    Ln <- do.call(rbind,lapply(result,'[[','Ln'))
    
    for(v in 1:length(tau))
    {
        Mn.sorted[[v]] <- sort(Mn[,v])
        Ln.sorted[[v]] <- sort(Ln[,v])
    }

    return(list(Mn.sorted=Mn.sorted,Ln.sorted=Ln.sorted,sigma2=sigma2,Sig=Sig))
}

# bootstrap
#' @importFrom stats rnorm var
bootstrap <- function(X,B,pairs,tau,Sig)
{
    
    if(is.matrix(X))
    {
        n <- nrow(X)
        p <- ncol(X)
        
        rtn <- sqrt(n)
        
        # sample Gaussian vectors
        
        if(p <= 100)
        {
            if(is.null(Sig)) Sig <- var(scale(X,scale=F))
            sigma2 <- diag(Sig)
            W <- MASS::mvrnorm(B,rep(0,p),Sig)
        }
        else
        {
            sigma2 <- apply(X,2,var)
            xi <- matrix(rnorm(n*B),B,n)
            W <- (xi %*% X - apply(xi,1,sum) %*% t(apply(X,2,mean))) / sqrt(n)
        }
        
        Mn.sorted <- list()
        Ln.sorted <- list()
        
        for(v in 1:length(tau))
        {
            sigma <- sqrt(sigma2)^tau[v]
            idx <- (sigma <= 0)
            S <- W / matrix(sigma,B,p,byrow=T)
            S[,idx] <- 0
            Mn <- apply(S,1,max)
            Ln <- apply(S,1,min)
            
            Mn.sorted[[v]] <- sort(Mn)
            Ln.sorted[[v]] <- sort(Ln)
        }
    }
    else
    {
        
        # size of each sample
        ns <- sapply(X,function(x){nrow(x)})
        
        # dimension
        p <- ncol(X[[1]])
        
        n <- sum(ns) # total sample size
        G <- length(X) # number of samples
        
        # bootstrap max statistic
        
        if(p <= 100)
        {
            # estimated covariance matrices of each sample if not provided
            if(is.null(Sig))
            {
                Sig <- lapply(X,function(x){var(scale(x,center=TRUE,scale=FALSE))})
            }
            
            # sample a Gaussian vector for each sample
            W <- list()
            for(g in 1:G)
            {
                W[[g]] <- MASS::mvrnorm(B,rep(0,p),Sig[[g]])
            }
            
            sigma2 <- lapply(Sig,diag)
            
        }
        else
        {
            sigma2 <- lapply(X,function(Z) apply(Z,2,var))
            
            xi <- lapply(ns,function(m) matrix(rnorm(m*B),B,m))
            
            W <- lapply(1:G,function(k){
                (xi[[k]] %*% X[[k]] - apply(xi[[k]],1,sum) %*% t(apply(X[[k]],2,mean))) / sqrt(ns[k])
            })
        }
        
        
        Mn.list <- list()
        Ln.list <- list()
        
        for(v in 1:length(tau))
        {
            Mn <- matrix(0,B,nrow(pairs))
            Ln <- matrix(0,B,nrow(pairs))
            
            for(q in 1:nrow(pairs))
            {
                j <- pairs[q,1]
                k <- pairs[q,2]
                
                lamj <- sqrt(ns[k]/(ns[j]+ns[k]))
                lamk <- sqrt(ns[j]/(ns[j]+ns[k]))
                
                sig2j <- sigma2[[j]] 
                sig2k <- sigma2[[k]]
                
                sigjk <- sqrt(lamj^2 * sig2j + lamk^2 * sig2k)^tau[v]
                S <- (lamj * W[[j]] - lamk * W[[k]]) / matrix(sigjk,B,p,byrow=T)
                idx <- (sigjk==0)
                S[,idx] <- 0
                Mn[,q] <- apply(S,1,max)
                Ln[,q] <- apply(S,1,min)
            }
            
            Mn.list[[v]] <- Mn
            Ln.list[[v]] <- Ln
        }
        
        Mn.sorted <- lapply(Mn.list, function(Mn) sort(apply(Mn,1,max)) )
        Ln.sorted <- lapply(Ln.list, function(Ln) sort(apply(Ln,1,min)) )
    }
    
    return(list(Mn.sorted=Mn.sorted,Ln.sorted=Ln.sorted,sigma2=sigma2,W=W))
}

# compute p-values for each candidate value of tau
pvalue <- function(X,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B=1000)
{
    if(is.null(Mn.sorted))
    {
        bs <- bootstrap(X,B,pairs,tau,Sig=NULL)
        sigma2 <- bs$sigma2
        Mn.sorted <- bs$Mn.sorted
        Ln.sorted <- bs$Ln.sorted
    }
    
    # only true if length(tau)==1
    if(!is.list(Mn.sorted)) Mn.sorted <- lapply(1:length(tau), function(v) Mn.sorted)
    if(!is.list(Ln.sorted)) Ln.sorted <- lapply(1:length(tau), function(v) Ln.sorted)
    
    # for one-sample
    if(is.matrix(X))
    {
        
        n <- nrow(X)
        
        pval <- sapply(1:length(tau), function(v){
            
            zl <- Inf
            zu <- -Inf
            
            sig <- sqrt(sigma2)^tau[v] 
            X.bar <- apply(X,2,mean)
            sqrt.n <- sqrt(n)
            
            idx <- (sig!=0)
            sig <- sig[idx]
            tmp <- X.bar[idx]
            
            zl <- min(zl,min(tmp*sqrt.n/sig))
            zu <- max(zu,max(tmp*sqrt.n/sig))
            
            eta.lower <- mean(Mn.sorted[[v]] >= zu)
            eta.upper <- mean(Ln.sorted[[v]] <= zl)
            eta.both <- min(1,2*min(eta.lower,eta.upper))
            
            switch(side,'both' = eta.both,
                   'lower' = eta.lower,
                   'upper' = eta.upper)
        })
    }
    else # for multiple-sample
    {

        ns <- sapply(X,function(x){nrow(x)})
        
        if(!is.list(Mn.sorted))
        {
            Mn.sorted <- list(Mn.sorted)
            Ln.sorted <- list(Ln.sorted)
        }
        
        pval <- sapply(1:length(tau), function(v){
            
            zl <- Inf
            zu <- -Inf
            
            for(q in 1:nrow(pairs))
            {
                j <- pairs[q,1]
                k <- pairs[q,2]
                
                lamj <- sqrt(ns[k]/(ns[j]+ns[k]))
                lamk <- sqrt(ns[j]/(ns[j]+ns[k]))
                
                sig2j <- sigma2[[j]]
                sig2k <- sigma2[[k]]
                
                sigjk <- sqrt(lamj^2 * sig2j + lamk^2 * sig2k)^tau[v] 
                
                X.bar <- apply(X[[j]],2,mean)
                Y.bar <- apply(X[[k]],2,mean)
                
                sqrt.harm.n <- sqrt(ns[j]*ns[k]/(ns[j]+ns[k]))
                
                idx <- (sigjk!=0)
                
                sigjk <- sigjk[idx]
                
                tmp <- X.bar-Y.bar
                tmp <- tmp[idx]
                
                zl <- min(zl,min(tmp*sqrt.harm.n/sigjk))
                zu <- max(zu,max(tmp*sqrt.harm.n/sigjk))
            }
            
            eta.lower <- mean(Mn.sorted[[v]] >= zu)
            eta.upper <- mean(Ln.sorted[[v]] <= zl)
            eta.both <- min(1,2*min(eta.lower,eta.upper))
            
            switch(side,'both' = eta.both,
                   'lower' = eta.lower,
                   'upper' = eta.upper)
        })
        
    }
    
    pval
}

# generate a zero-mean gaussian sample of size n with X's covariance matrix
mgauss <- function(X,n,Sig=NULL)
{
    p <- ncol(X)
    
    if(p <= 100)
    {
        if(is.null(Sig)) Sig <- var(scale(X,scale=F))
        W <- MASS::mvrnorm(n,rep(0,p),Sig)
    }
    else
    {
        xi <- matrix(rnorm(n*n),n,n)
        W <- (xi %*% X - apply(xi,1,sum) %*% t(apply(X,2,mean))) / sqrt(n)
    }
    W
}

hdsci.tau <- function(X,alpha,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method,verbose,R,ncore)
{
    D <- inspect.tau(X,tau,alpha=alpha,side=side,pairs=pairs,sigma2=sigma2,
                     Mn.sorted=Mn.sorted,Ln.sorted=Ln.sorted,
                     B=B,R=R,method=method,ncore=ncore)

    D$selected.tau <- choose.tau(alpha,tau,D$size,D$pval,ifelse(method %in% c('MGB','WB','RMGB'),yes='C',no='A'))
    
    if(verbose)
    {
        message(paste0('candidate values: ', paste0(tau,collapse=' ')))
        message(paste0('estimated sizes : ', paste0(D$size,collapse=' ')))
        message(paste0('p values        : ', paste0(D$pval,collapse=' ')))
        message(paste0('selected  tau   : ', D$selected.tau, ', using method=',method,', based on R=', R,' replicates.'))
    }
    
    return(D)
}

#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterExport
#' @import foreach
inspect.tau <- function(X,tau,alpha,side,pairs,sigma2,Mn.sorted,Ln.sorted,B,R,method,ncore)
{
    if(is.null(Mn.sorted) || is.null(Ln.sorted))
    {
        bs <- bootstrap(X,B,pairs,tau,Sig=NULL)
        sigma2 <- bs$sigma2
        Mn.sorted <- bs$Mn.sorted
        Ln.sorted <- bs$Ln.sorted
    }
    
    if(is.null(sigma2)){
        if(is.matrix(X)) sigma2 <- apply(X,2,var)
        else sigma2 <- lapply(X, function(x) apply(x,2,var))
    } 
    

    pval <- pvalue(X,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted)
    
    if(method %in% c('RMGB','RMGBA','WB','WBA'))
    {
        size <- size.tau(X,tau,alpha,side,B,pairs,FALSE,R=R,method=method,ncore=ncore)
    }
    else # in order to resuse Mn.sorted, etc, we do not use size.tau function
    {
        if(is.matrix(X))
        {
            X <- scale(X,scale=F)
            ns <- nrow(X)
        } 
        else
        {
            X <- lapply(X,function(z) scale(z,scale=F))
            ns <- sapply(X,function(x){nrow(x)})
        } 
        
        ncore_ <- ncore
        ncore <- min(parallel::detectCores(),ncore)
        
        if(ncore_ > 1 && ncore <= 1) warning('only one core is detected.')
        
        if(ncore==1)
        {
            test.result <- sapply(1:R, function(j)
            {
                if(is.matrix(X)) Y <- mgauss(X,ns,NULL)
                else Y <- lapply(1:length(ns), function(g) mgauss(X[[g]],ns[g],NULL))
                
                pvalue(Y,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B=B)
            })
        }
        else
        {

            cl <- parallel::makeCluster(ncore, type="SOCK")  
            parallel::clusterExport(cl, c('mgauss','pvalue'),envir = environment())
            
            doParallel::registerDoParallel(cl)  
            
            result <- foreach(i=1:R) %dopar% {
                
                        if(is.matrix(X)) Y <- mgauss(X,ns,NULL)
                        else Y <- lapply(1:length(ns), function(g) mgauss(X[[g]],ns[g],NULL))
                        
                        pvalue(Y,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B=B)
                        
                    }
            
            parallel::stopCluster(cl)
            test.result <- do.call(cbind,result)
        }

        rej <- test.result <= alpha
        size <- apply(rej,1,mean)
    }

    list(size=size,pval=pval,tau=tau)
}


# choose tau based on the estimated size and calculuated p-values
choose.tau <- function(alpha,tau,size,pv,mod='C',margin=0.01)
{
    if(mod == 'C')
    {
        s <- abs(size-alpha)
        smin <- min(s)
        
        # those with empirical size closest to alpha are "good" candidate tau values
        # allow a margin to account for computer finite-precision
        idx <- which( abs(s-smin) < margin*alpha) 
    }
    else # aggressive version, considering all candidate values with empirical size lower than alpha whenever possible
    {
        # no empirical sizes lower than alpha
        if(all(size > (1+margin)*alpha) == 0)
        {
            idx <- which(size==min(size))
        }
        else # those with empirical size lower than alpha are "good" candidate values
        {
            idx <- which(size <= (1+margin)*alpha)
        }
    }
    
    # consider only good candidate values
    pv[-idx] <- Inf 
    
    # now select tau to minimize p-value
    # when there are multiple such tau values, take the mean of them
    tau0 <- mean(tau[which(pv==min(pv))]) 
    
    # If tau0 is the mean of two or more candidate values,
    # then tau0 might not be in the prescribed list
    # In this case, we select a tau in the list that is closest to tau0
    tau[which.min(abs(tau-tau0))] 
}


#' Empirical Size Associated with Decay Parameter
#' @description Find the empirical size associated with the decay parameter conditional on a dataset
#' @param X a matrix (one-sample) or a list of matrices (multiple-samples), with each row representing an observation.
#' @param alpha significance level; default value: 0.05.
#' @param side either of \code{'lower','upper'} or \code{'both'}; default value: \code{'both'}.
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
size.tau <- function(X,tau,alpha=0.05,side='both',
                     B=ceiling(50/alpha),pairs=NULL,verbose=F,
                     method='MGB',R=ceiling(10/alpha),ncore=1,cuda=T,
                     nblock=32,tpb=64,seed=sample.int(2^30,1))
{
    
    if(ncore<=1 && cuda)
    {
        if('hdanova.cuda' %in% installed.packages()[,"Package"])
            return( eval(parse(text = 'hdanova.cuda::size.tau(X,tau,alpha,side,B,pairs,verbose,method,R,nblock,tpb,seed)')) )
        else
            message('Package hdanova.cuda is not detected. Automatically switch to the non-CUDA version.')
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
    
    test <- function(X,alpha,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method)
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
        
        
        pv <- pvalue(Y,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B=B)
        pv < alpha
    }
    
    ncore_ <- ncore
    ncore <- min(parallel::detectCores(),ncore_)
    
    if(ncore_ > 1 && ncore <= 1) warning('only one core is detected.')
    
    
    if(ncore==1)
        test.result <- sapply(1:R,function(j){
            test(X,alpha,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method)
        })
    else
    {
        
        cl <- parallel::makeCluster(ncore, type="SOCK")  
        parallel::clusterExport(cl, c('mgauss','pvalue'))
        
        doParallel::registerDoParallel(cl)  
        
        result <- foreach(i=1:R) %dopar% {
            test(X,alpha,side,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B,method)
        }
        parallel::stopCluster(cl)
        test.result <- do.call(cbind,result)
    }
    apply(test.result,1,mean)
}
