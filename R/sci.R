#' Simultaneous confidence interval (SCI) for high-dimensional data
#' @description Construct (1-\code{alpha}) simultaneous confidence interval (SCI)  for the mean or difference of means of high-dimensional vectors
#' @param X a matrix (one sample) or a list of matrices (two samples or more), with observations contained in rows
#' @param alpha significance level
#' @param side either of \code{c('lower','upper', or 'both')}. default value is 'both'
#' @param tau the decay parameter, automatically selected if set to \code{NULL}
#' @param B the number of bootstrap replicates
#' @param pairs a matrix with two columns, used when there are more than two populations, each row specifying a pair of populations for which the SCI are constructed, and if set to \code{NULL}, SCIs for all pairs are constructed.
#' @param Sig a matrix (one sample) or a list of matrices, each of them is the covariance matrix of a sample and automatically estimated if \code{NULL}
#' @return a list of the following objects: \code{sci.lower} (\code{sci.upper}) is a vector (one- or two-sample) or a list of vectors (three or more samples) specifying the lower (upper) bound of SCI for the mean (one sample) or the difference of means of each pair of samples; when number of samples is larger than two, \code{pairs} is a matrix of two columns, each row containing the a pair of indices of samples whose SCI is constructed
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lopes2019+}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples  
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0,10),diag((1:10)^(-0.5*g))))
#' 
#' # construct SCIs for the mean vectors with pairs={(1,3),(2,4)}
#' hdsci(X,alpha=0.05,pairs=matrix(1:4,2,2))
#' @export
hdsci <- function(X,alpha=0.05,side='both',tau=NULL,B=1000,pairs=NULL,Sig=NULL,verbose=F)
{
    if(is.matrix(X)) # one-sample
    {
        if(length(tau) != 1)
        {
            res <- hdsci.tau(X,alpha,tau,B,
                             method=c('avg'),
                             NULL,Sig,verbose)
            tau <- res$tau
        }
        return(hdsci1(X,alpha,side,tau,B,Sig,verbose))
        #return(hdsci1(X,alpha,side,tau,B,Sig,verbose))
    }
    else if(is.list(X))
    {
        K <- length(X)
        # now 2 or more samples
        if(is.null(pairs)) pairs <- t(combn(1:K,2))
        if(length(tau) != 1)
        {
            res <- hdsci.tau(X,alpha,
                             tau,B,method=c('avg'),
                             pairs,Sig,verbose)
            tau <- res$tau
        }
        
        sci <- hdsciK(X,alpha,side,tau,B,pairs,Sig)
        sci$tau <- tau
        return(sci)
    }
    else stop('X must be matrix or a list')
    
    
}

# for one sample
hdsci1 <- function(X,alpha,side,tau,B,Sig,verbose)
{
    n <- nrow(X)
    p <- ncol(X)
    
    rtn <- sqrt(n)
    
    # sample Gaussian vectors
    
    if(p <= 1000)
    {
        if(is.null(Sig)) Sig <- var(scale(X,scale=F))
        sigm2 <- diag(Sig)
        sigma <- sqrt(diag(Sig))^tau
        W <- MASS::mvrnorm(B,rep(0,p),Sig)
    }
    else
    {
        sigma2 <- apply(X,2,var)
        sigma <- sqrt(sigma2)^tau
        xi <- matrix(rnorm(m*B),B,n)
        W <- (xi %*% X - apply(xi,1,sum) %*% t(apply(X,2,mean))) / sqrt(n)
    }
    
    # bootstrap max statistic
    idx <- (sigma <= 0)
    S <- W / matrix(sigma,B,p,byrow=T)
    S[,idx] <- 0
    Mn <- apply(S,1,max)
    Ln <- apply(S,1,min)
    
    
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
    
    Ln.sorted <- sort(Ln)
    Mn.sorted <- sort(Mn)
    
    X.bar <- apply(X,2,mean)
    
    sci.lower <- X.bar - Mn.sorted[b2] * sigma / rtn
    sci.lower[idx] <- 0
    sci.upper <- X.bar - Ln.sorted[b1] * sigma / rtn
    sci.upper[idx] <- 0
    
    
    # output
    res <- list(tau=tau,
                Ln.sorted=Ln.sorted,
                Mn.sorted=Mn.sorted,
                Sig=Sig,
                sigma.tau=sigma,
                sid=side,
                sigma2=sigma2)
    
    if(side == 'both' || side == 'lower')
    {
        res$sci.lower <- sci.lower
    }
    
    if(side == 'both' || size == 'upper')
    {
        res$sci.upper <- sci.upper
    }
    
    return(res)
}

# for more than one sample
hdsciK <- function(X,alpha,side,tau,B,pairs,Sig,verbose)
{
    
    ns <- sapply(X,function(x){nrow(x)}) # size of each sample
    p <- ncol(X[[1]]) # dimension
    n <- sum(ns) # total sample size
    G <- length(X) # number of samples
    
    
    # bootstrap max statistic
    bres <- bootstrap(X,B,pairs,tau,Sig)
    Mn.sorted <- bres$Mn.sorted
    Ln.sorted <- bres$Ln.sorted
    sigma2 <- bres$sigma2
    
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
        
        sigjk <- sqrt(lamj^2 * sig2j + lamk^2 * sig2k)^tau 
        
        X.bar <- apply(X[[j]],2,mean)
        Y.bar <- apply(X[[k]],2,mean)
        sqrt.harm.n <- sqrt(ns[j]*ns[k]/(ns[j]+ns[k]))
        
        idx <- (sigjk==0)
        
        if(side == 'both' || side == 'lower')
        {
            tmp <- (X.bar-Y.bar) - Mn.sorted[b2] * sigjk / sqrt.harm.n
            tmp[idx] <- 0
            sci.lower[[q]] <- tmp
        }
        if(side == 'both' || size == 'upper')
        {
            tmp <- (X.bar-Y.bar) - Ln.sorted[b1] * sigjk / sqrt.harm.n
            tmp[idx] <- 0
            sci.upper[[q]] <- tmp
        }
    }
    
    # output
    res <- list(tau=tau,
                Ln.sorted=Ln.sorted,
                Mn.sorted=Mn.sorted,
                Sig=Sig,
                sid=side,
                pairs=pairs,
                sigma2=sigma2)
    
    if(side == 'both' || side == 'lower')
    {
        if(G <= 2) sci.lower <- sci.lower[[1]]
        res$sci.lower <- sci.lower
    }
    
    if(side == 'both' || size == 'upper')
    {
        if(G <= 2) sci.upper <- sci.upper[[1]]
        res$sci.upper <- sci.upper
    }
    
    return(res)
}

# bootstrap

bootstrap <- function(X,B,pairs,tau,Sig)
{
    
    # size of each sample
    ns <- sapply(X,function(x){nrow(x)})
    
    # dimension
    p <- ncol(X[[1]])
    
    n <- sum(ns) # total sample size
    G <- length(X) # number of samples
    
    # bootstrap max statistic
    Mn <- matrix(0,B,nrow(pairs))
    Ln <- matrix(0,B,nrow(pairs))
    
    if(p <= 1000)
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
    
    
    for(q in 1:nrow(pairs))
    {
        j <- pairs[q,1]
        k <- pairs[q,2]
        
        lamj <- sqrt(ns[k]/(ns[j]+ns[k]))
        lamk <- sqrt(ns[j]/(ns[j]+ns[k]))
        
        sig2j <- sigma2[[j]] 
        sig2k <- sigma2[[k]]
        
        sigjk <- sqrt(lamj^2 * sig2j + lamk^2 * sig2k)^tau 
        S <- (lamj * W[[j]] - lamk * W[[k]]) / matrix(sigjk,B,p,byrow=T)
        idx <- (sigjk==0)
        S[,idx] <- 0
        Mn[,q] <- apply(S,1,max)
        Ln[,q] <- apply(S,1,min)
    }
    
    Mn.sorted <- sort(apply(Mn,1,max))
    Ln.sorted <- sort(apply(Ln,1,min))
    
    return(list(Mn.sorted=Mn.sorted,Ln.sorted=Ln.sorted,sigma2=sigma2))
}

# select tau 
hdsci.tau <- function(X,alpha,tau,B,
                      method=c('min.max','max.min','avg'),
                      pairs,Sig,verbose)
{
    if(is.null(tau)) tau <- seq(from=0.1,to=0.9,by=0.1)
    
    res <- list()
    
    S <- data.frame(matrix(NA,length(tau),length(method)))
    colnames(S) <- method
    
    for(i in 1:length(tau))
    {
        #3. bootstraping to get SCI
        if(is.matrix(X)) sci <- hdsci1(X,alpha,'both',tau[i],B,Sig,verbose) # one-sample
        else sci <- hdsciK(X,alpha,'both',tau[i],B,pairs,Sig,verbose) # two or more
        
        #4. select tau for each method
        u <- unlist(sci$ci.upper)
        l <- unlist(sci$ci.lower)
        if('min.max' %in% method){
            S[i,'min.max'] <- max(u-l)
        }
        if('max.min' %in% method){
            S[i,'max.min'] <- min(u-l)
        }
        if('avg' %in% method){
            S[i,'avg'] <- sqrt(sum((u-l)^2))
        }
    }
    
    
    if('min.max' %in% method){
        res$min.max <- tau[which.min(S[,'min.max'])]
        if(verbose) print('min.max select tau=%f\n',res$min.max)
        res$tau <- res$min.max
    }
    if('max.min' %in% method){
        res$max.min <- tau[which.max(S[,'max.min'])]
        if(verbose) print('max.min select tau=%f\n',res$max.min)
        res$tau <- res$max.min
    }
    if('avg' %in% method){
        k <- which.min(S[,'avg'])
        res$avg <- tau[k]
        if(verbose) print(sprintf('avg select tau=%f',res$avg))
        res$tau <- res$avg
    }
    return(res)
}
