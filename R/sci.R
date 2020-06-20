#' Simultaneous confidence interval (SCI) for high-dimensional data
#' @description Construct (1-\code{alpha}) simultaneous confidence interval (SCI)  for the mean or difference of means of high-dimensional vectors
#' @param X a matrix (one sample) or a list of matrices (two samples or more), with observations contained in rows
#' @param alpha significance level
#' @param side either of \code{c('lower','upper', or 'both')}. default value is 'both'
#' @param tau the decay parameter, automatically selected if set to \code{NULL} or multiple values are provided
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
#' hdsci(X,alpha=0.05,pairs=matrix(1:4,2,2))$sci
#' @export
hdsci <- function(X,alpha=0.05,side='both',tau=NULL,B=1000,pairs=NULL,Sig=NULL,verbose=F)
{
    if(is.matrix(X)) # one-sample
    {
        sci <- hdsci1(X,alpha,side,tau,B,Sig,verbose)
        return(sci)
    }
    else if(is.list(X))
    {
        K <- length(X)
        # now 2 or more samples
        if(is.null(pairs)) pairs <- t(combn(1:K,2))
        sci <- hdsciK(X,alpha,side,tau,B,pairs,Sig,verbose)
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
    
    if(is.null(tau)) tau <- seq(0,1,by=0.1)
    
    
    # bootstrap max statistic
    bres <- bootstrap(X,B,pairs,tau,Sig)
    Mn.sorted <- bres$Mn.sorted
    Ln.sorted <- bres$Ln.sorted
    sigma2 <- bres$sigma2
    
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
        
        list(sci.lower=sci.lower,
             sci.upper=sci.upper,
             sigma2=sigma2,
             tau=tau[v],
             side=side,
             Mn.sorted=Mn.sorted[[v]],
             Ln.sorted=Ln.sorted[[v]])
    })
    
    
    # output
    res <- list(tau=tau,
                Sig=Sig,
                sid=side,
                sigma2=sigma2,
                pairs=NULL,
                sci.tau=sci.tau)
    
    if(length(tau) > 1){
        selected.tau <- hdsci.tau(X,alpha,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B)
        v <- which(tau==selected.tau)
        res$sci <- sci.tau[[v]]
        res$selected.tau <- selected.tau
    } 
    else
    {
        res$sci <- sci.tau[[1]]
        res$selected.tau <- res$tau
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
    
    if(is.null(tau)) tau <- seq(0,1,by=0.1)
    
    
    # bootstrap max statistic
    bres <- bootstrap(X,B,pairs,tau,Sig)
    Mn.sorted <- bres$Mn.sorted
    Ln.sorted <- bres$Ln.sorted
    sigma2 <- bres$sigma2
    
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
            
            if(side == 'both' || side == 'lower')
            {
                tmp <- (X.bar-Y.bar) - Mn.sorted[[v]][b2] * sigjk / sqrt.harm.n
                tmp[idx] <- 0
                sci.lower[[q]] <- tmp
            }
            if(side == 'both' || size == 'upper')
            {
                tmp <- (X.bar-Y.bar) - Ln.sorted[[v]][b1] * sigjk / sqrt.harm.n
                tmp[idx] <- 0
                sci.upper[[q]] <- tmp
            }
        }
        
        if(side == 'both' || side == 'lower')
        {
            if(G <= 2) sci.lower <- sci.lower[[1]]
        }
        
        if(side == 'both' || side == 'upper')
        {
            if(G <= 2) sci.upper <- sci.upper[[1]]
        }
        
        list(sci.lower=sci.lower,
             sci.upper=sci.upper,
             sigma2=sigma2,
             tau=tau[v],
             side=side,
             pairs=pairs,
             Mn.sorted=Mn.sorted[[v]],
             Ln.sorted=Ln.sorted[[v]])
    })
    
    # output
    res <- list(tau=tau,
                Sig=Sig,
                sid=side,
                pairs=pairs,
                sigma2=sigma2,
                sci.tau=sci.tau)
    
    if(length(tau) > 1){
        selected.tau <- hdsci.tau(X,alpha,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B)
        v <- which(tau==selected.tau)
        res$sci <- sci.tau[[v]]
        res$selected.tau <- selected.tau
    } 
    else
    {
        res$sci <- sci.tau[[1]]
        res$selected.tau <- res$tau
    }
    
    return(res)
}

# bootstrap

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

pvalue <- function(X,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B=1000)
{
    if(is.null(Mn.sorted))
    {
        bs <- bootstrap(X,B,pairs,tau,Sig=NULL)
        sigma2 <- bs$sigma2
        Mn.sorted <- bs$Mn.sorted
        Ln.sorted <- bs$Ln.sorted
    }
    
    if(is.matrix(X))
    {
        zl <- Inf
        zu <- -Inf
        n <- nrow(X)
        
        pval <- sapply(1:length(tau), function(v){
            
            sig <- sqrt(sigma2)^tau[v] 
            X.bar <- apply(X,2,mean)
            sqrt.n <- sqrt(n)
            
            idx <- (sig!=0)
            sig <- sig[idx]
            tmp <- X.bar[idx]
            
            zl <- min(zl,min(tmp*sqrt.n/sig))
            zu <- max(zu,max(tmp*sqrt.n/sig))
            
            
            eta <- mean(Mn.sorted[[v]] >= zu)
            eta <- min(eta,mean(Ln.sorted[[v]] <= zl))
            
            pval.v <- min(1,2*eta)
            pval.v
        })
    }
    else
    {
        zl <- Inf
        zu <- -Inf
        ns <- sapply(X,function(x){nrow(x)})
        
        if(!is.list(Mn.sorted))
        {
            Mn.sorted <- list(Mn.sorted)
            Ln.sorted <- list(Ln.sorted)
        }
        
        pval <- sapply(1:length(tau), function(v){
            
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
            
            eta <- mean(Mn.sorted[[v]] >= zu)
            eta <- min(eta,mean(Ln.sorted[[v]] <= zl))
            
            pval.v <- min(1,2*eta)
            pval.v
        })
        
    }
    
    pval
}


# select tau
hdsci.tau <- function(X,alpha,pairs,sigma2,tau,Mn.sorted,Ln.sorted,B)
{
    pv <- pvalue(X,pairs,sigma2,tau,Mn.sorted,Ln.sorted)
    
    if(is.matrix(X))
    {
        X <- scale(X,scale=F)
        n <- nrow(X)
        test <- sapply(1:100, function(j)
        {
            Y <- X[sample(1:n,n,replace=T),]
            pvalue(Y,pairs,NULL,tau,NULL,NULL,B=100)
        })
    }
    else
    {
        X <- lapply(X,function(z) scale(z,scale=F))
        ns <- sapply(X,function(x){nrow(x)})
        N <- sum(ns)
        test <- sapply(1:100, function(j)
        {
            Y <- lapply(X, function(z) z[sample(1:nrow(z),nrow(z),replace=T),])
            pvalue(Y,pairs,NULL,tau,NULL,NULL,B=100)
        })
    }
    
    rej <- test <= alpha
    size <- apply(rej,1,mean)
    
    s <- NULL
    
    if(all(size >alpha))
        s <- mean(tau[which(size==min(size))])
    else
    {
        idx <- size <= alpha
        pv[!idx] <- Inf
        s <- mean(tau[which(pv==min(pv))])
    }
    
    tau[which.min(abs(tau-s))]
}

