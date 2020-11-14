context("check routines in anova.R")

test_that("check routines in anova.R", {
    
    p <- 10
    alpha <- 0.05
    
    # one-sample, size of test
    n <- 30
    mu <- rep(0,p) # null is true
    Sig <- diag((1:p)^(-0.5))
    set.seed(353)
    X <- MASS::mvrnorm(n,mu,Sig)
    res <- hdtest(X,alpha,tau.method='MGB')
    expect_true(!res$reject)
    expect_true(res$pvalue >= alpha)
    
    
    
    # one-sample, power of test
    n <- 50
    mu <- 0.27*(1:p)^(-2) # null is false
    Sig <- diag((1:p)^(-0.5))
    
    set.seed(4353)
    X <- MASS::mvrnorm(n,mu,Sig)
    res <- hdtest(X,alpha,tau.method='MGBA')
    expect_true(res$reject)
    expect_true(res$pvalue < alpha)

    
    # 3-sample, size of test
    set.seed(41343)
    G <- 3
    n <- rep(30,G)
    mu <- lapply(1:G,function(g) rep(0,p))
    Sig <- lapply(1:G,function(g) diag((1:p)^(-0.5*g)))
    X <- lapply(1:G,function(g) MASS::mvrnorm(n[g],mu[[g]],Sig[[g]]))
    res <- hdtest(X,alpha,tau.method='MGB',B=200)
    expect_true(!res$reject)
    expect_true(res$pvalue >= alpha)
    
    
    # 3-sample, power of test
    set.seed(41343)
    G <- 3
    n <- rep(30,G)
    mu <- lapply(1:G,function(g) 0.45*g*(1:p)^(-2))
    Sig <- lapply(1:G,function(g) diag((1:p)^(-0.5*g)))
    X <- lapply(1:G,function(g) MASS::mvrnorm(n[g],mu[[g]],Sig[[g]]))
    res <- hdtest(X,alpha,B=200,R=100)
    expect_true(res$reject)
    expect_true(res$pvalue < alpha)
    
    
    # 3-sample, functional anova
    set.seed(41343)
    G <- 3
    n <- rep(30,G)
    m <- 100
    mu <- lapply(1:G,function(g) 0.45*g*(1:m)^(-2))
    Sig <- lapply(1:G,function(g) diag((1:m)^(-0.5*g)))
    X <- lapply(1:G,function(g) MASS::mvrnorm(n[g],mu[[g]],Sig[[g]]))
    res <- fdtest(X,alpha,B=200,R=100)
    expect_true(res$reject)
    expect_true(res$pvalue < alpha)
    
    
    # 3-sample, multicore
    if(all(c('doParallel','parallel') %in% installed.packages()[,"Package"]))
    {
        ncore <- parallel::detectCores()
        if(ncore > 1)
        {
            set.seed(41343)
            G <- 3
            n <- rep(30,G)
            mu <- lapply(1:G,function(g) 0.45*g*(1:p)^(-2))
            Sig <- lapply(1:G,function(g) diag((1:p)^(-0.5*g)))
            X <- lapply(1:G,function(g) MASS::mvrnorm(n[g],mu[[g]],Sig[[g]]))
            res <- hdtest(X,alpha,B=200,R=100,ncore=2)
            expect_true(res$reject)
            expect_true(res$pvalue < alpha)
        }
    }
})