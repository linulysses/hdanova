context("check routines in anova.R")

test_that("check routines in anova.R", {
    
    p <- 10
    alpha <- 0.05
    R <- 50
    
    # one-sample, size of test
    n <- 30
    mu <- rep(0,p) # null is true
    Sig <- diag((1:p)^(-0.5))
   
    
    res <- matrix(0,1,R)
    for(r in 1:R){
        X <- MASS::mvrnorm(n,mu,Sig)
        res[r] <- hdtest(X,alpha)$reject
    }   
    expect_true(mean(res) < 2*alpha)
    
    
    
    # one-sample, power of test
    n <- 50
    mu <- 0.5*(1:p)^(-2) # null is false
    Sig <- diag((1:p)^(-0.5))
    
    res <- matrix(0,1,R)
    for(r in 1:R){
        X <- MASS::mvrnorm(n,mu,Sig)
        res[r] <- hdtest(X,alpha)$reject
    }   
    expect_true(mean(res) > 0.1)
    
    
    
    # 3-sample, size of test
    G <- 3
    n <- rep(30,G)
    mu <- list()
    Sig <- list()
    for(g in 1:G)
    {
        mu[[g]] <- rep(0,p)
        Sig[[g]] <- diag((1:p)^(-0.5*g))
    }
    
    res <- matrix(0,1,R)
    for(r in 1:R){
        X <- list()
        for(g in 1:G)
            X[[g]] <- MASS::mvrnorm(n[g],mu[[g]],Sig[[g]])
        res[r] <- hdtest(X,alpha)$reject
    }   
    mean(res)
    expect_true(mean(res) < 2*alpha)
    
    
    # 3-sample, power of test
    G <- 3
    n <- rep(30,G)
    mu <- list()
    Sig <- list()
    for(g in 1:G)
    {
        mu[[g]] <- 0.5*g*(1:p)^(-2)
        Sig[[g]] <- diag((1:p)^(-0.5*g))
    }
    
    res <- matrix(0,1,R)
    for(r in 1:R){
        X <- list()
        for(g in 1:G)
            X[[g]] <- MASS::mvrnorm(n[g],mu[[g]],Sig[[g]])
        res[r] <- hdtest(X,alpha)$reject
    }   
    expect_true(mean(res) > 0.1)
})