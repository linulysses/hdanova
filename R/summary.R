#' Produce summary of HD ANOVA
#' @param object the object returned from \code{hdtest}, \code{hdsci} or \code{fdtest}
#' @param all.sci TRUE/FALSE, indicating whether to display the whole SCIs.
#' @param ... other arguments passed to summary or print.
#' @return An object with information on p-value and/or SCI 
#' @export
#' @examples 
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0.3*g,10),diag((1:10)^(-0.5*g))))
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' summary(hdtest(X,alpha=0.05,pairs=matrix(1:4,2,2),tau=c(0.4,0.5,0.6)))

summary.hdaov <- function(object, ..., all.sci = FALSE) 
{
    if(!inherits(object, "hdaov")) 
        stop(gettextf("object must be of class %s", 
                      dQuote("hdaov")), domain = NA)
    
    output <- list()
    
    # p value and rejected pairs
    if('reject' %in% names(object))
    {
        output$alpha <- object$alpha
        output$pvalue <- object$pvalue
        output$selected.tau <- object$selected.tau
        if(object$reject)
        {
            
            if(!is.null(object$pairs))
            {
                rej.pairs <- object$rej.pairs
                if(!is.matrix(rej.pairs)) rej.pairs <- matrix(rej.pairs,1,2)
                colnames(rej.pairs) <- c('G1','G2')
                output$rej.pairs <- rej.pairs
            }
        }
    }
    
    #SCI
    if('sci' %in% names(object))
    {
        if(is.null(object$pairs))
        {
            sci <- as.data.frame(cbind(object$sci$sci.lower,object$sci$sci.upper))
            sci <- cbind(attr(object,'vnames'),sci)
            colnames(sci) <- c(' ', 'lower','upper')
        }
        else
        {
            p <- length(object$sci$sci.lower[[1]])
            if(is.list(object$sci$sci.lower)) pair.no <- rep(c(1:length(object$sci$sci.lower)), c(rep(p,length(object$sci$sci.lower))))
            else pair.no <- rep(1, p)
            sci <- as.data.frame(cbind(unlist(object$sci$sci.lower),unlist(object$sci$sci.upper)))
            sci <- cbind(do.call(rbind,lapply(1:nrow(object$sci$pairs),function(k) matrix(object$sci$pairs[k,],p,2,byrow = T))),sci)
            sci <- cbind(pair.no,rep(attr(object,'vnames'),nrow(object$sci$pairs)),sci)
            colnames(sci) <- c( 'pair no.', 'variable', 'G1','G2', 'lower','upper')
        }
        
        rownames(sci) <- NULL
        output$sci <- sci
        
        z <- sci$lower > 0 | sci$upper < 0
        if(any(z)){
            sci <- cbind(sci,ifelse(z,'*',''))
            colnames(sci)[length(colnames(sci))] <- ' '
            output$sci <- sci
            output$sci.rej <- sci[z,]
        } 
        
        if(all.sci == FALSE)
        {
            output$sci <- NULL
        }
        
    }
    
    class(output) <- "summary.hdaov"
    output
}

#' @export
print.summary.hdaov <- function (x, ...) 
{
    # print p value and rejected pairs
    
    s <- sprintf('%f \n (p-value: %f, selected tau: %f)', x$alpha, x$pvalue, x$selected.tau)
    if(isTRUE(x$alpha>=x$pvalue))
    {
        cat(paste0('The null hypothesis is rejected at level ', s, '\n\n'))
        if('rej.pairs' %in% names(x))
        {
            cat('Pairs of samples with significant differences in mean:\n')
            print(x$rej.pairs, row.names = F)
            cat('\n')
        }
    }
    else cat(paste0('The null hypothesis cannot be rejected at level ', s, '\n\n'))
    
    
    
    # print SCI
    if('sci' %in% names(x))
    {
        cat('Simultaneous confidence intervals (all):\n')
        print(x$sci, row.names = F)
        cat('\n')
    }
    if('sci.rej' %in% names(x))
    {
        cat('Simultaneous confidence intervals (not containing 0):\n')
        print(x$sci.rej, row.names = F)
        cat('\n')
        
    }
    invisible(x)
}