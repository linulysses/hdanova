#' Produce summary of HD ANOVA
#' @param object the object returned from \code{hdtest}, \code{hdsci} or \code{fdtest}
#' @return An object with information on p-value and/or SCI 
#' @export
#' @examples 
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0.3*g,10),diag((1:10)^(-0.5*g))))
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' summary(hdtest(X,alpha=0.05,pairs=matrix(1:4,2,2),tau=c(0.4,0.5,0.6)))

summary.hdaov <- function (object) 
{
    if (!inherits(object, "hdaov")) 
        stop(gettextf("object must be of class %s", 
                      dQuote("hdaov")), domain = NA)
    
    x <- list()
    
    # print p value and rejected pairs
    if('reject' %in% names(object))
    {
        s <- sprintf('%f (p-value: %f, selected tau: %f)', object$alpha, object$pvalue, object$selected.tau)
        if(object$reject)
        {
            cat(paste0('The null hypothesis is rejected at level ', s, '\n\n'))
            cat('Pairs of samples with significant differences in mean:\n')
            if(!is.null(object$pairs))
            {
                rej.pairs <- object$rej.pairs
                if(!is.matrix(rej.pairs)) rej.pairs <- matrix(rej.pairs,1,2)
                colnames(rej.pairs) <- c('G1','G2')
                x$rej.pairs <- rej.pairs
                print(rej.pairs,row.names=F)
                cat('\n')
            }
        }
        else cat(paste0('The null hypothesis cannot be rejected at level ', s, '\n\n'))
        x$pvalue <- object$pvalue
    }
    
    # print SCI
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
            sci <- as.data.frame(cbind(unlist(object$sci$sci.lower),unlist(object$sci$sci.upper)))
            sci <- cbind(do.call(rbind,lapply(1:nrow(object$sci$pairs),function(k) matrix(object$sci$pairs[k,],p,2,byrow = T))),sci)
            sci <- cbind(rep(attr(object,'vnames'),nrow(object$sci$pairs)),sci)
            colnames(sci) <- c( 'variable', 'G1','G2', 'lower','upper')
        }
        
        rownames(sci) <- NULL
        x$sci <- sci
        
        z <- sci$lower > 0 | sci$upper < 0
        if(any(z)){
            sci <- cbind(sci,ifelse(z,'*',''))
            colnames(sci)[length(colnames(sci))] <- ' '
        } 
        
        cat('Simultaneous confidence intervals (all):\n')
        print(sci,row.names = F)
        cat('\n')
        
        if(any(z))
        {
            cat('Simultaneous confidence intervals (not containing 0):\n')
            print(sci[z,],row.names = FALSE)
        }
    }
    
    class(x) <- "summary.hdanova"
    invisible(x)
}