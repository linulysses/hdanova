#' SCI Plot for pairs that lead to rejection
#' @description Plot the confidence interval of there pairs that lead to rejection after high-dimensional data analysis.
#' @param x the object returned from \code{hdtest}, \code{hdsci} or \code{fdtest}
#' @param pair.no a vector of integers where the integers denote the orders in whole pairs in analysis; only these confidence intervals of pairs specified in \code{pair.no} or \code{pairs} will be plotted; default value: \code{NULL}, which is selected from \code{pairs}.
#' @param pairs a matrix with two columns where each row specifies a pair of populations in analysis; only these confidence intervals of pairs specified in \code{pair.no} or \code{pairs} will be plotted; default value: \code{NULL}, so that the specified pairs for plotting SCI are all pairs that lead to rejection.
#' @param sci.coordinate a vector of integers where each integer specifies the coordinate of SCI; only these specified coordinates of SCI will be plotted; default value: \code{NULL}, which is equivalent to all coordinates.
#' @param only.plot.reject TRUE/FALSE, indicating whether only plot these SCIs uncontaing 0; default value: \code{TRUE}, so that the plot only contains these SCIs which doesn't contain 0.
#' @param title the text for the title for plot; the title setting is the same as that in labels of \pkg{gglot2}; default value: \code{'Simultaneous Confidence Interval Plot'}.
#' @param subtitle the text for the subtitle for plot; the subtitle setting is the same as that in labels of \pkg{gglot2}; default value: \code{NULL} which indicates the pair number and group numbers.
#' @param xlabel the text for X coordinate for plot; the label setting is the same as that in labels of \pkg{gglot2}; default value: \code{'Coordinates'}.
#' @param ylabel the text for Y coordinate for plot; the label setting is the same as that in labels of \pkg{gglot2}; default value: \code{'Confidence Interval'}.
#' @param ... other arguments passed to plot.
#' @return plots of the SCIs of specified pairs and coordinates
#' @import ggplot2
#' @references 
#' \insertRef{Lopes2020}{hdanova}
#' 
#' \insertRef{Lin2020}{hdanova}
#' @examples
#' # simulate a dataset of 4 samples
#' X <- lapply(1:4, function(g) MASS::mvrnorm(30,rep(0.3*g,10),diag((1:10)^(-0.5*g))))
#' 
#' # test for the equality of mean vectors with pairs={(1,3),(2,4)}
#' plot(hdtest(X,alpha=0.05,pairs=matrix(1:4,2,2),tau=c(0.4,0.5,0.6), return.sci=TRUE))
#' @export
plot.hdaov <- function(x, pair.no = NULL, pairs = NULL, sci.coordinate = NULL, only.plot.reject = TRUE,
                       title = 'Simultaneous Confidence Interval Plot', subtitle = NULL, 
                       xlabel = 'Coordinates', ylabel = 'Confidence Interval', ...)
{
    if (!inherits(x, "hdaov")) 
        stop(gettextf("object must be of class %s", 
                      dQuote("hdaov")), domain = NA)
    
    if(!('sci' %in% names(x))) stop('\nObject must contain SCI.\nPlease run \'hdtest\' with \'return.sci = T\' again.')
    
    if(!('reject' %in% names(x))) cat('no rejection in object.\n') 
    
    if(isTRUE(ncol(pairs) != 2 && !is.null(pairs))) stop('wrong format in pairs')
    
    pair.no <- unique(pair.no)
    sci.coordinate <- unique(sci.coordinate)
    if(!is.null(pairs)) pairs <- pairs[!duplicated(pairs),]
    if(!is.null(pairs) && !is.matrix(pairs)) pairs <- matrix(pairs,1,2)
    

    original.pairs <- x$pairs[,c(1,2)]
    if(!is.matrix(original.pairs) && !is.null(original.pairs)) original.pairs <- matrix(original.pairs,nrow(x$pairs),2)
    
    #deal with pair number and pairs 
    if(is.null(pair.no) && is.null(pairs))
    {
        if(!is.null(x$pairs)) pair.no <- which(x$pairs[,'reject']==1)

    } else if(is.null(pair.no) && !is.null(pairs))
    {
        for (i in 1:nrow(pairs)) 
        {
            if(is.null(original.pairs))
            {
                warning('No need to set \'pairs\'. \nThe procedure does not involve the \'pairs\' setting in this case.')
                pairs = NULL
            }else {
                pair.location <- which(apply(original.pairs, 1, function(y) sum(y == pairs[i,]))==2) #check whether input pairs are involved in original pairs
                if(length(pair.location) == 0) stop('input wrong pairs')
                pair.no[i] <- pair.location
            }
        }
    } else if(!is.null(pair.no) && !is.null(pairs))
    {
        tmp.pair.no <- NULL
        for (i in 1:nrow(pairs)) 
        {
            if(is.null(original.pairs))
            {
                warning('No need to set \'pair.no\' and \'pairs\'. \nThe procedure does not involve neither \'pair.no\' nor \'pairs\' setting in this case.')
                pairs = NULL
                pair.no = NULL
            }else {
                pair.location <- which(apply(original.pairs, 1, function(y) sum(y == pairs[i,]))==2)
                if(length(pair.location) == 0) stop('input wrong pairs')
                tmp.pair.no[i] <- pair.location
            }
        }
        if(!is.null(pair.no) && sum(tmp.pair.no == pair.no[order(pair.no)]) != length(pair.no)) stop('pair number does not match pairs') #check whether the pair number corresponding to the input pairs are equal to the input pair.no
    }else
    {
        if(is.null(original.pairs)) 
        {
            warning('No need to set \'pair.no\'. \nThe procedure does not involve the \'pair.no\' setting in this case.')
            pair.no <- NULL
        }else if(max(pair.no)>nrow(original.pairs) | min(pair.no)<=0) stop('input wrong pair number')
    }
    
    
    # detect the coordinates to be plotted
    if(is.list(x$sci$sci.lower)) 
    {
        p <- length(x$sci$sci.lower[[1]])
    } else
    {
        p <- length(x$sci$sci.lower)
    }
    if(is.null(sci.coordinate)) sci.coordinate <- c(1:p)
    if(max(sci.coordinate) > p) stop('input wrong sci.coordinate')

    
    
    # plot sci
    if(is.list(x$sci$sci.lower)) #more than two samples
    {

            for (i in 1:length(pair.no)) 
            {
                sci.upper <- x$sci$sci.upper[[pair.no[i]]]
                sci.lower <- x$sci$sci.lower[[pair.no[i]]]
                if(only.plot.reject == TRUE)
                {
                    reject.coordinate <- which(sci.upper<0|sci.lower>0)
                    sci.coordinate.tmp <- sci.coordinate[sci.coordinate %in% reject.coordinate]
                    if(length(sci.coordinate.tmp)==0){
                        warning(paste0("Pair ",pair.no[i],": No rejection coordinate in user's input SCI coordinates.\n  The procedure will plot the input SCI coordinates.\n"))
                        sci.coordinate.tmp <- sci.coordinate
                    }
                }else{
                    sci.coordinate.tmp <- sci.coordinate
                }
                #control the maximum number of coordinates in one plot.
                if(length(sci.coordinate.tmp)>20){
                    cat(paste0("Pair ",pair.no[i],": Only plot first 20 coordinates.\n"))
                    sci.coordinate.tmp <- sci.coordinate.tmp[1:20]
                }
                # need control the coordinate number in output?
                plotsci(sci.coordinate.tmp, sci.lower, sci.upper, pair.no[i], 
                        original.pairs[pair.no[i],], title, subtitle, xlabel, ylabel)
            }
            
    }else if(is.vector(x$sci$sci.lower)) # two samples
    {
        sci.upper <- x$sci$sci.upper
        sci.lower <- x$sci$sci.lower
        if(only.plot.reject == TRUE)
        {
            reject.coordinate <- which(sci.upper<0|sci.lower>0)
            sci.coordinate.tmp <- sci.coordinate[sci.coordinate %in% reject.coordinate]
        }else{
            sci.coordinate.tmp <- sci.coordinate
        }
        #control the maximum number of coordinates in one plot.
        if(length(sci.coordinate.tmp)>20){
            cat("Only plot first 20 coordinates.\n")
            sci.coordinate.tmp <- sci.coordinate.tmp[1:20]
        }
        plotsci(sci.coordinate.tmp, sci.lower, sci.upper, pair.no, 
                original.pairs, title, subtitle, xlabel, ylabel)
    }
}


# sci plot function
plotsci <- function(sci.coordinate, sci.lower, sci.upper, pair.no, 
                    original.pairs, title, subtitle, xlabel, ylabel)
{ 
    if(is.null(subtitle)) 
    {
        if(is.null(original.pairs)) subtitle.default <- 'One sample SCI'
        else subtitle.default <- paste0('Pair ', pair.no, ': Group ', original.pairs[1], ' VS Group ', original.pairs[2])
    }else subtitle.default <- subtitle
    
    yl <- max(abs(cbind(sci.lower, sci.upper)))
    yl <- c(-yl, yl)
    print(ggplot(mapping = aes(x=as.character(sci.coordinate))) + 
             geom_hline(yintercept = 0, color = "black") +
             ylim(yl) +
             geom_errorbar(mapping=aes(ymin=sci.lower[sci.coordinate], ymax=sci.upper[sci.coordinate]), 
                           width=0.3, size=1, color = 4) + 
             geom_point(mapping = aes(x=as.character(sci.coordinate[which(sci.upper[sci.coordinate]<0|sci.lower[sci.coordinate]>0)]), 
                                      y=rep(0,length(which(sci.upper[sci.coordinate]<0|sci.lower[sci.coordinate]>0)))), 
                        color = "red", size = 3, fill = "white") + 
             scale_x_discrete(limits=c(as.character(sci.coordinate))) +
             labs(title = title, subtitle = subtitle.default, y = ylabel, x = xlabel) + 
             theme_bw())
}

