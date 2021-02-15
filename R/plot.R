#' SCI Plot for pairs that lead to rejection
#' @description Plot the confidence interval of there pairs that lead to rejection after high-dimensional data analysis.
#' @param x a list of matrices generated from \code{hdtest}.
#' @param pair.no a vector of integers where the integers denote the orders in whole pairs in analysis; only these confidence intervals of pairs specified in \code{pair.no} or \code{pairs} will be plotted; default value: \code{NULL}, which is selected from \code{pairs}.
#' @param pairs a matrix with two columns where each row specifies a pair of populations in analysis; only these confidence intervals of pairs specified in \code{pair.no} or \code{pairs} will be plotted; default value: \code{NULL}, so that the specified pairs for plotting SCI are all pairs that lead to rejection.
#' @param sci.coordinate a vector of integers where each integer specifies the coordinate of SCI; only these specified coordinates of SCI will be plotted; default value: \code{NULL}, which is equivalent to all coordinates.
#' @param hline.color the color of horizontal line in plot; the color setting is the same as that in aesthetic specifications of \pkg{gglot2}; default value: \code{'black'}.
#' @param sci.width the width of whisker in confidence interval of each coordinate; the width setting is the same as that in aesthetic specifications of \pkg{gglot2}; default value: \code{0.5}.
#' @param sci.size the thickness of line in confidence interval of each coordinate in plot; the size setting is the same as that in aesthetic specifications of \pkg{gglot2}; default value: \code{1}.
#' @param sci.color the color of line in confidence interval of each coordinate in plot; the color setting is the same as that in aesthetic specifications of \pkg{gglot2}; default value: \code{'blue'}.
#' @param rej.point.color the color of the points used to mark the coordinates of which confidence interval not containing 0; the color setting is the same as that in aesthetic specifications of \pkg{gglot2}; default value: \code{'red'}.
#' @param rej.point.size the diameter of the points used to mark the coordinates of which confidence interval not containing 0; the color setting is the same as that in aesthetic specifications of \pkg{gglot2}; default value: \code{4}.
#' @param title the text for the title for plot; the title setting is the same as that in labels of \pkg{gglot2}; default value: \code{'Simultaneous Confidence Interval Plot'}.
#' @param subtitle the text for the subtitle for plot; the subtitle setting is the same as that in labels of \pkg{gglot2}; default value: \code{NULL} which indicates the pair number and group numbers.
#' @param xlabel the text for X coordinate for plot; the label setting is the same as that in labels of \pkg{gglot2}; default value: \code{'Coordinates'}.
#' @param ylabel the text for Y coordinate for plot; the label setting is the same as that in labels of \pkg{gglot2}; default value: \code{'Confidence Interval'}.
#' @param ... other arguments passed to plot.
#' @return plots of the SCIs of specified pairs
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
plot.hdaov <- function(x, pair.no = NULL, pairs = NULL, sci.coordinate = NULL, 
                       hline.color = 'black', sci.width = 0.5, sci.size = 1, 
                       sci.color = 'blue', rej.point.color = 'red', rej.point.size = 4, 
                       title = 'Simultaneous Confidence Interval Plot', subtitle = NULL, 
                       xlabel = 'Coordinates', ylabel = 'Confidence Interval', ...)
{
    if (!inherits(x, "hdaov")) 
        stop(gettextf("object must be of class %s", 
                      dQuote("hdaov")), domain = NA)
    
    if(!('sci' %in% names(x))) stop('\nObject must contain SCI.\nPlease run \'hdtest\' with \'return.sci = T\' again.')
    
    if(isTRUE(ncol(pairs) != 2 && !is.null(pairs))) stop('wrong format in pairs')
    
    pair.no <- unique(pair.no)
    sci.coordinate <- unique(sci.coordinate)
    if(!is.null(pairs)) pairs <- pairs[!duplicated(pairs),]
    if(!is.null(pairs) && !is.matrix(pairs)) pairs <- matrix(pairs,1,2)
    
    ori.pairs <- x$pairs
    only.pairs <- ori.pairs[,c(1,2)]
    if(!is.matrix(only.pairs) && !is.null(only.pairs)) only.pairs <- matrix(only.pairs,nrow(ori.pairs),2)
    
    #deal with pair number and pairs 
    if(is.null(pair.no) && is.null(pairs))
    {
        if(!is.null(ori.pairs)) pair.no <- which(ori.pairs[,'reject']==1)
        if(!('reject' %in% names(x))) stop('no rejection in object')
    } else if(is.null(pair.no) && !is.null(pairs))
    {
        for (i in 1:nrow(pairs)) 
        {
            if(is.null(only.pairs))
            {
                warning('No need to set \'pairs\'. \nThe procedure does not involve the \'pairs\' setting in this case.')
                pairs = NULL
            }else {
            pair.location <- which(apply(only.pairs, 1, function(y) sum(y == pairs[i,]))==2)
            if(length(pair.location) == 0) stop('input wrong pairs')
            pair.no[i] <- pair.location
            }
        }
    } else if(!is.null(pair.no) && !is.null(pairs))
    {
        tmp.pair.no <- NULL
        for (i in 1:nrow(pairs)) 
        {
            if(is.null(only.pairs))
            {
                warning('No need to set \'pair.no\' and \'pairs\'. \nThe procedure does not involve neither \'pair.no\' nor \'pairs\' setting in this case.')
                pairs = NULL
                pair.no = NULL
            }else {
                pair.location <- which(apply(only.pairs, 1, function(y) sum(y == pairs[i,]))==2)
                if(length(pair.location) == 0) stop('input wrong pairs')
                tmp.pair.no[i] <- pair.location
            }
        }
        if(!is.null(pair.no) && sum(tmp.pair.no == pair.no[order(pair.no)]) != length(pair.no)) stop('pair number does not match pairs')
    }else
    {
        if(is.null(only.pairs)) 
            {
            warning('No need to set \'pair.no\'. \nThe procedure does not involve the \'pair.no\' setting in this case.')
            pair.no <- NULL
        }else if(max(pair.no)>nrow(only.pairs) | min(pair.no)<=0) stop('input wrong pair number')
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
    if(is.list(x$sci$sci.lower)) #more than three samples
    {
        for (i in 1:length(pair.no)) 
        {
            sci.upper <- x$sci$sci.upper[[pair.no[i]]]
            sci.lower <- x$sci$sci.lower[[pair.no[i]]]
            plotsci(sci.coordinate, hline.color, sci.lower, sci.upper, sci.width, 
                    sci.size, sci.color, rej.point.color, rej.point.size, pair.no[i], 
                    only.pairs[pair.no[i],], title, subtitle, xlabel, ylabel)
        }
    }else if(is.vector(x$sci$sci.lower)) # two samples
    {
        sci.upper <- x$sci$sci.upper
        sci.lower <- x$sci$sci.lower
        plotsci(sci.coordinate, hline.color, sci.lower, sci.upper, sci.width, 
                sci.size, sci.color, rej.point.color, rej.point.size, pair.no, 
                only.pairs, title, subtitle, xlabel, ylabel)
    }
}

# sci plot function
plotsci <- function(sci.coordinate, hline.color, sci.lower, sci.upper, sci.width, 
                    sci.size, sci.color, rej.point.color, rej.point.size, pair.no, 
                    only.pairs, title, subtitle, xlabel, ylabel)
{ 
    if(is.null(subtitle)) 
    {
        if(is.null(only.pairs)) subtitle.default <- 'One sample SCI'
        else subtitle.default <- paste0('Pair ', pair.no, ': Group ', only.pairs[1], ' VS Group ', only.pairs[2])
    }else subtitle.default <- subtitle
    plot(ggplot(mapping = aes(x=sci.coordinate)) + 
             geom_hline(yintercept = 0, color = hline.color) + 
             geom_errorbar(mapping=aes(ymin=sci.lower[sci.coordinate], ymax=sci.upper[sci.coordinate]), 
                           width=sci.width, size=sci.size, color = sci.color) + 
             geom_point(mapping = aes(x=(sci.coordinate[which(sci.upper[sci.coordinate]<0|sci.lower[sci.coordinate]>0)]), 
                                      y=rep(0,length(which(sci.upper[sci.coordinate]<0|sci.lower[sci.coordinate]>0)))), 
                        color = rej.point.color, size = rej.point.size) + 
             labs(title = title, subtitle = subtitle.default, y = ylabel, x = xlabel) + 
             theme_bw())
}

