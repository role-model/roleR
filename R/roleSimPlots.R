#' @title Plotting RoLE model simulations
#'
#' @description Plot the dynamic distributions of RoLE model simulations and the
#' timeseries of RoLE model summary statitistics
#'
#'
#' @param x a \code{roleSim} object
#' @param type a string specifying which type of data to plot
#'
#' @details Stub
#'
# @examples
#'
#' @return A \code{plotly} object
#'
#' @rdname roleSimPlots
#' @export


roleDistPlot <- function(x, type = c('Abundance', 'Trait', 'pi', 'phylo')) {
    type <- match.arg(type, c('Abundance', 'Trait', 'pi', 'phylo'))

    # browser()
    x <- x$local_comm[[type]][x$local_comm$Abundance > 0]

    dat <- data.frame(rank = 1:length(x), y = sort(as.numeric(x), TRUE))

    fig <- plotly::plot_ly(dat, x = ~rank, y = ~y,
                           type = 'scatter', mode = 'markers',
                           marker = list(color = 'transparent', size = 8,
                                         line = list(color = 'black', width = 1.5)))
    fig <- plotly::layout(fig,
                          xaxis = list(zeroline = FALSE, title = 'Species Rank'),
                          yaxis = list(zeroline = FALSE, title = type,
                                       type = ifelse(type != 'Trait', 'log', 'linear')))
    return(fig)
}


#' @rdname roleSimPlots
#' @export


roleTSPlot <- function(x, type = c('Abundance', 'Trait', 'Pi', 'phylo')) {
    'stub'
}

