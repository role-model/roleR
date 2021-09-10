#' @title Immigration
#'
#' @description Generic and methods for computing the immigration process
#'
#' @param x the object which determines method dispatch
#' @param i the index of the species undergoing immigration
#'
#' @rdname immigration
#' @export

setGeneric('immigration',
           function(x, i, ...) standardGeneric('immigration'),
           signature = 'x')

#' function to implement immigration for \code{localComm} class objects
#' @param x an object of class \code{localComm}
#' @param i the vector of indices of the species undergoing immigration
#' @param params a \code{roleParams} object

.immigrationLocal <- function(x, i, params) {

    #immigration moving from metacommunity to local?

    #increment abundance of species in local community by 1

    x@abundance[i] <- x@abundance[i] + 1

    return(x)
}

# set the method
setMethod('immigration', 'localComm', .immigrationLocal)

# TEST
# source("R/comm.R")
# source("R/immigration.R")
# foo <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
# foo@abundance[1]
# foo <- immigration(foo,1)
# foo@abundance[1]
