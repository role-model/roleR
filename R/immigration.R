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

.immigrationLocal <- function(x, i) {

    #increment abundance of species in local community by 1
    x@abundance[i] <- x@abundance[i] + 1

    return(x)
}
setMethod('immigration', 'localComm', .immigrationLocal)


#' function to implement immigration for \code{roleModel} class objects
#' @param x an object of class \code{roleModel}

.immigrationRole <- function(x) {

    i <- sample(x@metaComm@Smax, 1,
                prob = x@metaComm@abundance[1:x@metaComm@Smax])
    x@localComm <- immigration(x@localComm, i)
    return(x)
}
setMethod('immigration', 'roleModel', .immigrationRole)


# TEST
# source("R/comm.R")
# source("R/immigration.R")
# foo <- localComm(1:10,matrix(1:100,nrow = 10, ncol = 10),1:10,10)
# foo@abundance[1]
# foo <- immigration(foo,1)
# foo@abundance[1]
