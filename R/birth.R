#' @title Birth
#'
#' @description Generic and methods for computing the birth process
#'
#' @param x the object which determines method dispatch
#' @param i the index of the species undergoing birth
#'
#' @rdname birth
#' @export

setGeneric('birth',
           function(x, i, ...) standardGeneric('birth'),
           signature = 'x')

#' function to implement birth for \code{localComm} class objects
#' @param x an object of class \code{localComm}
#' @param i the vector of indices of the species undergoing birth
#' @param params a \code{roleParams} object

.birthLocal <- function(x, i, params) {

    #increment abundance of species in local community by 1
    x@abundance[i] <- x@abundance[i] + 1

    #intraspecific models will need code for trait inheritance etc later

    return(x)
}

# set the method
setMethod('birth', 'localComm', .birthLocal)
