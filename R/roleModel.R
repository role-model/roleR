#' @title An S4 class to specify the entire RoLE model
#'
#' @slot localComm an object of class \code{localComm}
#' @slot metaComm an object of class \code{metaComm}
#' @slot phylo an object of class \code{rolePhylo}
#' @slot params an object of calss \code{roleParams}
#'
#' @export
#' @include comm.R roleParams.R rolePhylo.R


setClass('roleModel',
         slots = c(localComm = 'localComm',
                   metaComm = 'metaComm',
                   phylo = 'rolePhylo',
                   params = 'roleParams'))


#' @title Create a RoLE model object
#'
#' @param localComm an object of class \code{localComm}
#' @param metaComm an object of class \code{metaComm}
#' @param phylo an object of class \code{rolePhylo}
#' @param params an object of class \code{roleParams}
#'
#' @return an object of class \code{roleModel}
#'
#' @export


roleModel <- function(local, meta, phy, p) {
    new('roleModel',
        localComm = local,
        metaComm = meta,
        phylo = phy,
        params = p)
}

# checker function for validation
#' @param r an object of class roleModel

checkModel <- function(r) {
    checks <- c()

    if(!validObject(r@localComm)){
        checks <- c(checks, 'model localComm invalid')
    }
    if(!validObject(r@metaComm)){
        checks <- c(checks, 'model metaComm invalid')
    }
    if(!validObject(r@phylo)){
        checks <- c(checks, 'model rolePhylo invalid')
    }
    if(!validObject(r@params)){
        checks <- c(checks, 'model params invalid')
    }


    #make sure local traits values match trait values in metacommunity
    localTraits = na.omit(r@localComm@traits)
    #if all non-NA local trait values match respective meta traits, we are good
    if(all(localTraits != r@metaComm@traits[localTraits[,1],])){
       checks <- c(checks, 'local trait values must match metacommunity
                   trait values')
    }

    # if any issues, return them, otherwise all OK
    if(length(checks) > 0) {
        return(checks)
    } else {
        return(TRUE)
    }
}

# validate
setValidity('roleModel', checkModel)
