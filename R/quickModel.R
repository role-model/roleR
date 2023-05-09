#' @title Internal functions to facilitate vignette examples.
#' @include roleParams.R
#' @export

quickParams <- function(){
    # its for some reason using the value of the last supplied non function in every fun 
    p <- roleParams(niter = 200, niterTimestep = 20, speciation_local = .2)
    return(p)
}

# function to create a model using a set of reasonable params
# used in testthat tests
# will eventually no longer be exported
#' @title quickModel
#'
#' @return a roleModel
#' @export
#'
quickModel <- function(){
    p <- quickParams()
    m <- roleModel(p)
    return(runRole(m)) 
}

# function to create a model using a set of reasonable params
# used in testthat tests
# will eventually no longer be exported
#' quickExp
#'
#' @return a run roleExperiment
#' @export
quickExp <- function(){
    
    p <- quickParams()
    expr <- roleExperiment(list(p,p,p))
    
    return(runRole(expr))
}

#' quickModelNonRun
#'
#' @return a setup but not run roleModel
#' @export
quickModelNonRun <- function(){
    p <- quickParams()
    m <- roleModel(p)
    return(m)
}