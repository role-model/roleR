#' @title An S4 class to specify parameters of the RoLE model
#'
#' @slot params a named list of numerical values of the parameters; if priors
#' are specified, \code{params} are used as initial values
#' @slot runType a character vector of length one specifying the type of run
#' (either "sim" or "fit")
#' @slot priors a names list of functions that generate samples from the
#' desired prior distributions
#'
#'
#' @export

setClass('roleParams',
         slots = c(params = 'list',
                   runType = 'character',
                   priors = 'list'))


#' @title Specify parameters for a RoLE model
#'
#' @param params a named list of parameter values, used as initial values if
#' priors are specified
#' @param runType can be either "sim" or "fit", specifying whether to run
#' simulations only, or also fit the model to data
#' @param priors a named list of functions to be used to generate samples from
#' the prior distributions of parameters
#'
#' @return an object of class \code{roleParams}
#'
#' @seealso \code{\link{roleParams-class}}
#'
#' @export

roleParams <- function(params, runType, priors) {
    new('roleParams',
        params = params, runType = runType, priors = priors)
}


#' checker function for validation
#' @param object an object of class roleParams

checkParams <- function(object) {
    checks <- c()

    # check parameter names
    allowedParams <- c('foo', 'bar')
    namedParams <- names(object@params)

    if(is.null(namedParams)) {
        checks <- c(checks,
                    'parameters must be named')
    } else {
        badParams <- namedParams[!(namedParams %in% allowedParams)]

        if(length(badParams) > 0) {
            checks <- c(checks,
                        paste('the following parameters are not allowed:',
                              paste(badParams, collapse = ', ')))
        }
    }

    # check that parameter values are numeric
    paramTypes <- sapply(object@params, class)
    if(!all(paramTypes %in% c('numeric', 'integer'))) {
        checks <- c(checks,
                    'all parameters must be coercible to numeric')
    }

    # check runtype
    runType <- object@runType

    if(length(runType) != 1) {
        checks <- c(checks,
                    'only one runtype can be specified')
    } else {
        if(!(runType %in% c('sim', 'fit'))) {
            checks <- c(checks,
                        '`runType` must be either "sim" or "fit"')
        }
    }


    # check that priors are functions or NULL
    priors <- object@priors

    if(!is.null(priors)) {
        if(!all(sapply(priors, class) == 'function')) {
            checks <- c(checks,
                        'all priors must be functions or NULL')
        }

        # check that param and prior names match
        if(is.null(names(priors)) | !all(names(priors) == namedParams)) {
            checks <- c(checks,
                        'parameter and prior names do not match')
        }
    } else {
        # if priors are null, we can't fit
        if(object@runType == 'fit') {
            checks <- c(checks,
                        'priors must be specified for `runType = "fit"`')
        }
    }


    # if any issues, return them, otherwise all OK
    if(length(checks) > 0) {
        return(checks)
    } else {
        return(TRUE)
    }
}

#' validate
setValidity('roleParams', checkParams)



# test
# bad <- roleParams(params = list(doo = 2, boo = 1), runType = 'sim',
#                   priors = list(foo = function(x) {x}))
# good <- roleParams(params = list(foo = 2), runType = 'sim',
#                    priors = list(foo = function(x) {x}))
