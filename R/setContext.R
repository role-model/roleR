
#' @title Set context metadata to an experiment
#' @param x The `roleExperiment` object set context metadata to. 
#' @param author string containing the author name of the RoLE experiment.
#' @param datestring string containing the date the RoLE experiment was run
#' @param description string containing a short description of the RoLE experiment's purpose or intent 
#' @details Each argument besides `x` should be a single string.

#' @rdname setContext
#' @export

setGeneric('setContext', 
           def = function(x, author, datestring, description) standardGeneric('setContext'), 
           signature = 'x')

#' @rdname setContext
#' @export

setMethod('setContext', 
          signature = 'roleExperiment', 
          definition = function(x, author, datestring, description) { 
              x@context <- setNames(c(author, datestring, description), c("author", "date","description"))
              return(x)
          }
)
