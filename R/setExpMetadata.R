
#' @title Set metadata from an experiment.
#' @param x The `roleExperiment` object set metadata to. 
#' @param author The author of the RoLE experiment.
#' @param datestring The date the RoLE experiment was run.
#' @param description A short description of the RoLE experiment.
#' @details Each argument besides `x` should be a single string.

#' @rdname setExpMetadata
#' @export

setGeneric('setExpMetadata', 
           def = function(x, author, datestring, description) standardGeneric('setExpMetadata'), 
           signature = 'x')

#' @rdname setExpMetadata
#' @export

setMethod('setExpMetadata', 
          signature = 'roleExperiment', 
          definition = function(x, author, datestring, description) { 
              x@auxMeta <- setNames(c(author, datestring, description), 
                                                        c("author", "date","description"))
              return(x)
          }
)
