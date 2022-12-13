
#' @title Set metadata from an experiment
#' @param x the object set metadata to

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
