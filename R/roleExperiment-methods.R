# selection method for `roleExperiment`
#' @export

setMethod("[", "roleExperiment", function(x, i, j, ..., drop = FALSE) {
    iMissing <- missing(i)
    jMissing <- missing(j)
    nargs <- nargs() # e.g., a[3,] gives 2 for nargs, a[3] gives 1.
    
    if(iMissing && jMissing) {
        i <- TRUE
        j <- TRUE
    } else if(jMissing && !iMissing) { 
        if (nargs == 2) {
            j <- i
            i <- TRUE
        } else {
            j <- TRUE
        }
    } else if(iMissing && !jMissing)
        i <- TRUE
    
    if(is.matrix(i)) {
        stop('matrix argument not supported')
    }
    
    if(any(is.na(i))) {
        stop('NAs not permitted in row index')
    }
    
    # vector of unique model run IDs before selecting
    origModID <- unique(x@experimentMeta$mod_id)
    
    # update metadata and runs
    x@experimentMeta <- x@experimentMeta[i, j, ..., drop = FALSE]
    x@modelRuns <- x@modelRuns[i]
    
    # check to see if any unique model runs have been dropped
    theseDropped <- !(origModID %in% x@experimentMeta$mod_id)
    
    # if so, drop them
    if(any(theseDropped)) {
        i <- origModID[!theseDropped]
        x@allParams <- x@allParams[i]
        
        # re-number IDs in metadata to remove any missing
        x@experimentMeta$mod_id <- order(x@experimentMeta$mod_id)
    }
    
    return(x)
})


# column selection method for `roleExperiment`
#' @export

setMethod("$", "roleExperiment", function(x, name) {
    if(name %in% names(x@experimentMeta)) {
        return(x@experimentMeta[, name])
    } else {
        stop('no $ method for object without metadata')
    }
})



# print method for `roleExperiment`
#' @export

setMethod('show', signature = signature(object = 'roleExperiment'),
          definition = function(object) {
              nmod <- length(unique(object@experimentMeta$mod_id))
              cat(sprintf('RoLE experiment with %s unique model%s', 
                          nmod, 
                          ifelse(nmod == 1, '', 's')), 
                  '\n')
          }
)
