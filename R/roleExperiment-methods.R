# 
# # Extract parts of a roleExperiment
# # @aliases [,roleExperiment-method
# # @docType methods
# # @param x a roleExperiment object
# # @param i i
# # @param j j
# # @param ... ignored
# # @param drop default F
# # @include roleExperiment.R
# 
# setMethod("[",
#          signature(x = "roleExperiment", i = "ANY", j = "ANY"),
#           function(x, i, j, ..., drop = FALSE) {
#     iMissing <- missing(i)
#     jMissing <- missing(j)
#     nargs <- nargs() # e.g., a[3,] gives 2 for nargs, a[3] gives 1.
# 
#     if(iMissing && jMissing) {
#         i <- TRUE
#         j <- TRUE
#     } else if(jMissing && !iMissing) {
#         if (nargs == 2) {
#             j <- i
#             i <- TRUE
#         } else {
#             j <- TRUE
#         }
#     } else if(iMissing && !jMissing)
#         i <- TRUE
# 
#     if(is.matrix(i)) {
#         stop('matrix argument not supported')
#     }
# 
#     if(any(is.na(i))) {
#         stop('NAs not permitted in row index')
#     }
# 
#     # vector of unique model run IDs before selecting
#     origModID <- unique(x@info$mod_id)
# 
#     # update metadata and runs
#     x@info <- x@info[i, j, ..., drop = FALSE]
#     x@modelRuns <- x@modelRuns[i]
# 
#     # check to see if any unique model runs have been dropped
#     theseDropped <- !(origModID %in% x@info$mod_id)
# 
#     # if so, drop them
#     if(any(theseDropped)) {
#         i <- origModID[!theseDropped]
#         x@allParams <- x@allParams[i]
# 
#         # re-number IDs in metadata to remove any missing
#         x@info$mod_id <- order(x@info$mod_id)
#     }
# 
#     return(x)
# })
# 
# 
# # Extract a name
# # @name $
# # @aliases $,roleExperiment-method
# # @param x a roleExperiment
# # @param name thing to extract
# # @docType methods
# # @rdname roleExperiment
# 
# setMethod("$", "roleExperiment", function(x, name) {
#     if(name %in% names(x@info)) {
#         return(x@info[, name])
#     } else {
#         stop('no $ method for object without metadata')
#     }
# })
# 
# 
# # print method for `roleExperiment`
# setMethod('show', signature = signature(object = 'roleExperiment'),
#           definition = function(object) {
#               nmod <- length(unique(object@info$mod_id))
# 
#               is_run <- !is.null(object@modelRuns[[2]]) # see if model has been run
# 
#               if(is_run){
#                   run_str <- "completed (run)"
#               }
# 
#               else{
#                   run_str <- "not-yet-run"
#               }
# 
#               cat(sprintf('%s RoLE experiment with %s unique model%s',
#                           run_str,
#                           nmod,
#                           ifelse(nmod == 1, '', 's')),
#                   '\n')
#           }
# )
# 
# 
# 
# # rbind method for `roleExperiment` class
# 
# setMethod('rbind2', signature = c('roleExperiment', 'roleExperiment'),
#           definition = function(x, y) {
#               out <- x
#               thisEx <- y
# 
#               # keep track of growing mod_id max index
#               j <- max(out@info$mod_id)
# 
#               # augment mod_id
#               thisEx@info$mod_id <- thisEx@info$mod_id + j
# 
#               # combine with `out`
#               out@info <- rbind(out@info, thisEx@info)
#               out@modelRuns <- c(out@modelRuns, thisEx@modelRuns)
#               out@allParams <- c(out@allParams, thisEx@allParams)
#               out@inits <- c(out@inits, thisEx@inits)
# 
#               return(out)
#           }
# )
# 
# 
# 
# 
# setMethod('rbind2', signature = c('roleExperiment', 'missing'),
#           definition = function(x, y) {
#               return(x)
#           })
# 
# 
