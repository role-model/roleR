# print method for `roleModel
#' @export

setMethod('show', signature = signature(object = 'roleModel'),
          definition = function(object) {
              niter <- 1
              #nmod <- length(unique(object@experimentMeta$mod_id))
              cat(sprintf('RoLE model with %s iterations', 
                          niter, 
                          ifelse(niter == 1, '', 's')), 
                  '\n')
          }
)