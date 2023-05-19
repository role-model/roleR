# print method for `roleModel

setMethod('show', signature = signature(object = 'roleModel'),
          definition = function(object) {
              is_run <- !is.null(object@modelSteps[[2]]) # see if model has been run
              
              if(is_run){
                  run_str <- "completed (run)"
              }
              
              else{
                  run_str <- "not-yet-run"
              }
              
              niter <- object@params@niter
              
              cat(sprintf('%s RoLE model with %s iterations',
                          run_str, niter), '\n')
          }
)
