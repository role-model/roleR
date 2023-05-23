# print method for `roleModel

setMethod('show', signature = signature(object = 'roleModel'),
          definition = function(object) {
              is_run <- !is.null(object@modelSteps[[2]]) 
              
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


# set coercion method from `roleModel` to `roleExperiment`
setAs(from = 'roleModel', to = 'roleExperiment',
      def = function(from) {
          # get metadata
          pout <- from@info
          
          # put an index column in that data.frame for model 
          pout <- cbind(mod_id = 1, pout)
          
          return(new('roleExperiment', 
                     info = pout,
                     modelRuns = from@modelSteps, 
                     allParams = list(from@params), 
                     inits = list()))
      }
)

