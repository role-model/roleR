#' @title Run a `roleModel` or `roleExperiment`.
#' @description Run a RoLE object to completion.
#' @param x `roleModel` or `roleExperiment` object to run.
#' @return `roleModel` or `roleExperiment` run to completion.
#' @details Prior to running the RoLE model(s), the parameters must be specified inside the RoLE object with `roleParams()`.
#' 
#' @examples 
#' # create and run a model
#' model <- roleModel(roleParams())
#' model <- runRole(model)
#' 
#' # create and run an experiment
#' p1 <- roleParams(speciation_local=0.2)
#' p2 <- roleParams(speciation_local=0.3)
#' p3 <- roleParams(speciation_local=0.4)
#' exp <- roleExperiment(list(p1,p2,p3))
#' exp <- runRole(exp)
#' 
#' @rdname runRole
#' @export

setGeneric('runRole', 
           def = function(x, cores=1) standardGeneric('runRole'), 
           signature = 'x')

setMethod('runRole', 
          signature = 'roleModel', 
          definition = function(x) {
                m <- x
                
              for(i in 2:length(m@modelSteps))
              {
                  if(!is.null(m@modelSteps[[i]]))
                  {
                      stop("Models can only be run once and this model has already been run")
                  }
              }
              
              pvals <- getValuesFromParams(m@params, 1:m@params@niter)
              
              # augment the data in the model based on the params
              m <- .bufferModelData(m)
                  
              # returns the new modelSteps (a list of roleData)
              m@modelSteps <- iterModelCpp(slot(m@modelSteps[[1]],"localComm"), 
                                           slot(m@modelSteps[[1]],"metaComm"),
                                           slot(m@modelSteps[[1]],"phylo"),
                                           pvals,print=F)
              # trim data, removing the unused buffer
              m <- .trimModelData(m)
              
              # increment species ids and edges to shift indexing from Cpp to R
              for(d in 1:length(m@modelSteps))
              {
                  m@modelSteps[[d]]@localComm@indSpecies <- m@modelSteps[[d]]@localComm@indSpecies + 1
                  m@modelSteps[[d]]@phylo@e <- m@modelSteps[[d]]@phylo@e + 1
              }
              return(m)
          }
)

setMethod('runRole', 
          signature = 'roleExperiment', 
          definition = function(x, cores = 1) {
             
              if(cores == 1){
                  x@modelRuns <- lapply(x@modelRuns, runRole)
              }
              else{
                  cl <- snow::makeCluster(cores,type="SOCK")
                  x@modelRuns <- snow::clusterApply(cl, x@modelRuns, runRole)
                  snow::stopCluster(cl)
              }
              return(x)
          }
)

# user-inaccessible helper to augment the data of the not-yet-run model based on
#   what is expected from the params 
# called right before the model is run in Cpp
.bufferModelData <- function(model){
    p <- model@params 
    
    # calculate expected number of new species using binom
    expec_n_spec <- qbinom(0.9,p@niter,prob = mean(p@speciation_local(1:niter)))
    el_add <- (expec_n_spec * 2 - 1) + 1
    at_add <- expec_n_spec + 1
    
    # buffer phylo 
    model@modelSteps[[1]]@phylo@e
    model@modelSteps[[1]]@phylo@e <- rbind(model@modelSteps[[1]]@phylo@e, matrix(-1, nrow = el_add, ncol = 2)) # edges get -1s
    model@modelSteps[[1]]@phylo@e
    model@modelSteps[[1]]@phylo@l
    model@modelSteps[[1]]@phylo@l <- c(model@modelSteps[[1]]@phylo@l, rep(0, el_add)) # lengths get 0s
    model@modelSteps[[1]]@phylo@l
    model@modelSteps[[1]]@phylo@alive
    model@modelSteps[[1]]@phylo@alive <- c(model@modelSteps[[1]]@phylo@alive, rep(FALSE, at_add)) # alives get FALSE
    model@modelSteps[[1]]@phylo@alive
    model@modelSteps[[1]]@phylo@tipNames
    model@modelSteps[[1]]@phylo@tipNames <- c(model@modelSteps[[1]]@phylo@tipNames, rep('', at_add)) # tipNames get ''
    model@modelSteps[[1]]@phylo@tipNames
    
    # calc buffer size for local species vects 
    # 1 is the expected number of new species plus a small add
    # 2 is the initial number of species plus the expected number of new species plus a small add
    local_add <- p@species_meta + expec_n_spec
    
    # buffer local species vectors with 0s
    model@modelSteps[[1]]@localComm@spAbund
    model@modelSteps[[1]]@localComm@spAbund <- c(model@modelSteps[[1]]@localComm@spAbund,rep(0,local_add))
    model@modelSteps[[1]]@localComm@spAbund
    model@modelSteps[[1]]@localComm@spTrait
    model@modelSteps[[1]]@localComm@spTrait <- c(model@modelSteps[[1]]@localComm@spTrait,rep(0,local_add))
    model@modelSteps[[1]]@localComm@spTrait
    model@modelSteps[[1]]@localComm@spAbundHarmMean
    model@modelSteps[[1]]@localComm@spAbundHarmMean <-  rep(0,length(model@modelSteps[[1]]@localComm@spAbund) + local_add)
    model@modelSteps[[1]]@localComm@spAbundHarmMean
    model@modelSteps[[1]]@localComm@spLastOriginStep
    model@modelSteps[[1]]@localComm@spLastOriginStep <-  rep(0,length(model@modelSteps[[1]]@localComm@spAbund) + local_add)
    model@modelSteps[[1]]@localComm@spLastOriginStep
    model@modelSteps[[1]]@localComm@spExtinctionStep
    model@modelSteps[[1]]@localComm@spExtinctionStep <-  rep(0,length(model@modelSteps[[1]]@localComm@spAbund) + local_add)
    model@modelSteps[[1]]@localComm@spExtinctionStep
    
    return(model)
}

# user-inaccessible helper to trim the data of unused indices after the model is run
.trimModelData <- function(model){
    
    for(i in 1:length(model@modelSteps)){
        # trim phylo
        model@modelSteps[[i]]@phylo@e <- model@modelSteps[[i]]@phylo@e[model@modelSteps[[i]]@phylo@e[,1] != -1,]
        model@modelSteps[[i]]@phylo@l <- model@modelSteps[[i]]@phylo@l[model@modelSteps[[i]]@phylo@l != 0] 
        last_alive_index <- tail(which(model@modelSteps[[i]]@phylo@alive == TRUE), n = 1)
        model@modelSteps[[i]]@phylo@alive <- model@modelSteps[[i]]@phylo@alive[1:last_alive_index]
        model@modelSteps[[i]]@phylo@tipNames <- model@modelSteps[[i]]@phylo@tipNames[model@modelSteps[[i]]@phylo@tipNames != ''] # this MAY cause errors
        
        # trim local
        # find place where augmented 0s started, which is the number of species to ever have existed 
        # I think this is the last alive index??
        model@modelSteps[[i]]@localComm@spAbund <- model@modelSteps[[i]]@localComm@spAbund[1:last_alive_index]
        model@modelSteps[[i]]@localComm@spTrait <- model@modelSteps[[i]]@localComm@spTrait[1:last_alive_index]
        model@modelSteps[[i]]@localComm@spAbundHarmMean  <- model@modelSteps[[i]]@localComm@spAbundHarmMean [1:last_alive_index]
        model@modelSteps[[i]]@localComm@spLastOriginStep <- model@modelSteps[[i]]@localComm@spLastOriginStep[1:last_alive_index]
        model@modelSteps[[i]]@localComm@spExtinctionStep <- model@modelSteps[[i]]@localComm@spExtinctionStep[1:last_alive_index]
    }
    
    return(model)
}


# run iter functions over params to generate a new object of class paramValues 
# that contains ONLY vectors of values and no functions
#' @param p an object of class `roleParams`
#' @param i the iterations over which to calculate parameter values
#' @return a list of parameter vectors, each equal in length to `length(i)`

getValuesFromParams <- function(p, i) {
    # slot names
    s <- slotNames(p)
    
    # apply over slots
    pvals <- lapply(s, function(x) {
        thisP <- slot(p, s)
        
        if(inherits(thisP, 'function')) {
            return(thisP(i))
        } else {
            return(rep(thisP, length(i)))
        }
    })
    
    return(pvals)
}


