#' @title roleExperiment - one or more models bundled collectively 
#' 
#' @description An S4 class to represent a self-enclosed modeling experiment or intentioned set of models.
#' It contains a list of `roleModel`objects, a list of `roleParams` to use for each model, a data.frame summary of the models, and some author metadata
#' 
#' @slot modelRuns a list of `roleModel`objects
#' @slot allParams a list of `roleParams` to use for each model
#' @slot authorMeta a named string vector that keep track of author metadata
#' It contains values for "author", "date", "description", "info", where each element is named by its respective string.
#' When the model is saved with `writeRole` a text file is generated using this metadata
#' @slot experimentMeta data.frame of summarizing metadata for all the models of the experiment 
#' (Need to chat with Andy about the exact intentions of this before writing more here)
#' 
# @examples 
# Create and run a roleExperiment that will contain three models with three different degrees of neutrality
# NOTE - this does not work yet as defaults have not yet been discussed and implemented
# p1 <- roleParams(neut_delta=0.3)
# p2 <- roleParams(neut_delta=0.6)
# p3 <- roleParams(neut_delta=1)
# expr <- roleExperiment(list(p1,p2,p3))
# 
# TEMPORARY WORKAROUND
# m <- quickModel()
# p1 <- m@params
# p2 <- p1
# p3 <- p1
# p1@neut_delta <- 0.3
# p2@neut_delta <- 0.6
# p3@neut_delta <- 1
# expr <- roleExperiment(list(p1,p2,p3))
# expr <- runRole(expr)
# 
# Attach author metadata to the experiment
# (in flux, could either make method to attach manually OR force attaching in constructor, which I like more)
#' 
#' @rdname roleExperiment
#' @export

setClass('roleExperiment',
         slots = c(authorMeta = 'character',
                   experimentMeta = 'data.frame',
                   modelRuns = 'list', 
                   allParams = 'list',
                   allFuns = 'list'))

#' @title Create a roleExperiment
#' @param allParams the list of model params to use for each model
#' @return a ready-to-run `roleExperiment`
#' 
#' @rdname roleExperiment
#' @export

roleExperiment <- function(allParams) {
    # if given priors or iterfuns, run functions on each model, param, and iter 
    # and add to new allParams object
    if(class(allParams) == 'rolePriors' | class(allParams) == 'roleIterFuns'){
        # save allFuns
        allFuns <- allParams
        # create new allParams
        allParams <- list()
        # for each model's set of priors
        for(p in 1:length(allFuns)){
            funs <- allFuns[[p]]
            params <- roleParams()
            # for each slot
            for(slot_name in slotNames(funs))
            {
                if(class(allFuns) == 'rolePriors'){
                    # get the list of funs, one per iteration
                    fun_list <- slot(priors,slot_name)
                    # for each iteration/functions
                    for(f in 1:length(fun_list)){
                        fun = fun_list[[f]]
                        slot(params,slot_name)[f] = fun()
                    }
                }
                else{
                    # get the list of funs, one per iteration
                    fun <- slot(iterfuns,slot_name)
                    # for each iteration/functions
                    for(i in 1:niter){
                        slot(params,slot_name)[i] = fun(i)
                    }
                }
                allParams <- list(allParams,params)
            }
        }
    }
    else{
        allFuns <- list(NULL)
    }
    
    # create models from params and add to experiment
    allModels <- list()
    for(i in 1:length(allParams)){
        allModels <- append(allModels,roleModel(allParams[[i]]))
    }
    
    return(new('roleExperiment', 
               experimentMeta = data.frame(), 
               modelRuns = allModels,
               allParams=allParams,
               allFuns=allFuns))
}


# set coercion method from `roleModel` to `roleExperiment`
setAs(from = 'roleModel', to = 'roleExperiment',
      def = function(from) {
          niter <- from@params@niter
          niterTimestep <- from@params@niterTimestep
          
          # index of when data were saved
          j <- c(0, which(1:niter %% niterTimestep == 0)) + 1
          
          # extract the param vals for those j indeces
          pnames <- slotNames(from@params)
          pnames <- pnames[!grepl('niter', pnames)]
          
          pout <- lapply(pnames, function(p) {
              p <- slot(from@params, p)
              
              if(length(p) > 1) {
                  p <- c(p[1], p)
                  return(p[j])
              } else {
                  return(rep(p, length(j)))
              }
          })
          
          # put param vals in a data.frame
          names(pout) <- pnames
          # browser()
          pout <- do.call(data.frame, pout)
          
          # add time information
          pout$iterations <- j - 1
          pout$generations <- pout$iterations / (pout$individuals_local / 2)
          
          # put an index column in that data.frame for model run
          pout <- cbind(mod_id = 1, pout)
          
          return(new('roleExperiment', 
                     experimentMeta = pout,
                     modelRuns = from@modelSteps, 
                     allParams = list(from@params), 
                     # iterFuns = list(from@iterFuns) # when iterFuns added to roleModel, uncomment this and delete line about `allFuns`
                     allFuns = list(NULL)))
      }
)


# rbind method for `roleExperiment` class

setMethod('rbind2', signature = c('roleExperiment', 'roleExperiment'), 
          definition = function(x, y) {
              out <- x
              thisEx <- y
              
              # keep track of growing mod_id max index
              j <- max(out@experimentMeta$mod_id)
              
              # augment mod_id
              thisEx@experimentMeta$mod_id <-
                  thisEx@experimentMeta$mod_id + j
              
              # combine with `out`
              out@experimentMeta <- rbind(out@experimentMeta,
                                          thisEx@experimentMeta)
              out@modelRuns <- c(out@modelRuns,
                                 thisEx@modelRuns)
              out@allParams <- c(out@allParams,
                                 thisEx@allParams)
              
              return(out)
          }
)

setMethod('rbind2', signature = c('roleExperiment', 'missing'), 
          definition = function(x, y) {
              return(x)
          })

# temporary helper to rep an S4 object
# used to create roleExperiments with many models of the same params
# i.e. roleExperiment(repS4(p,100)) does this 100 times using roleParams p
repS4 <- function(obj, n) {
    # Create an empty list to store replicated objects
    objs <- vector("list", n)
    
    # Replicate the object and store in the list
    for (i in 1:n) {
        objs[[i]] <- obj
    }
    
    # Return the list of replicated objects
    return(objs)
}
