#' @title A RoLE experiment: one or more models used collectively to ask a question or probe a hypothesis 
#' 
#' @description An S4 class to represent a complete self-enclosed experiment (a 
#'     model run or set of runs with associated meta data). Contains one or more 
#'     `roleModel` runs and the `roleParams` used to generate each of the runs
#' 
#' @slot auxMeta named string vector of author etc. 
#' @slot experimentMeta data.frame of model metadata
#' @slot modelRuns list of `roleData` objects containing outputs
#' @slot allParams list of `roleParams` objects containing input params for 
#'     each run
#' @examples 
#' Create a roleExperiment using a set of params
#' p1 <- roleParams(speciation_local=0.2)
#' p2 <- roleParams(speciation_local=0.3)
#' p3 <- roleParams(speciation_local=0.4)
#' exp <- roleExperiment(list(p1,p2,p3))
#' 
#' @rdname roleExperiment
#' @export

setClass('roleExperiment',
         slots = c(auxMeta = 'character',
                   experimentMeta = 'data.frame',
                   modelRuns = 'list', 
                   allParams = 'list',
                   allFuns = 'list'))

# constructor for roleExperiment
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