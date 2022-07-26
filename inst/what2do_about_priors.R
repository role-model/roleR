# option 1 ----

# make 1 model
params <- roleParams(....)
mod <- roleModel(params)
runMod <- run(mod) # doesn't need to return prior info...what about iter functions??

# make an experiment
exper <- roleExperiment(priors)
runExp <- run(exper) # `runExp` contains the data, the params, the priors

roleExperiment <- function(wtf) {
    if(class(wtf) == 'priors') {
        # sample from priors
        
        # return a list of params
        wtf <- listOfParams
    } else if(class(wtf) == 'iterFun') {
        # generate from iterFun
        
        # return a list of params
        wtf <- listOfParams
    }
    
    # now treat wtf as list of params
}

# in this option`niter`, `niterTimestep` and any "argument"-like params are still
# treated like things that we draw from priors, but we could draw them from a 
# "prior" like runif(1, 1000, 1000) (i.e. return the same value every time)

