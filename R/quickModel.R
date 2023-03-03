
quickParams <- function(){
    # its for some reason using the value of the last supplied non function in every fun 
    p <- roleParams()
    return(p)
}

# function to create a model using a set of reasonable params
# used in testthat tests
# will eventually no longer be exported
quickModel <- function(){
    p <- quickParams()
    m <- roleModel(p)
    return(runRole(m)) 
}

# function to create a model using a set of reasonable params
# used in testthat tests
# will eventually no longer be exported
quickExp <- function(){
    
    p <- quickParams()
    expr <- roleExperiment(list(p,p,p))
    
    return(runRole(expr))
}

quickModelNonRun <- function(){
    p <- quickParams()
    m <- roleModel(p)
    return(m)
}