cFun <- function(type="int",fun_name="sample_index_using_probs",
                   data=NULL,params=NULL,i=0, # used universally
                   probs=c(0), x = 0, # used in intFuns sample_zero_to_x and sample_index_from_probs
                   dead_index=0, # used universally
                   parent_indv=0, # used by call_birth and call_dispersal
                   dispersed_this_iter=TRUE, # used by call_speciation and update_speciation_local_meta
                   speciation_sp=0 # used in update_speciation_local_meta
                 ){
    
    if(!is.null(data)){
        l <- data@localComm
        m <- data@metaComm
        phy <- data@phylo
    }
    if(!is.null(params)){
        niter <- length(params@individuals_local)
    }
    
    if(type == "int"){
        return(intFunCpp(fun_name=fun_name,probs=probs,x=x))}
    if(type == "vector"){
        return(vectFunCpp(fun_name=fun_name,local=l,meta=m,phylo=phy
                            ,params=params,niter=niter,i=i))}
    if(type == "data"){
        print("calling dataFunCpp")
        return(dataFunCpp(fun_name=fun_name,local=l,meta=m,phylo=phy,
                          params=params,niter=niter,i=i,
                          dead_index=dead_index,
                          parent_indv=parent_indv,
                          probs=probs,
                          dispersed_this_iter=dispersed_this_iter,
                          speciation_sp=speciation_sp))}
}