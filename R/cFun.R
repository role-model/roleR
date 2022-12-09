cFun <- function(type="int",fun_name="sample_index_using_probs",
                   data=NULL,params=NULL,i=NULL,
                   probs=NULL,
                   dead_index=NULL,parent_indv=NULL){
    
    if(!is.null(data)){
        l <- data@localComm
        m <- data@metaComm
        phy <- data@phylo
    }
    if(!is.null(params)){
        niter <- length(params@individuals_local)
    }
    
    if(type == "int"){
        return(intFunCpp(fun_name=fun_name,data=data,params=params,probs=probs))}
    if(type == "vector"){
        return(vectorFunCpp(fun_name=fun_name,data=data,params=params,probs=probs))}
    if(type == "data"){
        return(dataFunCpp(fun_name=fun_name,local=l,meta=m,phylo=phy,
                          params=params,niter=niter,i=i,
                          dead_index=dead_index,parent_indv=parent_indv,probs=probs))}
}