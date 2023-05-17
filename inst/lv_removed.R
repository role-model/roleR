
# constructor
#' @rdname untbParams
#' @export

sp1_gr = 0.1
sp2_gr = 0.6
sp1_k = 50
sp2_k = 70
alpha12 = 0.1
alpha21 = 0.2
niter = 1000
lvParams <- function(sp1_n,sp2_n,sp1_gr,sp2_gr,sp1_k,sp2_k,alpha12,alpha21,niter) {
    max_gr_sp <- which.max(c(sp1_gr,sp2_gr))
    max_gr <- max(sp1_gr,sp2_gr)
    
    mu10 <- ifelse(max_gr_sp == 1, 0, max_gr - sp1_gr)
    mu20 <- ifelse(max_gr_sp == 2, 0, max_gr - sp2_gr)
    
    la10 <- sp1_gr + mu10
    la20 <- sp2_gr + mu20
    
    mu11 <- sp1_gr / sp1_k
    mu22 <- sp2_gr / sp2_k
    
    mu12 <- sp1_gr * alpha12 / sp1_k
    mu21 <- sp2_gr * alpha21 / sp2_k
    
    return(roleParams(
        individuals_local = sp1_k + sp2_k, # J (number of indv + rocks) = the total carrying capacity
        individuals_meta = individuals_meta,
        species_meta = species_meta,
        speciation_local = speciation,
        speciation_meta = 0.8,
        extinction_meta = 0.05,
        trait_sigma = 1,
        env_sigma = 1,
        comp_sigma = 1,
        neut_delta = 1, # makes the model neutral by ignoring env and comp sigmas
        env_comp_delta = 1,
        dispersal_prob = dispersal_prob,
        mutation_rate = 1e-7,
        equilib_escape = 1,
        alpha = 10,
        num_basepairs = 250,
        init_type = init_type, 
        niter = niter,
        niterTimestep = niterTimestep))
}
