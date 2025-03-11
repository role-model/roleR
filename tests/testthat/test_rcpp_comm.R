# test that you can give comm to rcpp and what it gives you back is the same
# roleR:::roleCommTester(m@modelSteps[[1]], p)


# test pure competition model
n <- 10
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 5, 
                individuals_meta = 100, 
                individuals_local = 10, 
                comp_sigma = 0.5,
                neut_delta = 0.5, env_comp_delta = 0)

m <- roleModel(p)
mr <- runRole(m)

# test pure filtering model

n <- 1000000
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 5, 
                individuals_meta = 100, 
                individuals_local = 10, 
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 0.5, env_comp_delta = 1)

m <- roleModel(p)

trt <- m@modelSteps[[1]]@localComm@indTrait
roleR:::envDistCalcTest(trt, p@env_optim, p@env_sigma)


mr <- runRole(m)


# test mixed model

n <- 1000000
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 5, 
                individuals_meta = 100, 
                individuals_local = 10, 
                comp_sigma = 0.5,
                env_sigma = 0.5,
                neut_delta = 0.5, env_comp_delta = 0.5)

m <- roleModel(p)
mr <- runRole(m)



# getSumStats(mr, funs = list("hill_abund" = hillAbund, "rich" = richness))


p <- roleParams(init_type = "oceanic_island",
                species_meta = 500,
                individuals_meta = 1e+06,
                individuals_local = 500,
                niter = 1, 
                niterTimestep = 1)
m <- roleModel(p)
rm <- runRole(m)

