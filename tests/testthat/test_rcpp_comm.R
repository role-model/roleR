# test that you can give comm to rcpp and what it gives you back is the same

p <- roleParams(niter = 10000000, 
                niterTimestep = 1000000,
                species_meta = 500, 
                individuals_meta = 100000, 
                individuals_local = 5000)

m <- roleModel(p)

# roleR:::roleCommTester(m@modelSteps[[1]], p)

mr <- runRole(m)

# getSumStats(mr, funs = list("hill_abund" = hillAbund, "rich" = richness))
