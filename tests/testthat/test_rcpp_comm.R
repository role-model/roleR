# test that you can give comm to rcpp and what it gives you back is the same

n <- 10000000
p <- roleParams(niter = n, 
                niterTimestep = n / 2,
                species_meta = 500, 
                individuals_meta = 100000, 
                individuals_local = 5000)

m <- roleModel(p)
mr <- runRole(m)



# roleR:::roleCommTester(m@modelSteps[[1]], p)
# getSumStats(mr, funs = list("hill_abund" = hillAbund, "rich" = richness))


p <- roleParams(init_type = "oceanic_island",
                species_meta = 500,
                individuals_meta = 1e+06,
                individuals_local = 500,
                niter = 1, 
                niterTimestep = 1)
m <- roleModel(p)
rm <- runRole(m)

