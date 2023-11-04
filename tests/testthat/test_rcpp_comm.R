# test that you can give comm to rcpp and what it gives you back is the same

p <- roleParams()

m <- roleModel(p)

roleR:::roleCommTester(m@modelSteps[[1]], p)
