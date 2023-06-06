library(roleR)

foopars <- roleParams(niter = 1000)
foomod <- roleModel(foopars)
foomod <- runRole(foomod)

getSumStats(foomod)
