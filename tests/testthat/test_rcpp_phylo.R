# test that you can give something to rcpp and what it gives you back is the same

test_that("rcpp trimming gives you back the same thing you give it", {
    m <- roleParams(species_meta = 5) |> 
        roleModel()
    r <- runRole(m)
    
    # initialized phylo (never saw rcpp or trimming)
    tre0 <- as(m@modelSteps[[1]]@phylo, "phylo") 
    
    # initialized phylo but having gone through rcpp stuff
    tre1 <- as(r@modelSteps[[1]]@phylo, "phylo") 
    
    expect_equivalent(tre0, tre1)
})



# test updatePhylo with simple phylo

library(ape)

# make a tree
x <- read.tree(text = "((A:2,D:1):2,(B:1,C:1):3);")
xrole <- as(x, "rolePhylo")
fossilID <- 2 # keep track of which tip is fossil 

# run it through role phylo stuff
iparent <- 1 # parent tip for speciation
s <- 0.1 # scale of converting iterations to edge length

# run update phylo and convert back to ape::phylo
yrole <- roleR:::testUpdatePhylo(xrole, iparent)
y <- as(yrole, "phylo")


# test ----

test_that("`updatePhylo` adds a new tip at the correct location", {
    # ancestor and sister of iparnet in new tree
    iAnc <- y$edge[y$edge[, 2] == iparent, 1]
    iSis <- y$edge[y$edge[, 1] == iAnc, 2]
    
    # new sister should be `Ntip(x) + 1`
    expect_equal(sort(iSis), c(iparent, Ntip(x) + 1))
})

# values that will be used in multiple subsequent tests
nspp <- Ntip(x)
newName <- y$tip.label[!(y$tip.label %in% x$tip.label)]

test_that("name of new species is `s<nspp + 1>`", {
    expect_equal(newName, paste0("s", nspp + 1))
})

test_that("dropping new species results in same topology as original tree", {
    expect_equal(drop.tip(y, newName)$edge, 
                 x$edge)
})


test_that("internal branches are uneffected", {
    yLessNew <- drop.tip(y, paste0("s", nspp + 1))
    
    expect_equal(x$edge.length[x$edge[, 2] > nspp], 
                 yLessNew$edge.length[yLessNew$edge[, 2] > nspp])
})
