# load rcpp function
Rcpp::sourceCpp("inst/specPhyloRCpp.cpp")

# make a phylo object
tre <- ape::read.tree(text = '(A:2,(B:1,C:1):1);')

# convert to rolePhylo
phy2role <- function(from) {
  # extract number of times
  n <- ape::Ntip(from)
  
  # extract edge matrix and edge lengths
  e <- from$edge
  l <- from$edge.length
  
  # extract tip labels
  tipNames <- from$tip.label
  
  # calculate alive or not
  tipAge <- ape::node.depth.edgelength(from)[1:n]
  
  alive <- rep(TRUE, n)
  alive[tipAge < max(tipAge)] <- FALSE
  
  
  # set default scale
  scale <- 1
  
  
  # buffer objects so we can add new species without augmenting objects
  addOn <- n * 4
  e <- rbind(e, matrix(-1, nrow = addOn, ncol = 2))
  l <- c(l, rep(0, addOn))
  alive <- c(alive, rep(FALSE, addOn))
  
  
  return(list(n = n, e = e, l = l, alive = alive, scale = scale))
}

# make a role phylo object (old s4 version)
rolePhy <- phy2role(tre)

# have a look at it
rolePhy

# run the rcpp function on it
specPhyloRCpp(i = 2, n = rolePhy$n, e = rolePhy$e, l = rolePhy$l, 
              alive = rolePhy$alive)

# looks like it's working correctly!!! (s4 object is updating correctly)
rolePhy
