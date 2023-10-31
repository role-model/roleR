# test that you can give something to rcpp and what it gives you back is the same

# test updatePhylo with simeple phylo
# 
# library(ape)
# 
# x <- read.tree(text = c("(A:3,(B:1,C:1):2);"))
# y <- unclass(x)
# padd <- 3
# y$alive <- c(rep(TRUE, length(y$tip.label)), rep(FALSE, padd))
# y$edge <- rbind(y$edge, matrix(-1, nrow = padd, ncol = 2))
# y$edge.length <- c(y$edge.length, rep(-1, padd))
# y$n <- 3
# 
# roleR:::testUpdatePhylo(y, 1)
