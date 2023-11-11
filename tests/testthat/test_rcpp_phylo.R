# test that you can give something to rcpp and what it gives you back is the same

# test updatePhylo with simeple phylo

library(ape)

# make a tree
x <- read.tree(text = "((A:2,D:1):2,(B:1,C:1):3);")


par(mfcol = c(2, 1), mar = c(2, 0, 0, 0) + 0.1)
plot(x, x.lim = c(0, 6))
axis(1)



nspp <- length(x$tip.label)

# convert to list
y <- unclass(x)

# add padding to list elements
padd <- 3
y$alive <- c(rep(TRUE, length(y$tip.label)), rep(FALSE, padd))
y$alive[y$tip.label == "D"] <- FALSE # D is extinct
y$edge <- rbind(y$edge, matrix(-1, nrow = padd, ncol = 2))
y$edge.length <- c(y$edge.length, rep(-1, padd))
y$tip.label <- c(y$tip.label, rep("", padd))
y$n <- nspp

# run the function to be tested
inew <- 3
s <- 0.5
xnew <- roleR:::testUpdatePhylo(y, inew, s)
xnew <- as(xnew, "phylo")
plot(xnew, x.lim = c(0, 6))
axis(1)

# new sp should be called "s" + "nspp + 1"
xnew$tip.label[nspp + 1] == paste0("s", nspp + 1)

# new sp should be attached to tip with index inew
newNode <- xnew$edge[xnew$edge[, 2] == inew, 1]
newSisterIDs <- xnew$edge[, 1] == newNode
all(c(inew, nspp + 1) %in% xnew$edge[newSisterIDs, 2])

# new sp branch length should be 1 * s
all(xnew$edge.length[newSisterIDs] == 1 * s)

# removing new sp should result in same tree topology
xold <- drop.tip(xnew, paste0("s", nspp + 1))
all(x$edge == xold$edge)

# tip lengths of old, but still extant, spp should be 1 * scale longer
# NOTE: this test as coded is only reliable if topologies are the same
xold$alive <- y$alive[1:nspp]
x$alive <- y$alive[1:nspp]
all(xold$edge.length[xold$edge[, 2] %in% (1:nspp)[xold$alive]] - 
              x$edge.length[x$edge[, 2] %in% (1:nspp)[x$alive]] == 1 * s)



# branch lengths of extinct spp should be unaffected 
# NOTE: this test as coded is only reliable if topologies are the same
all(xold$edge.length[xold$edge[, 2] %in% (1:nspp)[!xold$alive]] == 
        x$edge.length[x$edge[, 2] %in% (1:nspp)[!x$alive]])

# internal branches should be unaffected
# NOTE: this test as coded is only reliable if topologies are the same
all(xold$edge.length[xold$edge[, 2] > nspp] - x$edge.length[x$edge[, 2] > nspp] 
    == 0)

