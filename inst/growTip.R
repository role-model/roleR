library(ape)

x <- rphylo(6, 1, 0.8)
x$edge.length <- x$edge.length * 10

plot(x)
axis(1)
nodelabels()

n <- Ntip(x)

i <- 6

e <- x$edge
j <- which(e[, 2] == i)

e[e > n] <- e[e > n] + 1

newNode <- 2 * n + 1
e <- rbind(e, matrix(c(newNode, newNode,
                       e[j, 2], n + 1),
                     ncol = 2))

e[j, 2] <- newNode

l <- x$edge.length
l <- c(l, 0, 0)
l[e[, 2] <= n + 1] <- l[e[, 2] <= n + 1] + 1

y <- list(edge = e, edge.length = l,
          tip.label = paste0('p', e[e[, 2] <= n + 1, 2]),
          Nnode = n)

class(y) <- 'phylo'

plot(y)
nodelabels()
axis(1)

plot(x)
