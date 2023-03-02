# fun for Isaac 
newick <- function(phylo){
    ape_tree <- as(model@modelSteps[[1]]@phylo,"phylo")
    write.tree(ape_tree, file = 'temp.nwk', append = FALSE, digits = 10, tree.names = FALSE)
    nwk <- readChar('temp.nwk', file.info('temp.nwk')$size)
    unlink('temp.nwk')
    return(nwk)
}
# example
nwk <- newick(model@modelSteps[[1]]@phylo)
nwk
# note that this contains "\r\n" at the end