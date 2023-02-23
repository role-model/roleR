#library(ape)

#apePhylo <- read.tree(text="(A:2,(B:1,C:1):1);")
#plot(apePhylo)
#tre <- as(apePhylo,'rolePhylo')

test_that("update_speciation_phylo adds a new internal node", {
})

test_that("update_speciation_phylo adds new tips", {
})

test_that("update_speciation_phylo updates ancestry of internal nodes", {
})

test_that("update_speciation_phylo sets new edge length and augments old edge lengths", {
})

test_that("update_speciation_phylo DOESN'T update internal edges", {
})

test_that("update_speciation_phylo DOESN'T update extinct tips", {
})

test_that("update_speciation_phylo results in a phylo where the branches of all alive species lead to the present", {
})

test_that("update_speciation_phylo updates the alive vector", {
})

test_that("update_speciation_phylo increments the number of tips", {
})

test_that("update_speciation_phylo creates a rolePhylo object that can be coerced to an ape object", {
})

test_that("update_speciation_phylo creates the expected phylogeny given a trivial 3 species tree and deterministic inputs", {
})