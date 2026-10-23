#' @title Get gene sequences for the individuals in a `roleModel`
#'
#' @description msprime does not simulate sequences in roleR. Instead the 
#'     simulation returns the genotype matrix of the individuals (the 
#'     `genoMat` slot of the local community), and sequences are built from 
#'     it on demand by this function.
#'
#' @param x A `roleModel` object that has been run.
#' @param snapshot Integer index of the desired snapshot. Use `Inf` (default)
#'     for the last saved snapshot.
#' @param seed Integer seed used to make the sequences reproducible, or `NULL`
#'     to draw fresh random sequences. The global random number generator state
#'     is restored after the call, so setting a `seed` does not change 
#'     downstream random draws.
#'
#' @details The genotype matrix holds allele indices (0 is the ancestral allele
#'     of the site) but not the nucleotides or the positions of the sites along
#'     the sequence. Under the default mutation model (JC69) nucleotides are 
#'     interchangeable, so this function assigns each allele index a distinct 
#'     nucleotide at random (independently for each site), places the variable 
#'     sites at random positions along the sequence, and fills the remaining 
#'     positions with random nucleotides. The sequences have the same 
#'     distribution as sequences simulated by msprime, and all pairwise 
#'     differences, diversity statistics, site frequency spectra and trees 
#'     computed from them are identical to those computed from the genotype
#'     matrix. They are not the same strings msprime would have produced. 
#'     Because the sites (rows of the genotype matrix) are the same at every 
#'     snapshot, the same `seed` gives the same positions, nucleotides, and 
#'     filler at every snapshot. 
#'     
#'     Memory use is proportional to the number of individuals times the 
#'     number of base pairs.
#'
#' @return A character vector with one sequence per individual, in the same 
#'     order as `indSpecies` in the local community of the requested snapshot.
#'
#' @export

rawSeqs <- function(x, snapshot = Inf, seed = 1) {
    if(!is(x, "roleModel")) {
        stop("`x` must be a roleModel object", call. = FALSE)
    }
    
    i <- .resolveSnapshot(x@modelSteps, snapshot)
    lc <- x@modelSteps[[i]]@localComm
    G <- lc@genoMat
    
    if(ncol(G) != length(lc@indSpecies)) {
        stop("this snapshot has no genotype matrix; has the model been run ",
             "with `runRole`?", call. = FALSE)
    }
    
    synthSeqs(G, L = x@params@num_basepairs, seed = seed)
}


# make synthetic sequences from a genotype matrix `G` (sites x individuals, 
# entries are allele indices 0-3); positions of sites, nucleotides of each
# allele, and the nucleotides at non-variable positions are drawn at random
# `L` is the number of base pairs
# draws happen in a fixed order that depends only on `L` and the number of
# sites so that the same `seed` gives the same site table for any `G` with the
# same rows
synthSeqs <- function(G, L, seed = NULL) {
    if(!is.null(seed)) {
        hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        if(hadSeed) {
            oldSeed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        }
        on.exit({
            if(hadSeed) {
                assign(".Random.seed", oldSeed, envir = globalenv())
            } else if(exists(".Random.seed", envir = globalenv(), 
                             inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        }, add = TRUE)
        set.seed(seed)
    }
    
    bases <- c("A", "C", "G", "T")
    S <- nrow(G)
    
    pos <- sort(sample.int(L, S))
    alleleLetters <- t(vapply(seq_len(S), function(i) sample(bases), 
                              character(4)))
    ref <- sample(bases, L, replace = TRUE)
    
    return(seqsFromG(G, pos = pos, alleleLetters = alleleLetters, ref = ref))
}


# build sequences from a genotype matrix and a site table
# `G` sites x individuals matrix (dense or sparse) of allele indices (0-3)
# `pos` position (1-based) along the sequence of each of the S sites
# `alleleLetters` S x 4 character matrix, column k + 1 is the nucleotide for
#     allele k at that site
# `ref` character vector of length L, nucleotides at non-variable positions
# returns a character vector with one sequence per column of `G`
seqsFromG <- function(G, pos, alleleLetters, ref) {
    G <- as.matrix(G)
    S <- nrow(G)
    n <- ncol(G)
    
    out <- matrix(ref, nrow = n, ncol = length(ref), byrow = TRUE)
    
    if(S > 0) {
        lt <- alleleLetters[cbind(rep(seq_len(S), times = n), 
                                  as.vector(G) + 1)]
        out[, pos] <- t(matrix(lt, nrow = S))
    }
    
    return(apply(out, 1, paste, collapse = ""))
}
