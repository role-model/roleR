#' @title Simulate genetic variation
#' 
#' @description Calls msprime API to generate genetic data
#' 
#' @param localPop
#' @param metaPop
#' @param localSamp number of local individuals
#' @param mutationRate
#' @param immRate
#' @param maxTime
#' @param bp number of basepairs
#' @param ploidy
#' @param msp an R object wrapping msprime module (as returned by reticualte::import)


simGen <- function(localPop, metaPop, localSamp, 
                   mutationRate, immRate, maxTime, 
                   bp, ploidy, 
                   msp) {
    if(missing(msp)) {
        msprime <- reticulate::import('msprime')
    } else {
        msprime <- msp
    }
    
    # set-up demography object ----
    demography <- msprime$Demography()
    demography$add_population(name = 'loc', initial_size = localPop, 
                              growth_rate = 0)
    demography$add_population(name = 'meta', initial_size = metaPop, 
                              growth_rate = 0)
    
    # loc is source because we're backward in time
    demography$set_migration_rate(source = 'loc', dest = 'meta', rate = immRate)
    
    # set-up sample object ----
    samp <- msprime$SampleSet(localSamp, population = 'loc')
    
    ts <- msprime$sim_ancestry(samples=list(samp), demography = demography, 
                               sequence_length = bp, ploidy = ploidy, 
                               end_time = maxTime)
    
    # simulate mutations ----
    x <- msprime$sim_mutations(ts, rate = mutationRate)
    
    
    return(list(div = x$diversity() / bp,
                seq = processFasta(x$as_fasta()), 
                raw = x))
    
}

# 
# processFasta <- function(x) {
#     x <- gsub('[0-9]|\n', '', x)
#     
#     return(strsplit(x, '>n')[[1]][-1])
# }
# 
# 
# foo <- simGen(localPop = 1000, metaPop = 10000, localSamp = 5, 
#               mutationRate = 1, immRate = 1, maxTime = 10000, 
#               bp = 500, ploidy = 1)
# 
# cat(foo$raw$draw_text())
# foo$div
