#' @param x a roleModel object

getTimeseries <- function(x, FUN, ...) {
  o <- lapply(x, FUN, ...)
  
  return(do.call(rbind, o))
}

l <- replicate(n = 5, rpois(10, 1), simplify = FALSE)

getAbunds <- function(x) {round(x)}


getTimeseries(l, getAbunds)


getTimeseries(l, function(x) mean(x))





tres <- replicate(5, ape::rphylo(10, 1, 0.8), simplify = FALSE)

df <- data.frame(iter = 1:5, fullPhylo = I(tres))

df$fullPhylo[1]
