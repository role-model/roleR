#' qlseries
#'
#' @param p p
#' @param beta  beta
#' @param lower.tail lt
#' @param log l
#'
#' @return smthg
#' @export
#'
#' @importFrom extraDistr qlgser
qlseries <- function (p, beta, lower.tail = TRUE, log = FALSE) 
{
    extraDistr::qlgser(p, theta = exp(-beta), lower.tail = lower.tail, 
                       log.p = log)
}

#' 
#' #' sad from s and n
#' #'
#' #' @param S s
#' #' @param N m
#' #'
#' #' @return ard
#' #' @export
#' #'
#' sad_from_SN <- function(S, N) {
#'     
#'     qfun <- function(p, beta) {
#'         
#'         qlseries(p, beta)
#'     }
#'     
#'     # solve for alpha paramter
#'     asol <- uniroot(interval = c(.Machine$double.eps^0.25,
#'                                  .Machine$integer.max),
#'                     f = function(a) {
#'                         a * log(1 + N / a) - S
#'                     })
#'     
#'     # calculate p parameter and beta (as used by pika)
#'     p <- 1 - exp(-S / asol$root)
#'     beta <- -log(p)
#'     
#'     
#'     rank <- qfun(seq(1, 1/S, length = S) - 1/(2 * S), beta = beta)
#'     
#'     return(rank / sum(rank))
#' }