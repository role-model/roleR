#' @title Getter functions for `roleModel` objects
#'
#' @description Extract components from a `roleModel` at a given snapshot.
#'
#' @param x A `roleModel` object.
#' @param snapshot Integer index of the desired snapshot. Use `Inf` (default)
#'   for the last saved snapshot.
#'
#' @details `getMetaComm` takes only the `roleModel` argument because the
#'   metacommunity does not vary across snapshots; it returns the `metaComm`
#'   from the initial state. `getLocalComm`, `getPhylo`, and `getState` accept
#'   a `snapshot` argument. `getState` with `snapshot = Inf` is equivalent to
#'   the existing `getFinalState` for a `roleModel`.
#'
#' @include roleModel.R
#' @rdname getters
#' @export
setGeneric('getMetaComm',
           def = function(x) standardGeneric('getMetaComm'),
           signature = 'x')

#' @rdname getters
#' @export
setGeneric('getLocalComm',
           def = function(x, snapshot = Inf) standardGeneric('getLocalComm'),
           signature = 'x')

#' @rdname getters
#' @export
setGeneric('getPhylo',
           def = function(x, snapshot = Inf) standardGeneric('getPhylo'),
           signature = 'x')

#' @rdname getters
#' @export
setGeneric('getState',
           def = function(x, snapshot = Inf) standardGeneric('getState'),
           signature = 'x')


# Internal helper: resolve snapshot index to a valid list position.
# modelSteps may be pre-allocated with trailing NULLs; this finds the last
# non-NULL entry and clamps/maps the requested snapshot to it.
.resolveSnapshot <- function(steps, snapshot) {
    nValid <- sum(!vapply(steps, is.null, logical(1)))
    if (nValid == 0L) stop('no snapshots available in this roleModel', call. = FALSE)
    if (is.infinite(snapshot)) return(nValid)
    i <- as.integer(snapshot)
    if (i < 1L) stop('snapshot must be >= 1', call. = FALSE)
    if (i > nValid) stop('snapshot ', i, ' exceeds number of saved snapshots (', nValid, ')', call. = FALSE)
    i
}


#' @rdname getters
#' @aliases getMetaComm,roleModel-method
#' @export
setMethod('getMetaComm',
          signature(x = 'roleModel'),
          definition = function(x) {
              x@modelSteps[[1]]@metaComm
          }
)


#' @rdname getters
#' @aliases getLocalComm,roleModel-method
#' @export
setMethod('getLocalComm',
          signature(x = 'roleModel'),
          definition = function(x, snapshot = Inf) {
              i <- .resolveSnapshot(x@modelSteps, snapshot)
              x@modelSteps[[i]]@localComm
          }
)


#' @rdname getters
#' @aliases getPhylo,roleModel-method
#' @export
setMethod('getPhylo',
          signature(x = 'roleModel'),
          definition = function(x, snapshot = Inf) {
              i <- .resolveSnapshot(x@modelSteps, snapshot)
              x@modelSteps[[i]]@phylo
          }
)


#' @rdname getters
#' @aliases getState,roleModel-method
#' @export
setMethod('getState',
          signature(x = 'roleModel'),
          definition = function(x, snapshot = Inf) {
              i <- .resolveSnapshot(x@modelSteps, snapshot)
              x@modelSteps[[i]]
          }
)
