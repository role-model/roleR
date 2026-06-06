#' @title Update elements of an initialized roleModel
#'
#' @description These functions update the `localComm`, `metaComm`, or `phylo`
#'   element of the initial state of a `roleModel` that has not yet been run.
#'   Each function validates that the replacement is consistent with the model
#'   parameters before writing the change.
#'
#'   When `slot = NULL` (default), the entire element is replaced with `value`,
#'   which must be an object of the corresponding S4 class. When `slot` is a
#'   slot name (character), only that slot is updated and `value` should be the
#'   appropriate replacement object.
#'
#' @param model An initialized (pre-run) `roleModel` object.
#' @param value Replacement: a full S4 object when `slot = NULL`, or a vector /
#'   matrix matching the target slot when `slot` is specified.
#' @param slot Character name of the slot to update, or `NULL` to replace the
#'   entire element.
#' @return The updated `roleModel` object.
#'
#' @examples
#' m <- roleModel(roleParams(species_meta = 2, individuals_local = 10))
#'
#' # update a single slot
#' m <- updateMetaComm(m, c(0.3, 0.7), slot = "spAbund")
#' m <- updateLocalComm(m, rep(1L, 10), slot = "indSpecies")
#'
#' @rdname updateModel
#' @export
updateLocalComm <- function(model, value, slot = NULL) {
    assertPreRun(model)
    sm <- model@params@species_meta
    J  <- model@params@individuals_local(1)

    if (is.null(slot)) {
        validateLocalComm(value, sm, J)
        model@modelSteps[[1]]@localComm <- value
    } else {
        validateLocalCommSlot(slot, value, sm, J)
        methods::slot(model@modelSteps[[1]]@localComm, slot) <- value
    }

    model
}


#' @rdname updateModel
#' @export
updateMetaComm <- function(model, value, slot = NULL) {
    assertPreRun(model)
    sm <- model@params@species_meta

    if (is.null(slot)) {
        validateMetaComm(value, sm)
        model@modelSteps[[1]]@metaComm <- value
    } else {
        validateMetaCommSlot(slot, value, sm)
        methods::slot(model@modelSteps[[1]]@metaComm, slot) <- value
    }

    model
}


#' @rdname updateModel
#' @export
updateRolePhylo <- function(model, value, slot = NULL) {
    assertPreRun(model)

    if (is.null(slot)) {
        sm <- model@params@species_meta
        validateRolePhylo(value, sm)
        model@modelSteps[[1]]@phylo <- value
    } else {
        validateRolePhyloSlot(slot, value, model@modelSteps[[1]]@phylo)
        methods::slot(model@modelSteps[[1]]@phylo, slot) <- value
    }

    model
}


# check to make sure model not run
assertPreRun <- function(model) {
    n <- length(model@modelSteps)
    if (n != 1)
        stop(
            "update functions only work on initialized (pre-run) roleModel objects; ",
            "this model has ", n, " model steps",
            call. = FALSE
        )
    invisible(NULL)
}

# localComm slot-level validator
validateLocalCommSlot <- function(slot, value, sm, J) {
    ind_slots <- c("indSpecies", "indTrait", "indSeqs")
    spp_slots <- c("spAbund", "spAbundHarmMean", "spPropHarmMean",
                   "spLastOriginStep", "spLastOriginGen", "spGenDiv", "spTajD")
    all_slots <- c(ind_slots, spp_slots)

    if (!slot %in% all_slots)
        stop("unknown localComm slot '", slot, "'; valid slots: ",
             paste(all_slots, collapse = ", "), call. = FALSE)

    if (slot == "indSpecies") {
        if (length(value) != J)
            stop("indSpecies must have length individuals_local (", J, "), ",
                 "got ", length(value), call. = FALSE)
        bad <- unique(value[value < 1 | value > sm])
        if (length(bad) > 0)
            stop("indSpecies contains IDs (", paste(bad, collapse = ", "),
                 ") outside [1, species_meta=", sm, "]", call. = FALSE)

    } else if (slot == "indTrait") {
        if (!is.matrix(value) || nrow(value) != J)
            stop("indTrait must be a matrix with nrow = individuals_local (", J,
                 "), got nrow = ", if (is.matrix(value)) nrow(value) else "non-matrix",
                 call. = FALSE)

    } else if (slot == "indSeqs") {
        if (!is.character(value) || length(value) != J)
            stop("indSeqs must be a character vector of length individuals_local (",
                 J, "), got length ", length(value), call. = FALSE)

    } else if (slot %in% spp_slots) {
        if (length(value) != sm)
            stop("'", slot, "' must have length species_meta (", sm, "), ",
                 "got ", length(value), call. = FALSE)
    }

    invisible(NULL)
}

# localComm object-level validator
validateLocalComm <- function(lc, sm, J) {
    if (!is(lc, "localComm"))
        stop("value must be a localComm object", call. = FALSE)
    validateLocalCommSlot("indSpecies", lc@indSpecies, sm, J)
    validateLocalCommSlot("indTrait",   lc@indTrait,   sm, J)
    invisible(NULL)
}

# metaComm slot-level validator
validateMetaCommSlot <- function(slot, value, sm) {
    valid_slots <- c("spAbund", "spTrait")
    if (!slot %in% valid_slots)
        stop("unknown metaComm slot '", slot, "'; valid slots: ",
             paste(valid_slots, collapse = ", "), call. = FALSE)

    if (slot == "spAbund") {
        if (length(value) != sm)
            stop("spAbund must have length species_meta (", sm, "), ",
                 "got ", length(value), call. = FALSE)
        if (any(value < 0))
            stop("spAbund values must be non-negative", call. = FALSE)
        if (abs(sum(value) - 1) > 0.01)
            stop("spAbund must sum to 1 (got ", round(sum(value), 4), ")",
                 call. = FALSE)

    } else if (slot == "spTrait") {
        if (!is.matrix(value) || nrow(value) != sm)
            stop("spTrait must be a matrix with nrow = species_meta (", sm, "), ",
                 "got nrow = ", if (is.matrix(value)) nrow(value) else "non-matrix",
                 call. = FALSE)
    }

    invisible(NULL)
}

# metaComm object-level validator
validateMetaComm <- function(mc, sm) {
    if (!is(mc, "metaComm"))
        stop("value must be a metaComm object", call. = FALSE)
    validateMetaCommSlot("spAbund", mc@spAbund, sm)
    validateMetaCommSlot("spTrait", mc@spTrait, sm)
    invisible(NULL)
}

# rolePhylo slot-level validator
validateRolePhyloSlot <- function(slot, value, current_phylo) {
    valid_slots <- c("n", "e", "l", "alive", "tipNames")
    if (!slot %in% valid_slots)
        stop("unknown rolePhylo slot '", slot, "'; valid slots: ",
             paste(valid_slots, collapse = ", "), call. = FALSE)

    if (slot == "e") {
        if (!is.matrix(value) || ncol(value) != 2)
            stop("e must be a 2-column integer matrix", call. = FALSE)
    } else if (slot == "l") {
        ne <- nrow(current_phylo@e)
        if (length(value) != ne)
            stop("l must have length nrow(e) (", ne, "), got ", length(value),
                 call. = FALSE)
    } else if (slot == "alive") {
        if (!is.logical(value))
            stop("alive must be a logical vector", call. = FALSE)
    } else if (slot == "tipNames") {
        if (!is.character(value))
            stop("tipNames must be a character vector", call. = FALSE)
    }

    invisible(NULL)
}

# rolePhylo object-level validator
validateRolePhylo <- function(phy, sm) {
    if (!is(phy, "rolePhylo"))
        stop("value must be a rolePhylo object", call. = FALSE)
    if (phy@n != sm)
        stop("rolePhylo@n (", phy@n, ") must equal species_meta (", sm, ")",
             call. = FALSE)
    invisible(NULL)
}
