# Package-level environment that caches the sourced Python module so that
# role_msprime.py is parsed only once per R session rather than on every
# runRole() call.
.role_py_env <- new.env(parent = emptyenv())


# Return the sim_role Python function, initialising the Python environment on
# the first call.  Emits a clear error pointing to installMsprime() if any
# required package is missing.
.get_sim_role <- function() {
    if (!exists("sim_role", envir = .role_py_env, inherits = FALSE)) {
        required <- c("msprime", "tskit", "newick", "numpy", "pandas")
        missing  <- required[!vapply(required, reticulate::py_module_available,
                                     logical(1))]
        if (length(missing) > 0) {
            stop(
                "Required Python package(s) not available: ",
                paste(missing, collapse = ", "), ".\n",
                "Run installMsprime() to install all Python dependencies, ",
                "then restart R.",
                call. = FALSE
            )
        }
        pyfi <- system.file("python", "role_msprime.py", package = "roleR")
        reticulate::source_python(pyfi, envir = .role_py_env)
    }
    .role_py_env$sim_role
}


.onLoad <- function(libname, pkgname) {
    # Python is not initialised here — that happens lazily on first use.
    # We only register the reticulate delay-load configuration so that
    # reticulate knows which environment to prefer when it does start.
    # Users can override this by calling reticulate::use_condaenv() before
    # runRole().
    options(roleR.python_initialized = FALSE)
}
