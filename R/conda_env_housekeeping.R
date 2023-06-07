.env_present <- function() {
    tryCatch({ROLER_PYTHON_ENV %in% reticulate::conda_list()$name}, 
             error = function(cond) FALSE
    )
}

#' Activate roleR own dedicated Python environment
#'
#' This function attempts to activate a dedicated roleR Miniconda Python
#' environment previously set up via \code{setup_env}.
#'
#' @param quiet Should informative messages be printed to the console? Default
#'   is \code{FALSE}.
#'
#' @return No return value, called for side effects
#'
#' @export

init_env <- function(quiet = FALSE) {
    if (!.env_present()) {
        stop("Could not activate roleR's Python environment because it is ",
             "not present on your system \n",
             "To set up a dedicated Python environment you first need to run ", 
             "`roleR::setup_env()`",
             call. = FALSE)
    } else {
        reticulate::use_condaenv(ROLER_PYTHON_ENV, required = TRUE)
        
        # this is an awful workaround around the reticulate/Python bug which prevents
        # import_from_path (see zzz.R) from working properly -- I'm getting nonsensical
        #   Error in py_call_impl(callable, dots$args, dots$keywords) :
        #     TypeError: integer argument expected, got float
        # in places with no integer/float conversion in sight
        #
        # at least it prevents having to do things like:
        # reticulate::py_run_string("def get_pedigree_ids(ts): return [ind.metadata['pedigree_id']
        #                                                              for ind in ts.individuals()]")
        # (moved from ts_load() here because this is a better place for loading our Python functions)
        reticulate::source_python(file = system.file("pylib/pylib.py", package = "slendr"))
        
        if (!reticulate::py_module_available("msprime") ||
            !reticulate::py_module_available("tskit") ||
            !reticulate::py_module_available("pyslim")) {
            stop("Python environment ", PYTHON_ENV, " has been found but it",
                 " does not appear to have msprime, tskit and pyslim modules all",
                 " installed. Perhaps the environment got corrupted somehow?",
                 " Running `clear_env()` and `setup_env()` to reset the slendr's Python",
                 " environment is recommended.", call. = FALSE)
        } else {
            if (!quiet)
                message("The interface to all required Python modules has been activated.")
        }
    }
}

#' Setup a dedicated Python virtual environment for slendr
#'
#' This function will automatically download a Python miniconda distribution
#' dedicated to an R-Python interface. It will also create a slendr-specific
#' Python environment with all the required Python dependencies.
#'
#' @param quiet Should informative messages be printed to the console? Default
#'   is \code{FALSE}.
#' @param agree Automatically agree to all questions?
#' @param pip Should pip be used instead of conda for installing slendr's Python
#'   dependencies? Note that this will still use the conda distribution to
#'   install Python itself, but will change the repository from which slendr
#'   will install its Python dependencies. Unless explicitly set to \code{TRUE},
#'   Python dependencies will be installed from conda repositories by default,
#'   expect for the case of osx-arm64 Mac architecture, for which conda
#'   dependencies are broken.
#'
#' @return No return value, called for side effects
#'
#' @export
setup_env <- function(quiet = FALSE, agree = FALSE, pip = NULL) {
    if (.env_present()) {
        message("A required slendr Python environment is already present. You can activate\n",
                "it by calling init_env().")
    } else {
        if (agree)
            answer <- 2
        else
            answer <- utils::menu(
                c("No", "Yes"),
                title = paste0(
                    "This function will install a completely isolated Miniconda Python distribution\n",
                    "just for slendr and create an environment with all required Python modules.\n",
                    "\nEverything will be installed into a completely separate location into an\n",
                    "isolated environment in an R library directory. This won't affect your other\n",
                    "Python installations at all. You can always wipe out the automatically created\n",
                    "environment by running clear_env().\n\n",
                    "Do you wish to proceed with the automated Python environment setup?")
            )
        if (answer == 2) {
            message("=======================================================================")
            message("Installing slendr's Python environment. Please wait until")
            message("the installation procedure finishes. Do NOT interrupt the")
            message("process while the installation is still running.")
            message("======================================================================\n")
            Sys.sleep(10)
            
            if (!dir.exists(reticulate::miniconda_path()))
                reticulate::install_miniconda()
            
            # parse the Python env name back to the list of dependencies
            # (the environment is defined in .onAttach(), and this makes sure the
            # dependencies are defined all in one place)
            versions <- PYTHON_ENV %>% gsub("-", "==", .) %>% strsplit("_") %>% .[[1]]
            python_version <- gsub("Python==", "", versions[1])
            package_versions <- c(versions[-1], "pandas")
            
            reticulate::conda_create(envname = PYTHON_ENV, python_version = python_version)
            reticulate::use_condaenv(PYTHON_ENV, required = TRUE)
            
            # msprime/tskit conda dependency is broken on M1 Mac architecture, fallback
            # to pip in cases like this (otherwise use conda to avoid any potential
            # compilation issues such as missing libgsl)
            if (is.null(pip))
                pip <- all(Sys.info()[c("sysname", "machine")] == c("Darwin", "arm64"))
            
            reticulate::conda_install(envname = PYTHON_ENV, packages = package_versions, pip = pip)
            
            if (!quiet) {
                message("======================================================================")
                message("Python environment for slendr has been successfuly created, and ",
                        "the R\ninterface to msprime, tskit, and pyslim modules has been activated.\n")
                message("In future sessions, activate this environment by calling init_env().")
                message("=======================================================================")
            }
        } else
            warning("Your Python environment is not set up correctly which means that the tree\n",
                    "sequence functionality of slendr will not work.", call. = FALSE)
    }
}

#' Remove the automatically created slendr Python environment
#'
#' @param force Ask before deleting the environment?
#'
#' @return No return value, called for side effects
#'
#' @export
clear_env <- function(force = FALSE) {
    if (.env_present()) {
        path <- reticulate::conda_list() %>%
            dplyr::filter(grepl(PYTHON_ENV, name)) %>%
            { gsub("bin\\/python", "", .$python) }
        
        answer <- utils::menu(
            c("No", "Yes"),
            title = paste0(
                "Are you sure you want to delete the automatically created slendr ",
                "Python\nenvironment? It is located in:\n\n", path, "\n\n",
                "If you remove it, you can create it again by running",
                "`setup_env()`\nwithout any arguments in a new R session."
            )
        )
        if (answer == 2) reticulate::conda_remove(PYTHON_ENV)
        message("The slendr Python environment has been sucessfully removed.")
    } else
        warning("No automatic slendr Python environment has been found so there is\n",
                "nothing to delete.", call. = FALSE)
}

#' Check that the active Python environment is setup for slendr
#'
#' This function inspects the Python environment which has been activated by the
#' reticulate package and prints the versions of all slendr Python dependencies
#' to the console.
#'
#' @param verbose Should a log message be printed? If \code{FALSE}, only a logical
#'   value is returned (invisibly).
#'
#' @return Either \code{TRUE} (slendr Python environment is present) or \code{FALSE}
#'   (slendr Python environment is not present).
#'
#' @examples
#' \dontshow{check_dependencies(python = TRUE) # make sure dependencies are present
#' }
#' init_env()
#' check_env()
#' @export
check_env <- function(verbose = TRUE) {
    # if there is no Python available on user's system, don't immediately
    # jump to installing miniconda (let's deal with that in setup_env())
    orig_env <- Sys.getenv("RETICULATE_MINICONDA_ENABLED")
    Sys.setenv(RETICULATE_MINICONDA_ENABLED = FALSE)
    on.exit(Sys.setenv(RETICULATE_MINICONDA_ENABLED = orig_env))
    
    py <- reticulate::py_discover_config()
    
    has_tskit <- reticulate::py_module_available("tskit")
    has_msprime <- reticulate::py_module_available("msprime")
    has_pyslim <- reticulate::py_module_available("pyslim")
    # has_pylib <- !is.null(pylib)
    
    if (has_tskit)
        tskit_version <- paste("version", tskit[["_version"]]$tskit_version, "\u2713")
    else
        tskit_version <- "MISSING \u274C"
    
    if (has_msprime)
        msprime_version <- paste("version", msp[["_version"]]$version, "\u2713")
    else
        msprime_version <- "MISSING \u274C"
    
    if (has_pyslim)
        pyslim_version <- paste("version", pyslim$pyslim_version, "\u2713")
    else
        pyslim_version <- "MISSING \u274C"
    
    # if (has_pylib)
    #   pylib_status <- "successfully loaded \u2713"
    # else
    #   pylib_status <- "NOT LOADED \u274C"
    
    if (verbose) {
        cat("Summary of the currently active Python environment:\n\n")
        cat("Python binary:", py$python, "\n")
        cat("Python version:", py$version_string, "\n")
        
        cat("\nslendr requirements:\n")
        cat(" - tskit:", tskit_version, "\n")
        cat(" - msprime:", msprime_version, "\n")
        cat(" - pyslim:", pyslim_version, "\n")
        # cat(" - slendr module:", pylib_status, "\n")
    }
    
    if (!all(c(has_tskit, has_pyslim, has_msprime))) {
        return_value <- FALSE
        if (verbose)
            cat("\nNote that due to the technical limitations of embedded Python,",
                "if you\nwant to switch to another Python environment you will need",
                "to restart\nyour R session first.\n")
        # reference: https://github.com/rstudio/reticulate/issues/27#issuecomment-512256949
    } else
        return_value <- TRUE
    
    if (verbose)
        return(invisible(return_value))
    else
        return(return_value)
}
