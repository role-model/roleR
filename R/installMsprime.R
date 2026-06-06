#' @title Install Python dependencies for roleR
#'
#' @description Installs the Python packages required by roleR (\code{msprime},
#'   \code{tskit}, \code{newick}, \code{numpy}, \code{pandas}) into a conda or
#'   virtual environment managed by reticulate.  Run this once after installing
#'   the package; then restart R before calling \code{runRole()}.
#'
#' @param envname Name of the conda or virtual environment to install into.
#'   Defaults to \code{"r-roleR"}.
#' @param method Installation method passed to
#'   \code{\link[reticulate]{py_install}}: \code{"conda"} (default) or
#'   \code{"virtualenv"}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[reticulate]{py_install}}.
#'
#' @return Invisible \code{NULL}.
#'
#' @examples
#' \dontrun{
#' installMsprime()
#' }
#'
#' @export
installMsprime <- function(envname = "r-roleR", method = "conda", ...) {
    reqs_file <- system.file("python", "requirements.txt", package = "roleR")
    pkgs <- readLines(reqs_file)

    message("Installing Python dependencies into environment '", envname, "'...")
    reticulate::py_install(pkgs, envname = envname, method = method, ...)

    message(
        "Done. Restart R, then activate the environment with:\n",
        "  reticulate::use_condaenv(\"", envname, "\")\n",
        "before calling runRole()."
    )
    invisible(NULL)
}
