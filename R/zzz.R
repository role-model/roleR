# # Rcpp::loadModule("iterModelCpp", TRUE)
# #library(reticulate)
# out <- tryCatch(install_miniconda(),error=function(e) {})
# out <- tryCatch(conda_install("r-reticulate",
#               packages="msprime",
#               channel="conda-forge"),error=function(e) {})
# #install_miniconda(force=TRUE,update=FALSE)

# .onLoad <- function(libname, pkgname) {
#     packageStartupMessage('linking to python dependencies...')
#     
#     browser()
#     
#     reticulate::use_condaenv("r-reticulate")
#     pypack <- reticulate::py_list_packages()
#     pyreqs <- c('msprime', 'newick')
#     
#     booboo <- pyreqs[!(pyreqs %in% pypack$package)]
#     
#     if(length(booboo) > 0) {
#         s <- paste('missing python dependencies:', 
#                    paste(booboo, collapse = ', '))
#         packageStartupMessage(s)
#     } else {
#         packageStartupMessage('linking successful')
#     }
#     
# }
# 



# # global references to required Python packages - inspired by:
# # https://cran.r-project.org/web/packages/reticulate/vignettes/package.html
# # (Python environment initialization is now being done using init_env())
# tskit <- NULL
# pyslim <- NULL
# msp <- NULL
# pylib <- NULL

# define roleR's required Python dependencies and compose an environment name
# that will be used specifically for them
ROLER_PYTHON_ENV <- paste(c("Python", "msprime", "newick", "collections"), 
                          collapse = "_")
