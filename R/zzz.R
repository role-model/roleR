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
