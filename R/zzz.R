Rcpp::loadModule("iterModelCpp", TRUE)
#library(reticulate)
out <- tryCatch(install_miniconda(),error=function(e) {})
out <- tryCatch(conda_install("r-reticulate",
              packages="msprime",
              channel="conda-forge"),error=function(e) {})
#install_miniconda(force=TRUE,update=FALSE)
