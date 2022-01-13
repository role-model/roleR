# uninstall previous version
remove.packages("roleR")
#fix env 
Sys.setenv(R_INSTALL_STAGED = FALSE)
# restart R 
invisible(.rs.restartR())

Rcpp::compileAttributes()

# below doesn't work, just go Build -> Clean & Rebuild 
system("Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source ...roleR")

# rebuild documentation
devtools::document() 
library(roxygen2)
roxygen2::roxygenise()
