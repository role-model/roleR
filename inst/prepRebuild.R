# uninstall previous version
remove.packages("roleR")
#fix env 
Sys.setenv(R_INSTALL_STAGED = FALSE)
# restart R 
invisible(.rs.restartR())

Rcpp::compileAttributes()
system("Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source ...roleR")
