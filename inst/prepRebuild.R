# uninstall previous version
remove.packages("roleR")

#fix env 
Sys.setenv(R_INSTALL_STAGED = FALSE)

# restart R 
invisible(.rs.restartR())

# delete version in R library
unlink("C:/Users/hiest/Documents/R/win-library/4.1/roleR", recursive = TRUE)
utils::browseURL("C:/Users/hiest/Documents/R/win-library/4.1/roleR")

# restart R 
invisible(.rs.restartR())

# rebuild documentation
devtools::document() 
library(roxygen2)
roxygen2::roxygenise()

install.packages("RcppDeepState")
