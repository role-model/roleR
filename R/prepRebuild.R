

prepRebuild <- function()
{
  # uninstall previous version
  remove.packages("roleR")
  
  #fix env 
  Sys.setenv(R_INSTALL_STAGED = FALSE)
  
  # restart R 
  invisible(.rs.restartR())
  
  # delete version in R library
  unlink("C:/Users/hiest/Documents/R/win-library/4.1/roleR", recursive = TRUE)
  utils::browseURL("C:/OWNER/hiest/Documents/R/win-library/4.1/roleR")
}

