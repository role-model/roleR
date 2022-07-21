
prepRebuild <- function()
{
  # uninstall previous version
  remove.packages("roleR")
  
  #fix env 
  Sys.setenv(R_INSTALL_STAGED = FALSE)
  
  # restart R 
  invisible(.rs.restartR())
  
  # delete version in R library
  unlink("C:/OWNER/hiest/AppData/local/R/win-library/4.2/roleR", recursive = TRUE)
  utils::browseURL("C:/OWNER/hiest/AppData/local/R/win-library/4.2/roleR")
}

prepRebuild()
