
# used ONLY on Jacob's computer before rebuilding package to resolve a bug
.prepRebuild <- function()
{
  # uninstall previous version
  remove.packages("roleR")
  
  #fix env 
  Sys.setenv(R_INSTALL_STAGED = FALSE)
  
  # restart R 
  invisible(.rs.restartR())
  
  # delete version in R library
  #unlink("C:/Users/hiest/Documents/R/win-library/4.1/roleR", recursive = TRUE)
  #utils::browseURL("C:/OWNER/hiest/Documents/R/win-library/4.1/roleR")
  unlink("C:/Users/hiest/AppData/Local/R/win-library/4.2/00LOCK-roleR", recursive = TRUE)
  utils::browseURL("C:/OWNER/hiest/Documents/R/win-library/4.1/roleR")
}
