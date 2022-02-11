
#' @description returns a dataframe containing metadata extracted from sim, a roleModelCpp object
#' @param sim a roleModelCpp object representing a simulation 
#' @export
#' 

extractMetadata <- function(sim)
{
  # get param values from sim 
  values <- c(sim$params$values$dispersal_prob, sim$params$values$extinction_meta, 
              sim$params$values$individuals_local, sim$params$values$individuals_meta,
              sim$params$values$speciation_local, sim$params$values$speciation_meta,
              sim$params$values$species_meta, sim$params$values$trait_sigma,
              1000,"",R.Version()$version.string, packageVersion("roleR")) 
  values <- as.character(values)
  View(values)
  names <- c("dispersal prob","metacomm extinction prob",
             "num individuals local", "num individuals meta", 
             "localcomm speciation prob", "metacomm speciation prob",
             "num species meta prob", "trait sigma param","num iterations", 
             "Project Description", "R version", "roleR version")
  print(length(names))
  print(length(values))
  out <- data.frame(values)
  colnames(out) <- "Value"
  View(out)
  rownames(out) <- names
  return(out)
}

#' @description extracts a paramValuesCpp object from a metadata object - use to
#' easily transfer parameters from saved metadata to a new model simulation 
#' @param metadata a dataframe containing metadata extracted from a model using 
#' extractMetadata - TODO make this a strict object later
#' @export
#' 
# extracts a paramValuesCpp object from a metadata object
extractParams <- function(metadata) 
{
  out <- new(paramValuesCpp) 
  m <- as.numeric(metadata[,1])
  # ultimately these need to either be distributions or doubles
  # probably convert params to store as NumericVectors of distributions in paramValuesCpp 
  
  out$dispersal_prob = m[1]
  out$extinction_meta = m[2]
  out$individuals_local = m[3]
  out$individuals_meta = m[4]
  out$speciation_local = m[5]
  out$speciation_meta = m[6]
  out$species_meta = m[7]
  out$trait_sigma = m[8]
  
  return(out) 
}

# notes
# allow saving expression, text, executable code to parse to distributions of params etc
# how to capture priors in metadata
# optional field - user specified description of project 
# keyword field - based on description 
# version information R, roleR 