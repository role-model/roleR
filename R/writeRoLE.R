#' @title Write a roleExperiment to disk
#' @description writes to a .roleExperiment file that can be read into R using readRDs
#' @param x roleExperiment object to write.
#' @param directory directory (path) to write the object to.
#' @param filename name of the file to save.
#' @param save_txt T/F
#' File will be called *filename*.roleExperiment. 
#' @rdname writeRole
#' @export

setGeneric('writeRole', 
           def = function(x, directory=NULL, filename=NULL, save_txt=TRUE) standardGeneric('writeRole'), 
           signature = 'x')

# method for roleData
#' @rdname writeRole
#' @export

setMethod('writeRole', 
          signature = 'roleExperiment', 
          definition = function(x, directory=NULL, filename=NULL, save_txt=TRUE) { 
              if(is.null(directory)){
                  directory = getwd()
              }
              
              extension <- ".roleExperiment"
              info_tag <- ".roleExperimentInfo"
              title_line <- "Metadata for RoLE experiment R object - load into R using readRDS(*file_location_and_name*.roleExperiment)"
              info_line <- paste("Contains", length(x@modelRuns), "model runs")
              author_line <- paste("Author:", x@context["author"])
              date_line <- paste("Date:", x@context["date"])
              desc_line <- paste("Description:", x@context["description"])
              
              if(is.null(filename)){
                  filename = "unnamed_experiment"
              }
              
              saveRDS(x,paste0(directory,"/",filename,extension))
              if(save_txt){
                  con <- file(paste0(directory, "/", filename, info_tag, ".txt"))
                  writeLines(c(title_line,info_line,author_line,date_line,desc_line), con)
                  close(con)
              }
          }
)
