#' @title Write a roleModel, roleExperiment, or roleParams to disk.
#' @param x The object to write.
#' @param dir Directory (path) to write the object to.
#' @param filename The name of the file to save.
#' File will be called *filename*.roleexperiment. 
#' @rdname writeRoLE
#' @export

setGeneric('writeRole', 
           def = function(x, dir, filename, save_txt) standardGeneric('writeRole'), 
           signature = 'x')


# method for roleData
#' @rdname writeRole
#' @export

setMethod('writeRole', 
          signature = 'roleExperiment', 
          definition = function(x, dir, filename, save_txt) { 
              if(is.null(dir)){
                  dir = getwd()
              }
              
              if(class(x)[1] == "roleExperiment"){
                  extension <- ".roleexperiment"
                  info_tag <- ".roleexperimentinfo"
                  title_line <- "Metadata for RoLE experiment R object - load into R using readRDS(file_location_and_name)"
                  info_line <- paste("Contains", length(x@modelRuns), "model runs")
                  author_line <- paste("Author:", x@authorMeta["author"])
                  date_line <- paste("Date:", x@authorMeta["date"])
                  desc_line <- paste("Description:", x@authorMeta["description"])
              }
              else{
              }
              
              saveRDS(x,paste0(dir,"/",filename,extension))
              if(save_txt){
                  con <- file(paste0(dir, "/", filename, info_tag, ".txt"))
                  writeLines(c(title_line,info_line,author_line,date_line,desc_line), con)
                  close(con)
              }
          }
)
