#' @title Write a RoLE model, experiment, or params
#' @param x the object to write

#' @rdname writeRoLE
#' @export

setGeneric('writeRoLE', 
           def = function(x, dir, filename, save_txt) standardGeneric('writeRoLE'), 
           signature = 'x')


# method for roleData
#' @rdname writeRoLE
#' @export

setMethod('writeRoLE', 
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
                  author_line <- paste("Author:", x@auxMeta["author"])
                  date_line <- paste("Date:", x@auxMeta["date"])
                  desc_line <- paste("Description:", x@auxMeta["description"])
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
