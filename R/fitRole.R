
fit <- function(input, simrun)
{
  simrun <- sim
  sim$params$values
  
  ts_values <- data.frame(matrix(ncol = 40, nrow = length(simrun$timeseries)))
  
  #each row is a timestep of a sim
  for(i in 1:length(simrun$timeseries))
  {
    s <- simrun$timeseries[i]
    # make row of params and summary stats
    row <- c(s$params$values[1],s$params$values[2])
    ts_values[i] <- row #add to data
  }
  
  #rf <- randomForest(param ~ stat1 + stat2 + stat3 ... , data=ts_values, ..., subset, na.action=na.fail, test = input$param)
  #rf <- out 
}

