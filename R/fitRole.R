
fitRole <- function(input, role)
{
  # each row of all_stats_params will be a timestep of a model
  all_stats_params <- data.frame(matrix(ncol = 40, nrow = length(role@timeseries)))

  counter <- 1
  # for each model run
  for(r in 1:length(role@modelRuns))
  {
    model <- role@modelRuns[[r]]
    for(i in 1:length(model@timeseries))
    {
      # get summary stats df row
      stats <- model@timeseries[[1]]@stats[i,]
      
      # get params
      params <- model@paramValues
      
      # make row of needed params
      params <- c(params$speciation_local[i],
                  params$speciation_meta[i],
                  params$extinction_meta[i],
                  params$trait_sigma[i],
                  params$env_sigma[i],
                  params$comp_sigma[i],
                  params$dispersal_prob[i])
      params <- t(data.frame(params))
      rownames(params) <- NULL
      colnames(params) <- c("speciation_local","speciation_meta","extinction_meta","trait_sigma","env_sigma","comp_sigma","dispersal_prob")
      
      # make row of params and summary stats
      row <- cbind(stats,params)
      all_stats_params[counter] <- row #add to data
      counter <- counter + 1
    }
  }
  
  #rf <- randomForest(comp_sigma ~ stat1 + stat2 + stat3, data=all_stats_params, na.action=na.fail, test = input$)
  #return(rf)
}

