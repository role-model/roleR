
# returns a random forest to predict a parameter value using summary stats 
createPredModel <- function(expr,paramName)
{
  # each row of all_stats_params will be a timestep iteration of a model, each column is a param
  param_all_stats <- getSumStats(experiment, funs = list(hillAbund = hillAbund, 
                                                          hillTrait = hillTrait))
  param_vect <- c()
  expr@modelRuns
  for(m in 1:length(expr@modelRuns))
  {
      param_vect <- c(param_vect, expr@allParams[[m]]@comp_sigma)
  }
  param_all_stats <- cbind(param_all_stats,param_vect)
  
  rf <- randomForest(param ~ hillAbund_1 + hillAbund_2 + hillAbund_3 + hillAbund_4 + hillTrait_1 + hillTrait_2 + 
                         hillTrait_3 + hillTrait_4, data=param_all_stats, na.action=na.omit)
  return(rf)
}