#' @title createRolePredModel
#' 
#' @description create a random forest model to predict a parameter value using summary statistics

#' @param exp a `roleExperiment` to use models from for training and validation
#' @param predParam the parameter to predict 
#'     
#' @rdname createPredModel
#' @export

# returns a random forest to predict a parameter value using summary stats 
createRolePredModel <- function(exp,predParamName)
{
  # each row of all_stats_params will be a timestep iteration of a model, each column is a param
  param_all_stats <- getSumStats(exp, funs = list(hillAbund = hillAbund, 
                                                          hillTrait = hillTrait))
  param_vect <- c()
  # exp@modelRuns
  for(m in 1:length(exp@modelRuns))
  {
      param_vect <- c(param_vect, exp@allParams[[m]]@comp_sigma)
  }
  param_all_stats <- cbind(param_all_stats,param_vect)
  
  rf <- randomForest(param ~ hillAbund_1 + hillAbund_2 + hillAbund_3 + hillAbund_4 + hillTrait_1 + hillTrait_2 + 
                         hillTrait_3 + hillTrait_4, data=param_all_stats, na.action=na.omit)
  
  # gradient boosting
  
  # babys first machine learning
  
  # model assembly
  
  # model selection 
  
  # generate 10000 sims for neutral comand filtering 
  
  # classify with trained model 
  
  # model 3 model
  
  # continuous and discrete 
  
  # caret 
  
  # parameter estimation 
  
  return(rf)
}
