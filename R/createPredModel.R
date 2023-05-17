#' @title createRolePredModel
#' 
#' @description Creates a random forest model to predict a parameter value using summary statistics of simulation experiments.

#' @param expr A `roleExperiment` that contains models with simulation data to use for training and validation,
#' @param pred_param_name The name of the parameter to predict. See the `roleParams` documentation for a list of possible parameter choices. 
#' @include roleParams.R
#' @examples 
#' m <- quickExp()
#' p <- quickParams()
#' expr <- roleExperiment(repS4(p,100))
#' expr <- runRole(expr)
#' rf <- createRolePredModel(expr,pred_param_name="env_sigma")

# returns a random forest to predict a parameter value using summary stats 
createRolePredModel <- function(expr,pred_param_name)
{
  # each row of all_stats_params will be a timestep iteration of a model, each column is a param
  all_stats <- getSumStats(expr)[,1:16]
  is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
  all_stats[is.nan(all_stats)] <- 0

  # initialize empty vector of param to predict (could lapply this)
  pred_param_vect <- c()
  # for each model in expr
  for(m in 1:length(expr@modelRuns))
  {
      pred_param_vect <- c(pred_param_vect, slot(expr@allParams[[m]],pred_param_name))
  }
    
  all_stats_and_param <- cbind(all_stats,pred_param_vect)
  names(all_stats_and_param)[17] <- "pred_param"

  rf <- ranger::ranger(pred_param ~ hill_abund_1 + hill_abund_2 + hill_abund_3 + hill_abund_4 + 
                         hill_gen_1 + hill_gen_2 + hill_gen_3 + hill_gen_4 + 
                         hill_trait_1 + hill_trait_2 + hill_trait_3 + hill_trait_4 + 
                         hill_phy_1 + hill_phy_2 + hill_phy_3 + hill_phy_4,
                         data=all_stats_and_param)
  return(rf)
}
