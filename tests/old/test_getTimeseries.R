# works
test_that('getTimeseries works', {
  # new get timeseries
  library(roleR)
  role <- readRDS("data/example_out_in/test.roleexperiment")
  model <- role@modelRuns[[1]]
  ts <- getTimeseries(model,getAbundances,TRUE)
  
  entrScales <- c(3,4,5)
  test_exp <- addHillStats(test_exp,entrScales=entrScales)
  test_exp@modelRuns[[1]]@timeseries[[1]]@stats
  ts <- getTimeseries(test_exp,runNum=1,type="summary_stat",valueName="abundance_5")
  print(ts)
  ts <- getTimeseries(test_exp,runNum=1,type="summary_stat")
  print(ts)
  ts <- getTimeseries(test_exp,runNum=1,type="model_value")
  test_exp@modelRuns[[1]]@timeseries[[1]]@localComm@abundanceSp
  test_exp@modelRuns[[1]]@timeseries[[10]]@localComm@abundanceSp
})
