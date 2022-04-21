library(testthat)

context('several model runs work')

# works
test_that('experiment of 1 run and 100 iters works', {
  experiment <- roleR::dummyModel(R=TRUE,run=TRUE,fill_ts=TRUE,niters=100,return_experiment=TRUE)
  getTimeseries(experiment@modelRuns[[1]],valueName="abundanceSp")
  #writeRoleExperiment(experiment,fileName="test")
  #test <- readRDS("test.roleexperiment")
  #experiment@modelRuns[[1]]@timeseries[[10]]
})

#ecolottery contains list of dataframes