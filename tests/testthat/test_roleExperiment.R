library(testthat)

context('several model runs work')

# works
test_that('experiment of 1 run and 100 iters works', {
  experiment <- roleR::dummyModel(R=TRUE,run=TRUE,fill_ts=FALSE,niters=100,return_experiment=TRUE,print=TRUE)
  experiment@modelRuns[[1]]@timeseries[[1]]@localComm@abundanceSp
  experiment@modelRuns[[1]]@timeseries[[10]]@localComm@abundanceSp
})

#ecolottery contains list of dataframes