library(testthat)

context('several model runs work')

# works
test_that('experiment of 1 run and 100 iters works', {
  experiment <- dummyModel(R=TRUE,run=TRUE,fill_ts=TRUE,niters=100,return_experiment=TRUE)
  experiment@modelRuns[[1]]@timeseries[[10]]
})