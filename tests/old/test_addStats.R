library(testthat)

context('add hill stats works')

# works
test_that('compute hill works', {
  # load test run role experiment
  test_exp <- readRDS("data/example_out_in/test.roleexperiment")
  data <- test_exp@modelRuns[[1]]@timeseries[[1]]
  x  <- data@localComm@abundanceSp
  entropies <- c(3,4,5)
  traits <- data@localComm@traitsSp
  print(computeHill(x,type="abundance",entropies=entropies))
  print(computeHill(x,traits=traits,type="trait",entropies=entropies))
})

# works
test_that('addHillStats to exp works', {
  # load test run role experiment
  test_exp <- readRDS("data/example_out_in/test.roleexperiment")
  entrScales <- c(3,4,5)
  test_exp <- addHillStats(test_exp,entrScales=entrScales)
  print(test_exp@modelRuns[[1]]@timeseries[[1]]@stats)
})