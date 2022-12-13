
library(testthat)
library(roleR)

context('R helpers work')

# works
test_that('dummyModel (simple model creation) works', {
  
  model <- dummyModel(R=FALSE,run=TRUE)
  
  # fails sometimes with probs must be finite and non negative in iterSim
})

# works
test_that('roleDataFromCpp works', {
  
  model <- dummyModel(R=FALSE,run=TRUE)
  data <- roleDataFromCpp(model$timeseries[[1]])
  
  # no errors
})

# works
test_that('roleModelFromCpp works', {
  
  model <- dummyModel(R=FALSE,run=TRUE)
  model <- roleModelFromCpp(model)
  
  # no errors
})

#todo - continue testing helper functions, this will get the core roleSim working
# get timeseries working for gene sim