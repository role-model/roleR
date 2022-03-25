
library(testthat)
library(roleR)

context('R helpers work')

# works
test_that('dummyModel (simple model creation) works', {
  
  model <- dummyModel(R=FALSE,run=TRUE)
  
  # fails sometimes with probs must be finite and non negative in iterSim
})

# works
test_that('getTimeseries works', {
  
  model <- dummyModel(R=FALSE,run=TRUE)
  
  # the individual at the dead index should now be alive
  expect_equal(l$abundance_indv[1], 1)
  
  # the species id at the dead index should now be the species of the birth
  expect_equal(l$species_ids[1], 1)
  
  # the trait at the dead index should now be different
  expect_equal(l$traits[1], 2)
  
  # the species of the individual that gave birth (1) should be incremented by 1
  expect_equal(l$abundance_sp[2], 2)
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