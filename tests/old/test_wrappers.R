library(testthat)

context('R wrapper constructors and conversions work')

# works
test_that('local works', {
  
  local <- localComm(dataCpp$local$abundance_indv,dataCpp$local$species_ids,
                     dataCpp$local$traits,dataCpp$local$abundance_sp,
                     dataCpp$local$traits_sp,dataCpp$local$pi_sp)
})